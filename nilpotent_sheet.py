"""Two claims behind 'the sheet is where the local map is silent and the pressure is the whole dynamics', measured.

Along a path dA/dt = -A^2 - H. On a pure shear (a vortex sheet, locally) A is nilpotent: A^2 = 0, so the local term
vanishes and dA/dt = -H_dev entirely. Registered before running, on the high-vorticity set |w| > 0.5 max:
    nil  = <|A^2|> / <|A|^2>      nilpotency: 0 for shear, ~0.5 for a rotation core     prediction: sheet <= 0.2, classical >= 0.4
    dom  = <|H_dev|> / <|A^2|>    pressure dominance over the local term                 prediction: sheet >> 1, classical ~ 1 or below
plus the Q-R location <Q>/|A|^2 (0 = shear, > 0 rotation, < 0 strain) and the sign of xi.H_dev.xi (negative = the
pressure pushes with the stretching). usage: N=48 T=1.0 python nilpotent_sheet.py     (CPU, minutes)"""
import os, glob, time, numpy as np

N = int(os.environ.get("N", 48)); T = float(os.environ.get("T", 1.0))
fft, ifft = np.fft.fftn, np.fft.ifftn
k = np.fft.fftfreq(N, d=1.0 / N)
kx, ky, kz = np.meshgrid(k, k, k, indexing="ij"); K = [kx, ky, kz]
k2 = kx**2 + ky**2 + kz**2; k2s = k2.copy(); k2s[0, 0, 0] = 1.0
deal = (np.abs(kx) < N / 3) & (np.abs(ky) < N / 3) & (np.abs(kz) < N / 3)
x = np.linspace(0, 2 * np.pi, N, endpoint=False); X, Y, Z_ = np.meshgrid(x, x, x, indexing="ij")


def project(F):
    kd = sum(K[i] * F[i] for i in range(3)) / k2s
    return [F[i] - K[i] * kd for i in range(3)]


def transport(U):
    Ud = [Ui * deal for Ui in U]; u = [ifft(Ui).real for Ui in Ud]; out = []
    for i in range(3):
        adv = sum(u[j] * ifft(1j * K[j] * Ud[i]).real for j in range(3))
        div = sum(ifft(1j * K[j] * fft(u[j] * u[i])).real for j in range(3))
        out.append(-0.5 * fft(adv + div) * deal)
    return project(out)


def step(U, dt):
    a = transport(U); b = transport([U[i] + dt / 2 * a[i] for i in range(3)])
    c = transport([U[i] + dt / 2 * b[i] for i in range(3)]); d = transport([U[i] + dt * c[i] for i in range(3)])
    return [U[i] + dt / 6 * (a[i] + 2 * b[i] + 2 * c[i] + d[i]) for i in range(3)]


def diag(U):
    Ud = [Ui * deal for Ui in U]
    A = np.stack([np.stack([ifft(1j * K[i] * Ud[j]).real for j in range(3)], -1) for i in range(3)], -1)
    src = fft(np.einsum("...ij,...ji->...", A, A)); ph = src / k2s * deal
    H = np.stack([np.stack([ifft(-K[i] * K[j] * ph).real for j in range(3)], -1) for i in range(3)], -1)
    Hd = H - np.einsum("...ii->...", H)[..., None, None] / 3 * np.eye(3)
    w = np.stack([A[..., 1, 2] - A[..., 2, 1], A[..., 2, 0] - A[..., 0, 2], A[..., 0, 1] - A[..., 1, 0]], -1)
    wm = np.linalg.norm(w, axis=-1); high = wm > 0.5 * wm.max()
    Ah, Hdh, xi = A[high], Hd[high], (w[high] / (wm[high][:, None] + 1e-30))
    A2 = np.einsum("...ij,...jk->...ik", Ah, Ah)
    nA, nA2, nH = (np.linalg.norm(Ah, axis=(1, 2)), np.linalg.norm(A2, axis=(1, 2)), np.linalg.norm(Hdh, axis=(1, 2)))
    Q = -0.5 * np.einsum("...ij,...ji->...", Ah, Ah)
    xiHxi = np.einsum("...i,...ij,...j->...", xi, Hdh, xi)
    Z = 0.5 * np.mean(wm**2)
    return Z, nA2.mean() / (nA**2).mean(), nH.mean() / nA2.mean(), (Q / nA**2).mean(), (xiHxi < 0).mean(), xiHxi.mean() / (nA**2).mean()


flows = {}
p = os.environ.get("FOUND", "results/found/leashed64_dmin030.npz")
uf = np.load(p)["u"].astype(float); n0 = uf.shape[1]; U = []
for c in range(3):
    uh = fft(uf[c]) * (N / n0) ** 3; big = np.zeros((N, N, N), complex); h = n0 // 2
    for a in (slice(0, h), slice(-h, None)):
        for b in (slice(0, h), slice(-h, None)):
            for cc in (slice(0, h), slice(-h, None)):
                big[a, b, cc] = uh[a, b, cc]
    U.append(big)
flows["sheet (searcher's field)"] = project(U)
flows["kida-pelz"] = [fft(np.sin(X) * (np.cos(3 * Y) * np.cos(Z_) - np.cos(Y) * np.cos(3 * Z_))), fft(np.sin(Y) * (np.cos(3 * Z_) * np.cos(X) - np.cos(Z_) * np.cos(3 * X))), fft(np.sin(Z_) * (np.cos(3 * X) * np.cos(Y) - np.cos(X) * np.cos(3 * Y)))]
flows["taylor-green"] = [fft(np.sin(X) * np.cos(Y) * np.cos(Z_)), fft(-np.cos(X) * np.sin(Y) * np.cos(Z_)), fft(np.zeros_like(X))]
flows["pure shear (control: A^2 = 0)"] = [fft(np.sin(Y)), fft(np.zeros_like(X)), fft(np.zeros_like(X))]
print("N=%d^3, nu=0, Z0=0.375, high set |w| > 0.5 max.  nil = <|A^2|>/<|A|^2>, dom = <|H_dev|>/<|A^2|>" % N)
print("%-30s %5s   %7s   %7s   %7s   %9s   %12s   %s" % ("flow", "t", "Z/Z0", "nil", "dom", "<Q>/|A|2", "frac xiHxi<0", "<xiHxi>/|A|2"))
for name, U in flows.items():
    Zi = 0.5 * np.mean(sum(vi**2 for vi in [ifft(1j * K[a] * U[b] - 1j * K[b] * U[a]).real for a, b in ((1, 2), (2, 0), (0, 1))]))
    U = [Ui * np.sqrt(0.375 / Zi) for Ui in U]; Z0 = 0.375
    t, mark, t0 = 0.0, 0.0, time.time()
    while t <= T + 1e-9:
        if t >= mark - 1e-9:
            Z, nil, dom, q, fneg, hh = diag(U)
            print("%-30s %5.2f   %7.3f   %7.3f   %7.2f   %+9.3f   %12.3f   %+.3f   (%.0fs)" % (name, t, Z / Z0, nil, dom, q, fneg, hh, time.time() - t0), flush=True)
            mark += 0.5
            if t >= T - 1e-9 or "shear" in name: break
        umax = max(np.abs(ifft(Ui).real).max() for Ui in U)
        dt = min(2.0 / N, 0.5 * (2 * np.pi / N) / max(umax, 1e-9), mark - t + 1e-12)
        U = step(U, dt); t += dt
