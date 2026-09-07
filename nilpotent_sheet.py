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
    # the corrected picture's prediction: the sheet thins at the external compressive strain rate, across the sheet.
    S = 0.5 * (Ah + np.swapaxes(Ah, 1, 2)); lam, vec = np.linalg.eigh(S)                       # ascending: lam[:,0] most compressive
    gw = np.stack([ifft(1j * K[i] * fft(wm)).real for i in range(3)], -1)[high]                  # grad|w|: the sheet normal
    nrm = gw / (np.linalg.norm(gw, axis=-1, keepdims=True) + 1e-30)
    cos_n = np.abs(np.einsum("...i,...i->...", vec[..., 0], nrm))                               # |cos| between the compressive direction and the normal
    Snn = np.einsum("...i,...ij,...j->...", nrm, S, nrm)                                        # strain ACROSS the sheet: zero for the sheet's own shear, so this is the external compression that thins it
    alpha = np.einsum("...i,...ij,...j->...", xi, S, xi)                                        # exact growth rate of |w|: D|w|/Dt = |w| xi.S.xi
    tt = np.cross(xi, nrm); tt /= (np.linalg.norm(tt, axis=-1, keepdims=True) + 1e-30)          # in-plane transverse direction
    Stt = np.einsum("...i,...ij,...j->...", tt, S, tt)                                          # compression ALONG the sheet (narrowing); alpha + Snn + Stt = 0 exactly
    return Z, nA2.mean() / (nA**2).mean(), nH.mean() / nA2.mean(), (Q / nA**2).mean(), (xiHxi < 0).mean(), xiHxi.mean() / (nA**2).mean(), wm.max(), lam[:, 0].mean(), cos_n.mean(), Snn.mean(), alpha.mean(), Stt.mean()


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
print("N=%d^3, nu=0, Z0=0.375, high set |w| > 0.5 max.  nil = <|A^2|>/<|A|^2>, dom = <|H_dev|>/<|A^2|>" % N)
print("%-30s %5s   %7s   %7s   %7s   %9s   %8s   %10s   %10s   %s" % ("flow", "t", "Z/Z0", "nil", "dom", "<Q>/|A|2", "max|w|", "growth rate", "<lam_min>", "|cos(comp,normal)|"))
print("prediction (corrected): the sheet's OWN shear gives lam_min ~ -|w|/4 at 45 deg (|cos| ~ 0.7) and does not thin it; the EXTERNAL compression across the sheet, n.S.n, sets the growth rate: -<n.S.n> ~ d log max|w|/dt")
for name, U in flows.items():
    Zi = 0.5 * np.mean(sum(vi**2 for vi in [ifft(1j * K[a] * U[b] - 1j * K[b] * U[a]).real for a, b in ((1, 2), (2, 0), (0, 1))]))
    U = [Ui * np.sqrt(0.375 / Zi) for Ui in U]; Z0 = 0.375
    t, mark, t0 = 0.0, 0.0, time.time(); prev = None
    while t <= T + 1e-9:
        if t >= mark - 1e-9:
            Z, nil, dom, q, fneg, hh, wmax, lmin, cosn, snn, al, stt = diag(U)
            rate = (np.log(wmax / prev[0]) / (t - prev[1])) if prev else float("nan")
            print("%-24s t=%.2f  Z/Z0 %.3f  nil %.3f  max|w| %6.2f  growth %6s | stretch xi.S.xi %+.3f = thinning(-n.S.n) %+.3f + narrowing(-t.S.t) %+.3f   (%.0fs)" % (name[:24], t, Z / Z0, nil, wmax, ("%.3f" % rate) if prev else "-", al, -snn, -stt, time.time() - t0), flush=True)
            prev = (wmax, t)
            mark += 0.25
            if t >= T - 1e-9 or "shear" in name: break
        umax = max(np.abs(ifft(Ui).real).max() for Ui in U)
        dt = min(2.0 / N, 0.5 * (2 * np.pi / N) / max(umax, 1e-9), mark - t + 1e-12)
        U = step(U, dt); t += dt
