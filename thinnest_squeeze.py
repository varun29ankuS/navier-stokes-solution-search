"""The thinnest squeeze: the smallest singular value of the Lagrangian deformation gradient, against the jitter scale.

Every collapse - sheet, tube or point - is a material element pressed to zero thickness. Along a particle path the
deformation gradient F = grad_a X obeys dF/dt = A(X(t)) F with det F = 1 (Liouville), and its smallest singular value
sigma_min(F) is the thinnest that parcel has been pressed. Blow-up <=> sigma_min -> 0 somewhere in finite time. In the
Constantin-Iyer picture viscosity is jitter smearing labels over sqrt(2 nu t); the natural (unproved) inequality is
        sigma_min(F)  >~  sqrt(2 nu t)          "nothing can be squeezed thinner than the jitter smears it"
which would give |grad u| <~ 1/sqrt(nu t), the Type I scaling. Measured here along the searcher's sheet and Kida-Pelz
under Navier-Stokes (nu as given): particles seeded on the high-vorticity set, F integrated with A interpolated from
the spectral field. Reported: min and 1st percentile of sigma_min over particles, the jitter scale, their ratio, and
max|w| (for a sheet 1/sigma_min ~ max|w| x const). usage: IC=found|kp N=48 NU=2e-3 T=1.0 python thinnest_squeeze.py"""
import os, time, numpy as np

IC = os.environ.get("IC", "found"); N = int(os.environ.get("N", 48)); NU = float(os.environ.get("NU", 2e-3)); T = float(os.environ.get("T", 1.0))
NP = int(os.environ.get("NP", 4000)); DT = float(os.environ.get("DT", 0.01))
fft, ifft = np.fft.fftn, np.fft.ifftn
k = np.fft.fftfreq(N, d=1.0 / N); kx, ky, kz = np.meshgrid(k, k, k, indexing="ij"); K = [kx, ky, kz]
k2 = kx**2 + ky**2 + kz**2; k2s = k2.copy(); k2s[0, 0, 0] = 1.0
deal = (np.abs(kx) < N / 3) & (np.abs(ky) < N / 3) & (np.abs(kz) < N / 3)
L = 2 * np.pi; x = np.linspace(0, L, N, endpoint=False); X, Y, Z_ = np.meshgrid(x, x, x, indexing="ij"); dx = L / N


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
    f, f2 = np.exp(-NU * k2 * dt), np.exp(-NU * k2 * dt / 2)
    a = transport(U); b = transport([f2 * (U[i] + dt / 2 * a[i]) for i in range(3)])
    c = transport([f2 * U[i] + dt / 2 * b[i] for i in range(3)]); d = transport([f * U[i] + dt * f2 * c[i] for i in range(3)])
    return [f * U[i] + dt / 6 * (f * a[i] + 2 * f2 * b[i] + 2 * f2 * c[i] + d[i]) for i in range(3)]


def fields(U):
    Ud = [Ui * deal for Ui in U]
    u = [ifft(Ui).real for Ui in Ud]
    A = [[ifft(1j * K[i] * Ud[j]).real for j in range(3)] for i in range(3)]      # A[i][j] = d_i u_j
    return u, A


def interp(F, px, py, pz):
    g = np.stack([px, py, pz]) / dx; i0 = np.floor(g).astype(int); f = g - i0; i0 %= N; i1 = (i0 + 1) % N
    out = 0.0
    for a in (0, 1):
        for b in (0, 1):
            for c in (0, 1):
                w = (f[0] if a else 1 - f[0]) * (f[1] if b else 1 - f[1]) * (f[2] if c else 1 - f[2])
                out = out + w * F[(i1 if a else i0)[0], (i1 if b else i0)[1], (i1 if c else i0)[2]]
    return out


if IC == "kp":
    U = [fft(np.sin(X) * (np.cos(3 * Y) * np.cos(Z_) - np.cos(Y) * np.cos(3 * Z_))), fft(np.sin(Y) * (np.cos(3 * Z_) * np.cos(X) - np.cos(Z_) * np.cos(3 * X))), fft(np.sin(Z_) * (np.cos(3 * X) * np.cos(Y) - np.cos(X) * np.cos(3 * Y)))]
else:
    p = os.environ.get("FOUND", "results/found/leashed64_dmin030.npz"); uf = np.load(p)["u"].astype(float); n0 = uf.shape[1]; U = []
    for c in range(3):
        uh = fft(uf[c]) * (N / n0) ** 3; big = np.zeros((N, N, N), complex); h = n0 // 2
        for a in (slice(0, h), slice(-h, None)):
            for b in (slice(0, h), slice(-h, None)):
                for cc in (slice(0, h), slice(-h, None)):
                    big[a, b, cc] = uh[a, b, cc]
        U.append(big)
    U = project(U)
w = lambda U: [ifft(1j * ky * U[2] - 1j * kz * U[1]).real, ifft(1j * kz * U[0] - 1j * kx * U[2]).real, ifft(1j * kx * U[1] - 1j * ky * U[0]).real]
Zi = 0.5 * np.mean(sum(wi**2 for wi in w(U))); U = [Ui * np.sqrt(0.375 / Zi) for Ui in U]
wm = np.sqrt(sum(wi**2 for wi in w(U))); high = np.argwhere(wm > 0.5 * wm.max())
rng = np.random.default_rng(0); sel = high[rng.choice(len(high), min(NP, len(high)), replace=False)]
px, py, pz = (sel[:, 0] * dx + dx / 2, sel[:, 1] * dx + dx / 2, sel[:, 2] * dx + dx / 2)
Fp = np.tile(np.eye(3), (len(px), 1, 1))
print("IC=%s  N=%d^3  nu=%g  %d particles seeded on the high-vorticity set; F integrated along paths (dF/dt = A F, det F = 1)" % (IC, N, NU, len(px)))
print("   t     max|w|    det F (mean, Liouville=1)    sigma_min: min      1st pct     median     jitter sqrt(2 nu t)    min/jitter    1/sigma_min_min")
t, mark, t0 = 0.0, 0.0, time.time()
while t <= T + 1e-9:
    if t >= mark - 1e-9:
        sv = np.linalg.svd(Fp, compute_uv=False); smin = sv[:, -1]; det = np.linalg.det(Fp)
        jit = np.sqrt(2 * NU * t) if t > 0 else 0.0
        print("%5.2f   %7.2f   %.4f                     %.4f    %.4f     %.4f     %.4f                %s        %.1f   (%.0fs)" % (
            t, wm.max(), det.mean(), smin.min(), np.percentile(smin, 1), np.median(smin), jit, ("%.2f" % (smin.min() / jit)) if t > 0 else "  -  ", 1 / smin.min(), time.time() - t0), flush=True)
        mark += 0.25
        if t >= T - 1e-9: break
    u, A = fields(U)
    # RK2 for positions and F, with A interpolated at the particles
    def rates(qx, qy, qz, Fq):
        vx, vy, vz = interp(u[0], qx, qy, qz), interp(u[1], qx, qy, qz), interp(u[2], qx, qy, qz)
        Ap = np.stack([np.stack([interp(A[i][j], qx, qy, qz) for j in range(3)], -1) for i in range(3)], -2)   # Ap[...,i,j] = d_i u_j
        return vx, vy, vz, np.einsum("nji,njk->nik", Ap, Fq)                                                # dF/dt = (grad u)^T F, F_ik = dX_i/da_k
    vx, vy, vz, dF = rates(px, py, pz, Fp)
    mx, my, mz, Fm = (px + DT / 2 * vx) % L, (py + DT / 2 * vy) % L, (pz + DT / 2 * vz) % L, Fp + DT / 2 * dF
    vx, vy, vz, dF = rates(mx, my, mz, Fm)
    px, py, pz, Fp = (px + DT * vx) % L, (py + DT * vy) % L, (pz + DT * vz) % L, Fp + DT * dF
    U = step(U, DT); wm = np.sqrt(sum(wi**2 for wi in w(U))); t += DT
