"""The quadratic map inside Navier-Stokes, and the fraction of the flow the pressure holds inside the set.

Along a particle path the velocity gradient A = grad u obeys  dA/dt = -A^2 - H  (plus viscosity), H the pressure Hessian:
z -> z^2 + c in matrix form. Restricted Euler (Vieillefosse 1982; Cantwell 1992) keeps only the local part of H (its
trace, lap p = |w|^2/2 - |S|^2) and escapes to a finite-time singularity - the quadratic map running away. The full
equation has the nonlocal part of H, the global pressure, and stays bounded: the pressure Hessian is the "c" that
keeps the local orbit inside the set. This measures, on the high-vorticity set of a flow, the fraction of points whose
local map ESCAPES (|A| exceeds ESC x its initial value within a turnover) when H is frozen at
    (a) its local part only     (restricted Euler)
    (b) the full Hessian        (frozen at the current time)
and reports the RESCUED fraction, (a escapes) and not (b escapes): the part of the flow living on the global pressure.
Also the Q-R invariants (Q = -1/2 tr A^2, R = -1/3 tr A^3) and where the high set sits relative to the Vieillefosse tail
27 R^2 + 4 Q^3 = 0. usage: IC=found|tg|kp N=64 T=1.0 python mandelbrot_fraction.py"""
import os, glob, time, numpy as np

IC = os.environ.get("IC", "found")
N = int(os.environ.get("N", 64))
T = float(os.environ.get("T", 1.0))
ESC = float(os.environ.get("ESC", 10.0))
HORIZON = float(os.environ.get("HORIZON", 1.0))       # iterate the local map for one turnover of the flow's own time scale
fft, ifft = np.fft.fftn, np.fft.ifftn
k = np.fft.fftfreq(N, d=1.0 / N)
kx, ky, kz = np.meshgrid(k, k, k, indexing="ij")
K = [kx, ky, kz]
k2 = kx**2 + ky**2 + kz**2
k2s = k2.copy(); k2s[0, 0, 0] = 1.0
deal = (np.abs(kx) < N / 3) & (np.abs(ky) < N / 3) & (np.abs(kz) < N / 3)
x = np.linspace(0, 2 * np.pi, N, endpoint=False)
X, Y, Z_ = np.meshgrid(x, x, x, indexing="ij")


def project(F):
    kd = sum(K[i] * F[i] for i in range(3)) / k2s
    return [F[i] - K[i] * kd for i in range(3)]


def transport(U):
    Ud = [Ui * deal for Ui in U]
    u = [ifft(Ui).real for Ui in Ud]
    out = []
    for i in range(3):
        adv = sum(u[j] * ifft(1j * K[j] * Ud[i]).real for j in range(3))
        div = sum(ifft(1j * K[j] * fft(u[j] * u[i])).real for j in range(3))
        out.append(-0.5 * fft(adv + div) * deal)
    return project(out)


def step(U, dt):
    a = transport(U); b = transport([U[i] + dt / 2 * a[i] for i in range(3)])
    c = transport([U[i] + dt / 2 * b[i] for i in range(3)]); d = transport([U[i] + dt * c[i] for i in range(3)])
    return [U[i] + dt / 6 * (a[i] + 2 * b[i] + 2 * c[i] + d[i]) for i in range(3)]


def gradient_and_hessian(U):
    Ud = [Ui * deal for Ui in U]
    A = np.stack([np.stack([ifft(1j * K[i] * Ud[j]).real for j in range(3)], -1) for i in range(3)], -1)   # A[...,i,j] = d_i u_j
    src = fft(np.einsum("...ij,...ji->...", A, A))
    ph = src / k2s * deal
    H = np.stack([np.stack([ifft(-K[i] * K[j] * ph).real for j in range(3)], -1) for i in range(3)], -1)
    return A, H


def escapes(A0, Hf, dt, nsteps):
    """iterate the local map dA/dt = -A^2 - H with H frozen; escape = |A| > ESC |A0| at any step"""
    A = A0.copy(); n0 = np.sqrt(np.einsum("...ij,...ij->...", A0, A0)) + 1e-30
    esc = np.zeros(A.shape[:-2], bool)
    for _ in range(nsteps):
        A = A - dt * (np.einsum("...ij,...jk->...ik", A, A) + Hf)
        n = np.sqrt(np.einsum("...ij,...ij->...", A, A))
        esc |= ~np.isfinite(n) | (n > ESC * n0)
        A = np.where(esc[..., None, None], 0.0, A)
    return esc


if IC == "found":
    path = os.environ.get("FOUND", "results/found/leashed64_dmin030.npz")
    uf = np.load(path)["u"].astype(float); n0 = uf.shape[1]
    U = []
    for c in range(3):
        uh = fft(uf[c]) * (N / n0) ** 3; big = np.zeros((N, N, N), complex); h = n0 // 2
        for a in (slice(0, h), slice(-h, None)):
            for b in (slice(0, h), slice(-h, None)):
                for cc in (slice(0, h), slice(-h, None)):
                    big[a, b, cc] = uh[a, b, cc]
        U.append(big)
    U = project(U)
elif IC == "kp":
    U = [fft(np.sin(X) * (np.cos(3 * Y) * np.cos(Z_) - np.cos(Y) * np.cos(3 * Z_))),
         fft(np.sin(Y) * (np.cos(3 * Z_) * np.cos(X) - np.cos(Z_) * np.cos(3 * X))),
         fft(np.sin(Z_) * (np.cos(3 * X) * np.cos(Y) - np.cos(X) * np.cos(3 * Y)))]
else:
    U = [fft(np.sin(X) * np.cos(Y) * np.cos(Z_)), fft(-np.cos(X) * np.sin(Y) * np.cos(Z_)), fft(np.zeros_like(X))]
w = lambda U: [ifft(1j * ky * U[2] - 1j * kz * U[1]).real, ifft(1j * kz * U[0] - 1j * kx * U[2]).real, ifft(1j * kx * U[1] - 1j * ky * U[0]).real]
Zi = 0.5 * np.mean(sum(wi**2 for wi in w(U)))
U = [Ui * np.sqrt(0.375 / Zi) for Ui in U]
print("IC=%s  N=%d^3  nu=0  Z0=0.375; local map dA/dt = -A^2 - H iterated for %.1f x (1/max|A|) with H frozen; escape = |A| > %g |A0|" % (IC, N, HORIZON, ESC))
print("   t     Z/Z0     high-set size   escapes: local H only   full H   RESCUED by the global pressure   <Q>/|A|^2  <R>/|A|^3  frac past Vieillefosse tail")
t, mark, t0 = 0.0, 0.0, time.time()
while t <= T + 1e-9:
    if t >= mark - 1e-9:
        A, H = gradient_and_hessian(U)
        wm = np.sqrt(sum(wi**2 for wi in w(U))); high = wm > 0.5 * wm.max()
        Ah, Hh = A[high], H[high]
        trH = np.einsum("...ii->...", Hh)
        Hloc = np.eye(3)[None] * (trH / 3)[:, None, None]
        amax = np.sqrt(np.einsum("...ij,...ij->...", Ah, Ah)).max()
        dt = 0.02 / amax; nsteps = int(HORIZON / amax / dt)
        e_loc = escapes(Ah, Hloc, dt, nsteps); e_full = escapes(Ah, Hh, dt, nsteps)
        Q = -0.5 * np.einsum("...ij,...ji->...", Ah, Ah); R = -np.einsum("...ij,...jk,...ki->...", Ah, Ah, Ah) / 3
        a2 = np.einsum("...ij,...ij->...", Ah, Ah)
        disc = 27 * R**2 + 4 * Q**3
        print("%5.2f   %.3f     %6d          %.3f                  %.3f    %.3f                             %+.3f     %+.3f     %.3f   (%.0fs)" % (
            t, 0.5 * np.mean(wm**2) / 0.375, high.sum(), e_loc.mean(), e_full.mean(), (e_loc & ~e_full).mean(), (Q / a2).mean(), (R / a2**1.5).mean(), (disc > 0).mean(), time.time() - t0), flush=True)
        mark += 0.25
        if t >= T - 1e-9:
            break
    umax = max(np.abs(ifft(Ui).real).max() for Ui in U)
    dtf = min(2.0 / N, 0.5 * (2 * np.pi / N) / max(umax, 1e-9), mark - t + 1e-12)
    U = step(U, dtf); t += dtf
