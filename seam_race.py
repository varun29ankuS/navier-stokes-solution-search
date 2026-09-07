"""The seam race: does viscosity cut the twist before the roll-up completes?

The searcher's fastest field is built of antiparallel vortex sheets - the vorticity direction reverses across a seam
where |w| passes through zero. Squared coherence measures cannot see the reversal (a line field reads a Mobius strip as
a cylinder); the dynamics can: the seam is where Biot-Savart cancels (depletion), where Kelvin-Helmholtz rolls the
pair up, and - under Navier-Stokes - where RECONNECTION cuts and rejoins the lines. Hussain's group found reconnection
relieves the collapse in Kerr's antiparallel-tube setup. Measured here on the strong-vorticity set (|w| > 0.5 max):
    anti    fraction with an antiparallel neighbour two cells away (beta = 1 - xi.xi' > 1): the twist
    anti1   fraction with an antiparallel neighbour ONE cell away: the seam at the grid limit
    annih   fraction of twisted points where viscous cancellation exceeds stretching, -nu w.lap w > w.S.w: reconnection winning locally
    Z/Z0, max|w|, the strip delta with the clock.
REGISTERED: Euler - anti rises to ~0.9 and stays to the clock. Navier-Stokes - annih rises to O(1) before anti
saturates, anti then falls and Z peaks (the cut wins). KILL for the reconnection story: at nu = 2e-3 anti keeps rising
to the clock with annih still small (the roll-up outruns the cut).
usage: NU=2e-3 N=64 T=2.0 python seam_race.py"""
import os, time, numpy as np

N = int(os.environ.get("N", 64)); NU = float(os.environ.get("NU", 2e-3)); T = float(os.environ.get("T", 2.0)); EVERY = float(os.environ.get("EVERY", 0.1))
fft, ifft = np.fft.fftn, np.fft.ifftn
k = np.fft.fftfreq(N, d=1.0 / N); kx, ky, kz = np.meshgrid(k, k, k, indexing="ij"); K = [kx, ky, kz]
k2 = kx**2 + ky**2 + kz**2; k2s = k2.copy(); k2s[0, 0, 0] = 1.0
deal = (np.abs(kx) < N / 3) & (np.abs(ky) < N / 3) & (np.abs(kz) < N / 3)
KMAG = np.sqrt(k2); NB = int(N / 3); shells = [(KMAG >= n - 0.5) & (KMAG < n + 0.5) for n in range(1, NB)]


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


def strip(U):
    e = 0.5 * sum(np.abs(Ui) ** 2 for Ui in U) / N**6
    spec = np.array([e[s].sum() for s in shells]); ks = np.arange(1, NB); sel = (ks >= NB // 2) & (spec > 1e-300)
    return -np.polyfit(ks[sel], np.log(spec[sel]), 1)[0] / 2 if sel.sum() > 4 else np.nan


def diag(U):
    Ud = [Ui * deal for Ui in U]
    G = [[ifft(1j * K[i] * Ud[j]).real for j in range(3)] for i in range(3)]
    w = [G[1][2] - G[2][1], G[2][0] - G[0][2], G[0][1] - G[1][0]]
    wm = np.sqrt(sum(wi**2 for wi in w)) + 1e-30; xi = [wi / wm for wi in w]
    high = wm > 0.5 * wm.max()
    anti2 = np.zeros_like(high); anti1 = np.zeros_like(high); w2 = wm**2; soft = 0.0
    for ax in range(3):
        for sgn in (1, -1):
            b2 = 1 - sum(xi[c] * np.roll(xi[c], sgn * 2, axis=ax) for c in range(3)); b1 = 1 - sum(xi[c] * np.roll(xi[c], sgn * 1, axis=ax) for c in range(3))
            anti2 |= b2 > 1.0; anti1 |= b1 > 1.0
            if sgn == 1:                                   # strong twist: enstrophy-weighted sharp reversal at 1 and 2 cells (the searcher's penalised measure)
                soft += (w2 * np.roll(w2, 1, axis=ax) * np.maximum(b1 - 1, 0) ** 2).mean() + (w2 * np.roll(w2, 2, axis=ax) * np.maximum(b2 - 1, 0) ** 2).mean()
    soft /= (w2**2).mean()
    S = [[0.5 * (G[i][j] + G[j][i]) for j in range(3)] for i in range(3)]
    stretch = sum(w[i] * S[i][j] * w[j] for i in range(3) for j in range(3))
    wh = [fft(wi) for wi in w]
    lapw = [ifft(-k2 * whi).real for whi in wh]
    annih = -NU * sum(w[i] * lapw[i] for i in range(3))                     # viscous cancellation rate of |w|^2/2 (positive = destroying)
    tw = high & anti2
    return (0.5 * np.mean(wm**2), wm.max(), anti2[high].mean(), anti1[high].mean(),
            (annih[tw] > stretch[tw]).mean() if tw.sum() > 0 else 0.0, strip(U), soft)


p = os.environ.get("FOUND", "results/found/leashed64_dmin030.npz"); uf = np.load(p)["u"].astype(float); n0 = uf.shape[1]; U = []
for c in range(3):
    uh = fft(uf[c]) * (N / n0) ** 3; big = np.zeros((N, N, N), complex); h = n0 // 2
    for a in (slice(0, h), slice(-h, None)):
        for b in (slice(0, h), slice(-h, None)):
            for cc in (slice(0, h), slice(-h, None)):
                big[a, b, cc] = uh[a, b, cc]
    U.append(big)
U = project(U)
Z0 = diag(U)[0]; U = [Ui * np.sqrt(0.375 / Z0) for Ui in U]; Z0 = 0.375
print("seam race: searcher's sheet field, N=%d^3, nu=%g, T=%.1f; 2dx = %.4f" % (N, NU, T, 2 * 2 * np.pi / N))
print("   t     Z/Z0     max|w|    anti(2 cells)   anti(1 cell)   annih>stretch   strong twist   dlogZ/dt   delta")
t, mark, t0 = 0.0, 0.0, time.time(); prev = None
while t <= T + 1e-9:
    if t >= mark - 1e-9:
        Z, wmax, a2, a1, an, d, soft = diag(U)
        rate = (np.log(Z / prev[0]) / (t - prev[1])) if prev else float("nan")
        print("%5.2f   %6.3f   %7.2f      %.3f           %.3f          %.3f        %.5f      %+6.3f     %.3f%s   (%.0fs)" % (t, Z / Z0, wmax, a2, a1, an, soft, rate, d, "" if d > 2 * 2 * np.pi / N else "  <-- past the clock", time.time() - t0), flush=True)
        prev = (Z, t); mark += EVERY
        if t >= T - 1e-9: break
    umax = max(np.abs(ifft(Ui).real).max() for Ui in U)
    dt = min(2.0 / N, 0.5 * (2 * np.pi / N) / max(umax, 1e-9), mark - t + 1e-12)
    U = step(U, dt); t += dt
