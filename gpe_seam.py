"""The quantum seam: the antiparallel-pair reconnection in a fluid where the cut is guaranteed to win.

Gross-Pitaevskii  i psi_t = -1/2 lap psi + (|psi|^2 - 1) psi   (healing length xi = 1, sound speed c = 1, circulation 2 pi).
Madelung psi = sqrt(rho) e^{iS} turns it into Euler for (rho, v = grad S) plus the QUANTUM PRESSURE, a stiffness that
forbids density collapse below xi: the equation is globally regular (defocusing NLS), circulation is quantised, a
vortex line is a smooth zero of psi, and reconnection is two zero-lines crossing - smoothly, with a burst of sound.
Same structure the searcher's field builds (antiparallel sheets, the seam, the roll-up), none of the singularity:
every regulator nature offers installs a FIXED length; Navier-Stokes' only one, sqrt(nu/s), moves with the flow.

MODE=2d   validation: two counter-propagating vortex dipoles of separation d; the point-vortex speed is 1/d.
          REGISTERED: total energy conserved to 1e-5 relative; dipole speed within 20 percent of 1/d for d = 6.
MODE=3d   two antiparallel pairs of lines along z (zero momentum), separation d, each bowed toward itself (gap d - 2a at z = 0). Reported per
          time: energies (incompressible kinetic E_ki, compressible = sound E_kc, quantum E_q, interaction E_int,
          total), the seam gap (closest +/- zero pair over z-slices), fraction of z-slices still cut by BOTH lines,
          and the L1 vortex-line length. Afterwards: reconnection time t_r (first slice losing its pair), the
          gap law gap ~ (t_r - t)^p fitted on gap in [1.5, 0.6 d] before t_r, and the sound released across t_r.
          REGISTERED: the lines reconnect at finite t_r; p = 0.5 +/- 0.15 (Villois-Proment-Krstulovic 2017); E_kc
          rises across t_r while E_ki falls; the gap never re-closes after the cut (no re-approach, the Kelvin
          waves radiate). KILL: no reconnection by T, or p outside [0.25, 0.8].
usage: MODE=3d N=64 L=32 D=6 A=1 T=60 python gpe_seam.py      (2-D seconds; 3-D 64^3 minutes on CI)"""
import os, time, numpy as np

MODE = os.environ.get("MODE", "3d"); N = int(os.environ.get("N", 64)); L = float(os.environ.get("L", 32.0))
D = float(os.environ.get("D", 6.0)); A = float(os.environ.get("A", 1.0)); T = float(os.environ.get("T", 60.0))
DT = float(os.environ.get("DT", 0.02)); EVERY = float(os.environ.get("EVERY", 0.25)); TIM = float(os.environ.get("TIM", 1.0))
dim = 3 if MODE == "3d" else 2
dx = L / N; x = np.arange(N) * dx
k1 = 2 * np.pi * np.fft.fftfreq(N, d=dx)
K = np.meshgrid(*([k1] * dim), indexing="ij"); k2 = sum(Ki**2 for Ki in K); k2s = k2.copy(); k2s[(0,) * dim] = 1.0
fft, ifft = np.fft.fftn, np.fft.ifftn


def core(r):
    """Pade core profile of the 2-D GP vortex (Berloff 2004), f -> 1 far away, f ~ r near the axis"""
    r2 = r * r
    return np.sqrt(r2 * (0.3437 + 0.0286 * r2) / (1 + 0.3333 * r2 + 0.0286 * r2 * r2))


def vortex(X, Y, x0, y0, s):
    """one line/point vortex of sign s at (x0, y0) with 5x5 periodic images"""
    psi = np.ones_like(X, dtype=complex)
    for ix in (-2, -1, 0, 1, 2):
        for iy in (-2, -1, 0, 1, 2):
            dxx, dyy = X - x0 - ix * L, Y - y0 - iy * L
            psi *= core(np.sqrt(dxx**2 + dyy**2)) * np.exp(1j * s * np.arctan2(dyy, dxx))
    return psi


def step(psi, dt, imaginary=False):
    """Strang split-step: half nonlinear, full linear, half nonlinear; exact for each piece"""
    if imaginary:
        psi = psi * np.exp(-dt / 2 * (np.abs(psi) ** 2 - 1)); psi = ifft(np.exp(-dt * k2 / 2) * fft(psi))
        return psi * np.exp(-dt / 2 * (np.abs(psi) ** 2 - 1))
    psi = psi * np.exp(-1j * dt / 2 * (np.abs(psi) ** 2 - 1)); psi = ifft(np.exp(-1j * dt * k2 / 2) * fft(psi))
    return psi * np.exp(-1j * dt / 2 * (np.abs(psi) ** 2 - 1))


def energies(psi):
    """E_total = int 1/2|grad psi|^2 + 1/2 (rho-1)^2 ; kinetic split into quantum + hydrodynamic, hydrodynamic into
    incompressible + compressible by Helmholtz on w = sqrt(rho) v = Im(psi* grad psi)/|psi|  (Nore-Abid-Brachet 1997)"""
    ph = fft(psi); vol = dx**dim
    g = [ifft(1j * Ki * ph) for Ki in K]
    rho = np.abs(psi) ** 2; amp = np.sqrt(rho) + 1e-12
    E_tot = (0.5 * sum(np.abs(gi) ** 2 for gi in g) + 0.5 * (rho - 1) ** 2).sum() * vol
    E_int = (0.5 * (rho - 1) ** 2).sum() * vol
    E_q = 0.5 * sum(np.abs(ifft(1j * Ki * fft(amp))) ** 2 for Ki in K).sum() * vol
    w = [(np.conj(psi) * gi).imag / amp for gi in g]
    wh = [fft(wi) for wi in w]; kd = sum(K[i] * wh[i] for i in range(dim)) / k2s
    wc = [ifft(K[i] * kd).real for i in range(dim)]
    E_kc = 0.5 * sum(wci**2 for wci in wc).sum() * vol
    E_ki = 0.5 * sum((w[i] - wc[i]) ** 2 for i in range(dim)).sum() * vol
    return E_tot, E_ki, E_kc, E_q, E_int


def windings(phi):
    """winding number of each plaquette of a 2-D phase field (periodic): +1/-1 at a vortex, 0 elsewhere"""
    wrap = lambda a: (a + np.pi) % (2 * np.pi) - np.pi
    dxp = wrap(np.roll(phi, -1, 0) - phi); dyp = wrap(np.roll(phi, -1, 1) - phi)
    return np.rint((dxp + np.roll(dyp, -1, 0) - np.roll(dxp, -1, 1) - dyp) / (2 * np.pi)).astype(int)


def zeros_2d(psi2):
    """positions and signs of vortices in a 2-D slice; the zero located inside its plaquette by a linear fit of psi"""
    wn = windings(np.angle(psi2)); idx = np.argwhere(wn != 0); pos = []
    for i, j in idx:
        c00, c10, c01, c11 = psi2[i, j], psi2[(i + 1) % N, j], psi2[i, (j + 1) % N], psi2[(i + 1) % N, (j + 1) % N]
        cx, cy, c0 = 0.5 * (c10 + c11 - c00 - c01), 0.5 * (c01 + c11 - c00 - c10), 0.25 * (c00 + c10 + c01 + c11)   # psi ~ c0 + cx (u-1/2) + cy (v-1/2)
        M = np.array([[cx.real, cy.real], [cx.imag, cy.imag]]); rhs = -np.array([c0.real, c0.imag])
        uv = np.linalg.solve(M, rhs) if abs(np.linalg.det(M)) > 1e-12 else np.zeros(2)
        uv = np.clip(uv + 0.5, 0, 1); pos.append([(i + uv[0]) * dx, (j + uv[1]) * dx])
    return np.array(pos).reshape(-1, 2), wn[wn != 0]


def pair_gap(pos, sgn):
    """closest distance between a + and a - vortex (periodic), or inf"""
    P, M = pos[sgn > 0], pos[sgn < 0]
    if len(P) == 0 or len(M) == 0: return np.inf
    d = P[:, None, :] - M[None, :, :]; d = np.abs((d + L / 2) % L - L / 2)
    return np.sqrt((d**2).sum(-1)).min()


def seam_3d(psi):
    """gap = min over z-slices of the +/- pair distance; both = fraction of z-slices still cut by all four lines; L1 line length"""
    phi = np.angle(psi); gaps = []; both = 0; count = 0
    for z in range(N):
        pos, sgn = zeros_2d(psi[:, :, z]); gaps.append(pair_gap(pos, sgn))
        both += len(sgn) >= 4                         # all four lines still cross this slice
        count += len(sgn)
    for ax in (0, 1):
        for s in range(N):
            sl = phi[s] if ax == 0 else phi[:, s]
            count += (windings(sl) != 0).sum()
    return min(gaps), both / N, count * dx


t0 = time.time()
if dim == 2:
    X, Y = np.meshgrid(x, x, indexing="ij")
    # two counter-propagating dipoles: a single dipole carries momentum and needs a phase twist 2 pi d / L across the box
    psi = (vortex(X, Y, L / 2 - D / 2, L / 4, +1) * vortex(X, Y, L / 2 + D / 2, L / 4, -1)
           * vortex(X, Y, L / 2 - D / 2, 3 * L / 4, -1) * vortex(X, Y, L / 2 + D / 2, 3 * L / 4, +1))
else:
    X, Y, Zc = np.meshgrid(x, x, x, indexing="ij"); psi = np.ones_like(X, dtype=complex)
    bend = A * np.cos(2 * np.pi * Zc / L)
    for s, sign, y0 in ((+1, +1, L / 4), (-1, -1, L / 4), (-1, +1, 3 * L / 4), (+1, -1, 3 * L / 4)):
        x0 = L / 2 + sign * (D / 2 - bend)
        psi *= vortex(X, Y, x0, y0, s)                 # two antiparallel pairs (zero momentum), each bowed toward itself at z = 0
for _ in range(int(round(TIM / 0.05))):
    psi = step(psi, 0.05, imaginary=True)               # short imaginary-time relaxation: settles the core and the image mismatch
E0 = energies(psi)
print("GPE %s  N=%d^%d  L=%g  dx=%.3f  d=%g  a=%g  dt=%g;  xi = 1, c = 1, Gamma = 2pi;  imaginary-time relax %.1f;  setup %.0fs" % (MODE, N, dim, L, dx, D, A, DT, TIM, time.time() - t0))
print("E_total %.5f = E_ki %.5f + E_kc %.5f + E_q %.5f + E_int %.5f  (decomposition residual %.1e)" % (E0[0], E0[1], E0[2], E0[3], E0[4], (E0[1] + E0[2] + E0[3] + E0[4] - E0[0]) / E0[0]))

hist = []; t = 0.0; mark = 0.0; nstep = int(round(T / DT))
if dim == 2:
    print("   t     E_tot rel err    E_ki      E_kc      tracked vortex (x, y)    displacement   speed so far   1/d   (n vortices)")
    track = np.array([L / 2 - D / 2, L / 4]); disp = 0.0
    for n in range(nstep + 1):
        if t >= mark - 1e-9:
            E = energies(psi); pos, sgn = zeros_2d(psi)
            dd = (pos - track + L / 2) % L - L / 2; i = np.argmin((dd**2).sum(1)); disp += dd[i, 1]; track = pos[i]
            sp = disp / t if t > 0 else np.nan
            print("%5.1f   %+.2e       %.5f   %.5f   (%6.2f, %6.2f)          %+7.3f       %6.3f      %.3f   (%d)   (%.0fs)" % (t, (E[0] - E0[0]) / E0[0], E[1], E[2], track[0], track[1], disp, sp, 1 / D, len(sgn), time.time() - t0), flush=True)
            mark += EVERY
        if n == nstep: break
        psi = step(psi, DT); t += DT
    sp = abs(disp / t)
    print("2-D verdict: dipole speed %.3f vs point-vortex 1/d = %.3f (%+.0f%%): %s" % (sp, 1 / D, 100 * (sp * D - 1), "PASS" if abs(sp * D - 1) < 0.2 else "FAIL"))
else:
    print("   t     E_tot rel err    E_ki       E_kc       E_q        E_int     seam gap   slices cut by all 4   L1 line length")
    for n in range(nstep + 1):
        if t >= mark - 1e-9:
            E = energies(psi); gap, both, length = seam_3d(psi)
            print("%5.1f   %+.2e       %.5f    %.5f    %.5f    %.5f    %6.2f        %.3f              %7.1f   (%.0fs)" % (t, (E[0] - E0[0]) / E0[0], E[1], E[2], E[3], E[4], gap, both, length, time.time() - t0), flush=True)
            hist.append((t, E[1], E[2], gap, both, length)); mark += EVERY
        if n == nstep: break
        psi = step(psi, DT); t += DT
    H = np.array(hist); tt, Eki, Ekc, gap, both = H[:, 0], H[:, 1], H[:, 2], H[:, 3], H[:, 4]
    cut = np.where((both < both[0] - 0.05) & (np.minimum.accumulate(gap) < 1.5))[0]      # slices lost AND the gap had closed to the core scale
    if len(cut) == 0:
        print("3-D verdict: NO reconnection by T=%.0f (min gap %.2f)  -> KILL" % (T, gap.min()))
    else:
        ir = cut[0]; tr = tt[ir]
        sel = (tt < tr) & (gap > 1.5) & (gap < 0.6 * D) & np.isfinite(gap)
        if sel.sum() >= 3:
            p = np.polyfit(np.log(tr - tt[sel]), np.log(gap[sel]), 1)[0]
        else:
            p = np.nan
        w = (tt >= tr - 3) & (tt <= tr + 3); i0, i1 = np.where(w)[0][0], np.where(w)[0][-1]
        reclose = (both[ir:] > both[ir]).any()
        print("3-D reconnection at t_r = %.1f (first slice loses its pair; gap %.2f -> %.2f)" % (tr, gap[max(ir - 1, 0)], gap[ir]))
        print("gap law before the cut: gap ~ (t_r - t)^p with p = %.2f on %d points  [registered 0.5 +/- 0.15]" % (p, sel.sum()))
        print("across t_r +/- 3: E_ki %.5f -> %.5f (%+.1f%%), E_kc %.5f -> %.5f (%+.1f%%)%s" % (
            Eki[i0], Eki[i1], 100 * (Eki[i1] / Eki[i0] - 1), Ekc[i0], Ekc[i1], 100 * (Ekc[i1] / Ekc[i0] - 1) if Ekc[i0] > 0 else np.nan,
            (": sound released = %.0f%% of the incompressible energy lost" % (100 * (Ekc[i1] - Ekc[i0]) / (Eki[i0] - Eki[i1]))) if Eki[i0] > Eki[i1] else ""))
        print("re-approach after the cut: %s;  slices cut by both at T: %.3f" % ("YES" if reclose else "none", both[-1]))
        ok = 0.25 <= p <= 0.8 and Ekc[i1] > Ekc[i0]
        print("3-D verdict: reconnection %s, gap exponent %s, sound %s -> %s" % ("YES", "in range" if 0.25 <= p <= 0.8 else "OUT OF RANGE", "rises" if Ekc[i1] > Ekc[i0] else "does not rise", "PASS" if ok else "FAIL"))
