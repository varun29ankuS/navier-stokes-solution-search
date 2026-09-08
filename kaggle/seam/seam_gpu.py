"""The seam race at the viscous thickness: GPU spectral Navier-Stokes, 256^3, tracking the strong twist until it turns.

The searcher's fastest fields grow by pressing two intense antiparallel sheets together (the strong twist); forbid it
and growth collapses (C19). Locally the twist rises under every viscosity inside a 32^3 window (clock ~1.1) - the
pressing is inviscid, and the seam has not yet thinned to the viscous scale sqrt(nu/s) where reconnection cuts it.
Here the window is long enough to see the seam reach that scale. Per output time, on the high-vorticity set:
    Z/Z0, max|w|, strip delta with the clock (2dx)
    twist@sep  strong twist at FIXED physical separations 0.05, 0.1, 0.2, 0.4 (enstrophy-weighted sharp reversal;
               the earlier 1-2 cell version carried the grid scale); the verdict uses sep = 0.1
    anti       fraction with an antiparallel neighbour at separation 0.1
    Re_seam    |w| ell^2 / nu on the twisted set (median): the seam's own Reynolds number, Delta_u ell / nu
    ell        sheet thickness |w| / |grad |w||, median on the high set
    ell_nu     viscous thickness sqrt(nu / s), s = mean compression across the sheet (-n.S.n) on the high set
    race       ell / ell_nu: the seam reaches its viscous thickness at ~1
    cut        on the twisted high set, fraction where viscous cancellation -nu w.lap w exceeds POSITIVE stretching
REGISTERED (C20, 2026-09-07): at nu = 2e-3 the race variable reaches <= 1.5 before the clock, the strong twist peaks
there and falls by >= 20% before the clock, with max|w| growth < 3x; at nu = 0 the twist rises to the clock.
REGISTERED (C21, 2026-09-08, the feedback law): the twist obeys dtau/dt ~ (s - nu/ell^2) tau, so it peaks when the
seam Reynolds number crosses a fixed threshold: Re_seam at the twist peak is the SAME number (within x2) at
nu = 2e-3, 1e-3, 5e-4, while the peak arrives later and higher as nu falls. KILL: Re_seam at the peak rises by more
than x2 per halving of nu (the seam needs ever more Reynolds number to be cut: Delta_u growing with the collapse).
REGISTERED (C22, 2026-09-08, the descent law): the twist at fixed separations peaks in sequence from large to small
separation - one wave descending in scale - and above the viscous thickness its descent is nu-INDEPENDENT (peak
times at 0.4 ... 0.1 agree across nu = 2e-3, 1e-3, 5e-4 within 0.1) and ACCELERATES: the time per halving of the
separation shrinks as the wave descends between 0.4 and 0.1. Refuted by: constant or growing halving time (the
wave never arrives: phase 1 all the way down), or peak times that depend on nu above ell_nu.
REGISTERED (C23, 2026-09-08, generality): (a) Kida-Pelz, coherent (anti = 0): twist <= 0.01 at every separation
through its clock - no wave, the instrument is blind to symmetric focusing; (b) the concentration-rewarded field
(ckn64): a descending wave, peak times nu-independent within 0.05 at 2e-3 and 1e-3, accelerating octave time, its
own T*; (c) the tube pair: a descending wave, nu-independent, with CONSTANT octave time early (mutual induction at
fixed circulation: exponential approach) and acceleration only in the last octaves.
The V: gap arm fitted to the resolved peaks (gap = a (T* - t)), thickness arm to delta(t) (exponential), crossing
reported with sqrt(nu/s) there.
FORCED MODE (2026-09-08, the day the field moved): Fefferman's breakdown statements (C)/(D) allow a smooth force
f(x, t). FORCE=eps adds f = eps * (low-k part of the initial field, |k| <= FKMAX), divergence-free, smooth, periodic,
time-independent: the external squeeze held on forever. Energy is no longer fixed - the budget that makes the unforced
collapse "free but unpaid for" is replenished from outside. REGISTERED (C24): with FORCE > 0 the V bottom is no
longer the end of the descent - twist@0.05 keeps rising past the unforced turnover time, max|w| accelerates
instead of peaking, and the wave passes the viscous thickness (race < 1 with the twist still rising); at FORCE = 0
(the runs above) it turns over. Refuted by: a forced run whose twist still turns over and whose max|w| still peaks
inside the clock.
LAGRANGIAN MODE (2026-09-08, Kelvin's frame): LAGR=1 seeds NP particles on the high-vorticity set at t = TSEED (when the
pair is visible), splits them into the two sheets by the sign of omega . xi_ref (xi_ref the principal vorticity
direction of the high set), and advects them (RK2, trilinear on the GPU). Reported along the paths: the material
|w| (median over particles; Helmholtz: it should follow the stretching), the GAP (median distance from a particle of
one sheet to the nearest of the other, periodic) and the velocity JUMP across the pair (median |u_A - u_B| over those
nearest pairs). C25 says the jump stays bounded by the data; lambda is the gap's exponent. REGISTERED (C26): on the
sheet field at nu = 2e-3 the Lagrangian gap closes linearly (lambda = 1 +- 0.2 on the resolved window) and the jump
stays within 2x of its value at TSEED until the merge; the material |w| grows while the gap closes and turns over at
the merge. Refuted by: a jump that grows by more than 2x while the gap closes (C25 fails on this field), or lambda < 0.7.
TRACKING FORCE (2026-09-08): FMODE=track makes the force follow the collapse. Every output step the force is set to
eps x P[u_H], where u_H is the velocity induced (Biot-Savart) by the vorticity of the high set alone, |w| > 0.5 max,
smoothly masked: the pair's own self-induction, amplified - phase 2 fed directly - and nothing else. Smooth in x,
piecewise-steady in t; a legitimate f(x, t) for Fefferman (C)/(D) once recorded along the trajectory. REGISTERED (C27):
with FMODE=track and eps >= 1 at nu = 2e-3 the seam passes the floor inside the clock - twist@0.05 keeps rising past the
unforced turnover (t ~ 1.7), max|w| accelerates over the last quarter of the window, and the race variable drops below
1 with the twist still rising; the steady force (C24) did none of these. Refuted by: the twist turning over as
unforced under the tracking force too - which would say a force that merely amplifies the pair's self-induction is
still not the force the proofs need.
THE FLIP (2026-09-08): does the seam flip again below the cut? Reconnection leaves threads that are antiparallel to
each other at a smaller scale (Hussain's bridging; Yao-Hussain 2020's reconnection cascade). In the Lagrangian mode
the FLIP column is the fraction of tagged particles whose sign of omega . xi_ref has reversed since seeding: zero for
a clean cut, finite if the sheets have bridged; and the gap re-closing after the merge is a second seam. REGISTERED
(C28): at nu = 2e-3 the flip fraction rises from ~0 to 0.10-0.30 across the merge, and the material gap, after
re-opening, closes a second time before the clock. Refuted by: a flip fraction staying below 0.05, or no second
closing inside the clock.
WHO CLOSES THE GAP (2026-09-08): in the Lagrangian mode, for every nearest A-B pair, the closing rate from the full
velocity, (u_A - u_B) . d_hat, and the SELF-INDUCED closing rate from the pair's own Biot-Savart field u_H (the velocity
induced by the smoothly-masked high set), (u_H,A - u_H,B) . d_hat. Their ratio says whether the sheets close under
their own induction (Helmholtz) or under external strain. Also, for flipped particles, their distance to the other
sheet relative to the gap: near 0 = bridges between the sheets, ~1 = threads beside them. REGISTERED (C29): on the
descent (t = 1.0 to the merge) the self-induced rate accounts for >= 70% of the measured closing rate - the gap law
g = 0.51 (T* - t) is Helmholtz, not a fit. Refuted by: a self-induced share below 40% (external strain closes the gap).
usage: IC=found|kp|pair N=256 NU=2e-3 T=3.0 [FORCE=0.5 FKMAX=4 FMODE=steady|track] [LAGR=1 TSEED=1.0 NP=4000] python seam_gpu.py"""
import os, sys, time, math, subprocess, numpy as np, torch

# On Kaggle (one code file per kernel) this file schedules itself: each configuration runs as a subprocess of this
# script with SCHEDULE cleared, its log written to /kaggle/working.
SCHEDULE = os.environ.get("SCHEDULE")
if SCHEDULE is None and os.path.isdir("/kaggle/working"):
    SCHEDULE = "IC=found NU=2e-3 T=2.4 N=320 LAGR=1 TSEED=1.0;IC=found NU=1e-3 T=1.8 N=320 LAGR=1 TSEED=1.0"   # v11: who closes the gap
if SCHEDULE:
    for cfg in [c for c in SCHEDULE.split(";") if c.strip()]:
        env = dict(os.environ); env["SCHEDULE"] = ""; env.update(dict(kv.split("=") for kv in cfg.split()))
        tag = "_".join(kv.replace("=", "") for kv in cfg.split() if not kv.startswith("FOUND=")); out = os.path.join("/kaggle/working" if os.path.isdir("/kaggle/working") else ".", "seam_%s.txt" % tag)
        print("=== %s -> %s" % (cfg, out), flush=True)
        with open(out, "w") as f:
            pr = subprocess.Popen([sys.executable, __file__], env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
            for line in pr.stdout: print(line, end="", flush=True); f.write(line)
            pr.wait()
    raise SystemExit

IC = os.environ.get("IC", "found"); N = int(os.environ.get("N", 256)); NU = float(os.environ.get("NU", 2e-3)); T = float(os.environ.get("T", 3.0))
EVERY = float(os.environ.get("EVERY", 0.05)); DEV = "cuda" if torch.cuda.is_available() else "cpu"
SEPS = [float(v) for v in os.environ.get("SEPS", "0.05,0.07,0.1,0.14,0.2,0.28,0.4,0.56").split(",")]; VSEP = float(os.environ.get("VSEP", 0.1))
FOUND = os.environ.get("FOUND", "/kaggle/input/zef-found/leashed64_dmin030.npz")
FORCE = float(os.environ.get("FORCE", 0.0)); FKMAX = float(os.environ.get("FKMAX", 4.0)); FMODE = os.environ.get("FMODE", "steady")
SNAP = int(os.environ.get("SNAP", 0)); SNAPS = []
LAGR = int(os.environ.get("LAGR", 0)); TSEED = float(os.environ.get("TSEED", 1.0)); NP = int(os.environ.get("NP", 4000)); LAG = {}                                   # SNAP=1: save a slice of |w| and of the signed twist through the |w| maximum at every output
CD = torch.complex64; RD = torch.float32
k1 = torch.fft.fftfreq(N, d=1.0 / N).to(DEV).to(RD)
KX, KY, KZ = torch.meshgrid(k1, k1, k1, indexing="ij"); K = [KX, KY, KZ]
K2 = KX**2 + KY**2 + KZ**2; K2S = K2.clone(); K2S[0, 0, 0] = 1.0
deal = ((KX.abs() < N / 3) & (KY.abs() < N / 3) & (KZ.abs() < N / 3)).to(RD)
KMAG = torch.sqrt(K2); NB = int(N / 3)
fft, ifft = torch.fft.fftn, torch.fft.ifftn
dx = 2 * math.pi / N


def project(F):
    kd = sum(K[i] * F[i] for i in range(3)) / K2S
    return [F[i] - K[i] * kd for i in range(3)]


FHAT = None                                                                         # the smooth force, set after the initial field is built


def transport(U):
    Ud = [Ui * deal for Ui in U]; u = [ifft(Ui).real for Ui in Ud]; out = []
    for i in range(3):
        adv = sum(u[j] * ifft(1j * K[j] * Ud[i]).real for j in range(3))
        div = sum(ifft(1j * K[j] * fft(u[j] * u[i])).real for j in range(3))
        out.append(-0.5 * fft(adv + div) * deal + (FHAT[i] if FHAT is not None else 0.0))
        del adv, div
    return project(out)


def step(U, dt):
    f, f2 = torch.exp(-NU * K2 * dt), torch.exp(-NU * K2 * dt / 2)
    a = transport(U); b = transport([f2 * (U[i] + dt / 2 * a[i]) for i in range(3)])
    c = transport([f2 * U[i] + dt / 2 * b[i] for i in range(3)]); d = transport([f * U[i] + dt * f2 * c[i] for i in range(3)])
    return [f * U[i] + dt / 6 * (f * a[i] + 2 * f2 * b[i] + 2 * f2 * c[i] + d[i]) for i in range(3)]


def strip(U):
    """analyticity-strip width from the tail of the energy spectrum; nan while the tail is still empty (below 1e-12 of
    the peak shell: an upsampled low-mode field has nothing there but float32 round-off until the cascade arrives)"""
    e = 0.5 * sum(Ui.abs() ** 2 for Ui in U) / N**6
    spec = torch.stack([e[(KMAG >= n - 0.5) & (KMAG < n + 0.5)].sum() for n in range(1, NB)]).cpu().numpy().astype(np.float64)
    ks = np.arange(1, NB); sel = (ks >= NB // 2) & (spec > 1e-12 * spec.max())
    return -np.polyfit(ks[sel], np.log(spec[sel]), 1)[0] / 2 if sel.sum() > 4 else float("nan")


def energy(U):
    return 0.5 * sum((ifft(Ui * deal).real ** 2).mean() for Ui in U).item()


@torch.no_grad()
def interp3(F, P):
    """trilinear interpolation of a real field F [N,N,N] at particle positions P [n,3] in [0, 2pi)^3 (periodic)"""
    g = (P / dx) % N; i0 = torch.floor(g).long(); f = (g - i0.float()); i1 = (i0 + 1) % N; out = torch.zeros(P.shape[0], device=DEV, dtype=RD)
    for a in (0, 1):
        for b in (0, 1):
            for c in (0, 1):
                w = (f[:, 0] if a else 1 - f[:, 0]) * (f[:, 1] if b else 1 - f[:, 1]) * (f[:, 2] if c else 1 - f[:, 2])
                out = out + w * F[(i1 if a else i0)[:, 0], (i1 if b else i0)[:, 1], (i1 if c else i0)[:, 2]]
    return out


def lagr_seed(U):
    """particles on the high set, split into the two sheets by the sign of omega . xi_ref"""
    Ud = [Ui * deal for Ui in U]
    w = [ifft(1j * K[(i + 1) % 3] * Ud[(i + 2) % 3] - 1j * K[(i + 2) % 3] * Ud[(i + 1) % 3]).real for i in range(3)]
    wm = torch.sqrt(sum(wi**2 for wi in w)) + 1e-30; high = wm > 0.5 * wm.max()
    idx = torch.nonzero(high); sel = idx[torch.randperm(len(idx), device=DEV)[:NP]]
    P = (sel.float() + 0.5) * dx
    xi = torch.stack([wi[high] / wm[high] for wi in w], 1); M = xi.T @ xi; xi_ref = torch.linalg.eigh(M)[1][:, -1]
    sgn = torch.sign(sum(w[c][sel[:, 0], sel[:, 1], sel[:, 2]] * xi_ref[c] for c in range(3)))
    LAG["xi_ref"] = xi_ref
    return P, sgn


def lagr_advect(U, P, dt):
    u = [ifft(Ui * deal).real for Ui in U]
    v1 = torch.stack([interp3(ui, P) for ui in u], 1); Pm = (P + 0.5 * dt * v1) % (2 * math.pi)
    v2 = torch.stack([interp3(ui, Pm) for ui in u], 1)
    return (P + dt * v2) % (2 * math.pi)


def lagr_diag(U, P, sgn):
    Ud = [Ui * deal for Ui in U]; u = [ifft(Ui).real for Ui in Ud]
    w = [ifft(1j * K[(i + 1) % 3] * Ud[(i + 2) % 3] - 1j * K[(i + 2) % 3] * Ud[(i + 1) % 3]).real for i in range(3)]
    wm = torch.sqrt(sum(wi**2 for wi in w))
    wmp = interp3(wm, P); A = P[sgn > 0]; B = P[sgn < 0]
    if len(A) < 10 or len(B) < 10: return wmp.median().item(), float("nan"), float("nan"), float("nan")
    dvec = A[:, None, :] - B[None, :, :]; dvec = (dvec + math.pi) % (2 * math.pi) - math.pi
    dist = torch.sqrt((dvec**2).sum(-1)); j = torch.argmin(dist, 1); gap = dist[torch.arange(len(A), device=DEV), j]
    uA = torch.stack([interp3(ui, A) for ui in u], 1); uB = torch.stack([interp3(ui, B[j]) for ui in u], 1)
    jump = torch.sqrt(((uA - uB) ** 2).sum(1))
    wsign = torch.sign(sum(interp3(w[c], P) * LAG["xi_ref"][c] for c in range(3)))                # the sign of omega . xi_ref carried by each particle now
    flipped = wsign != sgn; flip = flipped.float().mean().item()
    # who closes the gap: closing rate from the full velocity and from the pair's own induced velocity, along the separation
    dhat = dvec[torch.arange(len(A), device=DEV), j] / (gap[:, None] + 1e-30)                      # unit vector from B to A (periodic)
    rate_full = ((uA - uB) * dhat).sum(1).median().item()                                        # d gap / dt from u (negative = closing)
    uH = induced_velocity(U); hA = torch.stack([interp3(hi, A) for hi in uH], 1); hB = torch.stack([interp3(hi, B[j]) for hi in uH], 1)
    rate_self = ((hA - hB) * dhat).sum(1).median().item()
    # where the flipped fluid sits: distance to the other sheet over the gap (0 = between the sheets, ~1 = beside them)
    fA = flipped[sgn > 0]
    where = (gap[fA] / gap.median()).median().item() if fA.sum() >= 5 else float("nan")
    return wmp.median().item(), gap.median().item(), jump.median().item(), flip, rate_full, rate_self, where


def induced_velocity(U):
    """u_H: the velocity induced (Biot-Savart) by the smoothly-masked high-vorticity set alone, in physical space"""
    Ud = [Ui * deal for Ui in U]
    w = [ifft(1j * K[(i + 1) % 3] * Ud[(i + 2) % 3] - 1j * K[(i + 2) % 3] * Ud[(i + 1) % 3]).real for i in range(3)]
    wm = torch.sqrt(sum(wi**2 for wi in w)); mask = torch.sigmoid((wm / wm.max() - 0.5) / 0.1)
    W = project([fft(wi * mask) for wi in w])
    uH = project([(1j * (K[(i + 1) % 3] * W[(i + 2) % 3] - K[(i + 2) % 3] * W[(i + 1) % 3]) / K2S) * deal for i in range(3)])
    return [ifft(Ui).real for Ui in uH]


def tracking_force(U):
    """eps x P[u_H]: the velocity induced by the smoothly-masked high-vorticity set alone"""
    Ud = [Ui * deal for Ui in U]
    w = [ifft(1j * K[(i + 1) % 3] * Ud[(i + 2) % 3] - 1j * K[(i + 2) % 3] * Ud[(i + 1) % 3]).real for i in range(3)]
    wm = torch.sqrt(sum(wi**2 for wi in w)); mask = torch.sigmoid((wm / wm.max() - 0.5) / 0.1)
    W = project([fft(wi * mask) for wi in w])
    uH = project([(1j * (K[(i + 1) % 3] * W[(i + 2) % 3] - K[(i + 2) % 3] * W[(i + 1) % 3]) / K2S) * deal for i in range(3)])
    return [FORCE * Ui for Ui in uH]


def diag(U):
    Ud = [Ui * deal for Ui in U]
    G = [[ifft(1j * K[i] * Ud[j]).real for j in range(3)] for i in range(3)]          # G[i][j] = d_i u_j
    w = [G[1][2] - G[2][1], G[2][0] - G[0][2], G[0][1] - G[1][0]]
    wm = torch.sqrt(sum(wi**2 for wi in w)) + 1e-30; xi = [wi / wm for wi in w]; w2 = wm**2
    high = wm > 0.5 * wm.max()
    anti = torch.zeros_like(high); softs = {}; norm = (w2**2).mean()
    for sep in SEPS:
        sh = max(1, int(round(sep / dx))); soft = 0.0
        for ax in range(3):
            beta = 1 - sum(xi[c] * torch.roll(xi[c], sh, dims=ax) for c in range(3))
            soft = soft + (w2 * torch.roll(w2, sh, dims=ax) * torch.relu(beta - 1) ** 2).mean()
            if sep == VSEP: anti |= (beta > 1) | (torch.roll(beta, -sh, dims=ax) > 1)
        softs[sep] = (soft / norm).item()
    soft = softs[VSEP]
    # sheet thickness and the viscous thickness
    wh = fft(wm); gw = [ifft(1j * K[i] * wh).real for i in range(3)]; gwm = torch.sqrt(sum(g**2 for g in gw)) + 1e-30
    ell = (wm / gwm)[high].median().item()
    n = [g / gwm for g in gw]
    S = [[0.5 * (G[i][j] + G[j][i]) for j in range(3)] for i in range(3)]
    snn = sum(n[i] * S[i][j] * n[j] for i in range(3) for j in range(3))
    s_comp = (-snn[high]).mean().clamp(min=1e-12).item()
    ell_nu = math.sqrt(NU / s_comp) if NU > 0 else float("nan")
    # reconnection: viscous cancellation vs POSITIVE stretching on the twisted high set
    stretch = sum(w[i] * S[i][j] * w[j] for i in range(3) for j in range(3))
    lapw = [ifft(-K2 * fft(wi)).real for wi in w]
    annih = -NU * sum(w[i] * lapw[i] for i in range(3))
    tw = high & anti
    cut = ((annih[tw] > torch.relu(stretch[tw])) & (stretch[tw] > 0)).float().mean().item() if tw.sum() > 0 else 0.0
    sel = tw if tw.sum() >= 0.01 * high.sum() else high                                  # the twisted set only once it is 1% of the high set
    re_seam = ((wm * (wm / gwm) ** 2)[sel].median() / NU).item() if NU > 0 else float("nan")     # |w| ell^2 / nu = Delta_u ell / nu on the seam
    valid = s_comp > 1e-6
    if SNAP:
        iz = int(torch.argmax(wm).item() % N)                                                  # the z-plane through the vorticity maximum (index of the last axis)
        sh = max(1, int(round(VSEP / dx))); bsum = torch.zeros_like(wm[:, :, iz])
        for ax in range(2):                                                                  # in-plane signed reversal at separation VSEP: max over the two in-plane directions
            beta = 1 - sum(xi[c] * torch.roll(xi[c], sh, dims=ax) for c in range(3))
            bsum = torch.maximum(bsum, beta[:, :, iz])
        SNAPS.append((wm[:, :, iz].cpu().numpy().astype(np.float32), bsum.cpu().numpy().astype(np.float32), iz))
    return (0.5 * w2.mean().item(), wm.max().item(), soft, anti[high].float().mean().item(), ell, ell_nu, cut, strip(U), softs, re_seam, valid)


x = torch.arange(N, device=DEV, dtype=RD) * dx
if IC == "pair":
    # two antiparallel Gaussian vortex tubes along x, separation D, core sigma, bowed toward each other (Kerr-type)
    D = float(os.environ.get("D", 0.8)); SIG = float(os.environ.get("SIG", 0.2)); A = float(os.environ.get("A", 0.15))   # apex gap D - 2A = 0.5 = 2.5 sigma: apart, then approaching
    X, Y, Z_ = torch.meshgrid(x, x, x, indexing="ij"); wx = torch.zeros_like(X)
    for sgn in (+1, -1):
        yc = math.pi + sgn * (D / 2 - A * torch.cos(X - math.pi)); r2 = ((Y - yc) % (2 * math.pi) - math.pi) ** 2 + ((Z_ - math.pi) % (2 * math.pi) - math.pi) ** 2
        wx = wx + sgn * torch.exp(-r2 / (2 * SIG**2))
    W = [fft(wx), torch.zeros_like(fft(wx)), torch.zeros_like(fft(wx))]
    W = project(W)                                                                     # divergence-free vorticity
    U = [(1j * (K[(i + 1) % 3] * W[(i + 2) % 3] - K[(i + 2) % 3] * W[(i + 1) % 3]) / K2S) * deal for i in range(3)]   # u = curl^-1 w
    U = project(U)
elif IC == "kp":
    X, Y, Z_ = torch.meshgrid(x, x, x, indexing="ij")
    U = [fft(torch.sin(X) * (torch.cos(3 * Y) * torch.cos(Z_) - torch.cos(Y) * torch.cos(3 * Z_))).to(CD),
         fft(torch.sin(Y) * (torch.cos(3 * Z_) * torch.cos(X) - torch.cos(Z_) * torch.cos(3 * X))).to(CD),
         fft(torch.sin(Z_) * (torch.cos(3 * X) * torch.cos(Y) - torch.cos(X) * torch.cos(3 * Y))).to(CD)]
    U = [Ui * deal for Ui in U]; U = project(U)
else:
    if not os.path.exists(FOUND):                                                       # Kaggle mounts datasets under varying paths: find the file by name
        import glob as _g; hits = _g.glob("/kaggle/input/**/" + os.path.basename(FOUND), recursive=True)
        if not hits: raise SystemExit("found field %s not under /kaggle/input: %s" % (os.path.basename(FOUND), _g.glob("/kaggle/input/**/*", recursive=True)[:20]))
        FOUND = hits[0]; print("found field at", FOUND, flush=True)
    uf = np.load(FOUND)["u"].astype(np.float32); n0 = uf.shape[1]; U = []
    for c in range(3):
        uh = np.fft.fftn(uf[c]) * (N / n0) ** 3; big = np.zeros((N, N, N), np.complex64); h = n0 // 2
        for a in (slice(0, h), slice(-h, None)):
            for b in (slice(0, h), slice(-h, None)):
                for cc in (slice(0, h), slice(-h, None)):
                    big[a, b, cc] = uh[a, b, cc]
        U.append(torch.tensor(big, device=DEV, dtype=CD))
    U = project(U)
with torch.no_grad():
    Z0 = 0.5 * sum((ifft(1j * K[a] * U[b] - 1j * K[b] * U[a]).real ** 2).mean() for a, b in ((1, 2), (2, 0), (0, 1))).item()
    U = [Ui * math.sqrt(0.375 / Z0) for Ui in U]; Z0 = 0.375; E0 = energy(U)
    if FORCE > 0 and FMODE == "steady":
        lowk = (KMAG <= FKMAX).to(RD); FHAT = [FORCE * Ui * lowk for Ui in U]              # f = eps * u0 restricted to |k| <= FKMAX: smooth, periodic, divergence-free, steady
        fpow = math.sqrt(sum((ifft(Fi).real ** 2).mean().item() for Fi in FHAT))
        print("forced: f = %.2f x (initial field, |k| <= %g), |f|_rms = %.4f, energy injection rate at t=0 = %.4f (vs 2 nu Z0 = %.4f)" % (
            FORCE, FKMAX, fpow, sum((ifft(Fi).real * ifft(Ui * deal).real).mean().item() for Fi, Ui in zip(FHAT, U)), 2 * NU * Z0), flush=True)
print("seam race on GPU: IC=%s  N=%d^3  nu=%g  T=%.1f  Z0=%.3f  clock 2dx=%.4f  device=%s  FORCE=%g (%s)" % (IC, N, NU, T, Z0, 2 * dx, DEV, FORCE, FMODE), flush=True)
print("   t     Z/Z0    max|w|   twist@" + " @".join("%g" % v for v in SEPS) + "   anti    ell     ell_nu   race   cut    Re_seam   delta     E/E0" + ("   | material|w|  gap   jump   flip   dgap/dt(u) dgap/dt(self)  flipped@" if LAGR else ""))
t, mark, t0 = 0.0, 0.0, time.time(); hist = []
while t <= T + 1e-9:
    if t >= mark - 1e-9:
        if FORCE > 0 and FMODE == "track":
            with torch.no_grad():
                FHAT = tracking_force(U); inj = sum((ifft(Fi).real * ifft(Ui * deal).real).mean().item() for Fi, Ui in zip(FHAT, U))
            if abs(t - round(t / 0.5) * 0.5) < 1e-6: print("   tracking force refreshed: injection rate %.4f (2 nu Z = %.4f)" % (inj, 2 * NU * 0.5 * sum((ifft(1j * K[a] * U[b] - 1j * K[b] * U[a]).real ** 2).mean().item() for a, b in ((1, 2), (2, 0), (0, 1)))), flush=True)
        Z, wmax, soft, anti, ell, ell_nu, cut, d, softs, re_seam, valid = diag(U); Er = energy(U) / E0
        lag = ""
        if LAGR and t >= TSEED - 1e-9:
            if "P" not in LAG:
                with torch.no_grad(): LAG["P"], LAG["sgn"] = lagr_seed(U)
                print("lagrangian: %d particles seeded on the high set at t = %.2f; sheets A/B = %d/%d" % (len(LAG["P"]), t, int((LAG["sgn"] > 0).sum()), int((LAG["sgn"] < 0).sum())), flush=True)
            with torch.no_grad(): mw, gap, jump, flip, rf, rs, where = lagr_diag(U, LAG["P"], LAG["sgn"])
            LAG.setdefault("rows", []).append((t, mw, gap, jump, flip, rf, rs, where)); lag = "   | %8.2f   %.4f   %.4f   %.3f   %+.4f  %+.4f  %s" % (mw, gap, jump, flip, rf, rs, ("%.2f" % where) if np.isfinite(where) else "  -  ")
        hist.append((t, Z / Z0, wmax, soft, anti, ell, ell_nu, cut, d if np.isfinite(d) else 0.0, re_seam, float(valid)) + tuple(softs[v] for v in SEPS) + (Er,))
        print("%5.2f   %6.3f   %7.2f   %s   %.3f   %.4f   %s   %5.2f   %.3f   %8.1f   %s   %.6f%s%s   (%.0fs)" % (
            t, Z / Z0, wmax, " ".join("%.5f" % softs[v] for v in SEPS), anti, ell, ("%.4f" % ell_nu) if valid else "   -  ", ell / ell_nu if (NU > 0 and valid) else float("nan"), cut, re_seam,
            ("%.4f" % d) if np.isfinite(d) else "(tail empty)", Er, "" if (np.isfinite(d) and d > 2 * dx) else ("" if not np.isfinite(d) else "  <-- past the clock"), lag, time.time() - t0), flush=True)
        mark += EVERY
        if t >= T - 1e-9: break
    with torch.no_grad():
        umax = max(ifft(Ui).real.abs().max().item() for Ui in U)
        dt = min(2.0 / N, 0.5 * dx / max(umax, 1e-9), mark - t + 1e-12)
        if LAGR and "P" in LAG: LAG["P"] = lagr_advect(U, LAG["P"], dt)
        U = step(U, dt); t += dt
if LAGR and LAG.get("rows"):
    Lr = np.array(LAG["rows"]); tl, mw, gp, jp, fl, rf, rs, wh = Lr[:, 0], Lr[:, 1], Lr[:, 2], Lr[:, 3], Lr[:, 4], Lr[:, 5], Lr[:, 6], Lr[:, 7]
    dcl = np.array([h[8] for h in hist]); tcl = np.array([h[0] for h in hist]); tclock = tcl[dcl > 2 * dx][-1] if (dcl > 2 * dx).any() else tl[0]
    m = (tl <= tclock) & np.isfinite(gp) & (gp > 0)
    imin = int(np.argmin(gp[m])) if m.any() else 0; tmerge = tl[m][imin]
    mm = m & (tl < tmerge) & (gp > 1.2 * gp[m].min())
    lam = np.nan
    if mm.sum() >= 4:
        a, b = np.polyfit(tl[mm], gp[mm], 1); Tst = -b / a if a < 0 else np.nan
        if np.isfinite(Tst): lam = np.polyfit(np.log(Tst - tl[mm]), np.log(gp[mm]), 1)[0]
    j0 = jp[m][0]; jmax = np.nanmax(jp[m & (tl <= tmerge)]) if (m & (tl <= tmerge)).any() else np.nan
    print("\nC26 Lagrangian: gap %.4f -> %.4f (min at t = %.2f, inside the clock %.2f); linear fit T* = %s, lambda = %s; jump %.4f at seed -> max %.4f before the merge (x%.2f); material |w| %.1f -> %.1f" % (
        gp[m][0], gp[m].min(), tmerge, tclock, ("%.2f" % Tst) if mm.sum() >= 4 and np.isfinite(Tst) else "-", ("%.2f" % lam) if np.isfinite(lam) else "-", j0, jmax, jmax / j0, mw[m][0], mw[m].max()))
    ok = np.isfinite(lam) and abs(lam - 1) <= 0.2 and jmax / j0 <= 2.0
    # C29: who closes the gap - the self-induced share of the closing rate on the descent
    desc = m & (tl < tmerge) & (rf < 0)
    if desc.sum() >= 3:
        share = float(np.median(rs[desc] / rf[desc])); fd = float(np.polyfit(tl[desc], gp[desc], 1)[0])
        print("C29 who closes the gap: on the descent the measured closing rate (u) is %+.4f median (finite-difference of the gap %+.4f), the self-induced rate (u_H) %+.4f: self-induced share %.2f" % (
            float(np.median(rf[desc])), fd, float(np.median(rs[desc])), share))
        print("REGISTERED C29: %s" % ("PASS: the sheets close under their own induction" if share >= 0.7 else ("KILL: external strain closes the gap (share %.2f)" % share if share < 0.4 else "between (share %.2f)" % share)))
    wpost = m & (tl > tmerge) & np.isfinite(wh)
    if wpost.any(): print("flipped fluid sits at %.2f of the gap from the other sheet (median over the post-merge rows): %s" % (float(np.median(wh[wpost])), "bridges between the sheets" if np.median(wh[wpost]) < 0.5 else "threads beside them"))
    # C28: the flip, and a second closing of the gap after the merge
    post = m & (tl > tmerge); second = False
    if post.sum() >= 3:
        gpost = gp[post]; ip = int(np.argmax(gpost)); second = ip < len(gpost) - 1 and gpost[ip:].min() < 0.9 * gpost[ip]
    fmerge = fl[m][imin] if m.any() else float("nan"); fend = fl[m][-1] if m.any() else float("nan")
    print("C28 flip: fraction of tagged particles with reversed omega . xi_ref: %.3f at seed -> %.3f at the merge -> %.3f at the clock; gap after the merge %s" % (
        fl[m][0], fmerge, fend, ("re-opens to %.4f then closes again to %.4f: a SECOND seam" % (gpost[ip], gpost[ip:].min())) if second else ("re-opens, no second closing inside the clock" if post.sum() >= 3 else "no rows after the merge inside the clock")))
    print("REGISTERED C28: %s" % ("PASS: the seam flips again" if (fend >= 0.10 and second) else ("KILL: clean cut (flip %.3f)" % fend if fend < 0.05 else "between")))
    print("REGISTERED C26: %s" % ("PASS: linear closing, bounded jump" if ok else ("KILL: the jump grows more than 2x while the gap closes (C25 fails here)" if jmax / j0 > 2.0 else ("KILL: lambda = %.2f < 0.7" % lam if np.isfinite(lam) and lam < 0.7 else "between (see rows)"))))
if SNAP:
    out = os.path.join("/kaggle/working" if os.path.isdir("/kaggle/working") else ".", "snaps_%s_nu%g.npz" % (IC, NU))
    np.savez_compressed(out, t=np.array([h[0] for h in hist]), wm=np.stack([s_[0] for s_ in SNAPS]), beta=np.stack([s_[1] for s_ in SNAPS]), iz=np.array([s_[2] for s_ in SNAPS]), maxw=np.array([h[2] for h in hist]), delta=np.array([h[8] for h in hist]))
    print("snapshots written to", out, flush=True)
H = np.array(hist); tt, tw, race, d, res, valid = H[:, 0], H[:, 3], H[:, 5] / np.where(H[:, 6] > 0, H[:, 6], np.nan), H[:, 8], H[:, 9], H[:, 10] > 0.5
race = np.where(valid, race, np.nan)
# C22: the descent - time of the twist peak at each separation (inside the clock only), and the halving time between scales
inclock0 = d > 2 * dx; last0 = np.where(inclock0)[0][-1] if inclock0.any() else 0
print("\nC22 descent: separation -> peak twist, peak time (clock expires at t = %.2f)" % tt[last0])
peaks = {}
for j, sep in enumerate(SEPS):
    col = H[:last0 + 1, 11 + j]; ip = int(np.argmax(col)); at_clock = (ip == last0) or (ip == 0)
    peaks[sep] = (tt[ip], col[ip], at_clock)
    print("   sep %.2f   peak %.5f at t = %.2f%s" % (sep, col[ip], tt[ip], "   (at the clock or at t = 0: not a peak)" if at_clock else ""))
good = [(sep, peaks[sep][0]) for sep in sorted(SEPS, reverse=True) if not peaks[sep][2] and peaks[sep][1] > 1e-4]
if len(good) >= 3:
    print("   halving times (time for the wave to descend one octave), from resolved peaks only:")
    for (s1, t1), (s2, t2) in zip(good[:-1], good[1:]):
        print("      %.2f -> %.2f : %.2f  per octave %.2f" % (s1, s2, t2 - t1, (t2 - t1) / max(math.log2(s1 / s2), 1e-9)))
    try:
        # the V: gap arm (linear in t through the resolved peaks) against the thickness arm (delta, exponential, inside the clock)
        gfit = [g for g in good if g[0] >= 0.1]                                          # the gap arm above the wall only
        ga, gb = np.polyfit([g[1] for g in gfit], [g[0] for g in gfit], 1); Tstar = -gb / ga if ga < 0 else float("nan")
        tpeak = tt[int(np.argmax(H[:last0 + 1, 11 + SEPS.index(VSEP)]))]                  # the thickness arm before the merge: first real delta -> twist@VSEP peak
        m = (tt <= tpeak) & (d > 2 * dx)
        if m.sum() >= 4 and np.isfinite(Tstar):
            c1, c0 = np.polyfit(tt[m], np.log(d[m]), 1)
            tq = np.linspace(tt[m][0], Tstar - 1e-3, 3000); gapq = ga * tq + gb; thq = np.exp(c0 + c1 * tq); iv = int(np.argmin(np.abs(gapq - thq)))
            jv = int(np.argmin(np.abs(tt - tq[iv]))); eln = H[jv, 6]
            print("   V: gap = %.3f (%.2f - t) from sep >= 0.1, T* = %.2f; thickness e-fold %.2f fitted on t in [%.2f, %.2f]; arms meet at t = %.2f, scale %.3f, sqrt(nu/s) there %.3f (%s the clock)" % (
                -ga, Tstar, Tstar, -1 / c1 if c1 < 0 else float("inf"), tt[m][0], tt[m][-1], tq[iv], thq[iv], eln, "inside" if tq[iv] <= tt[last0] else "beyond"))
    except Exception as e:
        print("   V: analysis failed (%s)" % e)

inclock = d > 2 * dx; ic = np.where(inclock)[0]
if len(ic) == 0:
    print("VERDICT: never inside the clock"); raise SystemExit
last = ic[-1]; ipk = int(np.nanargmax(tw[:last + 1])); tpk = tt[ipk]
fall = 1 - tw[last] / tw[ipk]
print("\nclock expires at t = %.2f (delta %.4f); strong twist peak %.5f at t = %.2f, at the clock %.5f (fall %.0f%%); max|w| growth to the clock %.2fx" % (tt[last], d[last], tw[ipk], tpk, tw[last], 100 * fall, H[last, 2] / H[0, 2]))
if NU > 0:
    if not np.isfinite(race[:last + 1]).any() or tw[ipk] <= 0:
        print("VERDICT: no seam formed (twist never rose) or the race variable is undefined - this run says nothing"); raise SystemExit
    rmin = np.nanmin(race[:last + 1]); print("race variable ell/ell_nu min inside the clock = %.2f (at t = %.2f); cut fraction at the clock %.3f" % (rmin, tt[int(np.nanargmin(race[:last + 1]))], H[last, 7]))
    print("C21: Re_seam at the twist peak (t = %.2f) = %.1f; twist peak %.5f; max|w| at the peak %.1f" % (tpk, res[ipk], tw[ipk], H[ipk, 2]))
    ok = rmin <= 1.5 and fall >= 0.2 and H[last, 2] / H[0, 2] < 3.0
    # KILL only if the twist is still rising at the clock AND max|w| is accelerating over the last quarter of the window
    # AND Re_seam is rising there - a rising twist with falling max|w| (Kida-Pelz) or a core-thickness race variable (a
    # diffusing pair) is not the phenomenon. Corrected after two spurious fires (pair v2, Kida-Pelz v4).
    q = max(1, (last + 1) // 4); wl = np.log(H[:last + 1, 2]); accel = (wl[last] - wl[last - q]) > (wl[last - q] - wl[max(last - 2 * q, 0)]) and wl[last] > wl[last - q]
    kill = (ipk == last) and rmin < 1.0 and accel and res[last] > res[max(last - q, 0)]
    print("REGISTERED C20 at nu=%g: %s" % (NU, "PASS: the seam reaches its viscous thickness and the twist turns over" if ok else ("KILL: the twist is still rising at the clock with the viscous scale resolved" if kill else "between (see rows)")))
else:
    print("REGISTERED C20 at nu=0: twist %s to the clock" % ("RISES" if ipk == last else "peaks at t=%.2f and falls %.0f%%" % (tpk, 100 * fall)))
