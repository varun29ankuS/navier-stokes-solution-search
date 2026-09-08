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
usage: IC=found|pair N=256 NU=2e-3 T=3.0 python seam_gpu.py"""
import os, sys, time, math, subprocess, numpy as np, torch

# On Kaggle (one code file per kernel) this file schedules itself: each configuration runs as a subprocess of this
# script with SCHEDULE cleared, its log written to /kaggle/working.
SCHEDULE = os.environ.get("SCHEDULE")
if SCHEDULE is None and os.path.isdir("/kaggle/working"):
    SCHEDULE = "IC=found NU=2e-3;IC=found NU=1e-3;IC=found NU=5e-4;IC=pair NU=2e-3 T=6"
if SCHEDULE:
    for cfg in [c for c in SCHEDULE.split(";") if c.strip()]:
        env = dict(os.environ); env["SCHEDULE"] = ""; env.update(dict(kv.split("=") for kv in cfg.split()))
        tag = "_".join(kv.replace("=", "") for kv in cfg.split()); out = os.path.join("/kaggle/working" if os.path.isdir("/kaggle/working") else ".", "seam_%s.txt" % tag)
        print("=== %s -> %s" % (cfg, out), flush=True)
        with open(out, "w") as f:
            pr = subprocess.Popen([sys.executable, __file__], env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
            for line in pr.stdout: print(line, end="", flush=True); f.write(line)
            pr.wait()
    raise SystemExit

IC = os.environ.get("IC", "found"); N = int(os.environ.get("N", 256)); NU = float(os.environ.get("NU", 2e-3)); T = float(os.environ.get("T", 3.0))
EVERY = float(os.environ.get("EVERY", 0.1)); DEV = "cuda" if torch.cuda.is_available() else "cpu"
SEPS = [float(v) for v in os.environ.get("SEPS", "0.05,0.1,0.2,0.4").split(",")]; VSEP = float(os.environ.get("VSEP", 0.1))
FOUND = os.environ.get("FOUND", "/kaggle/input/zef-found/leashed64_dmin030.npz")
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


def transport(U):
    Ud = [Ui * deal for Ui in U]; u = [ifft(Ui).real for Ui in Ud]; out = []
    for i in range(3):
        adv = sum(u[j] * ifft(1j * K[j] * Ud[i]).real for j in range(3))
        div = sum(ifft(1j * K[j] * fft(u[j] * u[i])).real for j in range(3))
        out.append(-0.5 * fft(adv + div) * deal)
        del adv, div
    return project(out)


def step(U, dt):
    f, f2 = torch.exp(-NU * K2 * dt), torch.exp(-NU * K2 * dt / 2)
    a = transport(U); b = transport([f2 * (U[i] + dt / 2 * a[i]) for i in range(3)])
    c = transport([f2 * U[i] + dt / 2 * b[i] for i in range(3)]); d = transport([f * U[i] + dt * f2 * c[i] for i in range(3)])
    return [f * U[i] + dt / 6 * (f * a[i] + 2 * f2 * b[i] + 2 * f2 * c[i] + d[i]) for i in range(3)]


def strip(U):
    e = 0.5 * sum(Ui.abs() ** 2 for Ui in U) / N**6
    spec = torch.stack([e[(KMAG >= n - 0.5) & (KMAG < n + 0.5)].sum() for n in range(1, NB)]).cpu().numpy()
    ks = np.arange(1, NB); sel = (ks >= NB // 2) & (spec > 1e-300)
    return -np.polyfit(ks[sel], np.log(spec[sel]), 1)[0] / 2 if sel.sum() > 4 else float("nan")


@torch.no_grad()
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
    sel = tw if tw.sum() > 0 else high
    re_seam = ((wm * (wm / gwm) ** 2)[sel].median() / NU).item() if NU > 0 else float("nan")     # |w| ell^2 / nu = Delta_u ell / nu on the seam
    valid = s_comp > 1e-6
    return (0.5 * w2.mean().item(), wm.max().item(), soft, anti[high].float().mean().item(), ell, ell_nu, cut, strip(U), softs, re_seam, valid)


x = torch.arange(N, device=DEV, dtype=RD) * dx
if IC == "pair":
    # two antiparallel Gaussian vortex tubes along x, separation D, core sigma, bowed toward each other (Kerr-type)
    D = float(os.environ.get("D", 0.7)); SIG = float(os.environ.get("SIG", 0.22)); A = float(os.environ.get("A", 0.2))
    X, Y, Z_ = torch.meshgrid(x, x, x, indexing="ij"); wx = torch.zeros_like(X)
    for sgn in (+1, -1):
        yc = math.pi + sgn * (D / 2 - A * torch.cos(X - math.pi)); r2 = ((Y - yc) % (2 * math.pi) - math.pi) ** 2 + ((Z_ - math.pi) % (2 * math.pi) - math.pi) ** 2
        wx = wx + sgn * torch.exp(-r2 / (2 * SIG**2))
    W = [fft(wx), torch.zeros_like(fft(wx)), torch.zeros_like(fft(wx))]
    W = project(W)                                                                     # divergence-free vorticity
    U = [(1j * (K[(i + 1) % 3] * W[(i + 2) % 3] - K[(i + 2) % 3] * W[(i + 1) % 3]) / K2S) * deal for i in range(3)]   # u = curl^-1 w
    U = project(U)
else:
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
    U = [Ui * math.sqrt(0.375 / Z0) for Ui in U]; Z0 = 0.375
print("seam race on GPU: IC=%s  N=%d^3  nu=%g  T=%.1f  Z0=%.3f  clock 2dx=%.4f  device=%s" % (IC, N, NU, T, Z0, 2 * dx, DEV), flush=True)
print("   t     Z/Z0    max|w|   twist@" + " @".join("%g" % v for v in SEPS) + "   anti    ell     ell_nu   race   cut    Re_seam   delta")
t, mark, t0 = 0.0, 0.0, time.time(); hist = []
while t <= T + 1e-9:
    if t >= mark - 1e-9:
        Z, wmax, soft, anti, ell, ell_nu, cut, d, softs, re_seam, valid = diag(U); hist.append((t, Z / Z0, wmax, soft, anti, ell, ell_nu, cut, d, re_seam, float(valid)))
        print("%5.2f   %6.3f   %7.2f   %s   %.3f   %.4f   %.4f   %5.2f   %.3f   %8.1f   %.4f%s   (%.0fs)" % (
            t, Z / Z0, wmax, " ".join("%.5f" % softs[v] for v in SEPS), anti, ell, ell_nu, ell / ell_nu if (NU > 0 and valid) else float("nan"), cut, re_seam, d,
            "" if d > 2 * dx else "  <-- past the clock", time.time() - t0), flush=True)
        mark += EVERY
        if t >= T - 1e-9: break
    with torch.no_grad():
        umax = max(ifft(Ui).real.abs().max().item() for Ui in U)
        dt = min(2.0 / N, 0.5 * dx / max(umax, 1e-9), mark - t + 1e-12)
        U = step(U, dt); t += dt
H = np.array(hist); tt, tw, race, d, res, valid = H[:, 0], H[:, 3], H[:, 5] / np.where(H[:, 6] > 0, H[:, 6], np.nan), H[:, 8], H[:, 9], H[:, 10] > 0.5
race = np.where(valid, race, np.nan)
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
    kill = (ipk == last) and rmin < 1.0
    print("REGISTERED C20 at nu=%g: %s" % (NU, "PASS: the seam reaches its viscous thickness and the twist turns over" if ok else ("KILL: the twist is still rising at the clock with the viscous scale resolved" if kill else "between (see rows)")))
else:
    print("REGISTERED C20 at nu=0: twist %s to the clock" % ("RISES" if ipk == last else "peaks at t=%.2f and falls %.0f%%" % (tpk, 100 * fall)))
