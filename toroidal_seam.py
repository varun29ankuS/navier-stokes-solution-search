"""Does the merging seam close into a toroidal (non-radiating) pair? A CPU test on the found field.

Hypothesis (2026-09-08): a flat antiparallel pair, stabilised against roll-up, becomes a self-screened object - a
toroidal vorticity distribution whose induced field is confined (the anapole of Zel'dovich / Afanasiev-Stepanovsky;
Hill's vortex, 1894). The GPU rows already show the self-screening (the pair's induced rate -> 0 at the merge) and a
radiation BURST at the cut, not a sealing. This asks the topological half: does poloidal circulation - fluid
circulating AROUND the pair, in the plane normal to the sheets' vorticity - grow as the pair's induction dies?
Measured on the tagged fluid: at every output, the material velocity decomposed against the frame (xi_ref = the
pair's vorticity direction; n = the sheet normal from the centroid offset; p = xi_ref x n), and the two circulations
    C_along  = median over particles of u . xi_ref-tangent loop  (Kelvin's loop: the jump; conserved)
    C_around = mean of (r x u) . xi_ref / |r|^2 over the tagged set = the poloidal angular velocity about the pair's axis
and the anisotropy of the tagged cloud (sheet-like: one small eigenvalue; toroidal: two comparable).
REGISTERED (C37): C_around stays below 0.1 x |w|_m through the merge and the cloud stays sheet-like (smallest / middle
eigenvalue ratio < 0.3): the pair merges flat, no toroidal closure. Refuted by: C_around rising to >= 0.3 |w|_m as the
pair's induction dies, with the cloud rounding (ratio > 0.5) - a toroidal pair.
usage: N=96 NU=2e-3 T=2.2 TSEED=1.0 python toroidal_seam.py      (CPU, ~10 minutes at 96^3)"""
import os, time, math, numpy as np

N = int(os.environ.get("N", 96)); NU = float(os.environ.get("NU", 2e-3)); T = float(os.environ.get("T", 2.2)); TSEED = float(os.environ.get("TSEED", 1.0))
NP = int(os.environ.get("NP", 3000)); EVERY = 0.1
fft, ifft = np.fft.fftn, np.fft.ifftn
k = np.fft.fftfreq(N, d=1.0 / N); kx, ky, kz = np.meshgrid(k, k, k, indexing="ij"); K = [kx, ky, kz]
k2 = kx**2 + ky**2 + kz**2; k2s = k2.copy(); k2s[0, 0, 0] = 1.0
deal = (np.abs(kx) < N / 3) & (np.abs(ky) < N / 3) & (np.abs(kz) < N / 3)
L = 2 * np.pi; dx = L / N


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


def interp(F, P):
    g = (P / dx) % N; i0 = np.floor(g).astype(int); fr = g - i0; i1 = (i0 + 1) % N; out = 0.0
    for a in (0, 1):
        for b in (0, 1):
            for c in (0, 1):
                w = (fr[:, 0] if a else 1 - fr[:, 0]) * (fr[:, 1] if b else 1 - fr[:, 1]) * (fr[:, 2] if c else 1 - fr[:, 2])
                out = out + w * F[(i1 if a else i0)[:, 0], (i1 if b else i0)[:, 1], (i1 if c else i0)[:, 2]]
    return out


def fields(U):
    Ud = [Ui * deal for Ui in U]; u = [ifft(Ui).real for Ui in Ud]
    w = [ifft(1j * K[(i + 1) % 3] * Ud[(i + 2) % 3] - 1j * K[(i + 2) % 3] * Ud[(i + 1) % 3]).real for i in range(3)]
    return u, w


def pcentroid(P):
    return np.array([np.arctan2(np.sin(P[:, c]).mean(), np.cos(P[:, c]).mean()) % L for c in range(3)])


uf = np.load(os.environ.get("FOUND", "results/found/leashed64_dmin030.npz"))["u"].astype(float); n0 = uf.shape[1]; U = []
for c in range(3):
    uh = fft(uf[c]) * (N / n0) ** 3; big = np.zeros((N, N, N), complex); h = n0 // 2
    for a in (slice(0, h), slice(-h, None)):
        for b in (slice(0, h), slice(-h, None)):
            for cc in (slice(0, h), slice(-h, None)):
                big[a, b, cc] = uh[a, b, cc]
    U.append(big)
U = project(U); _, w0 = fields(U); Z0 = 0.5 * np.mean(sum(wi**2 for wi in w0)); U = [Ui * math.sqrt(0.375 / Z0) for Ui in U]
print("toroidal seam test: N=%d^3 nu=%g T=%.1f seed at %.1f; clock 2dx = %.4f" % (N, NU, T, TSEED, 2 * dx), flush=True)
print("   t     |w|_m    pair induction   C_along(jump)   C_around/|w|_m   cloud axes (small/mid/large)   ratio   delta")
t, mark, t0 = 0.0, 0.0, time.time(); P = None; hist = []
while t <= T + 1e-9:
    if t >= mark - 1e-9:
        u, w = fields(U); wm = np.sqrt(sum(wi**2 for wi in w))
        e = 0.5 * sum(np.abs(Ui) ** 2 for Ui in U) / N**6; kmag = np.sqrt(k2); nb = int(N / 3)
        spec = np.array([e[(kmag >= n - 0.5) & (kmag < n + 0.5)].sum() for n in range(1, nb)]); ks = np.arange(1, nb); sel = (ks >= nb // 2) & (spec > 1e-12 * spec.max())
        delta = -np.polyfit(ks[sel], np.log(spec[sel]), 1)[0] / 2 if sel.sum() > 4 else float("nan")
        if P is None and t >= TSEED - 1e-9:
            high = np.argwhere(wm > 0.5 * wm.max()); rng = np.random.default_rng(0); sel_ = high[rng.choice(len(high), min(NP, len(high)), replace=False)]
            P = (sel_ + 0.5) * dx
            xi = np.stack([w[c][high[:, 0], high[:, 1], high[:, 2]] for c in range(3)], 1); xi /= np.linalg.norm(xi, axis=1, keepdims=True) + 1e-30
            xi_ref = np.linalg.eigh(xi.T @ xi)[1][:, -1]
            sgn = np.sign(sum(w[c][sel_[:, 0], sel_[:, 1], sel_[:, 2]] * xi_ref[c] for c in range(3)))
            print("seeded %d particles (A/B %d/%d), xi_ref = (%.2f, %.2f, %.2f)" % (len(P), (sgn > 0).sum(), (sgn < 0).sum(), *xi_ref), flush=True)
        if P is not None:
            wmp = np.median(interp(wm, P)); cen = pcentroid(P); rel = (P - cen + np.pi) % L - np.pi
            up = np.stack([interp(ui, P) for ui in u], 1)
            # pair induction: the velocity induced by the high set alone, projected on the separation of the two sheets
            mask = 1 / (1 + np.exp(-(wm / wm.max() - 0.5) / 0.1)); W = project([fft(wi * mask) for wi in w])
            uH = [ifft((1j * (K[(i + 1) % 3] * W[(i + 2) % 3] - K[(i + 2) % 3] * W[(i + 1) % 3]) / k2s) * deal).real for i in range(3)]
            A, B = P[sgn > 0], P[sgn < 0]; dvec = (A[:, None] - B[None]); dvec = (dvec + np.pi) % L - np.pi; dist = np.sqrt((dvec**2).sum(-1)); j = dist.argmin(1)
            dhat = dvec[np.arange(len(A)), j] / (dist[np.arange(len(A)), j][:, None] + 1e-30)
            hA = np.stack([interp(hi, A) for hi in uH], 1); hB = np.stack([interp(hi, B[j]) for hi in uH], 1); induct = np.median(((hA - hB) * dhat).sum(1))
            # Kelvin's loop along the sheets: the jump of the tangential velocity across the pair
            uA = np.stack([interp(ui, A) for ui in u], 1); uB = np.stack([interp(ui, B[j]) for ui in u], 1); c_along = np.median(np.linalg.norm(uA - uB, axis=1))
            # poloidal: angular velocity of the tagged fluid about the pair's axis (xi_ref through the centroid)
            r_perp = rel - np.outer(rel @ xi_ref, xi_ref); r2 = (r_perp**2).sum(1) + 1e-12
            c_around = np.mean(np.cross(r_perp, up) @ xi_ref / r2)
            ev = np.sort(np.linalg.eigvalsh(np.cov(rel.T)))
            row = (t, wmp, induct, c_along, c_around / wmp, np.sqrt(ev), np.sqrt(ev[0] / ev[1]), delta); hist.append(row)
            print("  %.2f   %6.2f   %+.4f          %.4f          %+.4f           %.3f / %.3f / %.3f          %.2f    %.4f%s" % (t, wmp, induct, c_along, c_around / wmp, *np.sqrt(ev), np.sqrt(ev[0] / ev[1]), delta, "" if delta > 2 * dx else "  <-- past the clock"), flush=True)
        else:
            print("  %.2f   (not yet seeded)   delta %.4f%s   (%.0fs)" % (t, delta, "" if delta > 2 * dx else "  <-- past the clock", time.time() - t0), flush=True)
        mark += EVERY
        if t >= T - 1e-9: break
    umax = max(np.abs(ifft(Ui).real).max() for Ui in U); dt = min(2.0 / N, 0.5 * dx / max(umax, 1e-9), mark - t + 1e-12)
    if P is not None:
        u, _ = fields(U); v1 = np.stack([interp(ui, P) for ui in u], 1); Pm = (P + 0.5 * dt * v1) % L
        v2 = np.stack([interp(ui, Pm) for ui in u], 1); P = (P + dt * v2) % L
    U = step(U, dt); t += dt
H = np.array([(h[0], h[1], h[2], h[3], h[4], h[6], h[7]) for h in hist]); ok = H[:, 6] > 2 * dx
ind = H[:, 2]; ca = H[:, 4]; ratio = H[:, 5]
i_die = np.where(ok & (np.abs(ind) < 0.2 * np.abs(ind[0])))[0]
print("\nC37: pair induction %+.4f at seed -> %+.4f at the clock (dies below 20%% at t = %s); C_around/|w|_m max inside the clock %.3f; cloud ratio %.2f -> %.2f" % (
    ind[0], ind[ok][-1], ("%.2f" % H[i_die[0], 0]) if len(i_die) else "never", np.abs(ca[ok]).max(), ratio[0], ratio[ok][-1]))
tor = np.abs(ca[ok]).max() >= 0.3 and ratio[ok][-1] > 0.5; flat = np.abs(ca[ok]).max() < 0.1 and ratio[ok][-1] < 0.3
print("REGISTERED C37: %s" % ("KILL: a toroidal pair forms" if tor else ("PASS: the pair merges flat, no toroidal closure" if flat else "between")))
