"""Figures for the seam race from the GPU run logs (results/seam_gpu/v3, v4, v5): the twist wave descending through the
scales, its nu-independence, the V (gap arm against thickness arm), and the time per octave for three structures.
Writes figures/twist_wave.png, figures/seam_V.png, figures/octaves.png. usage: python plot_seam.py (seconds)"""
import re, glob, os, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

SEPS = [0.05, 0.07, 0.1, 0.14, 0.2, 0.28, 0.4, 0.56]


def load(f):
    rows = [l.split() for l in open(f) if re.match(r"^ *\d+\.\d+ +\d", l)]; rows = [r for r in rows if len(r) >= 18]
    t = np.array([float(r[0]) for r in rows]); tw = np.array([[float(r[3 + j]) for j in range(8)] for r in rows])
    d = np.array([float(r[17]) if r[17][0].isdigit() else np.nan for r in rows]); wmax = np.array([float(r[2]) for r in rows])
    return t, tw, d, wmax


def peaks(t, tw, d, dx=2 * np.pi / 256):
    ok = d > 2 * dx; last = np.where(ok)[0][-1]; out = {}
    for j, s in enumerate(SEPS):
        col = tw[:last + 1, j]; ip = int(np.argmax(col))
        if 0 < ip < last: out[s] = t[ip]
    return out, t[last]


os.makedirs("figures", exist_ok=True)
runs = {2e-3: "results/seam_gpu/v3/seam_ICfound_NU2e-3_T2.4.txt", 1e-3: "results/seam_gpu/v3/seam_ICfound_NU1e-3_T2.0.txt", 5e-4: "results/seam_gpu/v3/seam_ICfound_NU5e-4_T1.6.txt"}

# 1. the twist wave: each separation's twist against time, nu = 2e-3, with the other viscosities dashed
fig, ax = plt.subplots(figsize=(8.5, 5))
cols = plt.cm.viridis(np.linspace(0.05, 0.95, 8))
for nu, ls, lw in ((2e-3, "-", 2.0), (1e-3, "--", 1.0), (5e-4, ":", 1.0)):
    t, tw, d, _ = load(runs[nu]); ok = d > 2 * 2 * np.pi / 256; last = np.where(ok)[0][-1]
    for j, s in enumerate(SEPS):
        ax.plot(t[:last + 1], tw[:last + 1, j], ls, color=cols[j], lw=lw, label=("sep %.2f" % s) if nu == 2e-3 else None)
ax.set_xlabel("t"); ax.set_ylabel("strong twist at fixed separation"); ax.set_yscale("log"); ax.set_ylim(1e-4, 1)
ax.set_title("The twist wave: the reversal descends through the scales\nsolid nu = 2e-3, dashed 1e-3, dotted 5e-4 (inside each run's clock)")
ax.legend(ncol=2, fontsize=8, loc="lower right"); ax.grid(alpha=0.3)
fig.tight_layout(); fig.savefig("figures/twist_wave.png", dpi=140); plt.close(fig)

# 2. the V: peak time vs separation (the gap arm) with the fitted line, and delta(t) (the thickness arm), nu = 2e-3 and 1e-3
fig, ax = plt.subplots(figsize=(8, 5))
mk = {2e-3: "o", 1e-3: "s", 5e-4: "^"}
pool = []
for nu in (2e-3, 1e-3, 5e-4):
    t, tw, d, _ = load(runs[nu]); pk, tc = peaks(t, tw, d)
    ax.plot([pk[s] for s in pk], [s for s in pk], mk[nu], ms=7, alpha=0.8, label="twist peaks, nu = %g (clock %.1f)" % (nu, tc))
    pool += [(pk[s], s) for s in pk if s >= 0.1 and nu in (2e-3, 1e-3)]
    if nu in (2e-3, 1e-3):
        m = (t >= 0.3) & (d > 2 * 2 * np.pi / 256) & (t <= 1.7)
        ax.plot(t[m], d[m], "-" if nu == 2e-3 else "--", color="gray", lw=1.5, label="thickness arm delta(t), nu = %g" % nu)
P = np.array(pool); a, b = np.polyfit(P[:, 0], P[:, 1], 1); tt = np.linspace(0.6, -b / a, 50)
ax.plot(tt, a * tt + b, "k-", lw=1, label="gap arm: %.2f (%.2f - t)" % (-a, -b / a))
ax.axhline(0.053, color="C3", ls=":", lw=1, label="sqrt(nu/s) at the merge, nu = 2e-3")
ax.set_xlabel("t"); ax.set_ylabel("separation / thickness"); ax.set_ylim(0, 0.62); ax.set_xlim(0.5, 1.9)
ax.set_title("The V: the gap between the sheets (inviscid, linear)
meets their thickness (viscous, exponential)")
ax.legend(fontsize=8); ax.grid(alpha=0.3)
fig.tight_layout(); fig.savefig("figures/seam_V.png", dpi=140); plt.close(fig)

# 3. time per octave against scale for the three structures
fig, ax = plt.subplots(figsize=(8, 5))
series = {"sheet field (adversary), nu = 2e-3": ("results/seam_gpu/v3/seam_ICfound_NU2e-3_T2.4.txt", "C0"),
          "sheet field, nu = 1e-3": ("results/seam_gpu/v3/seam_ICfound_NU1e-3_T2.0.txt", "C0"),
          "second adversarial field, nu = 1e-3": ("results/seam_gpu/v4/seam_ICfound_TAGckn_NU1e-3_T2.0.txt", "C1"),
          "tube pair (approaching), nu = 1e-3": ("results/seam_gpu/v5/seam_ICpair_NU1e-3_T8.txt", "C2")}
for name, (f, c) in series.items():
    if not os.path.exists(f): continue
    t, tw, d, _ = load(f); pk, _ = peaks(t, tw, d)
    ss = sorted([s for s in pk if s >= 0.1], reverse=True); xs, ys = [], []
    for s1, s2 in zip(ss[:-1], ss[1:]):
        dt = pk[s2] - pk[s1]
        if dt > 0: xs.append(np.sqrt(s1 * s2)); ys.append(dt / np.log2(s1 / s2))
    ax.plot(xs, ys, "o-" if "2e-3" in name or "pair" in name or "second" in name else "s--", color=c, label=name)
ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("separation (geometric mean of the octave)"); ax.set_ylabel("time per octave of descent")
xx = np.array([0.1, 0.5]); ax.plot(xx, 1.1 * xx, "k:", lw=1, label="proportional to the separation (rate ~ 1/ell, self-induced)")
ax.set_title("How fast the wave descends
equal octaves = phase 1 (external strain); shrinking = phase 2 (self-induced)")
ax.legend(fontsize=8); ax.grid(alpha=0.3, which="both")
fig.tight_layout(); fig.savefig("figures/octaves.png", dpi=140); plt.close(fig)
print("wrote figures/twist_wave.png, figures/seam_V.png, figures/octaves.png")
