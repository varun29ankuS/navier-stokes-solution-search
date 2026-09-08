"""figures/seam_race.gif from results/seam_gpu/snaps/snaps_found_nu0.002.npz: a slice of |w| through the vorticity
maximum (left), the signed reversal beta = 1 - xi.xi' at separation 0.1 on the same slice (right; > 1 is antiparallel,
the seam), and under them the strong twist at four separations with the frame's time marked and the clock. The seam
forms, descends through the scales, meets the viscous thickness and is cut. usage: python seam_gif.py (a minute)"""
import re, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

NU = 0.002
d = np.load("results/seam_gpu/snaps/snaps_found_nu%g.npz" % NU)
t, wm, beta, maxw, delta = d["t"], d["wm"], d["beta"], d["maxw"], d["delta"]
rows = [l.split() for l in open("results/seam_gpu/snaps/seam_ICfound_NU2e-3_T2.4_SNAP1.txt") if re.match(r"^ *\d+\.\d+ +\d", l)]; rows = [r for r in rows if len(r) >= 18]
tt = np.array([float(r[0]) for r in rows]); tw = {s: np.array([float(r[3 + j]) for r in rows]) for j, s in enumerate([0.05, 0.07, 0.1, 0.14, 0.2, 0.28, 0.4, 0.56])}
clock = tt[np.array([float(r[17]) if r[17][0].isdigit() else 1.0 for r in rows]) > 2 * 2 * np.pi / 256][-1]
n = len(t); vmax = np.percentile(wm, 99.7)

fig = plt.figure(figsize=(10, 7.2)); gs = fig.add_gridspec(2, 2, height_ratios=[3, 1.3], hspace=0.28, wspace=0.12)
axw, axb, axt = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[1, :])
imw = axw.imshow(wm[0].T, origin="lower", cmap="magma", vmin=0, vmax=vmax, extent=[0, 2 * np.pi, 0, 2 * np.pi])
imb = axb.imshow(np.clip(beta[0].T, 0, 2), origin="lower", cmap="RdBu_r", vmin=0, vmax=2, extent=[0, 2 * np.pi, 0, 2 * np.pi])
axw.set_title("|w| on the plane through its maximum"); axb.set_title("signed reversal beta at separation 0.1  (red > 1: the seam)")
for a in (axw, axb): a.set_xticks([]); a.set_yticks([])
cols = {0.4: "C8", 0.2: "C2", 0.1: "C0", 0.05: "C4"}
for s, c in cols.items(): axt.plot(tt, tw[s], color=c, lw=1.5, label="twist@%.2f" % s)
axt.axvline(clock, color="k", ls=":", lw=1); axt.text(clock, 0.5, " clock", fontsize=8, va="top")
axt.set_yscale("log"); axt.set_ylim(1e-4, 1); axt.set_xlim(0, tt[-1]); axt.set_xlabel("t"); axt.set_ylabel("strong twist"); axt.legend(fontsize=8, ncol=4, loc="lower right"); axt.grid(alpha=0.3)
cursor = axt.axvline(t[0], color="C3", lw=1.5)
title = fig.suptitle("", fontsize=11)


def update(i):
    imw.set_data(wm[i].T); imb.set_data(np.clip(beta[i].T, 0, 2)); cursor.set_xdata([t[i], t[i]])
    title.set_text("The seam race, nu = %g, 256^3 (slice %d/256)   t = %.2f   max|w| = %.1f   strip delta = %.3f%s" % (
        NU, d["iz"][i], t[i], maxw[i], delta[i], "" if delta[i] > 2 * 2 * np.pi / 256 else "   PAST THE CLOCK"))
    return imw, imb, cursor, title


anim = FuncAnimation(fig, update, frames=n, blit=False)
anim.save("figures/seam_race.gif", writer=PillowWriter(fps=6), dpi=80)
print("wrote figures/seam_race.gif with %d frames" % n)
