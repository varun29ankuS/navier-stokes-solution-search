"""Regularity needs a quantity that is both MONOTONE and DOMINATING; do the two sets intersect?

The equation:  find M[u] >= 0 with dM/dt <= 0 along the flow and M >= c |grad u|_inf.
Along trajectories, for a family of candidates, measure
    monotone   :  max over time of (M(t+dt) - M(t)) / (M(t) dt), the worst relative growth rate  (<= ~0: monotone)
    dominating :  min over time of M(t) / max|grad u|(t)^q  with q chosen so the ratio is dimensionless in the
                  regime of interest; if the ratio collapses toward 0 while the gradient grows, M does not dominate.
REGISTERED: 1-D inviscid Burgers (a theorem: it shocks) - the monotone set (int|u|, max|u|, total variation, by
Kruzkov) and the dominating set (int u_x^2, max|u_x|) are DISJOINT. 2-D Navier-Stokes (regular) - the sets INTERSECT:
enstrophy is monotone and dominates. Then 3-D is the question "do they intersect", already known not to within
every family this repository tried. usage: python monotone_vs_dominating.py   (1-D instant; 2-D a minute)"""
import numpy as np

print("=== 1-D inviscid Burgers, u0 = sin x + 0.5 sin(2x+1), run to 0.92 t*  (exact shock time t* = -1/min u0')")
N = 2048; x = np.linspace(0, 2 * np.pi, N, endpoint=False); k = np.fft.fftfreq(N, d=1.0 / N); deal = np.abs(k) < N / 3
u0 = np.sin(x) + 0.5 * np.sin(2 * x + 1); tstar = -1.0 / (np.cos(x) + np.cos(2 * x + 1)).min()
rhs = lambda uh: -np.fft.fft((np.fft.ifft(uh).real * np.fft.ifft(1j * k * uh).real + np.fft.ifft(1j * k * np.fft.fft(np.fft.ifft(uh).real ** 2)).real) / 3.0) * deal
uh = np.fft.fft(u0); t = 0.0; dt = 2e-4; hist = []
while t < 0.92 * tstar:
    a = rhs(uh); b = rhs(uh + dt / 2 * a); c = rhs(uh + dt / 2 * b); d = rhs(uh + dt * c); uh = uh + dt / 6 * (a + 2 * b + 2 * c + d); t += dt
    if len(hist) == 0 or t - hist[-1][0] >= 0.02:
        u = np.fft.ifft(uh).real; ux = np.fft.ifft(1j * k * uh).real
        hist.append((t, {"int|u|": np.mean(np.abs(u)), "max|u|": np.abs(u).max(), "total variation": np.mean(np.abs(ux)), "int u_x^2": np.mean(ux**2), "max|u_x|": np.abs(ux).max()}, np.abs(ux).max()))
names = list(hist[0][1].keys())
print("%-18s  %-22s  %-26s  verdict" % ("candidate", "worst growth rate", "M / max|grad|: first -> last"))
for nme in names:
    vals = np.array([h[1][nme] for h in hist]); ts = np.array([h[0] for h in hist]); g = np.array([h[2] for h in hist])
    growth = (np.diff(vals) / (vals[:-1] * np.diff(ts))).max()
    dom = vals / g
    mono = growth < 0.05; dominates = dom[-1] / dom[0] > 0.3
    print("%-18s  %+22.3f  %12.4f -> %-10.4f  %s" % (nme, growth, dom[0], dom[-1], ("MONOTONE" if mono else "grows") + " / " + ("DOMINATES" if dominates else "does not dominate")))
print("max|u_x| went %.1f -> %.1f (exact 1/(1-t/t*): %.1f)" % (hist[0][2], hist[-1][2], 1.0 / (1 - 0.92) * hist[0][2] / hist[0][2] * (-1.0 / (np.cos(x) + np.cos(2 * x + 1)).min()) ** 0 * hist[0][2]))
print("1-D verdict: monotone set and dominating set %s" % ("DISJOINT (as registered: Burgers shocks)" if not any((np.diff([h[1][n] for h in hist]) / (np.array([h[1][n] for h in hist][:-1]) * np.diff([h[0] for h in hist]))).max() < 0.05 and (np.array([h[1][n] for h in hist]) / np.array([h[2] for h in hist]))[-1] / (np.array([h[1][n] for h in hist]) / np.array([h[2] for h in hist]))[0] > 0.3 for n in names) else "INTERSECT"))

print("\n=== 2-D Navier-Stokes, decaying turbulence 128^2, nu = 2e-4, t = 0..3")
M = 128; nu = 2e-4; dt = 2e-3
kk = np.fft.fftfreq(M, d=1.0 / M); kx, ky = np.meshgrid(kk, kk, indexing="ij"); k2 = kx**2 + ky**2; k2s = k2.copy(); k2s[0, 0] = 1.0
deal2 = (np.abs(kx) < M / 3) & (np.abs(ky) < M / 3)
def vel(wh):
    psih = wh / k2s; psih[0, 0] = 0
    return np.fft.ifft2(1j * ky * psih).real, -np.fft.ifft2(1j * kx * psih).real
def rhs2(wh):
    u, v = vel(wh); w = np.fft.ifft2(wh).real
    return -0.5 * (np.fft.fft2(u * np.fft.ifft2(1j * kx * wh).real + v * np.fft.ifft2(1j * ky * wh).real) + 1j * kx * np.fft.fft2(u * w) + 1j * ky * np.fft.fft2(v * w)) * deal2
rng = np.random.default_rng(1); wh = np.fft.fft2(rng.standard_normal((M, M))) * ((k2 >= 9) & (k2 <= 36)) * deal2
u, v = vel(wh); wh *= 1 / np.sqrt(np.mean(u * u + v * v))
f, f2 = np.exp(-nu * k2 * dt), np.exp(-nu * k2 * dt / 2); t = 0.0; hist = []
while t < 3.0:
    a = rhs2(wh); b = rhs2(f2 * (wh + dt / 2 * a)); c = rhs2(f2 * wh + dt / 2 * b); d = rhs2(f * wh + dt * f2 * c); wh = f * wh + dt / 6 * (f * a + 2 * f2 * b + 2 * f2 * c + d); t += dt
    if len(hist) == 0 or t - hist[-1][0] >= 0.1:
        u, v = vel(wh); w = np.fft.ifft2(wh).real; gx, gy = np.fft.ifft2(1j * kx * wh).real, np.fft.ifft2(1j * ky * wh).real
        G = np.sqrt(sum(np.fft.ifft2(1j * ka * np.fft.fft2(f_)).real ** 2 for ka in (kx, ky) for f_ in (u, v))).max()
        hist.append((t, {"energy int|u|^2": np.mean(u * u + v * v), "max|u|": np.sqrt(u * u + v * v).max(), "enstrophy int w^2": np.mean(w * w), "max|w|": np.abs(w).max(), "palinstrophy int|grad w|^2": np.mean(gx**2 + gy**2)}, G))
names = list(hist[0][1].keys())
print("%-26s  %-22s  %-26s  verdict" % ("candidate", "worst growth rate", "M / max|grad u|: first -> last"))
inter = False
for nme in names:
    vals = np.array([h[1][nme] for h in hist]); ts = np.array([h[0] for h in hist]); g = np.array([h[2] for h in hist])
    growth = (np.diff(vals) / (vals[:-1] * np.diff(ts))).max(); dom = vals / g
    mono = growth < 0.05; dominates = dom[-1] / dom[0] > 0.3; inter |= (mono and dominates)
    print("%-26s  %+22.3f  %12.4f -> %-10.4f  %s" % (nme, growth, dom[0], dom[-1], ("MONOTONE" if mono else "grows") + " / " + ("DOMINATES" if dominates else "does not dominate")))
print("max|grad u| went %.2f -> %.2f" % (hist[0][2], hist[-1][2]))
print("2-D verdict: monotone set and dominating set %s" % ("INTERSECT (as registered: 2-D is regular)" if inter else "DISJOINT"))
