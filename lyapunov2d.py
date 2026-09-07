"""Positive control for the Lyapunov search: 2-D, where a monotone quantity is KNOWN to exist.

In 2-D the enstrophy Z = 1/2 int w^2 is conserved by Euler and decreases under viscosity - a local functional of the
state that never increases. The 3-D searches found no candidate in the local class and the adversary broke every one;
that verdict means something only if the same machine FINDS the known quantity in 2-D. Same construction as
lyapunov_search.py: M = Z_ref exp(Phi), Phi an enstrophy-weighted average of a small network over dimensionless local
features; dM/dt the exact Lie derivative by autograd against the 2-D vorticity equation (skew transport + exact
viscosity); an adversary searching low-k initial data for the trajectory along which M grows fastest; a strong
attack at the end. The prefactor is the PALINSTROPHY P = 1/2 int |grad w|^2, which grows transiently in 2-D (so a constant Phi fails),
and Phi may use the global scalars log(Z/P), log(E/Z) besides the local features; the known monotone quantity
Z = P (Z/P) is therefore reachable only by LEARNING Phi = log(Z/P) - the answer is available, not built in.
REGISTERED: PASS if the final candidate survives ATTACK (5 restarts x 60 iterations) with violation < 1e-3 and its
feature sensitivity is led by log(Z/P) (i.e. it rediscovered M = Z) or |w|^2; FAIL otherwise. A FAIL here retracts the 3-D verdict.
usage: N=64 NU=1e-3 T=1.0 ROUNDS=6 python lyapunov2d.py"""
import os, time, math, numpy as np, torch

torch.set_default_dtype(torch.float64)
torch.manual_seed(0)
N = int(os.environ.get("N", 64)); NU = float(os.environ.get("NU", 1e-3)); T = float(os.environ.get("T", 1.0))
ROUNDS = int(os.environ.get("ROUNDS", 6)); TRAIN = int(os.environ.get("TRAIN", 300)); ADV_ITERS = int(os.environ.get("ADV_ITERS", 20))
B = float(os.environ.get("B", 4.0)); TOL = float(os.environ.get("TOL", 1e-3)); RESTARTS = int(os.environ.get("RESTARTS", 5)); ATTACK_ITERS = int(os.environ.get("ATTACK_ITERS", 60))
fft, ifft = torch.fft.fftn, torch.fft.ifftn
k1 = torch.fft.fftfreq(N, d=1.0 / N) * 1.0
KX, KY = torch.meshgrid(k1, k1, indexing="ij")
K2 = KX**2 + KY**2; K2S = K2.clone(); K2S[0, 0] = 1.0
DEAL = ((KX.abs() < N / 3) & (KY.abs() < N / 3)).to(torch.float64)
LOWK = ((K2 <= 16) & (K2 > 0)).to(torch.float64)


def vel(wh):
    psih = wh / K2S
    return ifft(1j * KY * psih).real, -ifft(1j * KX * psih).real


def rhs(wh):
    """2-D vorticity equation: skew transport + exact viscous term, spectral tendency"""
    wh = wh * DEAL
    u, v = vel(wh); w = ifft(wh).real
    adv = u * ifft(1j * KX * wh).real + v * ifft(1j * KY * wh).real
    div = ifft(1j * KX * fft(u * w) + 1j * KY * fft(v * w)).real
    return -0.5 * fft(adv + div) * DEAL - NU * K2 * wh


def step(wh, dt):
    f, f2 = torch.exp(-NU * K2 * dt), torch.exp(-NU * K2 * dt / 2)
    tr = lambda V: (lambda Vd: -0.5 * fft(vel(Vd)[0] * ifft(1j * KX * Vd).real + vel(Vd)[1] * ifft(1j * KY * Vd).real + ifft(1j * KX * fft(vel(Vd)[0] * ifft(Vd).real) + 1j * KY * fft(vel(Vd)[1] * ifft(Vd).real)).real) * DEAL)(V * DEAL)
    a = tr(wh); b = tr(f2 * (wh + dt / 2 * a)); c = tr(f2 * wh + dt / 2 * b); d = tr(f * wh + dt * f2 * c)
    return f * wh + dt / 6 * (f * a + 2 * f2 * b + 2 * f2 * c + d)


def rollout(wh, T, every=None):
    t, out, nxt = 0.0, [], 0.0
    while t < T - 1e-12:
        if every is not None and t >= nxt - 1e-12:
            out.append((t, wh.detach().clone())); nxt += every
        u, v = vel(wh * DEAL); umax = float(torch.sqrt(u**2 + v**2).max())
        dt = min(2.0 / N, 0.5 * (2 * math.pi / N) / max(umax, 1e-9), T - t)
        wh = step(wh, dt); t += dt
    if every is not None:
        out.append((t, wh.detach().clone())); return out
    return wh


def enstrophy(wh): return 0.5 * (ifft(wh * DEAL).real ** 2).mean()
def energy(wh):
    u, v = vel(wh * DEAL); return 0.5 * (u**2 + v**2).mean()


FEATURES = ["|w|^2 / 2Z", "w / sqrt(2Z) (signed)", "|grad w|^2 / 2P", "|u|^2 / 2E", "|S|^2 / Z", "(w^2 - 2|S|^2)/Z", "|grad w|.u / (P E)^1/2", "G: log(Z/P)", "G: log(E/Z)"]


def features(wh):
    whd = wh * DEAL
    w = ifft(whd).real; u, v = vel(whd)
    Z = 0.5 * (w**2).mean(); E = 0.5 * (u**2 + v**2).mean()
    gx, gy = ifft(1j * KX * whd).real, ifft(1j * KY * whd).real
    P = 0.5 * (gx**2 + gy**2).mean()
    ux, uy, vx, vy = ifft(1j * KX * fft(u)).real, ifft(1j * KY * fft(u)).real, ifft(1j * KX * fft(v)).real, ifft(1j * KY * fft(v)).real
    s2 = ux**2 + vy**2 + 0.5 * (uy + vx) ** 2
    f = torch.stack([w**2 / (2 * Z), w / torch.sqrt(2 * Z), (gx**2 + gy**2) / (2 * P), (u**2 + v**2) / (2 * E), s2 / Z, (w**2 - 2 * s2) / Z, (gx * u + gy * v) / torch.sqrt(P * E),
                     torch.log(Z / P).expand_as(w), torch.log(E / Z).expand_as(w)], 0)
    return f, w**2 / (w**2).sum(), P


class G(torch.nn.Module):
    def __init__(self, nf=9, h=32):
        super().__init__()
        self.net = torch.nn.Sequential(torch.nn.Linear(nf, h), torch.nn.Tanh(), torch.nn.Linear(h, h), torch.nn.Tanh(), torch.nn.Linear(h, 1))
    def forward(self, f): return B * torch.tanh(self.net(f.permute(1, 2, 0)).squeeze(-1))      # Phi in [-B, B]: log(Z/P) is negative


g = G(); opt = torch.optim.Adam(g.parameters(), lr=3e-3)


def M_of(wh):
    f, wgt, P = features(wh)
    Phi = (wgt * g(f)).sum()
    return P * torch.exp(Phi), Phi, P


def violation(wh):
    wr = ifft(wh).real.detach().clone().requires_grad_(True)
    whh = fft(wr)
    M, Phi, Z = M_of(whh)
    grad = torch.autograd.grad(M, wr, create_graph=True)[0]
    F = ifft(rhs(whh.detach())).real
    return (grad * F).sum() / M, M


x = torch.linspace(0, 2 * math.pi, N + 1)[:-1]; X, Y = torch.meshgrid(x, x, indexing="ij")
rng = torch.Generator().manual_seed(0)


def field_from_params(P):
    wh = fft(P) * LOWK; wh = wh * DEAL
    return wh * torch.sqrt(torch.tensor(0.375) / enstrophy(wh))


def random_ic(seed): return field_from_params(torch.randn(N, N, generator=torch.Generator().manual_seed(seed)))


tg = fft(2 * torch.sin(X) * torch.sin(Y)); tg = tg * torch.sqrt(torch.tensor(0.375) / enstrophy(tg))
train_ics = {"taylor-green": tg, **{"random-%d" % s: random_ic(s) for s in range(5)}}
heldout = {"random-9": random_ic(9), "random-11": random_ic(11)}
EVERY = T / 8
print("2-D Lyapunov search (positive control): N=%d^2 nu=%g T=%.1f rounds=%d train=%d adversary iters=%d; M = P exp(Phi), P = palinstrophy (grows), global scalars available" % (N, NU, T, ROUNDS, TRAIN, ADV_ITERS), flush=True)


def trajectories(ics):
    out = []
    with torch.no_grad():
        for name, wh in ics.items():
            for (t, S) in rollout(wh, T, every=EVERY): out.append((name, t, S))
    return out


data = trajectories(train_ics); t0 = time.time()


def report(states, label):
    worst = (-1e9, None, None)
    for name, t, S in states:
        v, M = violation(S); v = v.item()
        if v > worst[0]: worst = (v, name, t)
    print("  %-20s worst relative dM/dt = %+.4e  (%s at t=%.2f)" % (label, *worst), flush=True); return worst[0]


def adversary(iters, seed, lr=0.05):
    torch.manual_seed(seed); P = torch.randn(N, N).requires_grad_(True); aopt = torch.optim.Adam([P], lr=lr); best = (-1e9, None)
    for it in range(iters):
        wh = field_from_params(P); worst = None; t = 0.0
        while t < T - 1e-12:
            u, v = vel(wh * DEAL); umax = float(torch.sqrt(u**2 + v**2).max())
            dt = min(2.0 / N, 0.5 * (2 * math.pi / N) / max(umax, 1e-9), T - t)
            wh = step(wh, dt); t += dt
            if int(t / EVERY) != int((t - dt) / EVERY):
                wr = ifft(wh).real; whh = fft(wr); M, Phi, Z = M_of(whh)
                grad = torch.autograd.grad(M, wr, create_graph=True)[0]
                v_ = (grad * ifft(rhs(whh)).real).sum() / M
                worst = v_ if worst is None else torch.maximum(worst, v_)
        aopt.zero_grad(); (-worst).backward(); aopt.step()
        if worst.item() > best[0]: best = (worst.item(), P.detach().clone())
    return best


for rnd in range(ROUNDS):
    for it in range(TRAIN):
        batch = [data[i] for i in torch.randperm(len(data))[:6].tolist()]
        loss = sum(torch.relu(violation(S)[0]) ** 2 for _, _, S in batch) / len(batch)
        opt.zero_grad(); loss.backward(); opt.step()
    print("== round %d" % rnd, flush=True)
    report(data, "training states"); report(trajectories(heldout), "held-out flows")
    bv, bP = adversary(ADV_ITERS, 100 + rnd)
    print("  adversary best violation %+.4e   (%.0fs)" % (bv, time.time() - t0), flush=True)
    with torch.no_grad():
        for (t, S) in rollout(field_from_params(bP), T, every=EVERY): data.append(("adversary-%d" % rnd, t, S))

print("\n== ATTACK: %d restarts x %d iterations against the final candidate" % (RESTARTS, ATTACK_ITERS), flush=True)
worst_attack = max(adversary(ATTACK_ITERS, 500 + r)[0] for r in range(RESTARTS))
sens = torch.zeros(9); cnt = 0
for name, t, S in trajectories(heldout):
    f, wgt, Z = features(S); f = f.detach().requires_grad_(True)
    gr = torch.autograd.grad((wgt * g(f)).sum(), f)[0]; sens += (gr * f.detach()).abs().sum(dim=(1, 2)); cnt += 1
print("feature sensitivity of the learned Phi (held-out):")
order = sorted(zip(FEATURES, (sens / cnt).tolist()), key=lambda z: -z[1])
for name, s_ in order: print("   %-36s %.3e" % (name, s_))
led = order[0][0].startswith("G: log(Z/P)") or order[0][0].startswith("|w|^2")
print("\nREGISTERED  attack worst violation %+.3e (tol %.0e); leading feature: %s -> %s" % (worst_attack, TOL, order[0][0], "PASS: the machine finds the known monotone quantity (P x Z/P = Z) in 2-D" if worst_attack < TOL and led else ("FAIL: the machine does not find enstrophy where it exists - the 3-D verdict is retracted" if worst_attack >= TOL else "between: survives the attack but is not enstrophy-led")))
torch.save(g.state_dict(), "results/lyapunov2d_g.pt")
