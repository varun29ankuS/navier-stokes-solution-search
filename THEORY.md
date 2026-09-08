# Isometry plus contraction: a formal note

This note states precisely what "freeze the conservative part, let only the dissipation act" means, proves the
elementary consequences, records the classical theorems that bound what it can do, and states the one hypothesis
under which it would decide regularity. Status of each statement is marked: **[proved here]**, **[classical]**,
**[open]**.

---

## 1. Setting

Let $H$ be a real Hilbert space with inner product $\langle\cdot,\cdot\rangle$ and norm $\|\cdot\|$. Consider an
evolution

$$\dot x  =  T(x)  +  D(x), \qquad x(0)=x_0 ,$$

with the two parts characterised by their action on the norm:

* **Transport** $T$ is *skew* (norm-preserving): $\langle x, T(x)\rangle = 0$ for all $x$ in its domain.
* **Dissipation** $D$ is *contractive*: $\langle x, D(x)\rangle = - \mathcal{E}(x) \le 0$, with $\mathcal{E}\ge 0$
  called the dissipation rate.

Examples. Burgers: $T(u) = -u u_x$, $D(u)=\nu u_{xx}$ on the torus; $\mathcal{E}(u)=\nu\|u_x\|^2$. Navier-Stokes on
the torus $\mathbb{T}^d$, divergence-free velocity fields: $T(u) = -P[(u\cdot\nabla)u]$ with $P$ the Leray projector,
$D(u)=\nu\Delta u$, $\mathcal{E}(u)=\nu\|\nabla u\|^2$. A recurrent memory: $x_{t+1} = R_t x_t$ with $R_t$
orthogonal is the discrete transport; a gate $x \mapsto (1-g)\odot x$ with $g\in[0,1]$ is the discrete dissipation.

**Definition 1 (entropy production).** For a differentiable flow map $\Phi_t$ on a finite-dimensional state space,
the entropy production rate is the volume contraction rate

$$\sigma(x)  =  - \nabla\cdot F(x) \qquad\text{for } \dot x = F(x),$$

and for a discrete map $x_{t+1}=\Phi(x_t)$ it is $\sigma = -\log|\det J_\Phi(x_t)|$. (This is the phase-space
contraction rate of Ruelle; for a thermostatted system it equals the thermodynamic entropy production.)

**Definition 2 (numerical entropy production).** For a numerical scheme applied to a system whose exact dynamics
satisfies $\frac{d}{dt}\|x\|^2 = -2\mathcal{E}(x)$, the numerical entropy production over $[0,t]$ is

$$\sigma_{\text{num}}  =  -\frac{1}{t}\Big(\log\frac{\|x_t\|^2}{\|x_0\|^2}  +  \frac{2}{\|x_0\|^2}\int_0^t \mathcal{E} ds\Big),$$

the dissipation the scheme produces beyond the physical one. For $\mathcal{E}\equiv 0$ it reduces to
$-\tfrac1t\log(\|x_t\|^2/\|x_0\|^2)$, the quantity measured in `burgers_entropy.py`.

---

## 2. What isometry plus contraction gives

**Proposition 1 (energy inequality). [proved here; classical for NS, Leray 1934]**
For any solution of $\dot x = T(x)+D(x)$ on $[0,t]$ smooth enough for the identities to hold,

$$\|x(t)\|^2 + 2\int_0^t \mathcal{E}(x(s)) ds  =  \|x_0\|^2 .$$

*Proof.* $\frac{d}{dt}\|x\|^2 = 2\langle x,\dot x\rangle = 2\langle x,T(x)\rangle + 2\langle x,D(x)\rangle = 0 - 2\mathcal{E}(x)$. Integrate. $\square$

**Proposition 2 (discrete energy inequality; zero numerical entropy). [proved here]**
Let a scheme advance $x_n \mapsto x_{n+1}$ by a map that is the composition of an exact isometry $R_n$
($\|R_n x\|=\|x\|$) and a contraction $C_n$ ($\|C_n x\|\le\|x\|$), in either order. Then
$\|x_{n+1}\|\le\|x_n\|$ for every $n$; the discrete solution is bounded for all $n$ by $\|x_0\|$; and if $C_n$
is the identity the scheme has $\sigma_{\text{num}} = 0$ exactly.

*Proof.* $\|C_n R_n x\| \le \|R_n x\| = \|x\|$, and likewise in the other order. Induct. If $C_n=\mathrm{id}$, the
norm is constant and $\sigma_{\text{num}}=0$ by Definition 2. $\square$

*Remark.* The pseudo-spectral scheme in this repository realises Proposition 2 up to time-stepping error: the
skew-symmetric form makes the semi-discrete transport exactly skew ($\langle \hat u, \hat T(\hat u)\rangle = 0$ on
the dealiased modes), the viscous integrating factor exp(-\nu k^2\Delta t) is an exact contraction per mode, and
RK4 adds an $O(\Delta t^4)$ deviation. Measured: $\sigma_{\text{num}} = -3e-15 (Burgers), energy conserved
to 1e-6 (3-D Euler).

**Proposition 3 (a frozen isometric memory carries a conserved quantity exactly; a dissipative one forgets at
the entropy production rate). [proved here]**
(i) Let $x_{t+1} = R_t x_t$ with each $R_t$ orthogonal, and let $\delta_t$ be the difference between two
trajectories with the same $R_t$. Then $\|\delta_t\| = \|\delta_0\|$ for all $t$.
(ii) Let $x_{t+1} = (1-g_t)\odot R_t x_t + b_t$ with $g_t\in[0,1]^d$ (gated recurrence with input $b_t$). Then
$\|\delta_t\| \le \prod_{s<t}(1-\min_i g_{s,i}) \|\delta_0\|$, and the per-step entropy production is
$\sigma_t = -\sum_i \log(1-g_{t,i})$.
(iii) (Position clock.) For rotation planes with angular rates $w_k$ driven by increments $\Delta_t$,
$\theta_{k,t+1} = \theta_{k,t} + w_k\Delta_t$, the phases satisfy $\theta_{k,t} = \theta_{k,0} + w_k\sum_{s<t}\Delta_s$
exactly, i.e. the state is an exact function of the accumulated quantity, for every $t$, with no trained parameter.
If the increments carry i.i.d. noise of variance $\varsigma^2$, the decoded quantity has error of order
$\varsigma\sqrt{t}$ (a random walk); if the true value is re-observed every $G$ steps, the error is of order
$\varsigma\sqrt{G}$, uniformly in $t$.

*Proof.* (i) $\delta_{t+1} = R_t\delta_t$, orthogonal. (ii) $\delta_{t+1} = (1-g_t)\odot R_t\delta_t$, and
$\|(1-g)\odot y\|\le (1-\min_i g_i)\|y\|$; the Jacobian is $\mathrm{diag}(1-g_t)R_t$ with
$\log|\det| = \sum_i\log(1-g_{t,i})$. (iii) Telescoping sum; the noise statements are the variance of a sum of
$t$ (resp. at most $G$) independent increments. $\square$

*Remark.* Part (ii) is the "forgetting law" and the measured $\sim 250$ nats/step of the trained drawing hand; part
(iii) is the measured 1.8 px at 4000 moves (no noise), $\sigma\sqrt{n}$ drift under poison, and the flat 6 px with a
glance every 50 moves. Note that in the Householder update $x\mapsto x-\beta u u^{\top}x$ the map is not
orthogonal unless $\beta\in\{0,2\}$: $\det = 1-\beta = \cos\theta$, so the "transport" of that architecture is
itself a contraction along $u$ and contributes $\log|\cos\theta|$ to $\sigma$. The exact-transport claim applies to
the rotation planes, not to the Householder erasure.

---

## 3. What isometry plus contraction cannot give

**Theorem 4 (the abstraction is insufficient for regularity). [classical: Tao 2016]**
There exists an averaged Navier-Stokes system $\dot u = \tilde T(u) + \nu\Delta u$ on $\mathbb{T}^3$ with $\tilde T$
skew (so Proposition 1 holds verbatim, with the same energy identity as Navier-Stokes) and smooth divergence-free
initial data whose solution blows up in finite time.

*Consequence.* No argument that uses only the structure "skew transport + contractive dissipation" and the energy
inequality can prove global regularity for 3-D Navier-Stokes, since the same argument would apply to Tao's system.
Any proof must use properties of the specific nonlinearity $P[(u\cdot\nabla)u]$ that are not consequences of
skewness.

*Reference.* T. Tao, *Finite time blowup for an averaged three-dimensional Navier-Stokes equation*, J. Amer. Math.
Soc. 29 (2016), 601-674.

**Theorem 5 (what suffices: a second controlled norm). [classical; conditional form proved here]**
Let $u$ be a Leray-Hopf solution on $\mathbb{T}^3\times[0,T)$ that is smooth on $[0,T)$. If

$$\int_0^T \|\omega(s)\|_{L^\infty} ds < \infty, \qquad \omega=\nabla\times u,$$

then $u$ extends smoothly past $T$ (Beale-Kato-Majda 1984). Consequently: **if there is a functional $M(u)\ge 0$ with
(a) $M(u) \ge c \|\omega\|_{L^\infty}$ (or any bound implying the BKM integral is finite), and (b)
$\frac{d}{dt}M(u)\le 0$ along solutions - i.e. $M$ is a norm in which the transport is skew or contractive and the
dissipation contractive - then solutions are globally smooth.**

*Proof of the consequence.* (b) gives $M(u(t))\le M(u_0)$ on $[0,T)$; (a) gives $\|\omega\|_{L^\infty}\le M(u_0)/c$;
the BKM integral is at most $T M(u_0)/c<\infty$; apply BKM. $\square$

*Remark (why 2-D is solved and 3-D is not).* In two dimensions enstrophy $Z=\tfrac12\|\omega\|^2$ satisfies
$\frac{d}{dt}Z = -\nu\|\nabla\omega\|^2 \le 0$: the vortex-stretching term $\int \omega\cdot(\omega\cdot\nabla)u$
vanishes identically, so the transport is skew in a *second* norm that controls the gradient, and global regularity
follows (Ladyzhenskaya 1959; for the periodic setting see Temam's text). In three dimensions the same computation gives
$\frac{d}{dt}Z = \int\omega\cdot(\omega\cdot\nabla)u - \nu\|\nabla\omega\|^2$, and the first term has no sign. The
transport is an isometry of $L^2$ only. Whether any functional satisfying (a) and (b) exists is **[open]**; it is
equivalent in spirit to the Millennium problem, and Theorem 4 says it cannot be found by abstract arguments alone.

---

## 4. What this repository establishes, stated exactly

1. The scheme is an instance of Proposition 2: its transport is skew on the dealiased modes and its dissipation is
   an exact contraction, so numerical entropy production is zero up to RK4 error. **[proved here + measured]**
2. On 1-D Burgers it reproduces the exact approach to the known singularity, $\max|u_x| = 1/(1-t)$, until the grid
   limit, which appears as a visible shortfall rather than as damping. **[measured]**
3. On 2-D Navier-Stokes (Taylor-Green) it reproduces a closed-form solution to 1e-14; its dissipation is exactly
   the physical $\nu\|\nabla u\|^2$. **[measured]**
4. On 3-D Euler (Taylor-Green) energy is conserved to 1e-6 at $32^3, 48^3, 64^3$, and enstrophy agrees across
   resolutions until the cascade reaches the grid scale, after which the resolutions *disagree* rather than
   converge to a wrong value. **[measured]**
5. None of 1-4 bears on the regularity question except as an instrument free of the artefact (numerical
   dissipation) that limits numerical searches for singular solutions. **[Theorem 4]**

## 5. The one open hypothesis, stated so that it can be attacked or refuted

**Hypothesis H.** There exists a functional $M$ on divergence-free fields on $\mathbb{T}^3$, satisfying (a) and (b)
of Theorem 5 along Navier-Stokes solutions.

If H holds, 3-D regularity follows by Theorem 5. If a smooth solution blows up, H is false. Theorem 4 shows H cannot
be established from skewness and dissipation alone; a candidate $M$ must use the structure of
$P[(u\cdot\nabla)u]$. The numerical route to evidence is to compute, for candidate $M$, the sign of
$\frac{d}{dt}M$ along flows near the strongest known amplification events, with a scheme whose own dissipation is
zero so that the sign is not an artefact - which is what Proposition 2 provides.

## 6. What a proof must control: the criteria, and what is measured

Each of the following is a theorem giving a sufficient condition for smoothness; a blow-up must violate all of them.
The right-hand column is what `criteria3d.py` measures along a flow, by resolution, with the budgets closed.

| criterion | statement | measured |
|---|---|---|
| Beale-Kato-Majda (1984) | int_0^T max\|w\| dt < inf implies smooth on [0,T] | max\|w\|(t) and its time integral |
| Ladyzhenskaya-Prodi-Serrin; Escauriaza-Seregin-Sverak (2003) | sup_t \|\|u\|\|_L3 < inf implies smooth (the scale-critical norm) | \|\|u(t)\|\|_L3 |
| Constantin-Fefferman (1993) | if the vorticity direction xi = w/\|w\| is Lipschitz where \|w\| is large, stretching is depleted and the solution is smooth | direction coherence rho = <1 - (xi(x).xi(x+h))^2> over \|w\| > 0.5 max, h = dx; and the local stretching alpha = xi.S.xi there |
| Caffarelli-Kohn-Nirenberg (1982) | the singular set has one-dimensional parabolic Hausdorff measure zero | (a blow-up, if any, is at points, not sheets or lines) |
| Tao (2016) | the skew-plus-dissipation structure alone cannot decide regularity | (limits what any argument built on Sections 1-2 can prove) |

A regularity proof shows one of the first three always holds. A blow-up construction violates all of them at a point
and uses the exact nonlinearity. The measurements do neither; they report which criterion is tightest on the classical
candidate flows (Taylor-Green; the perturbed ABC flow; the Kida-Pelz high-symmetry flow, whose apparent blow-up
dissolved at high resolution, Hou and Li 2006) and whether the critical norm stays bounded on them.

## 7. Further viewpoints, with what each lets us measure

**Topology: helicity and the linking of vortex lines.** Helicity $H = \int u\cdot\omega$ is a topological invariant
of ideal flow: it measures the linking of vortex lines (Moffatt 1969) and is conserved by Euler (measured here to six
digits under a seven-fold enstrophy growth). It bounds energy from below: $E \ge c\lvert H\rvert$ for knotted
vorticity (Arnold 1974; Freedman and He 1991), so linked vortex lines cannot relax away. Moffatt's conjecture is that
helicity inhibits the cascade. Two of our runs point the same way: the helical ABC flow amplified enstrophy seven-fold
where the zero-helicity Kida-Pelz flow amplified it twenty-five-fold. **Measured** (`adversarial_ic.py`,
`OBJ=helicity`): the maximal enstrophy amplification attainable at fixed initial enstrophy when the relative helicity
$H / 2\sqrt{EZ}$ is constrained to 0, 0.5, 0.9 of its Beltrami maximum. The registered prediction was that the fastest
amplifier is helicity-free and that attainable growth falls as helicity is imposed. **Result** (hard constraint,
32^3 truncated system, one time unit): 4.54, 4.68, 4.47, 3.11, 1.27 at relative helicity 0, 0.25, 0.5, 0.75, 0.9.
The prediction held in direction and failed in shape: helicity is a threshold, not a slope. Half the maximal
helicity costs nothing; the inhibition arrives between 0.5 and 0.9. A Beltrami-like field ($H/H_{\max} = 0.9$)
still amplifies enstrophy 1.27-fold, resolved, more than any classical flow at the same $Z_0$.

**Riemannian geometry: Euler as a geodesic.** Arnold (1966): an ideal flow is a geodesic on the group of
volume-preserving diffeomorphisms with the kinetic-energy metric; Ebin and Marsden (1970) proved local existence
this way. **Global regularity is exactly geodesic completeness** - whether the geodesic can be continued for all
time. The sectional curvature of the group along the flow is mostly negative (Arnold; Lukatskii), so nearby geodesics
diverge exponentially: Lagrangian unpredictability is curvature. The vortex-stretching term is how that curvature acts
on the velocity field. **Measured** (`OBJ=jacobi`): the maximal growth of a perturbation along the flow - the Jacobi
field, found by maximising $\lVert\delta u(T)\rVert / \lVert\delta u(0)\rVert$ through the differentiable solver - by
resolution and by flow. A blow-up would be a geodesic reaching the boundary of the group in finite time; a rapidly
growing Jacobi field is necessary for it and far from sufficient.

**The search itself: feedback, not imitation.** Lu and Doering (2008) and Ayala and Protas (2017) found, by adjoint
optimisation, initial data whose enstrophy grows far faster than any classical candidate. `adversarial_ic.py`
reproduces that search with autograd through the conserving solver and re-verifies every candidate at higher
resolution with the analyticity-strip clock, so a field that only wins by exploiting the search grid is rejected.
First run (Kida-Pelz enstrophy normalisation, 32^3 search, 64^3 verification): the found field amplified 8.75x on its
own grid and 20x at 64^3, with the analyticity strip collapsed to zero - outside the reliable window, therefore not
yet a result; the searcher is now run at a milder amplitude and verified at 96^3 and 128^3. The same machinery
pointed at $\delta$ itself is the honest numerical form of "look for a singularity": ask the optimiser to shrink the
analyticity strip, and let a verifier that cannot be fooled by the grid say whether it did.

**Thermodynamics of the truncated system: Liouville.** Write the Galerkin-truncated equation as $\dot U = F(U)$ on
the finite-dimensional phase space $\mathcal{P}$ of retained divergence-free modes.

**Proposition 6 (Lee 1952; measured here).** For $\nu = 0$, $\operatorname{div}_{\mathcal P} F = 0$: the flow on
$\mathcal P$ preserves phase volume, and the Gibbs entropy $S = -\int \rho\log\rho$ of any ensemble of solutions is
constant. For $\nu > 0$, $\operatorname{div}_{\mathcal P} F = -\nu (d-1) \sum_{k \in \mathcal K} \lvert k\rvert^2$,
a constant, so $S(t) = S(0) - \nu (d-1)\big(\sum_k \lvert k\rvert^2\big)\, t$ for every ensemble.

*Proof.* The quadratic term in mode $k$ is $\sum_{p+q=k} B_k(\hat u_p, \hat u_q)$; its derivative with respect to
$\hat u_k$ itself keeps only the terms with $p = k$ or $q = k$, i.e. $q = 0$ or $p = 0$, and the mean mode is zero
for a mean-free flow (or is not advected). The Leray projection and the dealiasing mask are linear and idempotent
and commute with this. So the diagonal of the Jacobian of the quadratic term vanishes identically, and the diagonal
of the viscous term is $-\nu\lvert k\rvert^2$ on each of the $d-1$ solenoidal directions per retained wavevector.
$\square$

This holds for the skew form, the advective form and the divergence form alike: Liouville is a property of the
convolution structure, not of energy conservation, and the advective form (which loses energy) preserves phase
volume. Measured (`liouville.py`, exact Jacobian trace by autograd): $-2.9\cdot10^{-16}$ (2-D) and
$+3.3\cdot10^{-16}$ (3-D) inviscid; $-24.20$ and $-82.32$ viscous, equal to the formula to every digit.

The consequence worth stating: the truncated inviscid system does not create entropy in *any* sense, per solution
(Proposition 2) or per ensemble (Proposition 6). What it does is stretch phase volume along the least stable
directions and compress it along the others, at the rate the Jacobi field measures (`jacobi_ladder.py`: local
exponent rising from 0.03 to 0.28 along Taylor-Green as the sheets form, 0.47 for Kida-Pelz at $t = 0.5$). The
information is not lost; it is moved below the scale a finite observer can read, and the analyticity-strip clock is
exactly the statement of when that has happened. A proof of regularity would have to control where the phase volume
goes, not how much of it there is; that is Theorem 4 again in thermodynamic dress.

**Complex singularities: the strip as a Riemann-surface statement.** The analyticity strip $\delta(t)$ is the
distance from the real domain of the nearest singularity of the solution continued to complex space (Sulem, Sulem
and Frisch 1983). A real blow-up at $t^*$ is $\delta(t^*) = 0$. The two classical decay laws, exponential (no real
singularity; Taylor-Green, Brachet et al.) and linear (singularity at $t^*$), are distinguishable inside the
reliable window by the sign of the change of the local decay rate $-\,d\log\delta/dt$: constant for exponential,
rising as $1/(t^* - t)$ for linear. **Measured** (`strip_tracker.py`) by flow and by resolution.

Status of all of the above: **[measured or in progress]**; none bears on Theorem 4, which stands.


## 8. The pressure conjecture, stated precisely (from the searches of 2026-09-06/07)

**Setting.** $u$ a smooth solution of 3-D Euler on $\mathbb{T}^3$, $\omega = \nabla\times u$, $\xi = \omega/|\omega|$,
$S = \tfrac12(\nabla u + \nabla u^{T})$, and $p$ the pressure, $-\Delta p = \partial_i u_j\,\partial_j u_i$. Write the
pressure Hessian $P = \nabla\nabla p$ and its traceless part $P^{\circ} = P - \tfrac13(\operatorname{tr}P)\,I$; the trace
is local ($\operatorname{tr}P = \Delta p = \tfrac12|\omega|^2 - |S|^2$), the traceless part is the nonlocal image of the
whole field. For a set $A\subset\mathbb{T}^3$ define the **global share**

$$ s_A(u) \;=\; \frac{\int_A |P^{\circ}|^2}{\int_A |P|^2}\,, \qquad
   h_A(u) \;=\; \frac{\int_A \xi\cdot P^{\circ}\xi}{\int_A \xi\cdot S^2\xi}\,, $$

with $A = \{|\omega| > \tfrac12 \max|\omega|\}$ the high-vorticity set; $h_A < 0$ means the nonlocal part of the pressure
Hessian pushes with the self-stretching along the vorticity.

**Observed** (`results/quiet*_pressure_*.txt`; search over data with $\hat u$ supported in $|k|\le 4$, enstrophy
$Z_0$ fixed, on $32^3$ and $48^3$ grids, over $T = 0.5, 1, 1.5$; verifier at $64$-$96^3$): the maximal amplification
$Z(T)/Z_0$ attainable with $s_A(u(T)) \le \sigma$ imposed is

| $\sigma$ | 0.5 | 0.45 | 0.4 | 0.35 | 0.3 |
|---|---|---|---|---|---|
| $\sup Z(T)/Z_0$, $48^3$, $T=1$ | 5.6 | 4.2 | 2.1 | 1.11 | 1.10 |

with the same knee and collapse at $32^3$ and at $T = 0.5, 1.5$; the searcher cannot realise $s_A < 0.35$ at all in
this class; and $h_A < 0$ in every fast field found (eleven of eleven). Along the fastest field, $h_A$ grows from
$-0.5$ to $-6.7$ up to the last resolved time and does not fade at any resolved scale (`help_fades_*`).

**Conjecture C17 (necessary condition).** There is $\sigma_* \in (0.35, 0.45)$ and a function $f$ with $f(\sigma)\to 1$
as $\sigma \downarrow \sigma_*$ such that for smooth Euler solutions on $\mathbb{T}^3$ with $\hat u_0$ supported in
$|k| \le K$,
$$ \sup_{t \le T}\, s_A(u(t)) \;\le\; \sigma \quad\Longrightarrow\quad \frac{Z(T)}{Z_0} \;\le\; f(\sigma)\,, $$
and moreover any solution with $Z(T)/Z_0 \ge 2$ has $h_A(u(t)) < 0$ on a set of times of positive measure in $[0,T]$.
In words: enstrophy cannot grow substantially unless the nonlocal part of the pressure carries at least a fixed share of
the Hessian on the high-vorticity set and acts with the stretching.

**What it would and would not give.** C17 is a necessary condition for growth, not a bound: the searcher's fields
satisfy it and keep satisfying it. A regularity argument along this line would need a second statement,
$$ \text{(C17b)}\qquad \int_A \big(\xi\cdot P^{\circ}\xi\big)_{-} \;\le\; C\, E^{a} Z^{b}\,, $$
a bound on the nonlocal help in terms of quantities the energy inequality controls, with exponents making the pair
critical. `help_fades` says C17b is *not* visible at any resolved scale in the class searched: the help grows to the
clock. If C17b holds it holds in the limit. Refutation of C17: a smooth datum with $Z(T)/Z_0 \ge 3$ and
$\sup_t s_A < 0.35$, or a fast amplifier with $h_A \ge 0$ throughout. Both are one run of `adversarial_ic.py`
(`OBJ=quiet`).

Status: **[conjecture with a stable curve; not a theorem; the closing statement C17b is unsupported at resolved scales]**.

## 9. The race, stated as a theorem (2026-09-08), and the hypothesis that would finish it

The runs of 2026-09-08 (`docs/LOG.md`, "The descent law and the V") measured, on the adversary's sheet field at
$256^3$, that the gap between the two antiparallel sheets closes linearly, $g(t) \approx 0.51\,(T^\ast - t)$, with
the peak times independent of viscosity to one sampling step, and that the pair merges where the gap meets the sheets'
thickness - at the viscous scale $\sqrt{\nu/s}$ for $\nu = 2\times10^{-3}$, where enstrophy and $\max|\omega|$ peak and
decay inside the reliable window. The following is the part of that picture that is a theorem.

**Theorem 7 (the race between a collapse scale and the viscous scale). [elementary; proved here]**
Let $u$ be a smooth solution of Navier-Stokes on $\mathbb{T}^3\times[0,T)$ and set $s(t) = \sup_x |\nabla u(x,t)|$,
the viscous scale $\ell_\nu(t) = \sqrt{\nu/s(t)}$. Suppose

$$\text{(Type I)}\qquad s(t) \;\le\; \frac{A}{T-t}\,,$$

and let $\ell(t)$ be any length attached to the solution (a gap, a core radius, a sheet thickness) with

$$\ell(t) \;\le\; L\,(T-t)^{\lambda}\qquad\text{on } [t_0, T).$$

Then

$$\frac{\ell(t)}{\ell_\nu(t)} \;\le\; L\sqrt{\frac{A}{\nu}}\;(T-t)^{\lambda - \frac12}\,,$$

so that: if $\lambda > \tfrac12$ the structure is inside its own viscous scale before $T$ ($\ell/\ell_\nu \to 0$); if
$\lambda = \tfrac12$ the ratio is bounded by the constant $L\sqrt{A/\nu}$, a Reynolds number of the structure; if
$\lambda < \tfrac12$ no such conclusion holds.

*Proof.* $\ell_\nu(t) = \sqrt{\nu/s(t)} \ge \sqrt{\nu (T-t)/A}$ by Type I. Divide. $\square$

**What it says and what it does not.** With the measured $\lambda = 1$ the sheet pair's gap falls below the viscous
scale before $T^\ast$ whatever the constants: a seam that closes linearly cannot outrun viscosity - which is what the
resolved runs show, and why the fastest collapse in Euler is the one viscosity kills most surely. The theorem does
**not** conclude regularity: "inside the viscous scale" must be handed to a criterion that finishes the argument
(Caffarelli-Kohn-Nirenberg's $\varepsilon$-regularity on the local dissipation, or Constantin-Fefferman on the
direction), and that hand-off is the real work. Nor is Type I known in general (it is excluded only under axisymmetry:
Chen-Strain-Tsai-Yau 2008, Koch-Nadirashvili-Seregin-Šverák 2009). What the theorem fixes is the *shape* of the
dangerous case: $\lambda = \tfrac12$ exactly - gap and floor descending in lockstep - with the outcome decided by one
number, the structure's Reynolds number at the turnover, which `seam_gpu.py` reports as `Re_seam` (C21). The Caltech
travelling profile of 2026-09-07 sits at $\lambda = 0.5$ by construction.

**Conjecture C25 (the seam's velocity jump is bounded by the data).** For smooth solutions whose high-vorticity set is a
pair of antiparallel sheets with gap $g(t)$ and velocity jump $\Delta u(t)$ across the pair, the self-induced strain
across the seam satisfies $s(t) \le C\,\Delta u(t)/g(t)$ (Biot-Savart for a sheet pair), and

$$\Delta u(t) \;\le\; C_0 \quad\text{on } [0,T)\quad\text{with } C_0 = C_0(u_0).$$

*If C25 held it would finish the seam class:* $\dot g \ge -c\,\Delta u \ge -c\,C_0$ gives $g(t) \ge g(t_0) - c C_0 (t -
t_0)$, i.e. $\lambda \ge 1$ for the gap, and Theorem 7 puts the pair inside its viscous scale before any collapse; the
measured $\lambda = 1$ is this bound saturated. *Why it is a conjecture and not a lemma:* a bounded velocity jump is a
bounded velocity, and bounded velocity implies regularity (Serrin). C25 is the wall, located: it names the quantity that
stayed bounded in every resolved run ($\Delta u \sim |\omega|\,\ell$ on the seam, `Re_seam`$\cdot\nu/\ell$), and it says
that proving it bounded *for the seam class only* - a geometric assumption in exchange for a global one - is what a
structural proof would have to do. Refutation of C25 (numerical): a resolved run in which $\Delta u$ across the seam
grows without bound while the gap closes - $\lambda < 1$ with `Re_seam` rising through the turnover. The Re ladder
(384^3 at $\nu = 10^{-3}$) is that test.

Status: **[Theorem 7 proved (elementary); C25 stated; the hand-off from "inside the viscous scale" to regularity not
carried out]**.
