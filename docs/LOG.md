# zero-entropy-flow - the research log

*This is the chronological log, written forward as the work happened, retractions included. The front page is
[../README.md](../README.md); the view from above is [MAP.md](MAP.md); the registered claims are
[../CLAIMS.md](../CLAIMS.md).*

# zero-entropy-flow

**One pseudo-spectral solver for the incompressible Navier-Stokes family in 1-D, 2-D and 3-D, built so that its own
numerical dissipation is zero** - the transport is exactly energy-conserving on the grid, viscosity is an exact
per-mode contraction - and then used to look at the things numerical dissipation normally hides: a singularity
forming, vortex stretching, the terms of every balance law.

What is here, each with the command that reproduces it and the observation that would refute it (`CLAIMS.md`):

- **Verified against five closed-form solutions**, one per dimension and then some, to machine precision:
  inviscid Burgers approaching its blow-up exactly as 1/(1-t); the Cole-Hopf viscous shock; 2-D Taylor-Green;
  a *random* 2-D exact solution; the 3-D viscous ABC flow. Errors 1e-12 to 1e-15 except through the shock, where
  the error is the grid's and says so.
- **Every budget closes term by term** - energy, enstrophy, palinstrophy - including the unsigned production terms,
  at residuals 1e-6 to 1e-9. The unsigned term is measured driving the 1-D blow-up, vanishing identically in 2-D,
  and in 3-D **vortex stretching is measured directly at 64^3, 96^3, 128^3 and converges upward** with resolution;
  energy 1.000000 throughout. Helicity conserved to six digits under 3-D dynamics that amplify enstrophy seven-fold.
- **A learning result:** a coarse simulator with frozen exact transport plus a learned correction that may *move*
  energy between scales but cannot *create* it cuts the long-rollout spectrum error six-fold against a fully learned
  model that drifts to 1.9x the true energy (`kolmogorov_v2.py`, reproduced on two machines). A dissipative-only
  learned closure did nothing in 2-D; the missing physics there is backscatter.
- **The regularity criteria, measured with a clock that says when to stop believing them:** BKM, the L3 norm,
  the Constantin-Fefferman direction coherence at fixed physical scale, and the analyticity-strip width delta(t),
  on Taylor-Green and Kida-Pelz at 64/96/128^3. Numbers past the clock were withdrawn from this page, and it says so.
- **Feedback, not imitation:** a differentiable copy of the solver searches initial data for the fastest enstrophy
  growth (Lu-Doering; Ayala-Protas) and a verifier at higher resolution with the clock rejects what only wins on the
  search grid - it has rejected three so far. The same search with helicity held fixed orders attainable growth by
  helicity (Moffatt's direction), and along Arnold's geodesic it measures the Jacobi-field exponent by resolution.
- **Liouville, exactly:** the truncated inviscid system conserves phase volume (Lee 1952) - exact Jacobian trace
  1e-16 by autograd - and the viscous contraction is a state-independent constant, so the Gibbs entropy of any
  ensemble falls linearly. Zero entropy production per solution *and* per ensemble; what the flow does is stretch.
- **Theory** (`THEORY.md`): the elementary propositions proved; Tao's 2016 theorem that this structure alone cannot
  decide 3-D regularity; the conditional theorem for what would; the open hypothesis stated so it can be attacked.
- **Not claimed:** anything about the regularity of 3-D Navier-Stokes.

Everything runs on numpy, on a laptop or on this repository's GitHub Actions runners, which produced the results in
`results/`. Issues and Discussions are open; refutations welcome, bring the output.

## Two animations

![2-D turbulence at 256^2, vorticity, with energy and enstrophy on every frame](../figures/turbulence_2d.gif)

*2-D decaying turbulence, 256^2, nu = 1e-4, vorticity. Vortices merge, filaments stretch and roll up, enstrophy
cascades to small scales; the energy on each frame decays only by the physical viscous rate. `flow_gif.py`.*

![1-D Burgers: inviscid blow-up (stopped at the grid limit) and the viscous shock on the Cole-Hopf exact solution](../figures/burgers_shock.gif)

*Left: inviscid Burgers, energy 1.000000000000 on every frame, the gradient climbing toward 1/(1-t) until the grid
limit, where the animation stops and says so. Right: viscous Burgers on top of the Cole-Hopf closed form through the
shock, error printed per frame.*

## Three pictures

![1-D Burgers: the frozen-rotation scheme sits on the exact 1/(1-t) blow-up curve while upwind drifts and loses energy](../figures/burgers_blowup.png)

*1-D Burgers. The scheme with zero numerical entropy production (black) lies on the exact blow-up curve (grey) until the
grid can no longer resolve the gradient; the standard dissipative scheme (orange) departs from it and, inset, loses the
energy the equation says is conserved.*

![3-D Euler: vortex stretching and enstrophy at three resolutions, converging upward](../figures/stretching_3d.png)

*3-D Euler, Taylor-Green. The vortex-stretching term - the one that decides the 3-D question - measured directly with the
budget closed, at 64^3, 96^3, 128^3. Energy is 1.000000 throughout. Coarser grids under-estimate the amplification; the
curves converge upward.*

![3-D Taylor-Green: vorticity isosurfaces in the symmetry cell and planar slices at t = 1..4](../figures/vortex_sheets_3d.png)

*3-D Euler, Taylor-Green at 48^3. Top: isosurface of |ω| in the symmetry cell [0,π]^3; bottom: |ω| on a plane. The
initial tori become sheets (t = 2) that thin and stretch (t = 3) and fold to the grid scale (t = 4) - the last frame
visibly under-resolved, as the resolution table says. Energy 1.000000 throughout. `vortex_iso.py`.*

![2-D turbulence: enstrophy decreasing while the unsigned palinstrophy production is large and positive](../figures/ladder_2d.png)

*2-D decaying turbulence. The unsigned production term one level above enstrophy is large and positive throughout, and
enstrophy falls anyway: regularity is not the absence of the dangerous term but the presence of a controlled norm one
level below it. Figures are regenerated by `plots.py` from the same code as the tables.*

## Poke holes
Every claim, what would refute it, the command that tests it, and the list of our own weak points: [CLAIMS.md](CLAIMS.md).
Issues and Discussions are open. Refutations are welcome; bring the output.

## The equations

**Burgers (1-D):**

$$u_t + u u_x = \nu u_{xx}, \qquad x \in [0, 2\pi),\ \text{periodic},\ u(x,0)=\sin x .$$

For $\nu = 0$ the solution is given implicitly by characteristics, $u = \sin(x - u t)$, and its gradient blows up in
finite time:

$$\max_x |u_x(\cdot,t)| = \frac{1}{1-t}, \qquad t^* = 1 .$$

Energy $E(t) = \tfrac{1}{2}\langle u^2\rangle$ is exactly conserved for $\nu=0$ while the solution is smooth, and for
$\nu>0$ obeys

$$\frac{dE}{dt} = - \nu \langle u_x^2\rangle .$$

**Navier-Stokes / Euler (2-D vorticity form and 3-D velocity form):**

$$\omega_t + (\mathbf{u}\cdot\nabla) \omega = \nu \Delta\omega, \qquad \mathbf{u} = \nabla^\perp\psi,\ \ \Delta\psi = \omega \quad (\text{2-D})$$

$$\mathbf{u}_t + (\mathbf{u}\cdot\nabla) \mathbf{u} = -\nabla p + \nu \Delta\mathbf{u}, \qquad \nabla\cdot\mathbf{u}=0 \quad (\text{3-D})$$

with the Taylor-Green initial condition $\mathbf{u}_0 = (\sin x\cos y\cos z,\ -\cos x\sin y\cos z,\ 0)$. In 2-D the
viscous Taylor-Green flow is exact: $u = e^{-2\nu t}\sin x\cos y,\ v = -e^{-2\nu t}\cos x\sin y$, so $E(t)=E_0 e^{-4\nu t}$.
In 3-D with $\nu=0$ energy is conserved exactly and enstrophy $Z = \tfrac12\langle|\omega|^2\rangle$ grows by vortex
stretching; whether that growth stays finite for all time is the open question.

## The scheme, and what "numerical entropy production" means

Write the equation as **conservative transport + dissipation**:

$$\partial_t \mathbf{u} = \underbrace{\mathcal{T}(\mathbf{u})}_{\text{conserves } E}  +  \underbrace{\nu \Delta\mathbf{u}}_{\text{dissipates}} .$$

The transport term is discretised so that the semi-discrete system conserves energy **exactly**: in Fourier space each
linear mode advances by an exact rotation (a unitary factor exp(-i k c \Delta t), which is an isometry), and the
nonlinearity is written in skew-symmetric form so that $\langle \mathbf{u}, \mathcal{T}(\mathbf{u})\rangle = 0$
identically. For Burgers:

$$u u_x  =  \tfrac{1}{3}\Big(u u_x + (u^2)_x\Big) ,$$

which conserves $\langle u^2\rangle$ term by term on the dealiased grid (2/3 rule). Viscosity enters as an exact
integrating factor exp(-\nu k^2 \Delta t) per mode; time stepping of the nonlinear term is RK4.

For any integrator define the **numerical entropy production rate** as the drift of the energy the equation says must
be conserved, per unit time:

$$\sigma_{\text{num}}  =  -\frac{1}{t} \log\frac{E(t)}{E(0)}\Big|_{\nu=0} .$$

For an exact isometric scheme $\sigma_{\text{num}} = 0$ up to round-off; for a dissipative scheme (upwind
differencing) $\sigma_{\text{num}} > 0$, and that spurious dissipation acts like an extra viscosity that damps precisely
the small-scale growth a blow-up produces. With physical viscosity present, a fair test is whether the **measured**
dissipation equals the **physical** one, $-\tfrac{dE}{dt} = \nu\langle u_x^2\rangle$, with nothing added.

The reference scheme is first-order upwind on the conservative flux $f = u^2/2$:

$$u_j^{n+1} = u_j^n - \frac{\Delta t}{\Delta x}\big(F_{j+1/2} - F_{j-1/2}\big) + \nu \Delta t \frac{u_{j+1}-2u_j+u_{j-1}}{\Delta x^2},
\qquad F_{j+1/2} = \begin{cases} f_j & u_j>0 \ f_{j+1} & u_j\le 0 \end{cases}$$

## Results so far

### 1-D Burgers, inviscid, `u0 = sin x` (`burgers_entropy.py`)
The gradient provably blows up at $t^*=1$ with $\max|u_x| = 1/(1-t)$. Grid $N=512$, dt = 2e-4.

| t | frozen-rotation scheme: E/E0, max\|u_x\|, L2 error | truth max\|u_x\| | first-order upwind: E/E0, L2 error |
|---|---|---|---|
| 0.50 | 1.00000, 2.00, 0.0000 | 2.00 | 0.9973, 0.0026 |
| 0.80 | 1.00000, 5.00, 0.0000 | 5.00 | 0.9952, 0.0066 |
| 0.90 | 1.00000, 9.99, 0.0000 | 10.00 | 0.9941, 0.0097 |
| 0.95 | 1.00000, 18.72, 0.0004 | 20.00 | 0.9934, 0.0116 |

Numerical entropy production sigma_num over $[0, 0.95]$: upwind 7e-3 per unit time; frozen -3e-15 (zero to machine precision). The frozen scheme tracks the exact approach to the singularity until the 512-point grid can no longer
resolve a gradient of 20 - a visible resolution limit, not hidden dissipation.

### 2-D viscous Taylor-Green, exact solution known (`taylor_green.py`, part A)
nu = 0.02, 128^2 grid, dt = 2e-3, viscosity as an exact integrating factor per mode.

| t | E/E_0 scheme | E/E_0 exact = exp(-4nu t) | L_2 error vs the exact field |
|---|---|---|---|
| 0.5 | 0.96078944 | 0.96078944 | 5e-15 |
| 1.0 | 0.92311635 | 0.92311635 | 1e-14 |
| 1.5 | 0.88692044 | 0.88692044 | 1.5e-14 |
| 2.0 | 0.85214379 | 0.85214379 | 2e-14 |

The dissipation is exactly the physical $\nu\langle|\nabla u|^2\rangle$ and nothing more. (Taylor-Green in 2-D is a
single decaying mode, so this is a precision test, not a cascade test.)

### 3-D inviscid Taylor-Green / Euler, the Brachet (1983) benchmark (`taylor_green.py`, part B)
Skew-symmetric nonlinearity, Leray projection, RK4, 2/3 dealiasing, $\Delta t = 2/N$.

| | t = 1 | t = 2 | t = 3 | t = 4 |
|---|---|---|---|---|
| E/E_0, all three grids | 1.000000 | 1.000000 | 1.000000 | 1.000000 |
| enstrophy Z, 32^3 | 0.417 | 0.574 | 0.914 | 1.570 |
| enstrophy Z, 48^3 | 0.417 | 0.574 | 0.931 | 1.724 |
| enstrophy Z, 64^3 | 0.417 | 0.574 | 0.939 | 1.823 |

Energy is conserved to six decimals in 3-D: zero numerical entropy production on the benchmark flow. Enstrophy
agrees across resolutions to three decimals through $t = 2$ and then fans out as the cascade reaches scales the
coarser grids cannot hold; the spread shrinks with resolution (convergence). Because nothing dissipates
numerically, under-resolution appears as **disagreement between grids** rather than as a smooth, plausible and wrong
curve - which is the point of the instrument. The trustworthy window at $64^3$ is roughly $t \le 3$; the literature
uses far higher resolution beyond that.

### Three more closed-form solutions, one per dimension (`exact_solutions.py`)
| exact solution | grid | error against the closed form |
|---|---|---|
| **1-D** viscous Burgers, u0 = sin x, nu = 0.02: Cole-Hopf, evaluated in heat-kernel form | 512 | 5e-14 (t=0.5), 6e-12 (t=1), 2e-6 through the shock (t=1.5-2, gradient 22 on dx = 0.012: the resolution floor) |
| **2-D** a *random* vorticity field on the single shell k^2 = 25 (nonlinear term vanishes identically): u0 exp(-25 nu t) | 128^2 | 2e-14 to 3e-14 out to t = 4 |
| **3-D** ABC flow with viscosity (Beltrami, curl u = u): u0 exp(-nu t) | 32^3 | 1e-15 to 7e-15 out to t = 8 |

![Exact solutions in 1-D, 2-D, 3-D and the scheme's error against each](../figures/exact_solutions.png)

*The 2-D case is a random exact solution, not a hand-picked mode; the 3-D case is the closed form the Beltrami
property provides. The 1-D shock comparison is the demanding one, and its error is the grid's, visibly.*

### 3-D Euler: the vortex-stretching term measured, by resolution (`budgets3d.py`, run on GitHub's runners)
$\dot Z = S$ with $S=\langle\omega\cdot(\omega\cdot\nabla)u\rangle$ measured from the field; energy $E/E_0 = 1.000000$ at
every grid and time; budget residual 1e-5-1e-6.

| t | 64^3: Z / S | 96^3: Z / S | 128^3: Z / S | S/Z^3/2 at 128^3 |
|---|---|---|---|---|
| 1 | 0.4169 / 0.0876 | 0.4169 / 0.0882 | 0.4169 / 0.0885 | 0.329 |
| 2 | 0.5743 / 0.238 | 0.5743 / 0.239 | 0.5743 / 0.239 | 0.550 |
| 3 | 0.9386 / 0.511 | 0.9430 / 0.526 | 0.9447 / 0.533 | 0.581 |
| 4 | 1.823 / 1.502 | 1.918 / 1.801 | 1.975 / 2.029 | 0.731 |

Converged to four digits through $t=2$. At $t=4$ the stretching still rises with resolution (1.50, 1.80, 2.03) with
shrinking differences: under-resolved grids **under**-estimate the amplification - they hide growth by inability, not
by damping, which is the opposite of a dissipative scheme's failure mode. The normalised rate $S/Z^{3/2}$ rises with
resolution as well (0.61, 0.68, 0.73 at $t=4$). This is a measurement of the term that decides the 3-D question, at
a resolution where its trend can be trusted, and nothing more.

### Mechanism, second invariant, and the BKM quantity (`mechanism3d.py`, $64^3$, from CI)
| flow | helicity < u.w> | enstrophy Z | max\|w\| (BKM integrand) | alignment with the intermediate strain eigenvector |
|---|---|---|---|---|
| Taylor-Green | 1e18 (zero by symmetry, preserved) | 0.375 → 1.823 | 2.0 → 13.7 | 0.356 at t=2 vs 0.317 / 0.327 (random 0.333) |
| perturbed ABC | **2.979345 → 2.979345** through t=3; 2.979187 at t=4 with energy 0.999979 | 1.58 → 11.3 | 3.1 → 47.5 | 0.339 vs 0.322 / 0.339 |

Helicity is conserved to six digits under dynamics that amplify enstrophy seven-fold; at $t=4$ the two invariants drift
together (within a factor 2.5), which is the grid being reached ($\max\lvert\omega\rvert = 47$ at $\Delta x = 0.1$).
The classical preference of vorticity for the intermediate strain eigenvector (Ashurst et al. 1987) is present but
weak here - a few points above random - because neither flow is developed turbulence by $t=4$; the forced run in
`spectra.py` is where the textbook signal belongs, and that check is pending.

### The regularity criteria along the classical candidates (`criteria3d.py`, CI)
**Reliability first.** The analyticity-strip clock (Sulem-Sulem-Frisch 1983) falls below 2 dx for full-box Kida-Pelz by
t = 1.0 at 64^3, 96^3 AND 128^3 (delta 0.124, 0.102, 0.079 against 0.196, 0.131, 0.098): the flow is trustworthy only to
t ~ 0.5-0.75 without symmetry reduction. Numbers previously shown here for Kida-Pelz at t = 3 were outside that window and
are withdrawn. Inside the window the picture converges and is physical:

| Kida-Pelz at t = 0.5 | Z | S | max\|w\| | CF coherence rho at h = 2pi/32 | Lipschitz exponent of the direction field |
|---|---|---|---|---|---|
| 64^3 | 5.4197 | 5.816 | 6.08 | 0.128 | 1.49 |
| 96^3 | 5.4197 | 5.817 | 6.08 | 0.127 | 1.51 |
| 128^3 | 5.4197 | 5.817 | 6.12 | 0.124 | 1.53 |

Four-digit agreement, and the geometric quantity is converged: the early roughening of the vorticity direction
(exponent 1.5 against 2 for a smooth field, from 2.7 at t = 0) is the flow's, not the grid's. The 64^3 run's later
"collapse" of the exponent to 0.5 was a grid artefact (96^3 and 128^3 give 1.15 at t = 1, itself past the clock). On
Taylor-Green and ABC (inside their windows) the critical L3 norm is flat to three digits, the direction field stays
coherent, and the high-vorticity set localises to under 1% of the volume. With nu = 1e-3, Kida-Pelz's enstrophy is
capped near 30 where the inviscid run reached 69, and the L3 norm falls 0.97 -> 0.83. Time integration is converged
(dt/2 changes the state by 3e-5). `THEORY.md` section 6 lists the criteria.

### Feedback, not imitation: searching initial data (`adversarial_ic.py`, CI)
Gradient ascent through the differentiable solver on the enstrophy amplification at fixed initial enstrophy, from smooth
|k| <= 4 data; every candidate re-verified at higher resolution with the analyticity clock.

| objective | search grid (32^3) | verified | classical best (reliable) | verdict |
|---|---|---|---|---|
| enstrophy, Kida-Pelz amplitude | 8.75x | 20.1x at 64^3, delta = 0 | 3.0x | outside the reliable window |
| enstrophy, Taylor-Green amplitude | 4.96x | 8.03x at 96^3, delta 0.06 < 0.13 | 1.11x | outside the reliable window |
| enstrophy with helicity = 0 (penalty) | 1.21x | 1.213x at 64^3, delta 0.28 > 0.20 | 1.112x | **pass** |
| enstrophy with helicity = 0.71 max (penalty) | 1.13x | 1.097x at 64^3 | 1.112x | fail (below classical) |
| enstrophy with helicity ~ 0.99 max (penalty) | 1.13x | 1.006x: no amplification | 1.112x | fail (below classical) |
| enstrophy, Taylor-Green amplitude, corrected verdict | 4.96x | 7.89x at 96^3 (delta 0.050), 8.70x at 128^3 (delta 0.040) | 1.112x | outside the window at both |
| enstrophy, leashed: delta(T) >= 0.30 on the 32^3 search grid | 4.44x, delta 0.27 | 8.28x at 128^3, delta 0.042 | 1.112x | outside the window |
| enstrophy, leashed: delta(T) >= 0.45 on the 32^3 search grid | 3.82x, delta 0.41 | 5.85x at 128^3, delta 0.047 | 1.112x | outside the window |
| enstrophy, leashed on a 64^3 search grid (checkpointed autograd) | 3.09x, delta 0.26 | 3.175x at 128^3 (delta 0.071 vs 0.098), 3.182x at 192^3 (delta 0.055 vs 0.065) | 1.112x | outside the window at both |
| helicity held at 0 (hard constraint, corrected bound) | 4.54x | 7.38x at 96^3, delta 0.047 | 1.112x | outside the window |
| helicity held at 0.25 max | 4.68x | 7.11x at 96^3, delta 0.065 | 1.112x | outside the window |
| helicity held at 0.50 max | 4.47x | 6.35x at 96^3, delta 0.069 | 1.112x | outside the window |
| helicity held at 0.75 max | 3.11x | 3.99x at 96^3, delta 0.076 | 1.112x | outside the window |
| helicity held at 0.90 max | 1.27x | 1.266x at 96^3, delta 0.20 > 0.13 | 1.112x | **pass** |
| Jacobi growth along Taylor-Green | 1.90x spreading | - | - | measurement only |

Pointed at enstrophy, the searcher twice found smooth low-k data that cascades to the grid cutoff within one time unit
while every classical flow stays smooth - the Lu-Doering / Ayala-Protas phenomenon - and twice the verifier rejected
it as unresolved, which is the design working. The leashed search (a penalty when delta(T) on the search grid drops
below a threshold) satisfied its leash on the 32^3 grid and still failed at 128^3, and that failure is the finding:
the coarse grid's truncation blocks the cascade the field triggers, so the field's spectrum tail *looks* resolved at
32^3 (delta 0.27-0.41) while its true continuation has delta 0.04. Resolution measured on the search grid is not
resolution. The leash has to be measured where the cascade lives. On a 64^3 search grid (gradient checkpointing, 2.2 GB)
the amplification is nearly converged between grids (3.09 at 64^3, 3.175 at 128^3) and the strip at 128^3 is 0.071
against a threshold of 0.098: a three-fold amplifier from smooth |k| <= 4 data that is almost, not yet, resolved.
At 192^3 the amplification is converged (3.182) but the strip is 0.055 against 0.065: delta keeps *falling*
with resolution (0.26 at 64^3, 0.071 at 128^3, 0.055 at 192^3) as each finer grid sees more of a flattening tail.
Registered FAIL. The separation is the finding: the enstrophy amplification of this field is a robust number, and its
small scales are not resolved at any grid tried. The field is committed (`results/found/leashed64_dmin030.npz`, 64^3
float32, 3 MB) for anyone with a bigger box. Topology, with the
helicity held fixed by a hard constraint (Newton projection onto the level set after every step, `HELMODE=project`):
on the 32^3 truncated system the maximal amplification over one time unit is *flat* up to relative helicity 0.5
(4.54, 4.68, 4.47) and then collapses (3.11 at 0.75, 1.27 at 0.9). Helicity inhibits the cascade, as Moffatt's
conjecture wants, but as a threshold rather than a slope: half the maximal helicity costs nothing, ninety percent
costs a factor of four.

![Maximal enstrophy amplification versus relative helicity held fixed](../figures/helicity_threshold.png)

The four fast fields cascade to the cutoff at 96^3 and are registered as outside the window;
the 0.9 field stays resolved and beats every classical flow (1.266 vs 1.112, pass). The earlier penalty-method rows
above, which showed a gentle monotone decline, were the penalty failing to explore, not the physics. (Correction, 2026-09-06: the first version of this table labelled these
0.5 and 0.7; the code's helicity bound was 2 sqrt(2 E Z) instead of the Cauchy-Schwarz 2 sqrt(E Z), a factor
sqrt 2. Fixed in `adversarial_ic.py`; a hard-constraint ladder at 0, 0.25, 0.5, 0.75, 0.9 is the replacement.) Arnold's
geodesic spreading along mild Taylor-Green is 1.9x over one time unit; it needs a ladder before it means more.

### Learned coarse simulators of 2-D Kolmogorov flow (`kolmogorov_v2.py`, run on CI)
$32^2$ coarse models against a $128^2$ truth, 2000-step rollouts from held-out states, $Re\approx1250$.

| coarse simulator | energy / truth | spectrum error (2nd half, averaged) |
|---|---|---|
| fully learned CNN (width 64, 16-step unroll) | 0.95 - **1.91** | 0.50 |
| frozen exact transport, no closure | 0.96 - 1.07 | 0.53 |
| frozen + dissipative gate only (nu_tge0; v1) | 0.95 - 1.07 | 0.53 |
| **frozen + energy-neutral learned redistribution + gate** | 0.87 - 1.07 | **0.09** |

In 2-D the missing physics at coarse resolution is backscatter, which a dissipative closure cannot express: the gate
alone changes nothing (0.53 = no closure). A learned term projected to move energy between scales *without creating
it* cuts the long-rollout spectrum error six-fold, while the fully learned model drifts to 1.9x the true energy.
Registered: spectrum ratio 0.18 (bar 0.6), energy band within [0.8, 1.2], learned leaves it - all pass. The dip to 0.87
is NOT the first-order leak of the skew projection: v3 (`kolmogorov_v3.py`, exact energy-neutral rescaling) reproduces
v2 to the digit, so the dip comes from the gate, which trained beside the skew term learns to dissipate more than
viscosity alone - and that is what best matches the truth's spectrum. A prediction of ours refuted by its own test.

### Energy spectra (`spectra.py`, CI) and exact time integration (`midpoint.py`)
![Energy spectra: 2-D k^-3 and 3-D k^-5/3 with a pile-up at the grid cutoff](../figures/spectra.png)

*2-D decaying turbulence at $256^2$ follows Kraichnan's $k^{-3}$ over a decade (a little steeper, as decaying 2-D
turbulence does). 3-D forced turbulence at $48^3$ follows $k^{-5/3}$ for a few wavenumbers and then **piles up at
the cutoff**: energy arrives faster than viscosity removes it on a grid this small. A dissipative scheme would have
damped that into a plausible slope; the conserving scheme shows the under-resolution. $96^3$ is queued.*

Implicit midpoint on the skew transport makes the discrete energy conservation exact: 3-D drift
3.5e-8 (RK4) -> 5.7e-12 (midpoint), round-off. Weak point 5 of `CLAIMS.md`, closed.

### Budget closure in 1-D and 2-D, where Hypothesis H is known (`budgets.py`)
Every balance law checked term by term along the flow (measured $d/dt$ across one step vs the right-hand side from
the field), including the **unsigned production terms** whose 3-D cousin is the regularity problem.

| case | budget | frozen scheme residual | upwind residual | what it shows |
|---|---|---|---|---|
| 1-D inviscid (H false) | gradient, -(1/2)<u_x^3> production | 1e8 | 1-5 % | production +0.13, +0.47, +1.59 at t = 0.3, 0.6, 0.8: the blow-up, measured |
| 1-D viscous (H true) | energy | 1e9 | **10 %** | upwind's numerical dissipation is a tenth of the physical viscosity |
| 1-D viscous | maximum principle max\|u\| | 0.994, 0.988, 0.984 | - | never increases: the controlled norm that gives 1-D regularity |
| 2-D decaying (H true, M=Z) | energy / enstrophy / palinstrophy | 1e8 / 1e7 / 3-9e-6 | - | palinstrophy production **+1387, +2149, +1579, +576** (unsigned, large) while enstrophy falls 12.86 -> 6.09 |

The 2-D row is Theorem 5 of `THEORY.md` with its terms filled in: the dangerous unsigned term is present and large,
and regularity holds because a norm one level below it is controlled. The registered residual bar of 1e-6 is
narrowly missed by the 2-D palinstrophy budget (terms of order $1e3$, RK4 error); halving the time step would
clear it, and the bar is left where it was.

## Liouville: the ensemble version of zero entropy production

Lee (1952) showed that the Galerkin-truncated Euler equations are a divergence-free vector field on the phase space
of retained modes: phase-space volume, and with it the Gibbs entropy of any ensemble of solutions, is exactly
conserved. That is the ensemble form of the statement this repository is built on. `liouville.py` measures it with
autograd, taking the exact trace of the Jacobian of the transport operator (one reverse pass per coordinate):

```
                          div F = tr(dF/dU)         expected
2-D  16^2  nu = 0         -2.9e-16                  0
3-D  12^3  nu = 0         +3.3e-16                  0
2-D  16^2  nu = 0.01      -24.20                    -nu (d-1) sum_k k^2 = -24.20
3-D  12^3  nu = 0.01      -82.32                    -nu (d-1) sum_k k^2 = -82.32
```

Two things worth noticing. The inviscid divergence is zero for the skew form, the advective form and the divergence
form alike: Liouville is a weaker property than energy conservation (the advective form loses energy but not phase
volume). And the viscous contraction is a constant independent of the state, so along truncated Navier-Stokes the
Gibbs entropy of any ensemble decreases linearly, S(t) = S(0) - nu (d-1) (sum_k k^2) t, however turbulent the flow.
What the flow does is not create or destroy phase volume but stretch it: the Jacobi ladder below measures the rate.

## The Riemannian view, measured: Jacobi fields along the Euler geodesic

Arnold (1966): an Euler flow is a geodesic on the group of volume-preserving maps with the L2 metric, and the
separation of two nearby geodesics (a Jacobi field) is governed by the sectional curvature; negative curvature
means exponential separation. `jacobi_ladder.py` integrates a flow and a copy perturbed by 1e-5 of a random
solenoidal field and reports the separation growth and its local exponent lambda, with the reliability clock on
every row (rows past the clock are omitted here).

```
Taylor-Green            growth |du|/(eps|v|)              local exponent lambda
    t         32^3      48^3      64^3              32^3     48^3     64^3
   0.5      1.0088    1.0128    1.0143              0.018    0.025    0.028
   1.0      1.0466    1.0553    1.0572              0.073    0.082    0.083
   1.5      1.1134    1.1301    1.1348              0.124    0.137    0.142
   2.0      1.2193    1.2488    1.2593              0.182    0.200    0.208
   2.5         -         -      1.4485                -        -      0.280

Kida-Pelz
   0.5         -      1.2301    1.2660                -      0.414    0.472
```

![Jacobi-field growth along Taylor-Green and Kida-Pelz by resolution](../figures/jacobi_ladder.png)

Three readings. The growth converges with resolution inside the window (Taylor-Green at t = 2: 48^3 and 64^3
agree to 1%). The local exponent rises steadily with time, 0.03 to 0.28 for Taylor-Green: the geodesic moves into
more negatively curved parts of the group as the vortex sheets form, which is Arnold's prediction in the one
quantity a computer can measure. And Kida-Pelz is four times more unstable already at t = 0.5, before its enstrophy
has grown by more than 30%: it is not yet converged at 64^3 (3% change from 48^3), so the number is provisional.
Energy in both copies stays at 1 to 1e-8 throughout; nu = 0.

Together with Liouville above this is the whole thermodynamic picture of the truncated inviscid system in two
numbers: phase volume is conserved exactly (div F = 0), and it is stretched at rate lambda in the least stable
direction and therefore compressed elsewhere. The flow does not lose information; it moves it to where a finite
observer cannot read it. That is also, word for word, what the analyticity-strip clock measures.

## Complex singularities: the analyticity strip as a Riemann-surface question

The clock used throughout this page is itself the classical singularity detector (Sulem, Sulem and Frisch 1983;
Frisch, Matsumoto and Bec 2003): E(k) ~ exp(-2 delta(t) k) where delta(t) is the distance of the nearest
complex-space singularity of the solution from the real domain, and a real finite-time singularity is delta(t*) = 0.
The two classical hypotheses are distinguishable *inside the reliable window* by the trend of the local decay rate
-d log(delta)/dt: constant for exponential decay (no real singularity), rising as 1/(t* - t) for linear decay.
`strip_tracker.py` samples delta every 0.025-0.05 time units and fits both laws over the window only.

![Analyticity-strip width delta(t) by resolution, Taylor-Green and Kida-Pelz](../figures/strip_decay.png)

```
flow            N^3   window          exponential tau   (2nd half)   linear t*   (2nd half)   decay rate start -> end   better fit
Taylor-Green     48   0.15 - 2.15       0.93            1.20          2.09        2.84           2.3 -> 0.51             exponential
Taylor-Green     64   0.15 - 2.60       1.04            1.69          2.41        3.65            -  -> 0.37             exponential
Taylor-Green     96   0.20 - 2.30       1.03            0.95          2.39        2.72            -  -> 0.77             exponential
Kida-Pelz        64   0.15 - 0.65       0.40            0.44          0.81        0.97           3.3 -> 2.0              exponential
Kida-Pelz        96   0.15 - 0.80       0.44            0.46          0.92        1.11           3.9 -> 1.7              exponential
Kida-Pelz       128   0.15 - 0.85       0.42            0.45          0.91        1.14           4.1 -> 1.9              exponential
```

Extended to 160^3 and 192^3 (`strip_hi`): tau = 0.41 and 0.375 (second half 0.50, 0.58), decay rate flat at 2.0,
linear t* retreating to 1.24 and 1.34. Five resolutions agree. The front-arrival test - the time at which the
cascade reaches the grid's reliability line 2dx = 4pi/N, which converges to t* for a real singularity (the dyadic
model below) and grows without bound otherwise:

```
   N^3     96      128     160     192       slope per ln N     prediction for 256^3
   t_arr   0.830   0.880   0.920   0.950     0.17, 0.18, 0.17   0.999
```

Logarithmic growth with a constant slope over four resolutions: no saturation. (The slope is below tau, so a single
exponential delta is too simple a model of the tail; the free-prefactor fit is unstable even with 32 shells and is
not used.)

Kida-Pelz is the case the literature argued about (Boratav and Pelz 1994 for a singularity; Hou and Li 2006, 2008
against). Here the exponential time constant is converged across 64/96/128^3 (0.40, 0.44, 0.42; second half
0.44, 0.46, 0.45), the local decay rate falls by half inside every window instead of rising, and the linear
extrapolation's t* retreats as the window lengthens - all three signatures of no real singularity in the reach of
these runs. Taylor-Green likewise (Brachet et al. 1983). The absolute delta shifts down with resolution because the
prefactor exponent n in k^-n exp(-2 delta k) is held at 0 (the clock's definition) and the fit range moves to
higher k; the decay *law* is what is compared, and it is resolution-independent. This is not a proof of anything:
it says that within t < 0.85 for Kida-Pelz and t < 2.6 for Taylor-Green, at up to 128^3, the nearest complex
singularity is moving away from the real axis at a slowing rate, exactly as the regularity side of the argument
predicts, and that a claim of a singularity from either flow would have to come from beyond where this instrument
can see.

## Searching for the missing functional, with an adversary in the loop

Hypothesis H needs a functional M(u) >= Z that never increases along Navier-Stokes trajectories. Theorem 4 says it
must be built from what Tao averaged away - local structure in physical space - so `lyapunov_search.py` looks for
one there: M = Z exp(Phi), Phi an enstrophy-weighted average of a small network over eight dimensionless pointwise
features of the vorticity and strain fields (stretching rate along the vorticity, strain magnitude and determinant,
direction roughness |grad xi|, local enstrophy and energy density, two alignment measures). No grid quantity enters,
so the trivial truncated-system bound is unavailable. dM/dt is the exact Lie derivative by autograd. After each
training round the differentiable solver searches initial data for the trajectory along which M grows fastest, and
those states join the training set. Registered verdict: PASS only if the final adversary and the held-out classical
flows show no violation above 1e-3.

```
24^3, nu = 2e-3, T = 0.6, Z0 = 0.375, three rounds (results/lyapunov_24_local.txt)
round   training worst   held-out worst (Kida-Pelz, random)   adversary's best violation (relative dM/dt)
  0        -0.020           -0.081                              +0.40
  1        -0.051           -0.283                              +2.40
  2        -0.039           -0.076                              +0.71
feature sensitivity of the learned Phi:  |w|^2/2Z 1.08,  xi.S.xi/sqrt Z 0.33,  |S|^2/Z 0.22,  |grad xi|/k_rms 0.09, ...
REGISTERED -> FAIL
```

Reading: the learner finds, every round, a local M that decreases along every classical flow and every previous
counterexample - and the adversary finds, every round, a new smooth low-k field along which it grows at 40-240% per
unit time. No trend toward closure in three rounds. The candidate leans on exactly the Constantin-Fefferman
quantities (local enstrophy density, stretching rate along the vorticity, strain magnitude), which is where a proof
would have to live; the search says that within this class - bounded, local, first-derivative features, one
time unit - the adversary always has a move. A society of six candidates with M = Z exp(min_i Phi_i) (multiple Lyapunov functions, Branicky 1998: each candidate
only has to decrease where it is the one in force, and a switch can only lower M) does no better - adversary +0.47,
+1.37, +0.66 - and is worse out of sample (Kida-Pelz now violates, +0.07). So the adversary's move is not a region no
candidate covers; from any M in this class it finds a smooth field along which M grows. The class is the problem: every
feature is local in the vorticity and strain at a point, and the quantity that decides how the stretching will
*change* is the pressure Hessian, a nonlocal singular integral of the whole field - the restricted-Euler dynamics
without it blows up in finite time (Vieillefosse 1982; Cantwell 1992), and it is precisely the specific singular
integral Tao's averaging destroys. Runs on CI at 32^3 for eight rounds then showed the one encouraging series of the week - adversary violations
falling from 1.0 to 0.03-0.09 with local features only - and the register was built for exactly this moment: the
same candidate attacked afresh by a stronger adversary (four restarts, forty iterations; `results/lyapunov_32_p0_attack.txt`)
is broken at +0.40, +0.35, +0.55. The fall was the fixed-budget searcher tiring, not the candidate hardening. With
pressure-Hessian features the eight-round series plateaued at 0.3 and the attack gives +0.74-0.79. With the signed-beta features the eight-round series reached 0.01, the best seen; attacked, it breaks at +0.42
(`results/lyapunov_32_b1_attack.txt`). A twenty-round candidate with a 25-iteration training adversary fell to 0.01-0.04 over its last five rounds
and, attacked with five restarts of sixty iterations, breaks at +1.03 (`results/lyapunov_32_r20_attack.txt`). Three
candidates, three falling series, three breaks. A candidate is only as good as the strongest adversary that has
failed against it; none has.

**Retraction (2026-09-07).** The verdict drawn from those failures - that the class of local functionals is closed in
3-D - is withdrawn. A positive control (`lyapunov2d.py`, `results/lyapunov2d_64*.txt`) ran the same machine in 2-D,
where a monotone quantity is a theorem (the enstrophy), with the answer available as a feature (log(Z/P), so that
Phi = log(Z/P) gives M = Z exactly). It failed the same way: attack violations +0.45 (viscous) and +0.41 (Euler),
leaning on the right feature and never reaching it. The 3-D failures were evidence about the *learner*, not about
Navier-Stokes. What survives: the searcher is a good refuter (every candidate it broke was broken), and no candidate
has passed - but "no candidate exists in this class" was never established. The learner needs a positive control it
can pass before its 3-D silence means anything.

**Control, second attempt (v2 learner: the global scalars enter linearly so the answer is exactly representable;
least-squares with a margin).** Viscous 2-D: passes the attack (no violation found, best -0.004) with learned
coefficients a = [1.001, -3.4] on [log(Z/P), log(E/Z)] - i.e. it recovered the enstrophy coefficient to three digits
and found M ~ Z^4.4 / E^3.4, which is a genuine monotone quantity under 2-D Navier-Stokes (from Z^2 <= E P) that the
author had not thought of; the registered wording "led by log(Z/P)" was too narrow and it prints "between", but in
substance it passes. Inviscid 2-D: fails (a = 1.114 where only a = 1 exactly is monotone; attack +0.34): the learner
finds monotone quantities where dissipation gives slack and cannot land on an exact conservation law. The 3-D problem
is the viscous one, so the 3-D run with this learner is the first whose verdict counts.

**3-D with the passed learner** (`results/lyapunov3d_v2_g1*`, `lyapunov3d_v3_*`). First run: the learner slipped into
energy - a = [-0.60 log Z, -1.03 log k_rms, -0.94 helicity] is sqrt(E) e^{-0.94 h} Z^{-0.12}, monotone because energy
is and useless because M < Z; attack +0.077. The energy loophole closed (a dominance term M >= Z, min Phi 2.2-3.2 on
held-out states): the constrained candidate leans on local enstrophy density, strain and the stretching rate - the
Constantin-Fefferman quantities - and its training adversary weakened to +0.03 and once came up empty; the strong
attack breaks it at +0.6 to +1.3. The worst-case (soft-max) loss does worse: +4.1. So, with a control behind it: in
the class of bounded local vorticity/strain functionals plus global scalars, over one turnover on low-k data, the
monotone set and the dominating set do not intersect in 3-D. The learner rediscovered Tao's wall from the inside -
handed everything, it found energy - and then found nothing beyond it. This is feedback,
not imitation, applied to the proof itself: the machine cannot produce a theorem, but it produces the counterexamples
a human would need to see before trying to.

## The price of the projection: what the global does for the local

The nonlinearity of Navier-Stokes, (u.grad)u, is purely local. Euler is that operator followed by the Leray
projection onto divergence-free fields, and the pressure *is* the projection: it is the only nonlocal thing in the
equation, decided by a Poisson equation over the whole box. Drop the projection and the equation is 3-D Burgers,
which blows up in finite time by a theorem (along characteristics the gradient obeys DA/Dt = -A^2, so
A(t) = A0 (I + A0 t)^-1 exactly and the first shock is at t* = -1/min lambda(A0)). `projection_price.py` runs both
from the same initial field with the same instrument (results at 48^3, all flows at the same initial enstrophy):

```
                            Euler (projected)                        Burgers (unprojected, local)
flow            exact t*    strip decay rate, t=0.2 -> 0.8           strip decay rate      clock stops at    max|grad u| vs exact A0(I+A0 t)^-1
Taylor-Green    1.000       2.71 -> 1.21  (falling)                  2.77 -> 3.59 (rising)  80% of t*        1.999 / 2.000, 2.490 / 2.500 until the grid fails
Kida-Pelz       1.277       1.95 -> 0.86  (falling)                  2.65 -> 1.68           71% of t*        1.001 / 1.003, 1.261 / 1.290
searcher's field 0.422      2.71 -> 1.40  (falling, Z x3)            4.34 -> 2.99           95% of t*        5.449 / 5.415
```

The unprojected run reproduces the exact characteristics solution to four digits until the grid fails, and its
clock stops at 71-95% of the exact shock time: a sixth closed-form check, in 3-D, on a blow-up. With the projection
on, from the same data, the strip's decay rate falls instead of rising in every case. Nothing else differs. Two
further readings. Energy: Euler conserves it *because* of the projection (the unprojected energy drifts by the
physical amount 1/2 int |u|^2 div u, up to 0.8% here); the conservation this repository verifies to 1e-14 is a
consequence of the nonlocal term, not independent of it. And the search: even the adversary's fastest field, whose
unprojected copy shocks at t = 0.42, is held by the projection to a three-fold enstrophy growth with a falling
decay rate.

How much of the pressure Hessian is global? Its trace is local (lap p = |w|^2/2 - |S|^2, decided at the point);
the traceless part is decided by the whole field. Restricted Euler keeps only the local part and blows up
(Vieillefosse 1982; Cantwell 1992). `pressure_share.py` measures the global share <|P_dev|^2>/<|P|^2> on the
high-vorticity set:

```
                    t = 0     0.25    0.5     0.75    1.0       Z(1)/Z0
Taylor-Green        0.39      0.41    0.49    0.58    0.67      1.11
Kida-Pelz           0.34      0.34    0.38    0.48    0.57      1.10
ABC (steady)        0.725     0.725   0.725   0.725   0.725     1.00
searcher's field    0.52      0.47    0.47    0.49    0.50      3.04
```

In the classical flows the global share of the pressure Hessian *rises* as the sheets form: the whole box
increasingly informs each point. In the adversary's field it stays flat at one half while the enstrophy triples.
That is a measurement of how the searcher wins: it finds data along which the global part of the pressure does not
grow with the stretching. Read with the Lyapunov result above, the two say the same thing from opposite sides: a
local functional cannot see the term that decides the future, and the fields that defeat it are the ones where
that term is weakest.

## The signed angle: what the squared coherence cannot see

Every alignment quantity in the literature and on this page so far is squared - the Constantin-Fefferman coherence
is 1 - (xi.xi')^2 - which folds antiparallel onto parallel: two vortex lines pointing opposite ways register as
perfectly coherent. That is right for CF's depletion argument (the Biot-Savart kernel depends on |sin|) and blind for
a diagnostic, because antiparallel neighbouring vorticity is the geometry of sheets and of the classic singularity
candidates (Kerr 1993; Hou-Li 2006). `pressure_share.py` now also reports the signed two-point angle
beta = 1 - xi(x).xi(x+h) in [0, 2] and the fraction of the strong-vorticity set whose neighbour two cells away is
antiparallel (beta > 1), at 64^3:

```
                     mean beta(h=2)  t=0 / 0.5 / 1.0     antiparallel fraction  t=0 / 0.5 / 1.0     Z(1)/Z0
Taylor-Green         0.004 / 0.006 / 0.012                0.000 / 0.000 / 0.000                       1.11
Kida-Pelz            0.009 / 0.014 / 0.039                0.000 / 0.000 / 0.000                       1.10
searcher's field     0.052 / 0.148 / 0.542                0.000 / 0.184 / 0.906                       3.12
```

By t = 1, ninety-one percent of the strongest vorticity in the searcher's field has an antiparallel neighbour: it is
built almost entirely of antiparallel pairs, while the classical flows have none at all - and the squared coherence
reported all three at the same 0.13. This is the sharpest single picture of how the adversary wins, and it names
the geometry: antiparallel sheets. Read with the concentration exponent (alpha ~ 4, flat) it says the pairs thin
without sharpening, which is what the projection does to them.

## Water's memory: reversibility and the compression half

Euler is time-reversible. `jacobi_ladder.py` with `REVERSE=1` runs a flow to T, flips the velocity, runs the same
instrument back to 0 and reports the round-trip error - the amount the flow has forgotten, which for an exact
inviscid flow is zero and for a truncated one is the information that passed below the cutoff - and the Jacobi
exponent of the reversed leg, which is minus the most CONTRACTING exponent of the forward flow: the compression half
that Liouville pairs with every stretch and that nobody measures.

```
                        round-trip error        forward stretching lambda_f (t=0.5-1.0)   reversed lambda_b (same window)   lambda_b - lambda_f
Taylor-Green   64^3     1.5e-9                  0.083                                     0.016                             -0.07
Kida-Pelz      96^3     2.0e-6                  1.27 (past the clock)                     0.57                              -0.70
searcher's     96^3     5.0e-5                  0.76 (past the clock)                     0.04                              -0.72
```

The memory is perfect until the cascade carries information below the resolvable scale: five orders of magnitude
separate the round-trip error of Taylor-Green from the searcher's field, in the order of their cascades - the
truncated analogue of what viscosity does. And in every flow the strongest compression is far weaker than the
strongest stretching: Liouville balances the books with one concentrated stretch and many diffuse squeezes. A
collapse to a point would need the compression concentrated as well (lambda_b >= lambda_f); no flow here shows it,
the adversary's field least of all. Provisional: the reversed leg starts from the least-resolved state, so most of
its rows are past the clock; the ordering is the claim, not the sizes.

## Water's memory as an algorithm: Navier-Stokes from Cauchy's memory and jitter

Cauchy (1815): vorticity is the initial vorticity remembered along particle paths. Constantin and Iyer (2008):
Navier-Stokes is the same memory read along Brownian-jittered paths and averaged - viscosity is the jitter, not a
separate term. `memory_paths.py` runs it in 2-D (Chorin's random vortex method): particles carry w0, move with the
flow plus sqrt(2 nu) dW, are deposited to the grid, and the velocity comes from the deposited vorticity by
Biot-Savart. No viscous term is ever applied. Against the spectral instrument, enstrophy at t = 2:

```
                              spectral NS     particles + jitter (2 / 4 / 8 per cell)     no-noise control
nu = 2e-3, 128^2              0.938           0.939 / 0.933 / 0.932                        0.991
nu = 2e-3, 256^2              0.956           0.956                                        0.999
nu = 1e-2, 128^2              0.766           0.770                                        0.991
nu = 0                        1.000           0.991                                        0.991
field error vs spectral, 128^2:  9.1% -> 4.5% -> 2.4% as particles per cell double (1/sqrt N, Monte Carlo)
```

The viscous decay rate is reproduced to 0.5-1% at two viscosities by jitter alone; the no-noise particles keep their
memory (the 0.99 is deposition smoothing at t = 0, identical at nu = 0). `wave_particle.py` adds the exact 1-D
version: the spectral (wave) solution evaluated at the particle positions agrees with the memory the particles carry
to 7e-15 while the wave is resolved, and the discrepancy then grows by seven orders of magnitude exactly as the
analyticity strip collapses toward the grid - the two descriptions are one flow, and the difference between them is
the observer's. Two more phase facts from the same script: the same amplitude spectrum shocks anywhere between
t = 0.16 and 0.38 depending on its phases (2000 random phase sets, factor 2.37); and in 2-D a developed turbulent field
with every phase randomised - identical spectrum, energy and enstrophy - has its palinstrophy production cut by a
factor of 250 (ratio 0.004), then regrows its phase correlations within a quarter time unit and cascades at the
original rate. The energy is in the amplitudes; the cascade is interference. And the interference is *between* scales
(`results/wave_particle_2d.txt`, T4): scramble the phases only above a cutoff k_c and only below it, and the two
halves together retain 9% (k_c = 8), 20% (16), 40% (32) of the cascade - the rest was in phase relations spanning
the cut. Keeping the large scales alone coherent leaves 0.8%; the small alone, 30%. The cascade is not in the
large scales or the small; it is the agreement between them, tightest around k = 8-16.

## Frequency matching: does incommensurability protect a flow?

Quasicrystals form because atoms hold two incommensurate lengths no periodic arrangement can satisfy at once. Does
data whose two scales cannot lock cascade more slowly or shock later? `frequency_matching.py`, registered prediction
"no": Conway's local-isomorphism property means every local configuration recurs, so the worst case is not avoided,
only relocated. 1-D Burgers, exact, equal energy: the harmonic pair (5, 10) has shock time 0.097 at its best phase and
0.067 at its worst (spread 1.46x); the incommensurate pair (5, 7 ~ sqrt 2) has 0.086 and 0.083 (spread 1.03x). Both
reach the fully aligned slope k + q somewhere; the incommensurate pair reaches it *whatever the phase*. Aperiodic
completeness is the reason it cannot be a shield: it guarantees the worst configuration is present. In 2-D the
question turned out not to be cleanly posable with two modes: a first design (5,10) versus (5,7) was confounded by
wavenumber content and is withdrawn; the fair design - two waves of identical |k| at 90, 53 and 37 degrees - gave
nine identical runs with no cascade at all, because w = k^2 psi makes any such field an exact steady state of 2-D
Euler (equal-magnitude waves do not interfere in 2-D). Recorded as a null with its reason.

In 3-D the memory carries a vector that the deformation stretches - Cauchy's formula proper, D w/Dt = (w.grad)u along
the path - and stretching is the whole 3-D problem. `memory_paths3d.py` runs it on Taylor-Green (results
`memory3d_*.txt`), enstrophy at t = 1 relative to each method's own start:

```
                      spectral        particles (Cauchy + jitter, no viscous term)     field error
32^3  nu = 0          1.1118          1.1029 (2/cell)  1.1023 (3/cell)                 3.4%
48^3  nu = 0          1.1130          1.1092                                           1.5%
32^3  nu = 2e-3       1.0960          1.0875                                           3.3%
48^3  nu = 2e-3       1.0971          1.0942                                           1.9%
```

The enstrophy growth by vortex stretching is reproduced to 0.3% at 48^3, and the viscous reduction of it (1.6% in the
spectral run) is reproduced by the Brownian kick alone (1.5%). The error is set by the trilinear read and deposit,
not by the particle count (3.40% to 3.25% from 2 to 3 per cell, 3.4% to 1.5% from 32^3 to 48^3): the pass criterion
"error falls with particles" is met only weakly, "falls with resolution" cleanly.

**The coherence race, and where it lands.** Scramble the phases above k_c and time the recovery of the cross-scale
cascade (`results/wave_particle_2d.txt`, T5): the time to rebuild half of the lost coherence is about 0.1 whether
the cut is at k = 4, 8, 16 or 32 - a speed limit, set by the large-scale strain rate. Coherence propagates through
scale at the strain rate. A singularity would need coherence to reach k = infinity in finite time, i.e. the strain
rate's time integral to diverge: that is the Beale-Kato-Majda criterion (1984), verbatim. The wave view, made
precise, is BKM in phase language - which also explains the 1-D control (Burgers' strain diverges, its coherence
time goes to zero, it shocks) and why the analyticity-strip clock works (its decay rate is the coherence
propagation rate). It gives three instruments and a mechanism for a criterion from 1984; it does not give a new
inequality, and this page does not claim one.

## Must the local move the global? The pressure held quiet

The searcher's fields kept the global (traceless) share of the pressure Hessian flat at 0.50 while their enstrophy
tripled, where the classical flows' share rises. Conjecture, registered before the run: the global part of the
pressure cannot stay quiet while stretching grows without bound. Test (`OBJ=quiet`, `results/quiet_pressure_*.txt`):
maximise enstrophy growth over one time unit with the global share at T capped.

```
cap on the global pressure share     growth reached (32^3 search)     share used     at 96^3
<= 0.5                               4.23x                            0.484          7.7x  (unresolved, as fast fields are)
<= 0.4                               2.39x                            0.409          3.1x  (unresolved)
<= 0.3                               1.09x                            0.434          1.13x (resolved; Taylor-Green level)
```

Robustness (`quiet_pressure2`, 48^3 search and 32^3 at T = 0.5 / 1.5), growth reached under each cap:

```
cap on global share    48^3 T=1    32^3 T=1    32^3 T=0.5    32^3 T=1.5     sign <xi.Pdev.xi>/<xi.S2.xi>
<= 0.5                 5.64x       4.23x       2.30x         7.38x          negative in all eleven runs (-0.7 to -3.5)
<= 0.45                4.18x        -           -             -
<= 0.4                 2.14x       2.39x       1.19x         4.63x
<= 0.35                1.11x        -           -             -
<= 0.3                 1.10x       1.09x       1.03x         1.19x
```

Same knee at 0.4-0.45 and same collapse by 0.35 at both grids and all three windows; the searcher could not push
the share below ~0.35-0.42 at all - capped, it gave up growth rather than quiet the pressure. And the sign: in every
fast field the global part of the pressure Hessian acts *with* the stretching along the vorticity. The sheet
recruits the box. Stated as a conjecture with a stable curve behind it: amplification beyond ~2x over a turnover
requires the traceless pressure Hessian to carry at least ~0.42 of the total on the high-vorticity set and to act with
the stretching; cap it at 0.35 and no amplification is reachable in this class of data.

Does the help fade as the sheet thins - the second half a proof along this line would need? (`help_fades`, the sheet
field sampled every 0.1 at 96^3 and 128^3, `results/help_fades_*.txt`.) Inside the reliable window it does not: the
global Hessian's push along the vorticity grows from -0.5 to -6.7 at t = 0.8 (the last resolved row at 128^3) while
the global share stays 0.45-0.52. A fade appears only past the clock and is weaker at 128^3 than at 96^3 at equal
time (-4.8 vs -3.7 at t = 1.6): a truncation signature, not physics. Kida-Pelz, the flow that does not grow, is the
one whose help fades (-18 to -1.1 as its share rises to 0.6). So the necessary condition C17 is met and stays met for
as long as the instrument can see; if a proof exists along this line, the bound on the box's help must come from the
limit, where no grid reaches - the wall, seen from the pressure side.

Forbid the pressure's global part from answering and the growth dies: half of it at 0.4, all of it at 0.3, where
the searcher gave up growing rather than quiet the pressure and still could not reach the cap. The "flat 0.50" was
not the pressure staying quiet; it was the pressure already carrying as much as fast growth requires. The local
cannot stretch without the global responding. A curve on a 32^3 search grid over one time unit: a conjecture with
evidence, not a theorem - but it is about the pressure, which is the term Tao's theorem says a proof must use.

## The quadratic map inside the equation

Along a particle path the velocity gradient obeys dA/dt = -A^2 - H, H the pressure Hessian: z -> z^2 + c in matrix
form. Restricted Euler keeps only the local part of H and escapes to a singularity (Vieillefosse 1982; Cantwell
1992); the full equation has the nonlocal part, the global pressure, as the c that keeps the orbit bounded.
`mandelbrot_fraction.py` iterates the local map on the high-vorticity set with H frozen at its local part and at
the full Hessian, over ten local turnovers, and counts escapes (`results/mandelbrot_fraction_48.txt`):

```
                        <Q>/|A|^2      local map escapes    with the frozen global H    rescued by the box
Kida-Pelz  (t=0..1)     +0.41 -> +0.20   100% -> 67%          54% -> 78%                  46% -> 17%
searcher's sheet        +0.12 -> +0.08    62% -> 46%          58% -> 38%                  19% -> 24%
```

The classical flows are rotation-dominated, deep inside, where the local map spins and the global c holds nearly
half of it. The sheet sits at Q ~ 0.08, on the edge where rotation and strain balance and c decides - and there,
at t = 0.25, more points escape with the global Hessian than without it (48% vs 36%): the box's pressure pushes the
sheet's orbits outward, the sign result of the quiet-pressure searches seen from inside the map. The adversary's
strategy, in Mandelbrot's language: go to the boundary of the set and recruit c. A snapshot diagnostic (H frozen,
horizon a choice), not a theorem.

**The sheet under the map, measured** (`nilpotent_sheet.py`, 48^3, `results/nilpotent_sheet_48.txt`). Registered:
on a sheet the gradient is nilpotent (a pure shear has A^2 = 0), so the local term of dA/dt = -A^2 - H is silent.
Holds: on the high-vorticity set the nilpotency ratio |A^2|/|A|^2 falls from 0.44 to 0.17 as the searcher's sheets
form, against 0.46-0.47 for Kida-Pelz and Taylor-Green and exactly 0 for the pure-shear control. Fails as first
stated: the pressure does not then dominate in norm - |H_dev|/|A^2| is 0.67 on the sheet, lower than the classical
flows' 0.75-0.89 - because a flat sheet has no pressure either: A^2 = 0 and H_dev = 0, a plane shear is steady.
Everything that happens to a real sheet is external and weak: the large-scale strain thinning it (the residual 0.17)
and curvature-induced pressure rolling it, against viscosity thickening it - which is why it thins slowly and
exponentially rather than collapsing. Caveat added to the "help" sign above: on rotation cores xi.H_dev.xi < 0 is
forced by the trace subtraction (lap p = 2Q > 0 at a pressure minimum), so the negative sign on the classical flows
is partly local; on the sheet, where Q ~ 0, the sign is genuinely nonlocal and small (-0.04).

**The strain budget on the sheet, closed** (`results/nilpotent_sheet_48b.txt`). |w| grows at exactly the rate of
stretching along the vortex lines, xi.S.xi, and incompressibility splits that into compression across the sheet
(thinning, -n.S.n) plus compression along it (narrowing, -t.S.t):

```
sheet, t     stretch xi.S.xi     thinning (-n.S.n)     narrowing (-t.S.t)
0.25         1.07                0.57                  0.33
0.50         1.05                0.83                  0.05
0.75         0.81                0.30                  0.44
1.00         0.60                0.16                  0.38
```

Two phases. To t = 0.5 the sheet is squeezed thinner by the surrounding strain; from t = 0.75 the across-sheet
compression fades and in-plane compression takes over - the sheet folds toward a tube, which is where the
antiparallel fraction jumped from 0 to 0.9. The pressure's role is the second phase: the curvature pressure that
rolls a shear layer. Neither phase is the local map (nilpotency 0.17 throughout). The object a proof would have to
bound is the second phase - the roll-up of a thinning sheet under its own curvature pressure against viscosity:
Birkhoff-Rott with thickness, in three dimensions. (The triad xi, n, t is not exactly orthonormal, so the budget
closes to ~85%.)

**The thinnest squeeze, and a candidate inequality refuted** (`thinnest_squeeze.py`, `results/thinnest_squeeze_48.txt`).
Every collapse is a material element pressed to zero thickness; along a path the deformation gradient F obeys
dF/dt = A F with det F = 1 (measured 1.0000 on every path), and sigma_min(F) is the thinnest that parcel has been
pressed. The natural guess in the memory picture - nothing can be squeezed thinner than the jitter sqrt(2 nu t) smears
it, which would give the Type I bound - is false: along the searcher's sheet under Navier-Stokes (nu = 2e-3) sigma_min
falls by half every quarter time unit and crosses below the jitter scale between t = 0.75 and 1.0 (ratio 1.5 -> 0.68),
in a flow that is perfectly resolved and regular. The reason: material lines in any smooth flow thin exponentially
without bound (Batchelor's regime for a dye filament) while the velocity stays smooth; the label map's features below
the diffusive scale are averaged away by the jitter and say nothing about grad u. What does hold: 1/sigma_min tracks
max|w| (23 vs 21 on the sheet, 1.6 vs 1.8 on Kida-Pelz) - sheet thickness is 1/|w| - so the Lagrangian side returns
the same exponential thinning law as the Eulerian side.

## The constant that is really a function

Navier-Stokes assumes friction -nu Laplacian u with nu a constant. Make it a function and regularity is a theorem:
Lions (1969) for (-Laplacian)^alpha with alpha >= 5/4 (Tao 2009: 5/4 minus a log); Ladyzhenskaya (1967) for a
viscosity rising with the strain rate. Real water does both (kinetic regime below the molecular relaxation scale;
temperature dependence). The constant-nu Laplacian is the one open case. `friction_function.py` runs the searcher's
antiparallel-sheet field under each (48^3, nu = 2e-3, t = 1):

```
friction                                       Z(1)/Z0    max|w|(1)    strip delta(1)
none (Euler)                                   3.04       23.6         0.149
constant nu, Laplacian (Navier-Stokes)         2.47       20.7         0.177
hyperviscosity alpha = 1.10                    2.47       20.7         0.179
hyperviscosity alpha = 1.25 (Lions, regular)   2.46       20.6         0.182
Ladyzhenskaya nu(|grad u|) (regular)           2.03       18.0         0.213
```

Friction as a function of the local strain bites at once - it acts on the sheet itself and cuts its growth by a
third. Friction as a function of scale is invisible here: alpha = 1.25 is provably regular and alpha = 1 is the open
case, and at this resolution and window they cannot be told apart, because the extra damping lives at scales the
sheet has not reached. That is the shape of the wall from the side the theorems name: the tie between stretching
and constant-nu friction is a statement about the limit, which no finite window shows.

## The blow-up type, in the frame the theorems use

The combined CKN / Seregin-Sverak / KNSS program zooms into a would-be singularity and asks for its type: max|w| ~
(t* - t)^-gamma with gamma >= 1 (the BKM gate; Type I is gamma = 1), faster (Type II), or exponential (no
singularity). `strip_tracker.py` now fits both inside the reliable window. Kida-Pelz at 96^3: a power law with
gamma = 0.21, far below the gate. The searcher's field at 96^3 and 128^3: exponential growth of max|w| at a rate of
about 1.6 per time unit fits at least as well as any power law, and the local growth rate falls toward the end of
the 128^3 window (1.74 -> 0.91). Exponential thinning is what a strained sheet does; it is not a singularity type.

## Keep the state, track the dissipation: the CKN concentration exponent

Navier-Stokes carries its state by an exact transport and leaks energy only through nu |grad u|^2, whose time
integral is bounded a priori by E(0). Caffarelli, Kohn and Nirenberg (1982) turned that into the strongest partial
result there is: a point is regular unless the dissipation in the parabolic cylinder of radius r around it is at
least eps r, so the singular set has parabolic dimension at most one. `ckn_exponent.py` measures how the dissipation
actually concentrates around the strongest vortex along a flow, D(r) ~ r^alpha over a ladder of r (uniform smooth
dissipation: alpha = 5; a sheet: 4; a tube: 3; CKN-critical: 1), reporting a radius only once the cylinder fits
inside the elapsed time. 64^3, nu = 1e-3, same initial enstrophy:

```
                     fit alpha (r = 2-6 cells)     smallest resolved pair (3 -> 4 cells)    over t = 0.4 .. 1.0
Taylor-Green         5.4                           4.4                                      flat, Z x 1.1
searcher's field     4.5 - 4.6                     3.6 - 3.8                                flat, Z x 2.8, max|w| x 3
```

The adversary's field concentrates its dissipation like sheets, a full exponent below Taylor-Green at the smallest
resolved scale, and the exponent does not fall while the enstrophy nearly triples: three and a half exponents from
the critical value and not moving toward it. This is the sixth view of the same flow (strip, direction coherence,
projection, pressure share, helicity, and now dissipation concentration), and they say one thing between them: the
searcher makes enstrophy grow by building sheets in the helicity-free regime where the global part of the pressure
stays weak, and the projection keeps the sheets from sharpening into anything a singularity needs. The exponent is
also the first search objective that points at the theorem itself: minimise alpha rather than maximise enstrophy.
Done (`OBJ=ckn`, 64^3 search, 128^3 verification, `results/sheet_conjecture_64.txt`): rewarded for concentration the
searcher builds something sharper than a sheet - spatial exponent alpha_s = 1.25 on the search grid, 1.36 at 128^3
(a tube; uniform is 3, a sheet 2) - with only 1.8x enstrophy growth, and it is unresolved at 128^3 (delta 0.078 vs
0.098), registered FAIL. So growth buys sheets and concentration buys tubes, and neither is resolved. The sheet
conjecture stands for the enstrophy objective; the stronger claim, that nothing sharper can be built, is refuted.

## The equation, and the two solved cases as controls

Everything this page has tried to find reduces to one requirement: a quantity M[u] >= 0 that is non-increasing along
the flow and dominates the gradient, M >= c |grad u|. Two sets of candidates - the MONOTONE ones and the DOMINATING
ones - and regularity is the statement that they intersect. `monotone_vs_dominating.py` measures both properties for
a family of norms along trajectories where the answer is a theorem (`results/monotone_vs_dominating.txt`):

```
1-D Burgers (shocks), gradient x12:   int|u|, max|u|, total variation: monotone, ratio to the gradient falls 10x
                                      int u_x^2, max|u_x|: dominate, grow 10-18x        -> the sets are DISJOINT
2-D Navier-Stokes (regular):          enstrophy: non-increasing (-0.009) AND ratio to the gradient within 30%
                                                                                          -> the sets INTERSECT
```

(2-D caveat: the gradient grew only 27% in this window, so "dominates" is weakly tested there; enstrophy is the real
member of both sets, energy passes only because little happened.) In 3-D the question is whether the sets intersect;
this repository's searches - hand-picked norms, the learner, the learner with a passed 2-D control - have not found a
member of both in any family tried. The searcher's antiparallel-sheet field is the counterexample to every
topological candidate (helicity, linking, crossing number all vanish on it while it is pressed without bound), and
the pressure's cubic flux is what breaks every energy-like candidate. If a 3-D member exists it is neither.

## The quantum seam: the same collapse in a fluid where the cut is guaranteed to win

Every physical fluid stops being the idealised equation before a singularity could form - mean free path, thermal
noise, and for a superfluid the quantum pressure. The quantum case is the clean one, because it is *solved*. The
Gross-Pitaevskii equation i psi_t = -1/2 lap psi + (|psi|^2 - 1) psi becomes, under psi = sqrt(rho) e^{iS} (Madelung
1927), Euler for (rho, v = grad S) plus one extra term, the quantum pressure (1/2) grad(lap sqrt(rho)/sqrt(rho)): a
stiffness that forbids density collapse below the healing length xi. The equation is globally regular (defocusing NLS,
3-D cubic is H^1-subcritical), circulation is quantised (2 pi), a vortex line is a smooth zero of psi, and
reconnection - our seam - is two zero-lines crossing, smoothly, with a burst of sound (Koplik-Levine 1993; Kerr 2011;
Villois-Proment-Krstulovic 2017: the gap closes as (t_r - t)^{1/2}).

`gpe_seam.py` is a split-step spectral GP solver (exact in each half-step, energy to 1e-6) with the Nore-Abid-Brachet
energy decomposition (incompressible kinetic / sound / quantum / interaction) and vortex lines located as phase
windings, sub-cell. Registered before the CI run (CLAIMS C18): 2-D, two counter-propagating dipoles at separation
d = 6 travel within 20% of the point-vortex speed 1/d (local check: +6%, energy drift 1e-6); 3-D, two antiparallel
pairs bowed toward each other reconnect at a finite t_r, the gap law before the cut has exponent 0.5 +/- 0.15, the
sound energy rises across the cut while the incompressible energy falls, and the lines never re-approach. Result
(`results/gpe_seam_3d_64.txt`, `_96.txt`; the two resolutions agree on every number below):

```
                              64^3 (dx 0.50)      96^3 (dx 0.33)
reconnection t_r              15.8                15.5           (closest gap 3.5 xi at t = 0; all slices cut until then)
gap ~ (t_r - t)^p, gap 1.5-3.6   p = 0.39            0.37        (includes the Crow-instability phase, not a power law)
                  gap 1.0-2.5   p = 0.65            0.55         (the asymptotic window: consistent with 1/2)
sound E_kc across t_r +/- 3   +32%                +36%
sound E_kc, cut -> peak       15.5 -> 37.5 (t 36)  15.5 -> 36.6 (t 36)   2.4x, released over ~20 time units, not at the instant
E_ki across t_r +/- 3         +3%                 +3%            <- registered "falls": it does not, at the cut
E_ki, t_r -> t = 36           600 -> 576 (-4%)    601 -> 577     it falls over the release window
slices still cut by all 4     0.89 -> 0.58        0.89 -> 0.59   (t = 18 -> 60): the rings retract steadily
re-approach after the cut     none                none
```

C18 stands: the pair reconnects at a finite time in both resolutions, the gap law in the asymptotic window is the
literature's 1/2, sound rises across the cut and never stops rising for twenty time units, the lines never come back.
The one registered element that failed as written: the incompressible energy does not drop *at* the cut (+3%); it
drops over the release window that follows. The cut itself is cheap; the price is paid afterwards, by the retracting
cusps and the Kelvin waves radiating - which is the template: in the quantum fluid the seam is cut on contact, and
what would have been a collapse leaves as sound.

**What this says we are missing - exactly.** Not a term that produces energy (the ledger above is exact) but a
*regulator whose length does not move with the flow.* The quantum pressure acts at fourth order (lap^2 on the
density) with a fixed coefficient; viscosity acts at second order and its scale sqrt(nu/s) is set by the strain
the flow itself generates, so the roll-up can lower its own floor. The gap is quantifiable: Navier-Stokes with
(-lap)^alpha in place of -lap is *proven* regular for alpha >= 5/4 (Lions 1969; Tao 2009 shaves a logarithm off
5/4), because at 5/4 the energy budget becomes scale-invariant and the free collapse of the section above stops
being free - it costs a fixed, finite amount, and the budget forbids it. The quantum fluid regularises at alpha = 2,
comfortably past the line; real viscosity sits at alpha = 1, a quarter of a Laplacian short. Everything this
repository has measured on the sheet - phase 1 external squeeze stopped at sqrt(nu/s), phase 2 self-generated
roll-up that the budget cannot price - is what that quarter-power looks like from inside one flow.

## Growth needs the twist? The Möbius conjecture, adversarially

The literature has two theorems about the direction of vorticity and no theorem about the gap between them:
Constantin-Fefferman (1993) - if the direction is Lipschitz-coherent where vorticity is intense, no blow-up; and
Grujic (2025, on Moffatt-Kimura antiparallel tubes) - where the direction *oscillates*, a cancellation in the
stretching term. Partial coherence is covered by neither. The searcher's fastest fields sit at 90% antiparallel on
the high-vorticity set (`seam_race.py`), while every classical flow sits near 0. So: is the twist *necessary* for
fast growth? `adversarial_ic.py` now takes `TWISTW`, a differentiable enstrophy-weighted antiparallel penalty at T
(1- and 2-cell neighbours), with the hard anti fraction printed in the search and in the verification. Registered as
C19 before the run: at matched leash and iterations, forbidding the twist halves the attainable growth at both 32^3
and 64^3. Refuted by a fast, coherent amplifier - which would bypass both theorems and would be new.
Results in `results/twist*.txt` when the CI job lands. First rows (32^3 search, 128^3 verification): the three
classical flows have anti = 0.000 at 128^3 with growth 1.01-1.11; at weight 3 the searcher removes the twist, stalls
near 1.3, buys it back and ends at 3.15 (4.32 at 128^3, anti 0.32, unresolved); at weight 10 the *strong* twist (the
enstrophy-weighted measure) is held at 0.0035 - a tenth of the weight-3 value - and growth ends at 1.49 against the
8.28 baseline, while the plain count of antiparallel neighbours creeps to 0.40 by placing faint reversed vorticity
next to strong structures. Re-scoping of C19 recorded there: the count was the wrong instrument; growth tracks the
strong reversal. Denied the twist, the searcher first went helical (H/Hmax to -0.32, the regime where stretching
dies) and abandoned it as soon as any weak twist was available.

**Is the strong twist the monotone quantity?** If growth needs it, it dominates; the other half is whether
Navier-Stokes makes it non-increasing. Quick local test (`seam_race.py`, strong-twist column added, 32^3, T = 1.2,
the clock at ~1.1):

```
                          strong twist 0 -> 1.2     Z/Z0 at 1.2     dlogZ/dt at t = 1
nu = 0                    0.026 -> 0.298            1.25            +0.27
nu = 1e-3                 0.026 -> 0.267            1.17            +0.19
nu = 2e-3                 0.026 -> 0.239            1.09            +0.12
nu = 5e-3                 0.026 -> 0.174            0.90            -0.07
```

It rises at every viscosity, every step, including nu = 5e-3 where the enstrophy falls throughout: the sheets are
pressed together by the external strain, which viscosity does not touch at these thicknesses. Refuted as the
monotone quantity within the window; it sits in the dominating set, like every candidate before it. What would turn
it over is reconnection at the viscous seam thickness, which this window does not reach - the seam race again.

**The 64^3 run refutes C19 as registered (2026-09-08).** With the twist forbidden at 64^3 (weight 10, 45 iterations,
leash 0.30) the searcher reached 2.22x on the search grid and **2.24x at 128^3, resolved** (delta 0.113 > 0.098,
anti 0.000, alpha_s 2.24 - a coherent sheet), against the twisted baseline's 3.18x (which is *not* resolved at 128^3,
though its value is stable to 192^3). That is a 30% cut, not the registered 50%. At the same measurement scale the
baseline's high set is 91% antiparallel at t = 1 and the coherent field's is 9%. So at 64^3 the twist buys about
40% more growth; it is not necessary for growth. The 32^3 numbers (8.28 vs 1.2-1.5) overstated the effect because
the 32^3 baseline's 8.28 is unresolved growth. Full table, all at 128^3 verification:

```
                         search   growth   at 128^3   resolved   anti at 128^3
baseline (twisted)       32^3     8.28     (8.28)     no         -
w = 3                    32^3     3.15     4.32       no         0.32
w = 10                   32^3     1.49     1.52       YES        0.00
w = 30                   32^3     1.19     1.19       YES        0.00
baseline (twisted)       64^3     3.18     3.18       no (stable to 192^3)
w = 10                   64^3     2.22     2.24       YES        0.00
classical                -        -        1.01-1.11  YES        0.00
```

Two things survive. The trade-off is real and monotone at every resolution tried - forbidding the twist always
costs growth, and denied it the searcher walks to the helicity plateau (H/Hmax = -0.50 at w = 30, the ladder's
edge). And the twist penalty is the first objective that produced *resolved* found fields: 1.19, 1.52, 2.24 - the
coherent regime's growth is now measured rather than extrapolated. A defect in the instrument, recorded: the twist
is measured at one and two grid cells, so it carries the grid scale (0.39 at 32^3, 0.20 at 64^3, 0.05 at 256^3);
the 256^3 seam run below reads the same 64^3 field as twist 0.000 at t = 0 for exactly this reason. It has to be
defined at a fixed physical separation before the curve is compared across resolutions.

## The seam race at the viscous thickness: 256^3 on a GPU, resolved end to end

`kaggle/seam/seam_gpu.py` (PyTorch, same scheme, T4). The searcher's 64^3 sheet field at 256^3, T = 3, clock
2dx = 0.049. Race variable = sheet thickness |w|/|grad|w|| on the high set over the viscous thickness sqrt(nu/s).
Registered as C20 before the run. `results/seam_gpu/`.

```
nu = 2e-3, resolved throughout (delta >= 0.0585 > 0.049 to t = 3):
  t      Z/Z0    max|w|   twist(0.05)  anti    ell     ell_nu   race   cut(visc > +stretch on the seam)
 0.0    1.000     7.2     0.00000      0.000   0.636   0.261    2.4    0.00
 0.9    2.306    22.2     0.00001      0.000   0.296   0.092    3.2    0.00
 1.2    2.825    22.6     0.00051      0.013   0.260   0.065    4.0    0.82
 1.5    3.163    25.1     0.00762      0.081   0.205   0.051    4.0    0.63
 1.8    3.161    41.1     0.01329      0.083   0.223   0.086    2.6    0.81     <- Z peak, max|w| peak (5.7x)
 2.1    2.941    37.1     0.00558      0.049   0.188   0.062    3.1    0.53
 3.0    2.421    32.5     0.00376      0.016   0.156   0.104    1.5    0.65
twist peaks 0.0147 at t = 1.7 and falls 74%; the seam is cut (viscous cancellation beats positive stretching on
50-80% of it) from the moment it forms; enstrophy and max|w| peak and decay inside the clock.

nu = 0:      the clock expires at t = 0.9 (delta 0.2255 -> 0.050, e-folding ~0.6) with max|w| 3.7x and no twist
             at the 0.05 scale yet; everything after (Z "167x" at t = 3) is unresolved and is not a number.
nu = 1e-2:   Z peaks at 1.31 and decays; viscous from the start.
pair IC:     setup failed - the tubes (D 1.2, sigma 0.35) did not interact in T = 3 (max|w| 6.6 -> 7.5) and the
             compression across them is ~0, so the race variable is undefined; the printed PASS is void.
```

Registered clauses at nu = 2e-3: twist falls >= 20% - yes (74%); race variable <= 1.5 inside the clock - 1.49 at
the last row, barely; max|w| growth < 3x - **no**, 4.55x at the clock (peak 5.7x). Formally "between"; the physics
is not ambiguous: at this Reynolds number (~500, low-k data) **the cut wins**, and for the first time the whole
race is inside the reliable window. The Euler seam, by contrast, is out of reach of any spectral box: the
analyticity strip collapses exponentially and 256^3 buys t = 0.9. Two honest limits: the Reynolds number is low
and the data is the searcher's, not Kerr's; the tube-pair control has to be rebuilt (closer, longer) before it
says anything.

### The twist is a wave in scale, and the Re ladder ran out of clock (2026-09-08, `results/seam_gpu/v2/`)

Twist re-measured at fixed physical separations (0.05, 0.1, 0.2, 0.4) on the same 256^3 runs, plus nu = 1e-3 and
5e-4, plus a rebuilt tube pair. The new instrument shows something the grid-scale one could not:

```
nu = 2e-3 (resolved to t = 3):   twist@0.4 peaks 0.48 at t = 0.9
                                 twist@0.2 peaks 0.30 at t = 1.5
                                 twist@0.1 peaks 0.11 at t = 1.6   (falls 62% by t = 3)
                                 twist@0.05 peaks 0.013 at t = 1.8  (falls 72%)
```

The reversal migrates down in scale as the sheets are pressed together - a wave in scale, each separation's twist
rising and falling as the seam passes through it - and at nu = 2e-3 the wave is arrested at ~0.05-0.09, the viscous
thickness sqrt(nu/s), where the cut takes over (cut fraction 0.35-0.47). That is the seam race seen whole.

C21 (Re_seam at the twist peak the same across nu) is **inconclusive**: the resolved peak exists only at nu = 2e-3
(Re_seam = 224 at t = 1.6; my "order 10" was wrong by 20x). At nu = 1e-3 the clock expires at t = 1.6 and at
nu = 5e-4 at t = 1.1, both *before* the twist@0.1 turns over - the numbers there (497, 446) are values at the clock,
not peaks. Halving nu needs 1.68x the grid; 256^3 resolves the peak at 2e-3 only. Everything past the clock in
those runs (Z 11x, max|w| 115 at 5e-4) is unresolved and not a number. C20's max|w| clause at the clock: 4.55x,
6.3x, 3.8x - all above the registered 3x.

**A registered KILL fired, on the pair control**, and I think it is the criterion, not the physics - stated so it can
be checked: the rebuilt pair (D 0.7, sigma 0.22, T = 6, resolved throughout) has twist@0.1 still rising at t = 6
(0.014) with the race variable at 0.56, which is the KILL clause as written. But Z is flat (1.06), max|w| grew
2.1x over six time units, Re_seam is *falling* (83 -> 53), and the race variable is small because it is measuring
the diffusing Gaussian *core* (ell ~ sqrt(nu t)), not the gap between the tubes; twist@0.4 is falling while @0.1
rises - the pair is slowly approaching, not collapsing. The clause conflated core thickness with seam gap. Recorded
as a KILL that fired; the corrected clause for the next run: twist rising AND max|w| accelerating AND Re_seam
rising. Re_seam's definition also switches sets when the twisted set is tiny (the 47 at nu = 1e-3, t = 0.6);
it needs a minimum count.

### The descent law and the V (2026-09-08, `results/seam_gpu/v3/`)

Eight separations in half-octave steps, output every 0.05, three viscosities. Registered as C22 before the run.

```
peak time of twist@sep        nu = 2e-3   nu = 1e-3   nu = 5e-4
  0.56                          0.70        0.70        0.70
  0.40                          0.95        0.95        0.90
  0.28                          1.25        1.25        (clock)
  0.20                          1.40        1.35        (clock)
  0.14                          1.50        1.50        (clock)
time per octave, 0.56 -> 0.28   0.52, 0.58  0.52, 0.58
time per octave, 0.28 -> 0.14   0.31, 0.19  0.21, 0.29
```

**C22 passes both clauses.** Every resolved peak agrees across viscosities within one sampling step: the descent of
the reversal through the scales is inviscid, measured resolved at finite nu. And it accelerates - about 0.55 per
octave at the top of the ladder, about 0.25 at the bottom, identically at both rungs. Successive half-octave times
shrink by a factor 0.76 on average against 1/sqrt(2) = 0.71 for a rate proportional to 1/ell: the wave is rhythmic
in the collapse's own clock and accelerating in ours. Pooling the resolved peaks, the gap between the sheets follows
gap = 0.51 (1.78 - t): a straight line to zero at T* = 1.78 - if nothing intervened.

**The V.** Two arms descend at different speeds: the gap (inviscid, linear, above) and the sheets' own thickness,
the analyticity strip delta(t) (viscous, exponential, slower). Where they meet the pair merges into one doubled
structure, before the gap reaches zero. Registered before reading the lower rungs: the meeting time is nearly
nu-independent (1.65-1.75), its scale falls with nu, and at lower nu the sheets meet above the viscous thickness.

```
nu        V bottom t    V scale    ell_nu there    clock    reading
2e-3      1.65          0.068      0.053           2.40     merge at the viscous scale: the merge is the cut
1e-3      1.70          0.043      0.025           1.60     merge ABOVE the viscous scale: an Euler-like merge
5e-4      1.74          0.022      0.023           1.10     both arms extrapolated past the clock; not counted
```

At nu = 2e-3 the V bottom coincides with everything the run itself flagged - the twist@0.07 peak (1.65), the Z peak
(1.65), sqrt(nu/s) = 0.053 - three lengths meeting at one point. At nu = 1e-3 the sheets meet while viscosity still
cannot cut them: that is where the Euler seam race starts, located from resolved data at a place the clock cannot
reach. Caveat recorded: the thickness arm is not nu-independent (e-folding 3.4, 1.6, 0.8) and part of that is the
fit window; the gap arm is the clean one. C21 (Re_seam at the twist peak across nu) remains inconclusive: only the
2e-3 peak at separation 0.1 is resolved.

### Generality (C23), read from v4 (`results/seam_gpu/v4/`, 2026-09-08)

```
                              nu = 2e-3                        nu = 1e-3
Kida-Pelz (coherent)          no descending wave; twist at      same; max|w| 0.85x
                              all separations rises together
                              (0.56: 0.33, 0.4: 0.28 at t=2.4,
                              below 0.005 under 0.14); max|w|
                              0.90x, Z 1.37
second adversarial field      peaks 0.28@0.95 0.20@1.00 0.14@1.00   0.28@0.90 0.20@0.95 0.14@1.00
(ckn64, tube-builder)         per octave ~0.9 then ~0.1            0.82, 0.97, then 0.10, 0.10
                              wave stops at 0.10-0.14 near t=1;    same; delta ~0.05 flat
                              delta ~0.06-0.07 flat from t=0.7
touching tube pair            0.56@2.40, 0.40@3.75, rest at the    0.56@2.15 0.40@3.25 0.28@3.75 0.20@4.05
(apex gap 1 sigma, T=6)       clock (6.0); max|w| 1.7x             0.14@4.30, then 0.1/0.07/0.05 at the clock (4.5)
                                                                   per octave 2.27, 0.97, 0.62, ~0.5; max|w| 3.4x
```

(a) Kida-Pelz: no wave, as registered - but the numeric clause (twist <= 0.01 at every separation) is violated at
0.4 and 0.56 by the pattern's own periodicity growing in place; recorded as violated. The C20 KILL fired there under
the old clause (race variable from a row where the compression across the "sheet" was ~0, with max|w| falling):
spurious, and the reason the clause was corrected. (b) The second field: nu-independent within 0.05 and
accelerating, so C23(b) holds - but not the sheet field's smooth 1/ell law: the octave time drops from ~0.9 to ~0.1 in
one step near t = 0.9, and the wave stops where the gap meets this field's flat thickness (~0.06), at t ~ 1.0 and
scale ~0.1 - above sqrt(nu/s) at 1e-3 (0.04-0.05), at it at 2e-3. Same V, T* ~ 1.3 instead of 1.78. (The in-script
V line for this run is invalid: v4 still counted the t = 0 maximum at 0.56 as a peak.) (c) The touching pair: not the
registered approach test (audit finding F), but at nu = 1e-3 it shows the two phases in one run anyway - 2.3 per
octave at the top, then 0.97, 0.62, ~0.5 as the tubes close - and the wave reaches 0.1 exactly as the clock expires
(t = 4.5), with max|w| 12 -> 42 inside the clock, 62 past it. Under the corrected KILL clause it does not fire:
max|w| is decelerating at the clock (36.7, 40.5, 41.9). The clean approach (apex gap 2.5 sigma) is v5.

### The approaching pair (C23c), v5 (`results/seam_gpu/v5/`, 2026-09-08)

Apex gap 2.5 sigma, T = 8, 256^3.

```
peak time of twist@sep      nu = 2e-3     nu = 1e-3      per octave (1e-3)
  0.56                       4.10          3.55
  0.40                       6.75          5.35           3.7
  0.28                       (clock)       6.75           2.7
  0.20                       (clock)       ~7.9 plateau   2.4
  0.14                       (clock)       7.35
  0.10                       (clock)       7.70           0.6-0.7
max|w| at t = 8              0.78x         2.76x
```

A descending wave, and the cleanest two-phase record on the page: 3.7, 2.7, 2.4 per octave while the tubes close under
mutual induction, then 0.6-0.7 in the last two octaves - a four-fold jump into the self-induced phase, on a hand-built
structure with no adversary involved. The nu-independence clause is **refuted** for the pair: peaks shift by 0.5-1.4
between the viscosities, because at 2e-3 the Gaussian cores diffuse over an 8-unit run (max|w| falls to 0.78x) and the
approach slows. The descent law is inviscid when the descent is fast (the sheet field, one time unit); a slow descent
gives viscosity time to act on the cores. The V at 1e-3: gap = 0.099 (9.25 - t), five times slower than the sheet
field's wave, meeting the thickness arm at t ~ 8.8, scale 0.048, above sqrt(nu/s) = 0.029 - beyond the clock, not
counted. The old C20 KILL clause fired at 2e-3 with max|w| falling 22%; under the corrected clause it does not.

### The 320^3 rung (`results/seam_gpu/v7/`, 2026-09-08): C21 KILL, the V corrected

```
                        256^3 (biased thickness fit)        320^3 (pre-merge fit)
nu = 2e-3   peaks       0.70 0.95 1.25 1.40 1.50 1.60      0.70 0.95 1.25 1.40 1.50 1.60 1.65 1.70   (identical)
            V           t 1.65, scale 0.068, ell_nu 0.053  t 1.69, scale 0.052, ell_nu 0.051
nu = 1e-3   peaks       0.70 0.95 1.25 1.35 1.50 (0.1 at   0.70 0.95 1.25 1.35 1.50 1.60 1.65 1.70   (0.1 resolved:
                        the clock)                          clock 1.75; peak 0.182 at 1.60, falls 15%)  0.182 at 1.60)
            V           t 1.70, scale 0.043, ell_nu 0.025  t 1.72, scale 0.031, ell_nu 0.028
Re_seam at the peak     224 (2e-3)                          594 (1e-3):  x2.65 per halving of nu
```

Two corrections. (1) The V bottom sits **on** sqrt(nu/s) at both viscosities once the thickness arm is fitted before
the merge (audit D) - the 256^3 reading "above the viscous scale at 1e-3" was the biased fit and is retracted; the
merge is the cut at both rungs, and the merge scale tracks sqrt(nu/s). (2) C21's KILL clause fired: Re_seam at the
turnover rises x2.65 per halving of nu, i.e. the velocity jump across the seam at the merge grows (Delta_u 2.2 -> 3.3)
as viscosity falls. Recorded as the first registered number in the repository pointing toward a singularity - across
viscosities, at the merge; whether the jump grows along one solution is C26 (the Lagrangian frame), running. The
320^3 run at 2e-3 reproduces the 256^3 descent to the sampling step: the picture is resolution-converged.

### Instrument audit of `seam_gpu.py` (2026-09-08, before the generality run is read)

Six findings, all fixed in the next version; none touches C22 (its peak times are read straight off the twist
columns, which are correct). (A) The GPU port printed no energy - the numpy solver's 1e-16 was the correctness
evidence and the complex64 port had never been asked; E/E0 is now a column and the nu = 0 run is its certificate. **Closed (v5b):** at nu = 0, E/E0 = 1.000000 at every output to T = 1 at both 256^3 and 128^3; and Z(1)/Z0 = 3.175 at 128^3 / 3.184 at 256^3 against the numpy solver's 3.175 at 128^3 / 3.182 at 192^3 for the same field - two solvers, same growth to three digits.
(B) Every initial condition has an empty spectral tail at t = 0, so the strip fit read float32 round-off until the
cascade arrived: the "past the clock" flags before t ~ 0.25 were this, not the clock; now "(tail empty)". (C) The
in-script gap arm used the wall-affected peaks at 0.07 and 0.05; now sep >= 0.1 (the hand fit above already did).
(D) The thickness arm was fitted through the post-merge plateau, biasing its e-fold (3.35 vs 2.1 pre-merge); now
fitted from the first real delta to the twist peak. (E) A column maximal at t = 0 counted as a peak; excluded.
(F) The tube pair's apex gap was D - 2A = 0.2 = one core radius: the tubes overlapped at t = 0, so the v4 pair rows
test a touching pair, not an approach; the approach (apex gap 2.5 sigma, T = 8) is rerun as v5.

## The forced seam race (C24), v6 (`results/seam_gpu/v6/`, 2026-09-08)

A smooth, steady, divergence-free force f = eps x (the initial field, |k| <= 4) - Fefferman's (C)/(D) admit such an f -
held on for the whole run. Energy injection at t = 0: 0.022 (eps 0.5) and 0.066 (eps 1.5) against a viscous drain of
0.0015. Registered: the V bottom stops being the end of the descent (twist@0.05 keeps rising past t ~ 1.7, max|w|
accelerates, the wave passes sqrt(nu/s) with the twist still rising).

```
nu = 2e-3, eps = 0.5, resolved to T = 3 (delta 0.051-0.078):
  t      Z/Z0    max|w|   twist@0.1   twist@0.05   anti    race   cut    E/E0
 0.9     3.32     26.7     0.004       0.00002      0.10    2.7    0.21   1.58
 1.2     4.50     29.0     0.057       0.0027       0.44    5.0    0.29   1.74
 1.5     5.43     46.3     0.106       0.019        0.57    4.4    0.28   1.88     <- twist@0.1 peaks at 1.4 (0.110), as unforced
 1.8     5.78     52.1     0.068       0.010        0.42    2.5    0.30   2.01
 2.4     6.45     69.2     0.050       0.007        0.60    2.0    0.29   2.25
 2.7     6.89     73.9     0.070       0.012        0.73    2.4    0.30   2.37     <- a second seam forms
 3.0     7.21     70.6     0.045       0.006        0.73    1.2    0.49   2.47
eps = 1.5 at 2e-3: the clock expires at 1.25 (the pumped cascade fills the tail); 7.9x at the clock, twist@0.1 at its
peak there. nu = 1e-3: clocks at 1.10 (eps 0.5) and 0.90 (eps 1.5); nothing past them is a number.
```

**C24 refuted on its mechanism clauses.** With the force on, the twist@0.1 still peaks at t = 1.4 and falls 59%; the
twist@0.05 peaks at 1.5 and falls; the race variable never drops below 1.18; the cut fraction sits at 0.3-0.5. The seam
is cut exactly as it was unforced. What changes is what happens *after*: unforced, max|w| peaked at 41 and decayed and
Z peaked at 3.2; forced, the pumped large scales rebuild seams (a second twist rise at 2.4-2.7) and max|w| climbs to 74
at t = 2.7 (9.9x at the clock, still rising) with Z at 7.2 and the energy at 2.5x. A steady smooth force gives a driven
flow that keeps cutting and rebuilding - forced turbulence - not a collapse that passes the floor. A forced blow-up
needs a force that tracks the collapse in space and time, which is how the Cordoba-Martinez-Zoroa forces are built.
That is the refutation clause as registered, and it is the informative outcome: the forced route is not "push harder",
it is "push exactly where the self-similar solution needs it".

## Tao's wall, in pictures: an energy-conserving equation that provably blows up

Theorem 4 (Tao 2016) says that the exact structure this repository verifies - energy conservation, the scaling, the
quadratic skew nonlinearity - is not enough to rule out a singularity, because an equation with exactly that
structure blows up. The simplest such equation is the dyadic shell model of Katz and Pavlovic (2005), and its
inviscid blow-up is a theorem (Katz-Pavlovic 2005; Kiselev-Zlatos 2005). `dyadic.py` runs it through the same
instrument: RK4 with exact viscous integrating factor, energy drift reported, the analyticity width delta(t)
fitted from the shell spectrum with a reliability rule, and the same two decay laws fitted inside the window that
`strip_tracker.py` fits for Kida-Pelz.

![The dyadic model and Kida-Pelz Euler through the same diagnostic](../figures/taos_wall.png)

```
                              energy drift    delta(t)                      decay rate across window     verdict
dyadic, nu = 0, 34 shells     3e-8            straight to 0 at t* = 0.5585   7.8 -> 1.4e6  (x 178000)     singularity (theorem)
dyadic, nu = 1e-6 / 1e-4      (dissipates)    front stalls at shell 12 / 7   turns negative               no singularity
Kida-Pelz Euler 64/96/128^3   1e-8            exponential, tau 0.42          4.1 -> 1.9  (falls by half)  no singularity indicated
```

Three things this shows. Energy is conserved to 3e-8 while the enstrophy goes as (t* - t)^-1.7 to the truncation:
zero numerical entropy production is a property of the *instrument*, not evidence about the *equation*. The decay
rate of delta is the discriminator that separates the two cases by five orders of magnitude, so when it falls for
Kida-Pelz that is a statement, not noise. And the dyadic model does not satisfy Liouville (div F = -sum k_{j+1}
a_{j+1}, state-dependent, printed on every row), while the true truncated Euler system does exactly. That is *not* a
property that separates the true equation from Tao's class: a multiplier-averaged nonlinearity M B(Mu, Mu), a member
of Tao's averaged family, conserves phase volume to 1e-17 as well (`liouville.py`, form `tao-class`;
`results/liouville_taoclass.txt`), because Liouville follows from the convolution structure and the dyadic model is
not a convolution. So phase-volume conservation sits on Tao's side of the wall along with energy. What the true
equation has and the averaged one lacks remains the short list Tao averaged away: locality in physical space,
particle paths, and the geometry of stretching. A proof has to be made of those.

## The minimal blow-up datum: an objective the searcher taught us not to pose

Rusin and Sverak (2011) showed that if any critical-space datum blows up, one of minimal critical norm does.
`OBJ=minimal` tried to hunt it: minimise |u0|_{H^1/2} subject to the flow reaching the grid's cutoff by T. Three
versions, three exploits by the searcher, all caught by the verifier (`results/minimal_datum_*.txt`): the first
shrank the amplitude until the strip fit read the round-off floor as delta = 0; the second was trivially met because
the "tail" band overlapped the initial data; the third - energy fraction 1e-4 in the top sixth of the modes - was met
by a tiny field whose weak third-order interactions leak a trace into that band (32^3: norm 0.071; 48^3: 0.330,
"rising with resolution" as predicted, for the wrong reason; both completely smooth at 64-128^3, delta 0.2-0.75).
The result is the lesson: every finite-resolution criterion for "this datum becomes singular" has a cheap
non-singular satisfier, and a gradient searcher finds it. The minimal blow-up datum is defined by loss of
regularity, which is exactly what no finite grid certifies; a numerical version needs a criterion that is the open
question itself. This is CKN's untestability at finite resolution, in the searcher's mouth. Closed, recorded.

## A boundary of the method, from the recent literature

Wang, Lai, Gomez-Serrano and Buckmaster (arXiv 2509.14185, 2025; 2511.22819) found new families of *unstable*
self-similar singularities for Boussinesq, IPM and Euler with boundary by searching profile space directly with
physics-informed networks, and state the hypothesis that unstable singularities are the route for boundary-free
Euler and Navier-Stokes. An unstable singularity lives on a measure-zero set of initial data. A gradient search from
generic smooth data - which is what `adversarial_ic.py` is - cannot find one by construction; it finds the fastest
*growth*, and what it found is sheets and tubes. Every "no singularity indicated" on this page is a statement about
the flows run and the windows resolved, never about the existence of unstable singularities, which this instrument
does not look for. Barker (arXiv 2510.20757, 2025-26) gives quantitative bounds near a potential blow-up for Hou's
approximately axisymmetric Navier-Stokes candidate that are explicitly amenable to numerical testing; that candidate
lives in a cylinder with a boundary, outside this periodic box, and is the natural next target for an instrument
built like this one. Shahmurov (arXiv 2604.09949, 2026) proves for axisymmetric model equations that the sign of the
elliptic response decides blow-up - the physical sign is global, the reversed sign blows up - which is the analytic
face of the projection experiment above.

## The field on 2026-09-08: the forced route

Two developments this week bear directly on how to read everything above.

**Verified.** Buckmaster and Alpöge posted finite-time blow-up *with smooth forcing* for incompressible porous media,
Boussinesq and 3-D incompressible Euler, Lean-verified (Euler on 2026-08-22), built on the Córdoba-Martínez-Zoroa
program of forced blow-ups and pushed to smooth forcing with heavy use of language models. They report a
hypo-dissipative Navier-Stokes blow-up in preparation. (Buckmaster's statement: cims.nyu.edu/~tristanb/statement.pdf.)

**Claimed, unseen.** Per the same statement, OpenAI reported to Buckmaster an internal ~100-page proof of forced
blow-up for Navier-Stokes on R^3 and T^3 with smooth forcing. Nobody outside has seen it as of this writing.

**Why forcing is the Clay problem and not a loophole.** Fefferman's official statement offers four options: (A), (B)
existence and smoothness on R^3 / T^3 with f = 0; (C), (D) breakdown on R^3 / T^3 for *some* smooth divergence-free
u0 and *some* smooth f(x, t) with the stated decay. A correct forced blow-up on T^3 is statement (D).

**What it means for this page.** Everything here concerns the unforced dynamics - the physical question, (A)/(B). The
leak budget, the free collapse, the seam race, the descent law and the V are all statements about a fluid nobody is
pushing. With a smooth force the budget is replenished from outside: the external squeeze (phase 1 above) can be
held on indefinitely and never has to hand over to the self-induced roll-up; the thickness arm can be driven down
without the cut winning any race. The forced route goes around the wall by removing the one assumption every
energy-based argument rests on. If the forced proofs hold, the prize is resolved on the negative side and the
question "does water do this on its own" stays exactly where it was, with the measurements above as one answer's
worth of evidence about the mechanism.

## What this is and is not
- It is a measurement of an **instrument property**: no artefact dissipation. Numerical searches for self-similar
  blow-up (Hou; Gomez-Serrano, Buckmaster et al. 2022; and later neural-network-assisted searches) are limited by
  exactly this artefact, which damps the small-scale growth they are trying to detect.
- It is **not** a statement about 3-D Navier-Stokes regularity. Nothing here proves anything; a proof is a theorem
  or nothing.

## Further reading
[CLAIMS.md](CLAIMS.md) - every claim, its refutation condition, its command, and our weak points.

[EQUATIONS.md](EQUATIONS.md) - the one-family system in $d=1,2,3$, the skew/contractive split, the discretisation, and the budget line with its
production term per dimension: [EQUATIONS.md](EQUATIONS.md). One solver for all three: `ns_d.py`.

[THEORY.md](THEORY.md) - definitions, the elementary propositions with proofs, the classical limits (Tao 2016; Beale-Kato-Majda), and the one open
hypothesis stated attackably: [THEORY.md](THEORY.md).

## Run
```
python burgers_entropy.py      # ~1 minute, CPU
python taylor_green.py         # a few minutes, CPU
```
numpy only.
