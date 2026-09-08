# Glossary

Every symbol, measure and claim number used on these pages, one line each. Conventions: the domain is the periodic
box $[0,2\pi)^3$, $\nu$ the viscosity, $\omega = \nabla\times u$ the vorticity, $S$ the strain (symmetric part of
$\nabla u$), $\xi = \omega/|\omega|$ the vorticity direction. All 3-D runs are normalised to the same initial
enstrophy $Z_0 = 0.375$.

## Flow quantities

| symbol | meaning |
|---|---|
| $E$ | kinetic energy $\tfrac12\langle\lvert u\rvert^2\rangle$. Unforced Navier-Stokes: $dE/dt = -2\nu Z$ exactly; the solver reproduces this to 1e-6 and conserves $E$ to six digits at $\nu = 0$. |
| $Z$ | enstrophy $\tfrac12\langle\lvert\omega\rvert^2\rangle$; $Z/Z_0$ is the growth reported everywhere. In 2-D non-increasing (regular); in 3-D it can grow by vortex stretching. |
| max$\lvert\omega\rvert$ | the peak vorticity; blow-up needs $\int_0^T \max\lvert\omega\rvert\,dt = \infty$ (Beale-Kato-Majda). |
| $H$ | helicity $\langle u\cdot\omega\rangle$; $H/H_{\max}$ with $H_{\max} = 2\sqrt{E Z}$. Aligned (Beltrami-like) flow has $\lvert H/H_{\max}\rvert \to 1$ and weak stretching. |
| $\delta$, the strip | analyticity-strip width: the exponential decay rate of the energy spectrum's tail (Sulem-Sulem-Frisch). The distance of the nearest complex singularity from the real axis. |
| **the clock** | a 3-D number is trusted only while $\delta > 2\,dx = 4\pi/N$. Past that the grid cannot represent the flow and the number is discarded ("past the clock"). "(tail empty)" marks the start-up before the cascade fills the spectrum, when the fit has nothing to fit. |
| $\alpha_s$ (CKN exponent) | how dissipation concentrates around the vorticity maximum: $\int_{B_r}\lvert\nabla u\rvert^2 \sim r^{\alpha_s}$. Uniform 3, sheet 2, tube 1, point 0; a singularity needs $\le 1$. Measured on the adversary's field: ~4. |
| $s$ | strain rate; in the seam diagnostics the mean compression across the sheet, $-n\cdot S\,n$, on the high-vorticity set. |
| the high set | $\{\lvert\omega\rvert > \tfrac12\max\lvert\omega\rvert\}$: where the seam diagnostics are evaluated. |

## Seam diagnostics (`seam_race.py`, `kaggle/seam/seam_gpu.py`)

| symbol | meaning |
|---|---|
| $\beta$ | signed angle between the vorticity direction at a point and at a neighbour: $\beta = 1 - \xi\cdot\xi'$. 0 aligned, 1 orthogonal, 2 antiparallel. $\beta > 1$ is a reversal. Squared-coherence measures cannot see it. |
| **twist@sep** ("strong twist") | enstrophy-weighted reversal at a fixed physical separation: mean of $\lvert\omega\rvert^2\lvert\omega'\rvert^2\,(\beta-1)_+^2$ over neighbours at that separation, normalised by $\langle\lvert\omega\rvert^4\rangle$. Large only where two *intense* regions point opposite ways - a seam. Separations 0.05 ... 0.56 in half-octave steps; the verdict column uses 0.1. |
| **anti** | the hard count: fraction of the high set with an antiparallel neighbour ($\beta > 1$) at separation 0.1. Sees weak reversals the strong twist ignores. |
| $\ell$ | sheet thickness, $\lvert\omega\rvert/\lvert\nabla\lvert\omega\rvert\rvert$, median on the high set. |
| $\ell_\nu$ | viscous thickness $\sqrt{\nu/s}$: the scale at which viscosity cuts a sheet pressed at rate $s$ (the Burgers-layer scale). |
| **race** | $\ell/\ell_\nu$: the seam reaches its viscous thickness at ~1. |
| **cut** | on the twisted high set, the fraction where viscous cancellation $-\nu\,\omega\cdot\Delta\omega$ exceeds the *positive* stretching $\omega\cdot S\omega$: reconnection winning locally. |
| **Re_seam** | $\lvert\omega\rvert\,\ell^2/\nu = \Delta u\,\ell/\nu$ on the twisted set: the seam's own Reynolds number. The gain of the twist's feedback loop changes sign at a fixed value of it (C21). |
| **the seam** | where two sheets meet antiparallel: $\lvert\omega\rvert \to 0$ across a surface, $\beta > 1$ across it. Biot-Savart cancels there, Kelvin-Helmholtz rolls it, viscosity reconnects it. |
| **the twist wave** | the reversal measured at fixed separations peaking in sequence from large to small: one structure descending through the scales. |
| **gap arm** | the line through the peak times: $g(t) = a\,(T^\ast - t)$ on the sheet field; inviscid (same at all $\nu$). |
| **thickness arm** | $\delta(t)$, the sheets' own thickness; exponential, viscous. |
| **the V** | where the gap arm meets the thickness arm: the pair merges into one doubled structure before the gap reaches zero. |
| $\lambda$ | the collapse exponent: gap $\propto (T^\ast - t)^\lambda$. The sheet field's closing is exponential relaxation onto the thickness (C32), for which $\lambda = 1$ was a short-window tangent. The critical class is $\lambda = \tfrac12$ (Theorem 7). |
| $T^\ast$ | the zero of the straight line fitted to the peak times (1.78 for the sheet field). Not an arrival time: the true law is exponential and never reaches zero. |
| $k$ | the relaxation rate of the gap onto the thickness, $dg/dt = -k(g - \delta)$; $1.7 \pm 0.1$ on the sheet field at $
u = 2	imes10^{-3}, 10^{-3}, 5	imes10^{-4}$; $pprox 0.16$ and not clean on the tube pair (C32). |
| **time per octave** | time for the wave to descend one halving of separation. Shrinking in proportion to the separation = constant closing speed (lambda = 1) seen on a log axis; equal octaves = exponential thinning. |
| **phase 1 / phase 2** | the budget argument: an external squeeze at fixed rate thins exponentially and stops at $\sqrt{\nu/s}$; a collapse whose speed does not fall is unpriced by the budget. The sheet field's descent is constant-speed closing driven by the sheets' own (cancelling) induction - not self-accelerating; that earlier reading is withdrawn (C30). |
| **the floor** | the length below which a fluid cannot be squeezed: the healing length (quantum), $\sqrt{\nu/s}$ (Navier-Stokes, moving), none (Euler). |
| **FORCE** ($\varepsilon$) | a smooth, steady, divergence-free force $\varepsilon\times$(the initial field's modes with $\lvert k\rvert \le 4$), Fefferman's (C)/(D) admit such an $f$. |

## Searcher and verifier (`adversarial_ic.py`)

| term | meaning |
|---|---|
| **the adversary / searcher** | gradient ascent (autograd through the solver) on smooth low-$k$ initial data for the largest $Z(T)/Z_0$, or for other objectives (helicity, CKN exponent, quiet pressure, minimal norm). |
| **the leash** (DMIN) | a penalty keeping the strip $\delta(T)$ above a floor on the search grid, so the searcher cannot win by leaving the resolved regime. |
| **the verifier** | the same initial datum re-run at 128^3-192^3 with the clock; "resolved" means $\delta > 2dx$ there. |
| **found field** | an initial condition the searcher produced, stored in `results/found/*.npz`; `leashed64_dmin030` is "the sheet field" used throughout the seam work, `ckn64` "the second field". |
| **TWISTW** | the twist penalty: forbid strong reversals in the objective (C19). |

## The claims (`CLAIMS.md`)

Each was registered with the observation that would refute it before its run; status in brackets.

| # | one line | status |
|---|---|---|
| C1 | numerical entropy production is zero to round-off | stands |
| C2 | every balance law closes term by term | stands |
| C3 | 1-D Burgers follows the exact 1/(1-t) to the grid limit | stands |
| C4, C4b, C4c | closed-form solutions reproduced; time integration converged | stand |
| C5 | 3-D vortex stretching converges upward with resolution | stands |
| C6 | vorticity aligns with the intermediate strain eigenvector | stands |
| C7 | helicity conserved to the order of energy | stands |
| C8 | a learned simulator with frozen exact transport holds energy | stands |
| C9 | (negative) an eddy-viscosity closure learns nothing in 2-D | stands |
| C10 | skew-projected learned correction is energy-neutral to first order | stands |
| C11 | Kida-Pelz early evolution converged; direction roughens at a fixed scale | stands |
| C12 | the searcher finds data the verifier rejects as unresolved | stands (it has, repeatedly) |
| C13 | attainable growth flat in helicity to half the maximum, then collapses | stands |
| C14 | the strip separates a proven singularity from Kida-Pelz by its decay trend | stands |
| C15 | (negative, re-scoped) the Lyapunov learner cannot find a monotone dominating quantity | stands as a statement about the method |
| C16 | the fastest field concentrates dissipation like a sheet | stands |
| C17 | fast growth requires and recruits the global pressure | stands; closing statement C17b unsupported |
| C18 | the Gross-Pitaevskii seam reconnects with the known gap law and releases sound | stands; "E_ki falls at the cut" retracted |
| C19 | growth needs the twist | 32^3 met, **64^3 refuted** (30% cut, not 50%); re-scoped to the strong measure |
| C20 | the seam reaches its viscous thickness and the twist turns over | twist and race clauses met, max\|w\| clause missed; the cut wins at Re ~ 500 |
| C21 | Re_seam at the twist peak is a fixed number across $\nu$ | inconclusive (only 2e-3 resolved at 256^3); 320^3 rung running |
| C22 | the descent law: inviscid and accelerating | **passes** on every clause |
| C23 | the descent law generalises (Kida-Pelz none; second field; pair) | (a), (b) hold; (c) two phases hold, $\nu$-independence refuted for the slow pair |
| C24 | a steady smooth force carries the seam through the floor | **refuted**: it drives, rebuilding seams; growth by re-supply |
| C25 | (conjecture, `THEORY.md` §9-10) the seam's velocity jump is bounded by the data, i.e. finite total in-plane compression | holds along five resolved solutions (C26) |
| C26 | Kelvin's frame: the Lagrangian gap closes linearly, the jump stays bounded | **passes** x5 (lambda 0.99; jump within 9%) |
| C27 | a tracking force carries the seam through the floor | between: faster, same shape, still cut at 0.1; a reversal at 0.05 rising at the clock |
| C28 | the seam flips again below the cut | between: the cut is a sign reversal on material fluid (0.13 by the clock, threads beside the sheets); no second seam inside the clock |
| C29 | the sheets close under their own induction | between at both nu: 0.54 / 0.52 from the top half of the sheets; 0.84 from the whole sheets |
| C30 | the pair closes by its own curvature (Da Rios) | **refuted**: curvature rises x5 while the induced closing falls to zero - antiparallel induction cancels on approach |
| C31 | the far field is part of what closes the gap | between on the sheet field: 0.16 at all three nu, rising to 0.29 at the merge; not C17's 0.42. KILL on the pair: the far field opposes (-0.21) |
| C32 | the gap relaxes exponentially onto the thickness, k viscosity-independent | **passes**: k = 1.7 +- 0.1 at 2e-3, 1e-3, 5e-4; the pair k ~ 0.16 |
| C33 | the compression across the sheet is a dipole gradient, s ~ g^-3; the seam is one slow variable | closure withdrawn (estimator); the model gives the merge scale from the measured strain, not yet its time |
| C34 | the strain budget closes on the material sheets | **KILL** at 2e-3 and 1e-3, outcome (ii): viscous diffusion removes ~half the strain input; t.S.t has the opposite sign to the model |
| C35 | the diffusive loss fraction falls with viscosity and does not vanish | **passes** on the fall: 0.48, 0.33, 0.07, 0.00 - it does vanish; the seam has an Euler limit |
| C36 | the seam rolls (Lundgren's spiral) | **KILL**: 0.06-0.11 turn at one-twentieth of the solid-body rate with the vorticity growing - the antiparallel pair is stabilised against K-H and merges flat |

## Theorems and propositions (`THEORY.md`)

Propositions 1-3 (energy inequality, discrete energy inequality, frozen memory carries a conserved quantity),
Theorem 4 (Tao: the abstraction is insufficient for regularity), Theorem 5 (what suffices: a second controlled norm,
conditional), Proposition 6 (Liouville for the truncated system), **Theorem 7** (the race: a collapse with exponent
$\lambda > \tfrac12$ is inside its viscous scale before $T$ under Type I strain; $\lambda = \tfrac12$ is decided by
Re_seam).

## The problem statement

Fefferman's four options: **(A)**, **(B)** existence and smoothness on $\mathbb{R}^3$ / $\mathbb{T}^3$ with $f = 0$;
**(C)**, **(D)** breakdown on $\mathbb{R}^3$ / $\mathbb{T}^3$ for *some* smooth data and *some* smooth force with the
stated decay. Everything on these pages is about (A)/(B) unless marked "forced". **Type I**: $\sup\lvert\nabla u\rvert
\le A/(T-t)$, the self-similar rate; excluded under axisymmetry, open in general.
