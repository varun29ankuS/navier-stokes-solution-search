# Scripts by theme

Every script is standalone, run from the repository root, configured by environment variables named in its
docstring, and prints a registered verdict where one applies. Result files are the script's stdout under `results/`.

## The solver and its verification

- `ns_d.py` - One solver for the incompressible Navier-Stokes family in d = 1, 2, 3 dimensions on the periodic box [0,2pi)^d:
- `exact_solutions.py` - Three more exact solutions, one per dimension, compared with the scheme to machine precision where it can.
- `taylor_green.py` - Taylor-Green vortex with the frozen-rotation spectral scheme. Two parts.
- `budgets.py` - Budget closure: with zero numerical entropy production, every balance law should close to machine precision
- `budgets3d.py` - 3-D Euler enstrophy budget with the vortex-stretching term measured directly, Taylor-Green initial data.
- `midpoint.py` - Exact discrete energy conservation: the implicit midpoint rule on the skew transport.
- `verify_skew_closure.py` - Verification of the claim in EQUATIONS.md: a learned correction projected to be skew cannot change the energy
- `liouville.py` - Liouville's theorem for the truncated system (Lee 1952; Kraichnan 1973): the Galerkin-truncated Euler equations are a
- `spectra.py` - Energy spectra: the two turbulence laws everyone recognises, produced by the energy-conserving scheme.

## 1-D and 2-D controls

- `burgers_entropy.py` - Does zero numerical entropy production let an integrator SEE a blow-up that a dissipative one hides?
- `monotone_vs_dominating.py` - Regularity needs a quantity that is both MONOTONE and DOMINATING; do the two sets intersect?
- `lyapunov2d.py` - Positive control for the Lyapunov search: 2-D, where a monotone quantity is KNOWN to exist.
- `memory_paths.py` - Water's memory as an algorithm: Navier-Stokes solved by remembering Euler along noisy particle paths (2-D).
- `wave_particle.py` - Wave and particle: the same flow in both descriptions, and what the phases decide.
- `frequency_matching.py` - Does incommensurability protect a flow? Two-scale initial data, commensurate (k, 2k) versus incommensurate (k, ~sqrt2 k).
- `kolmogorov_closure.py` - Learned coarse simulators of 2-D Kolmogorov flow: does freezing the conservative transport and learning only the
- `kolmogorov_v2.py` - v2. Learned coarse simulators of 2-D Kolmogorov flow. Changes from v1 (which found the dissipative gate inert):
- `kolmogorov_v3.py` - v3: as v2, but the learned redistribution is made EXACTLY energy-neutral by rescaling the state after the correction

## The adversarial searcher

- `adversarial_ic.py` - Feedback, not imitation: search for the initial condition that amplifies enstrophy the most.
- `lyapunov_search.py` - Searching for the missing functional, with an adversary in the loop.

## 3-D diagnostics along a flow

- `criteria3d.py` - The regularity criteria, measured along 3-D Euler flows (nu = 0 unless NU is set).
- `mechanism3d.py` - 3-D Euler: the mechanism of vortex stretching, a second invariant, and the BKM quantity.
- `strip_tracker.py` - Complex-singularity tracking (Sulem, Sulem & Frisch 1983; Frisch, Matsumoto & Bec 2003): the energy spectrum of an
- `jacobi_ladder.py` - Arnold's picture, measured: Euler flow is a geodesic on the group of volume-preserving maps, and the growth of a
- `ckn_exponent.py` - Keep the state, track the dissipation: the Caffarelli-Kohn-Nirenberg concentration exponent.
- `pressure_share.py` - Does the global inform the local? The pressure Hessian P = grad grad p splits into a LOCAL part, its trace
- `projection_price.py` - The price of the projection: the same nonlinearity with and without the Leray projection, from the same data.
- `nilpotent_sheet.py` - Two claims behind 'the sheet is where the local map is silent and the pressure is the whole dynamics', measured.
- `mandelbrot_fraction.py` - The quadratic map inside Navier-Stokes, and the fraction of the flow the pressure holds inside the set.
- `thinnest_squeeze.py` - The thinnest squeeze: the smallest singular value of the Lagrangian deformation gradient, against the jitter scale.
- `memory_paths3d.py` - Water's memory with stretching: 3-D Navier-Stokes from Cauchy's formula along jittered particle paths.
- `friction_function.py` - The constant that is really a function: the same flow under four frictions.

## The seam

- `seam_race.py` - The seam race: does viscosity cut the twist before the roll-up completes?
- `kaggle/seam/seam_gpu.py` - The seam race at the viscous thickness: GPU spectral Navier-Stokes, 256^3, tracking the strong twist until it turns.
- `gpe_seam.py` - The quantum seam: the antiparallel-pair reconnection in a fluid where the cut is guaranteed to win.

## Tao's wall

- `dyadic.py` - What Tao's wall looks like: an energy-conserving equation that provably blows up, run through the same instrument.

## Figures

- `plots.py` - Figures for the README, regenerated from the same code that produced the tables. Writes figures/*.png.
- `plot_dyadic.py` - figures/taos_wall.png: the dyadic model (energy-conserving, proven blow-up) and Kida-Pelz Euler (open question), through
- `plot_helicity.py` - figures/helicity_threshold.png from results/helicity_proj_*.txt: maximal enstrophy amplification over one time unit
- `plot_jacobi.py` - figures/jacobi_ladder.png from results/jacobi_*.txt: separation growth and local exponent by resolution, hollow markers
- `plot_strip.py` - figures/strip_decay.png from results/strip_*.txt: the analyticity-strip width delta(t) on a log axis by resolution
- `flow_gif.py` - Animated flows from the same solver, as GIFs, with the conserved quantities printed on every frame.
- `vortex_iso.py` - 3-D Taylor-Green Euler: isosurfaces of |vorticity| at t = 1, 2, 3, 4 - the vortex sheets forming.

## Running on CI

`.github/workflows/experiments.yml`, dispatched with `which=<job>`; the job names are listed in the workflow's `description`. The commit step pulls with `-X theirs` so concurrent jobs do not race on `results/`.
