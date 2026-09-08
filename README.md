# navier-stokes-solution-search

A search, not a solution: an instrument for the mechanism of singularity formation in incompressible flow, and an
honest log of what it found. (Formerly `zero-entropy-flow`, after the solver's defining property.)

One pseudo-spectral solver for the Navier-Stokes family in 1-D, 2-D and 3-D whose own numerical dissipation is zero
(skew transport exactly energy-conserving on the grid, viscosity an exact per-mode contraction), a differentiable copy
of it that searches initial data adversarially, a higher-resolution verifier with a clock that says when to stop
believing a number, and a GPU port certified against the numpy original. Every claim on these pages was registered
before its run and is scored against that registration in `CLAIMS.md` - twenty-four so far, about half retracted in
public. Nothing here claims anything about the regularity of 3-D Navier-Stokes.

**Read in this order:** this page - [`docs/MAP.md`](docs/MAP.md) (the view from above) - [`CLAIMS.md`](CLAIMS.md)
(what was predicted, what happened) - [`docs/LOG.md`](docs/LOG.md) (the full chronological log, 1200 lines, with every
table) - [`THEORY.md`](THEORY.md) (the propositions proved and the theorem that bounds the method).

## The standard

- **Registered.** A prediction is written down, with the observation that would refute it, before the run.
- **Clocked.** A 3-D number is reported only while the analyticity-strip width exceeds two grid spacings; past that it
  is not a number, and the page says so. This discarded the most exciting rows this repository ever produced.
- **Retracted in the open.** Wrong predictions stay on the page with the word "wrong" next to them.
- **Cross-checked.** Two solvers (numpy, PyTorch/GPU), two precisions, same growth to three digits; energy 1.000000.

## What it found, unforced

1. **Growth builds sheets, and the sheets are sub-singular in every view.** The adversary's fastest fields are
   antiparallel vortex sheets; their dissipation concentrates like a sheet (CKN exponent ~4, a singularity needs ~1),
   their strip decays as a smooth flow's, their coherence is preserved, and the pressure they recruit is global.
2. **Every energy-class candidate for a monotone-dominating quantity fails, for a structural reason.** Learned
   functionals with a passed 2-D control, topology, helicity, memory functionals, the twist: each lands in the
   dominating set and never in the monotone one. In 1-D Burgers the two sets are disjoint and it shocks; in 2-D they
   intersect and it is regular; in 3-D no member of both was found in any family tried. Tao's theorem says why.
3. **The twist helps growth and is not necessary for it.** Forbidding strong antiparallel reversals cut attainable
   growth by 30% at 64^3 (a registered 50% was refuted) and produced the first resolved found fields (1.19, 1.52, 2.24x).
   Denied the twist the searcher walks to the helicity plateau, where stretching dies.
4. **The seam descends through the scales as a wave, and the descent is a measurable inviscid law.** Twist measured at
   fixed separations peaks in sequence 0.56, 0.40, ..., 0.05; the peak times agree across nu = 2e-3, 1e-3, 5e-4 within
   one sampling step; the time per octave shrinks from ~0.55 to ~0.25 - rhythmic in the collapse's own clock,
   accelerating in ours; gap = 0.51 (1.78 - t). A hand-built tube pair shows the two phases in one run (3.7, 2.7, 2.4
   then 0.6-0.7 per octave). A second adversarial field descends in one abrupt step. Kida-Pelz shows no wave.
5. **The V.** The gap arm (inviscid, linear) meets the thickness arm (viscous, exponential) at t = 1.65-1.74 - at the
   viscous scale for nu = 2e-3 (the merge is the cut: enstrophy and max|w| peak and decay, resolved to T = 3), above it
   for 1e-3. The Euler seam race starts there, located from data above the wall at a place the clock cannot reach.
6. **The quantum control.** In Gross-Pitaevskii the same antiparallel pair reconnects at a fixed floor (the healing
   length) with the literature's gap law, and never returns; every regulator nature offers is a fixed length, and
   Navier-Stokes' only one, sqrt(nu/s), moves with the flow. The gap between them is a quarter of a Laplacian
   (regularity is a theorem for (-Delta)^alpha, alpha >= 5/4).

## The field on 2026-09-08

Buckmaster and Alpöge posted Lean-verified finite-time blow-up *with smooth forcing* for 3-D Euler (and Boussinesq,
IPM), on the Córdoba-Martínez-Zoroa program. OpenAI reported to them, unseen, a forced Navier-Stokes blow-up on R^3
and T^3 - Fefferman's statements (C)/(D), which admit a smooth force and are the Clay problem's breakdown direction.
Ganeshram, Duruisseaux and Anandkumar posted evidence with a proof framework for an unforced, travelling, axisymmetric
self-similar Euler singularity on R^3 (swirl-driven, in Luo-Hou's variables; profile not yet public). None of these
answers the unforced physical question, (A)/(B); a smooth force replenishes the budget that makes the unforced
collapse "free but unpaid for". This repository's forced seam race (C24) watches how; its unforced results stand as
measurements of the mechanism the fluid has on its own. Details: `docs/LOG.md`, "The field on 2026-09-08".

## The structures

Sheet, seam, twist wave, thickness arm, viscous floor, V - and the two phases: external squeeze that the budget stops
at sqrt(nu/s), and self-induced roll-up the budget cannot price. `docs/MAP.md` puts them next to each other.

## Run one thing

```
pip install numpy torch matplotlib
python burgers_entropy.py                       # 1-D: the blow-up an integrator with zero numerical entropy can see (seconds)
python monotone_vs_dominating.py                # 1-D and 2-D: the two sets, disjoint then intersecting (a minute)
NU=2e-3 N=32 T=1.2 python seam_race.py          # the seam race at 32^3 with the strong-twist column (seconds)
DMIN=0.30 TWISTW=10 N=32 ITERS=60 NVER=128 python adversarial_ic.py     # forbid the twist (CI-scale, ~10 min)
IC=found N=64 NU=2e-3 T=3 python kaggle/seam/seam_gpu.py                 # the seam race, GPU or CPU, any N
```

Every script is standalone, reads the found fields from `results/found/`, and prints its registered verdict.
`docs/SCRIPTS.md` lists all 43 by theme with their result files; `.github/workflows/experiments.yml` runs them on CI
(`which=<job>`); `results/README.md` indexes the 170 result files.

## Boundaries

Resolution: 256^3 on a T4, Re ~ 500-1500 - small by modern standards and the Euler clock expires at t ~ 0.9 even
there. Data: mostly one adversarial field, one seed; the generality runs (Kida-Pelz, a second field, a tube pair) are
one day old. The instrument measures; it does not prove. Its contribution is the map and the standard.

Issues and Discussions are open. Refutations welcome; bring the output.
