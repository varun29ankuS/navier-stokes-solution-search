# The view from above

*Rewritten at milestones, not appended to. Last rewrite: 2026-09-08, after the forced run. The chronological record is [LOG.md](LOG.md);
the registered scores are [../CLAIMS.md](../CLAIMS.md).*

## The objects

**Sheet.** What the adversary builds when asked for growth: a thin layer of intense vorticity, nilpotent locally
(A^2 = 0, the local part of the quadratic map is silent), pressed by external strain across it and narrowed by strain
along it. Dissipation concentrates on it with CKN exponent ~4; a singularity needs ~1. The budget forbids a sheet from
collapsing (constant cost per halving); it permits a line or a point (cost -> 0). CKN says the singular set has
dimension <= 1. So a sheet must become something lower-dimensional first.

**Seam.** Where two sheets meet antiparallel: vorticity reverses across a surface on which |w| passes through zero.
Squared-coherence measures read it as a smooth cylinder; the signed angle beta = 1 - xi.xi' > 1 sees it. Biot-Savart
cancels there, Kelvin-Helmholtz rolls it, viscosity reconnects it. Every hand-built singularity candidate in the
literature (Kerr, Moffatt-Kimura, Luo-Hou's odd symmetry plane) sits on one.

**Twist wave.** The reversal measured at fixed physical separations peaks in sequence from large to small - one
structure descending through the scales as the sheets are pressed together. Above the viscous thickness the descent
is inviscid (peak times agree across viscosities to one sampling step) when it is fast; a slow descent (a tube pair
over eight time units) gives viscosity time to act on the cores and the inviscid clause fails.

**Two phases.** Phase 1: external strain, fixed rate s, exponential thinning with equal time per octave - the budget
runs out at sqrt(nu/s), the classical viscous sheet, safe. Phase 2: self-induced strain s ~ Delta_u / ell, time per
octave proportional to the separation, a geometric series - the budget is indifferent because the collapse's energy
shrinks with it. The sheet field is already in phase 2 when we meet it (0.55 -> 0.25 per octave); the tube pair shows
the handover in one run (3.7, 2.7, 2.4 then 0.6-0.7).

**Thickness arm.** The analyticity strip delta(t): the sharpest scale in the flow, the sheets' own thickness. It decays
exponentially and never reaches zero on its own. It is viscous (its e-folding changes with nu) where the gap arm is not.

**The V.** Gap arm (inviscid, linear, gap = 0.51 (1.78 - t) on the sheet field) meets thickness arm (viscous,
exponential) at t = 1.65-1.74: the pair merges into one doubled structure before the gap reaches zero. At nu = 2e-3 the
merge is at sqrt(nu/s) and merging is cutting - enstrophy and max|w| peak and decay, resolved to T = 3. At 1e-3 the
merge is above sqrt(nu/s): the Euler seam race starts there, and the clock cannot see it.

**The floor.** Every regulator nature offers installs a fixed length: the healing length in a quantum fluid (cut on
contact, guaranteed, the literature's gap law), the mean free path, thermal noise. Navier-Stokes has one length,
sqrt(nu/s), and s is the flow's own. Regularity is a theorem for (-Delta)^alpha with alpha >= 5/4; viscosity is
alpha = 1. The forced route (Fefferman (C)/(D)) removes the floor's foundation instead: a smooth force replenishes the
budget. Measured (C24, v6): a *steady* smooth large-scale force does not carry the seam through the floor - the twist
turns over and the seam is cut exactly as unforced - but growth continues by re-supply, the pumped large scales
rebuilding seams (max|w| 9.9x and rising at T = 3, against 5.7x peak-and-decay). A driven flow, not a collapse. The
forced proofs use forces that track the collapse; a force that merely pushes gives turbulence.

## The patterns, across the whole

1. **Every candidate fails the same way, and the way tracks the theorem.** Burgers: monotone and dominating sets
   disjoint, and it shocks. 2-D: they intersect, and it is regular. 3-D: every quantity built - energy-class,
   topological, learned, memory, twist - landed in the dominating set and never in the monotone one. The table has
   been drawing the regularity theory's boundary from outside the whole time.
2. **Three fluids, one shape: a self-similar descent meets a floor.** Quantum: fixed floor, cut on contact. NS: a
   floor that moves with the flow, but present at Re ~ 500. Euler: no floor. The GPE control, the seam race and the
   descent law are one experiment at three floor heights.
3. **The descent in scale is a descent in dimension.** The twist wave is a sheet pair becoming a tube, 2 -> 1; CKN
   allows only <= 1 to be singular; the budget forbids 2 and frees 1. Scale and dimension descend together.
4. **The rhythm is in the collapse's own clock.** Successive half-octave times shrink by ~0.76 against 1/sqrt(2) for a
   rate proportional to 1/ell: even octaves in -log(T* - t), acceleration in laboratory time. The constant that is
   really a function: the period is constant in the structure's time and a function in ours.
5. **Prediction errors have a direction.** Across the week's misses (attack size, placement, eye list, twist halving,
   Re_seam): systems are more capable, and physical thresholds less clean, than predicted from the mechanism in view.
   Registered predictions are now checked against this bias before they are written.

## Closed

Energy-class functionals (Tao's wall, in pictures: the dyadic model). Topology (the antiparallel counterexample).
Helicity as a slope. Liouville as a discriminator. The local Lyapunov class, with a passed 2-D control. The minimal
blow-up datum. Incommensurability. The thinnest-squeeze inequality. Huisken-type entropies. The discretely
self-similar orbit (Chae; Chae-Wolf). The averaged jitter bound (Mahithitarmmatorn 2026). The twist as the monotone
quantity (rises under every viscosity). Sheet-only concentration as a proven statement (Grujić's sparseness is the
rigorous cousin; the converse is open).

## Open, and reachable

- **The Re ladder for the V and for Re_seam at the turnover.** 256^3 resolves the twist peak at nu = 2e-3 only; 384^3
  reaches 1e-3 (a T4 at its limit); 5e-4 needs ~640^3.
- **The Lagrangian descent.** Kelvin's frame: tag the two sheets as material surfaces, measure gap and circulation
  along paths; the rate is then Helmholtz's law with circulation fixed, and "self-induced" is tested directly.
- **The forced seam race with a tracking force.** C24 showed a steady force drives rather than collapses; the next
  question is the smallest *time-dependent* force that carries the seam through sqrt(nu/s) - the shape of the force
  the proofs need, measured.
- **The travelling axisymmetric profile.** When public: evolve it in this solver, perturbed and periodised, with the
  clock; a stable singularity should show a finite-T* descent in the coherent regime.
- **Generality with more seeds.** Everything above is one adversarial field, one second field, one pair, Kida-Pelz.

## Open, and not reachable here

Whether the wave arrives at zero in Euler (the clock expires at t ~ 0.9 even at 256^3; the strip collapses
exponentially). Whether the cut keeps winning as Re -> infinity. Any theorem.
