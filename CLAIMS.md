# Claims, how to refute them, and where we are weakest

Every claim below is registered with the observation that would refute it and the command that produces the number.
If you find one of these fails on your machine, open an issue with the output. If you think a claim is weaker than
stated, say so; the weak points we already know about are listed at the bottom.

## Claims

**C1. The scheme's numerical entropy production is zero to round-off.**
Refuted by: energy drift for inviscid Burgers or 3-D Euler larger than ~1e-6 relative with the default settings.
Test: `python burgers_entropy.py` (expect sigma_num ~ 1e-15) and `DIM=3 N=48 NU=0 T=4 python ns_d.py` (expect E/E0 = 1.000000).

**C2. Every balance law closes term by term, including the unsigned production terms.**
Refuted by: a relative residual above 1e-6 for any budget except the 2-D palinstrophy budget (3-9e-6, RK4 error on
terms of size 1e3; halve dt to clear it).
Test: `python budgets.py`; `DIM=1|2|3 python ns_d.py`.

**C3. The 1-D Burgers gradient follows the exact 1/(1-t) until the grid limit.**
Refuted by: max|u_x| departing from 1/(1-t) by more than the resolution shortfall before t = 0.9 at N = 512.
Test: `python burgers_entropy.py`; figure `figures/burgers_blowup.png` from `python plots.py`.

**C4. 2-D viscous Taylor-Green is reproduced to machine precision.**
Refuted by: L2 error against e^{-2 nu t} sin x cos y above 1e-12 at t = 2.
Test: `python taylor_green.py` (part A) or `DIM=2 N=128 NU=0.02 T=2 DT=2e-3 IC=tg python ns_d.py`.

**C4b. Three further closed forms are reproduced: Cole-Hopf viscous Burgers to 6e-12 before the shock and 2e-6 through
it (grid floor); a random single-shell 2-D field to 3e-14; the viscous ABC flow in 3-D to 7e-15 over t = 8.**
Refuted by: errors materially above these at the stated grids and time steps.
Test: `python exact_solutions.py`.

**C4c. Time integration is converged: halving every step changes the state by ~3e-5 (Kida-Pelz 32^3, t = 1); the
printed budget residual is the centred-difference diagnostic's own second-order error, not the integrator's.**
Refuted by: a state difference between DT_DIV=1 and DT_DIV=2 runs above 1e-3 relative.
Test: `IC=kp N=32 T=1 python criteria3d.py` and the same with `DT_DIV=2`.

**C5. In 3-D Euler (Taylor-Green) the vortex-stretching term converges upward with resolution: 1.50, 1.80, 2.03 at
t = 4 for 64^3, 96^3, 128^3, and the normalised rate S/Z^{3/2} rises with resolution.**
Refuted by: a higher-resolution run (192^3 or above) giving S(t=4) below 2.03, or a Z(t=4) inconsistent with the
published Taylor-Green curves (Brachet et al. 1983 and later).
Test: `N=128 T=4 python budgets3d.py` (~40 min on 2 cores); results in `results/budgets3d_*.txt`.

**C6. The vorticity aligns preferentially with the intermediate strain eigenvector (Ashurst et al. 1987).**
Refuted by: the intermediate-eigenvector fraction not exceeding the other two at 64^3 by t = 3.
Test: `IC=tg N=64 T=4 python mechanism3d.py`. Status: met, weakly (0.341 vs 0.328/0.331 at t = 3; 0.356 at t = 2); note that by the analyticity clock the 64^3 Taylor-Green run is past its reliable window at t = 3, so the t = 2 number is the one that counts. Neither
flow is developed turbulence by t = 4; the strong textbook signal is expected in the forced run and is not yet shown.

**C7. Helicity is conserved under 3-D Euler dynamics to the same order as energy.**
Refuted by: helicity drift exceeding energy drift by more than a factor of ~10 for the perturbed ABC flow.
Test: `IC=abc N=64 T=4 python mechanism3d.py`. Status: pass (helicity 2.979345 through t = 3; drift at t = 4 within 2.5x of
energy drift, as the grid is reached).

**C8. A learned coarse simulator with frozen exact transport holds energy within a few percent over 2000 steps where
a fully learned one drifts by tens of percent.**
Refuted by: a fully learned baseline of comparable size that keeps energy within the same band at the same
resolution and Reynolds number.
Test: `python kolmogorov_closure.py` (v1, ~40 min); `python kolmogorov_v2.py` (v2).

**C9 (negative). A non-negative eddy-viscosity closure learns nothing in 2-D at 32^2 and Re ~ 1250.**
Refuted by: the same closure producing a measurably different rollout from no closure.
Test: `python kolmogorov_closure.py`.

**C10. The skew-projected learned correction is energy-neutral to first order in the step (not exactly).**
Refuted by: a per-step relative energy change from the projected correction comparable to the unprojected one
(3.5e-4), or exactly zero (which would mean the second-order argument is wrong).
Test: `python verify_skew_closure.py` (expect ~1e-5 vs 3.5e-4).

**C11. The Kida-Pelz early evolution is converged, and its vorticity direction roughens at a fixed physical scale.**
At t = 0.5 (inside the reliability window at every resolution) 64/96/128^3 agree on Z = 5.4197 and S = 5.817 to four
digits, and the Constantin-Fefferman direction coherence at fixed physical scale h = 2pi/32 is 0.128 / 0.127 / 0.124
with a local Lipschitz exponent 1.49 / 1.51 / 1.53 (2 would be a smooth direction field). Refuted by: any of these
moving with resolution beyond the third digit, or the exponent rising toward 2 at higher N.
Test: `IC=kp N=64 T=1 python criteria3d.py` (and 96, 128). Earlier t = 3 Kida-Pelz numbers on this page were
withdrawn: the analyticity strip is below 2 dx by t = 1.0 at all three resolutions, so nothing later is a statement
about Euler.

**C12. The adversarial searcher finds smooth low-k initial data that the verifier rejects as unresolved.**
Twice, a |k| <= 4 field found by gradient ascent on enstrophy amplification over t = 1 cascaded to the grid cutoff
(8.75x on the 32^3 search grid / 20x at 64^3 with delta(T) < 0; 4.96x / 8.0x at 96^3 with delta 0.06 < 0.13), and
twice the registered verdict is FAIL because the verifier's clock says the number is not a statement about Euler.
That is the design working; it is not a blow-up. Refuted by: a searcher result that stays resolved (delta > 2 dx on
the verification grid) and beats the best classical flow. Test: `python adversarial_ic.py` and the `leashed` CI job.

**C13. Attainable enstrophy growth is flat in helicity up to half the maximum, then collapses.** With the helicity
held fixed by a hard constraint on the 32^3 truncated system, the maximal amplification over t = 1 at fixed Z0 is
4.54, 4.68, 4.47, 3.11, 1.27 at relative helicity 0, 0.25, 0.5, 0.75, 0.9 (`results/helicity_proj_*.txt`). The
direction is Moffatt's; the shape (a threshold, not a slope) is the claim. The 0.9 field is resolved at 96^3 and
beats Taylor-Green (1.266 vs 1.112); the others cascade to the cutoff and are outside the window. Refuted by: a
64^3 search grid giving a different ordering, or a helical field at 0.75-0.9 that matches the helicity-free growth.
(The first version of this claim, from a penalty method with a sqrt 2 error in the bound, read 1.21 / 1.10 / 1.01
at 0 / 0.71 / 0.99; the penalty was failing to explore. Both are kept in the README table.)

**C14. The strip diagnostic separates a proven singularity from Kida-Pelz by the trend of its decay rate.** On the
inviscid dyadic model (blow-up is a theorem) the local decay rate -d log(delta)/dt rises 1.8e5-fold inside the
reliable window and delta reaches zero at t* = 0.5585 with energy conserved to 3e-8; on Kida-Pelz Euler at
64/96/128^3 it falls by half. Refuted by: a resolution at which the Kida-Pelz decay rate rises, or a choice of the
fit range that makes the dyadic rate flat. Test: `python dyadic.py`, `IC=kp N=128 python strip_tracker.py`.

**C15 (negative, about the METHOD - re-scoped 2026-09-07). The Lyapunov learner, as built, cannot find a monotone
functional even where one is a theorem: in 2-D with log(Z/P) supplied as a feature it fails its positive control
(attack +0.45 / +0.41, `results/lyapunov2d_64*.txt`). Its 3-D failures therefore say nothing about 3-D. Original
statement, withdrawn as a claim about the equation: Within the class of bounded local functionals M = Z exp(Phi), Phi a learned enstrophy-weighted
average of pointwise vorticity/strain features, an adversary finds a violating trajectory every round.** Three rounds
at 24^3: adversary violations +0.40, +2.40, +0.71 (relative dM/dt) with no closing trend; held-out classical flows
never violate. An eight-round 32^3 series that fell to 0.03-0.09 was re-attacked with a stronger adversary and broken at +0.40 to
+0.55 (`results/lyapunov_32_p0_attack.txt`): fixed-budget adversaries understate violations, so every future
candidate is judged by `ATTACK` mode, not by its own training adversary. A twenty-round candidate broke at +1.03 under a 5 x 60 attack. Refuted by: a candidate that survives an
ATTACK run with violations below 1e-3 (which would yield a candidate inequality, not a theorem).

**C16. The fastest-growing field found concentrates dissipation like sheets, not like a singularity.** CKN exponent
alpha = 4.5-4.6 (3.6-3.8 at the smallest resolved pair) at 64^3, flat over the window while Z x 2.8; Taylor-Green
5.4 (4.4). Refuted by: alpha falling toward 1 with time or with resolution for any field found by the searcher.
Test: `IC=found N=64 NU=1e-3 python ckn_exponent.py`.

**C17. Fast enstrophy growth requires the global part of the pressure, and recruits it.** With the traceless share of
the pressure Hessian on the high-vorticity set capped at T, the maximal amplification over a turnover falls from
5.6x (cap 0.5) through 2.1x (0.4) to 1.1x (0.35) at 48^3, with the same knee and collapse at 32^3 and at T = 0.5 and
1.5; the searcher cannot push the share below ~0.35-0.42; and the global Hessian's projection along the vorticity is
negative (aids the stretching) in all eleven runs. Refuted by: a fast amplifier (> 3x) found with share < 0.35, or
one whose global Hessian opposes the stretching. Test: `OBJ=quiet SMAX=0.35 N=48 python adversarial_ic.py`.

**C18 (registered 2026-09-07, before the CI run). The Gross-Pitaevskii seam reconnects, with the known gap law, and
releases sound.** 2-D: two counter-propagating dipoles at d = 6 in a 64^2 box travel within 20% of the point-vortex
speed 1/d with total energy conserved to 1e-5 (local pilot: +6%, 1e-6). 3-D: two antiparallel pairs of lines bowed to
a closest gap of 4 xi reconnect at a finite t_r; the gap closes as (t_r - t)^p with p = 0.5 +/- 0.15; the compressible
(sound) energy rises across t_r while the incompressible energy falls; the lines do not re-approach afterwards.
Refuted by: no reconnection by T = 60, or p outside [0.25, 0.8], or sound energy not rising across the cut. Test:
`MODE=3d N=64 L=32 D=6 A=1 T=60 python gpe_seam.py`. This is a control for the viscous seam race, not a statement
about Navier-Stokes. **Result (2026-09-07, 64^3 and 96^3 agree):** t_r = 15.8 / 15.5; p = 0.65 / 0.55 on gap in
[1, 2.5] (0.39 / 0.37 on the full window, which includes the Crow phase); E_kc +32% / +36% across the cut, 2.4x by
t = 36; no re-approach. Failed as written: E_ki rises 3% across t_r +/- 3 and falls only over the release window
after it. Stands on the parts that matter; the "E_ki falls at the cut" clause was wrong and is retracted.

**C19 (registered 2026-09-07, before the run). Growth needs the twist.** The searcher's fastest fields are ~90%
antiparallel on their high-vorticity set (the seam). Forbid the twist - penalise the enstrophy-weighted antiparallel
measure at T (`TWISTW`) - and the attainable amplification at matched leash (DMIN 0.30) and iterations falls below
HALF of the twisted baseline (32^3: 8.28 -> < 4.1; 64^3: 3.18 -> < 1.6) for every weight at which the final anti
fraction on the high set is below 0.2. Refuted by: a field with anti fraction < 0.1 reaching >= 0.8 of the baseline.
Why it matters: Constantin-Fefferman covers coherent direction, Grujic's cancellation covers incoherent direction; if
fast growth lives only in the twisted regime, the seam is the only place left. Test: `DMIN=0.30 TWISTW=10 N=32 T=1.0
ITERS=60 NVER=128 python adversarial_ic.py`.
**Result (2026-09-08):** 32^3 clauses met (1.19-1.49 vs 8.28) but the 32^3 baseline is unresolved growth; the 64^3
clause is **refuted** - coherent 2.24x resolved at 128^3 vs twisted 3.18x, a 30% cut. Re-scoped: the twist adds ~40%
at 64^3 and is not necessary. The measure carries the grid scale (1-2 cells) and must be redefined at a fixed physical
separation before any cross-resolution comparison. Kept: the monotone trade-off at each resolution, the helicity
plateau, and three resolved found fields (1.19, 1.52, 2.24).

**C20 (registered 2026-09-07, before the run). The seam reaches its viscous thickness and the twist turns over.**
GPU spectral run (`kaggle/seam/seam_gpu.py`), 256^3, clock 2dx = 0.049, the searcher's sheet field and a Kerr-type
antiparallel tube pair, T = 3. Race variable = sheet thickness |w|/|grad|w|| (median on the high set) over the viscous
thickness sqrt(nu/s), s the mean compression across the sheet. At nu = 2e-3 the race variable reaches <= 1.5 inside
the clock, the strong twist (enstrophy-weighted sharp reversal) peaks there and falls by >= 20% before the clock, and
max|w| grows < 3x. At nu = 0 the twist rises to the clock. Refuted by: the twist still rising at the clock with the
race variable < 1 (the roll-up outruns the cut at a resolved viscous scale) - which would be the first resolved
evidence for phase 2 winning, and the most important number this repository could produce.
**Result (2026-09-08):** nu = 2e-3, 256^3, resolved to T = 3: twist peaks at t = 1.7 and falls 74% (met); race
variable 1.49 at the end (met, barely); max|w| growth 4.55x at the clock, peak 5.7x (**not** met, registered < 3x);
the seam is cut from the moment it forms (50-80%); Z and max|w| peak and decay inside the clock. Formally between,
physically the cut wins at Re ~ 500. nu = 0: the clock expires at t = 0.9. Tube-pair control: setup failed, void.

**C21 (registered 2026-09-08, before the run). The twist is a feedback variable, and its turning point is a fixed
seam Reynolds number.** With tau the strong twist (now at a fixed physical separation, 0.1), ell the sheet thickness
and s the compression across it, the run data are consistent with dtau/dt ~ (s - nu/ell^2) tau: pressing (inviscid,
partly self-induced) against reconnection. The gain changes sign at a fixed value of Re_seam = |w| ell^2 / nu, so the
twist peaks at the SAME Re_seam (within x2) at nu = 2e-3, 1e-3, 5e-4 on the searcher's sheet field at 256^3, while the
peak arrives later and higher as nu falls. Refuted by: Re_seam at the peak rising by more than x2 per halving of nu
(the seam needing ever more Reynolds number to be cut - the velocity jump growing with the collapse, the direction
of a singularity). Also re-tests C20's max|w| clause at lower nu and rebuilds the tube-pair control (D 0.7, sigma
0.22, T 6). Test: `kaggle/seam/seam_gpu.py` (schedule inside).
**Result (2026-09-08, 320^3): the KILL clause fired.** The nu = 1e-3 twist peak is resolved at 320^3 (0.182 at t = 1.60,
clock 1.75). Re_seam at the turnover: 224 (2e-3) -> 594 (1e-3), x2.65 per halving against the registered x2. The seam
at the merge carries a larger velocity jump at lower viscosity (Delta_u 2.2 -> 3.3): it needs, and gets, more Reynolds
number before it is cut. Two rungs, one field; the 5e-4 rung needs ~450^3. This is the first registered number in the
repository pointing in the direction of a singularity, and it is recorded as such. The 320^3 run at 2e-3 reproduces
the 256^3 descent to three digits.
**Re-read after C26 (v8):** the direct Lagrangian velocity jump is the same at both viscosities (max x1.08-1.09 over
seeding, ~0.95 at the merge). The x2.65 was the estimator: thinner sheets at the same jump raise peak |w| and the
proxy ell does not shrink in step. The KILL stands as a fired criterion on Re_seam; its interpretation is withdrawn.
**Result (2026-09-08): inconclusive.** Only nu = 2e-3 resolves the twist peak (Re_seam 224, t = 1.6); at 1e-3 and
5e-4 the 256^3 clock expires (t = 1.6, 1.1) before the peak, so the ladder cannot be compared. The prediction
"order 10" was wrong by 20x. New, resolved at 2e-3: the twist is a wave in scale (peaks at 0.4, 0.2, 0.1, 0.05 in
sequence) arrested at the viscous thickness. The pair control's KILL clause fired (twist rising, race 0.56) with Z
flat and Re_seam falling; judged post hoc to be a defect of the clause (core thickness vs seam gap) - recorded as
fired, clause corrected for the next run, judgement flagged as post hoc.

**C22 (registered 2026-09-08, before the run). The descent law.** The twist measured at fixed separations peaks in
sequence from large to small separation - one wave descending in scale, arrested at the viscous thickness. Above
that thickness the descent is nu-independent (peak times at separations 0.56 ... 0.1 agree across nu = 2e-3, 1e-3,
5e-4 within 0.1) and it accelerates: the time per octave of descent shrinks between 0.4 and 0.1 (self-induced
strain, phase 2). Refuted by: constant or growing time per octave (external strain only - the wave never reaches
zero), or peak times that move with nu above ell_nu. Why it matters: the descent above ell_nu is an inviscid fact
measured resolved at finite nu; finite-time arrival at zero is the Euler question. Test: `kaggle/seam/seam_gpu.py`
with the eight-separation ladder (schedule inside).
**Result (2026-09-08): PASS.** Resolved peak times agree across nu = 2e-3, 1e-3 (and 5e-4 where resolved) within
0.05; time per octave falls from ~0.55 to ~0.25 between 0.56 and 0.14 at both rungs; gap = 0.51 (1.78 - t). The
V (gap arm vs thickness arm) meets at t = 1.65 / 1.70 / 1.74, scale 0.068 / 0.043 / 0.022; at 1e-3 the merge is
above sqrt(nu/s). Registered predictions for the V (time nearly nu-independent, scale falling, merge above the
viscous scale at lower nu) met at 1e-3; the 5e-4 row is extrapolated and not counted.

**C23 (registered 2026-09-08, before the run). The descent law is a property of antiparallel sheets, not of one
field.** Same eight-separation ladder, 256^3, nu = 2e-3 and 1e-3, three initial conditions. (a) Kida-Pelz (coherent,
anti = 0): twist <= 0.01 at every separation through its clock - no wave; the instrument is blind to symmetric
focusing, as it should be. (b) The concentration-rewarded adversarial field (`results/found/ckn64.npz`): a
descending wave, peak times nu-independent within 0.05, accelerating octave time, its own T*. (c) A Kerr-type
antiparallel tube pair (D 0.8, sigma 0.2, bow 0.3, T = 6): a descending wave, nu-independent, with CONSTANT time per
octave early (mutual induction at fixed circulation: exponential approach) and acceleration only in the last octaves.
Each run reports its V (gap arm, thickness arm, crossing, sqrt(nu/s) there). Refuted by: a wave on Kida-Pelz; peak
times that move with nu on (b) or (c); or a pair whose octave time never shrinks (no self-induced phase at all).
**Result (2026-09-08, v4):** (a) holds qualitatively (no descent), numeric clause violated at 0.4/0.56 (pattern
growth in place). (b) holds: nu-independent within 0.05, accelerating (0.9 -> 0.1 per octave in one step), V at
t ~ 1.0, scale ~ 0.1, above the viscous scale at 1e-3. (c) run with a touching pair (audit F); at 1e-3 it nonetheless
shows slow-then-fast octaves (2.3, 0.97, 0.62, ~0.5); the registered approach test is v5.
**(c) result (v5, apex gap 2.5 sigma, T = 8):** descending wave - yes; constant-then-shrinking octave time - yes (3.7,
2.7, 2.4 then 0.6-0.7 at 1e-3); nu-independence - **refuted** (peaks shift 0.5-1.4 between 2e-3 and 1e-3; the cores
diffuse over the slow approach). The two-phase structure holds; the inviscid clause holds only for fast descents.

**C24 (registered 2026-09-08, before the run). The forced route, watched.** Fefferman's (C)/(D) allow a smooth force.
`FORCE=eps` adds f = eps x (the initial field restricted to |k| <= 4): smooth, periodic, divergence-free, steady - the
external squeeze held on forever, the budget replenished from outside. With eps = 0.5 and 1.5 at nu = 2e-3 and 1e-3,
256^3: the V bottom is no longer the end of the descent - twist@0.05 keeps rising past the unforced turnover
(t ~ 1.7), max|w| accelerates instead of peaking, and the wave passes the viscous thickness (race variable < 1 with
the twist still rising) inside the clock. Refuted by: a forced run whose twist still turns over and whose max|w|
still peaks inside the clock - which would mean a steady large-scale force does not by itself defeat the cut, and
the forced blow-ups in the literature need a force that tracks the collapse.
**Result (2026-09-08, v6, resolved at nu = 2e-3, eps = 0.5, T = 3): REFUTED on the mechanism.** The twist still
turns over (peak 0.110 at t = 1.4, fall 59%), the race variable stays >= 1.18, the cut fraction 0.3-0.5 - the seam is
cut as unforced. Growth continues by re-supply instead: max|w| 9.9x at the clock and rising (unforced: 5.7x peak and
decay), Z 7.2, E/E0 2.5, a second seam at t ~ 2.6. A steady force drives; it does not carry the seam through the floor.
eps = 1.5 and nu = 1e-3 leave the clock before t = 1.3 and are not counted.

**C26 (registered 2026-09-08, before the run). Kelvin's frame: the Lagrangian gap closes linearly and the velocity
jump stays bounded.** Particles seeded on the high set at t = 1.0 and split into the two sheets by the sign of
omega . xi_ref, advected with the flow (`LAGR=1`). On the sheet field at nu = 2e-3 (256^3) and 1e-3 (320^3): the
median distance from one sheet's particles to the other's closes with lambda = 1 +- 0.2 on the resolved window, the
median velocity jump across the nearest pairs stays within 2x of its value at seeding until the merge, and the material
|w| grows while the gap closes and turns over at the merge. Refuted by: a jump that grows more than 2x while the gap
closes (C25 fails on this field along one solution), or lambda < 0.7. Context: the Eulerian rungs (C21) already show
the jump at the merge rising x1.5 per halving of nu across solutions.
**Result (2026-09-08, v8): PASS at both viscosities.** Lagrangian gap 0.27 -> 0.15 with lambda = 0.99 at nu = 2e-3
(256^3) and 0.99 at 1e-3 (320^3); the gap minimum (t = 1.65 / 1.70) coincides with the Eulerian V (1.69 / 1.72); the
velocity jump across the pair is 1.01 -> max 1.10 (x1.09) and 1.08 -> max 1.17 (x1.08), ~0.95 at the merge at both nu
- bounded along each solution and viscosity-independent. The material |w| turns over (t ~ 1.25-1.3) and falls: the cut
acting on the tagged sheets. Consequence for C21: its KILL fired on the Eulerian estimator (|w| ell^2 / nu on the
twisted set), whose growth across nu comes from thinner sheets at the same jump (peak |w| 45 vs 29), not from a growing
jump; the interpretation "the velocity jump grows with the collapse" is refuted by the direct measurement. C25 holds
on this field along two solutions.

**C27 (registered 2026-09-08, before the run). The tracking force carries the seam through the floor.** `FMODE=track`:
every output step the force is set to eps x P[u_H], u_H the velocity induced by the smoothly-masked high-vorticity set
alone - the pair's own self-induction amplified, phase 2 fed directly, nothing else; smooth in x, piecewise-steady in
t. At nu = 2e-3, 256^3, eps = 1 and 3: inside the clock the twist@0.05 keeps rising past the unforced turnover
(t ~ 1.7), max|w| accelerates over the last quarter of the window, and the race variable drops below 1 with the twist
still rising - the seam passes sqrt(nu/s), which the steady force (C24) never achieved. Refuted by: the twist turning
over as unforced under the tracking force too, which would say amplifying the pair's self-induction is still not the
force the forced proofs use.
**Result (2026-09-08, v9, 320^3, nu = 2e-3): between.** The tracking force speeds the collapse without changing its
shape - gap arm 0.51 (1.78 - t) -> 0.72 (1.42 - t) -> 0.88 (1.20 - t) for eps 0 / 0.5 / 1, the twist@0.1 peak the same
0.11 in all three, arriving earlier; max|w| 4.6x -> 8.7x -> 15.3x at the clock. The seam at 0.1 still turns over
(fall 11%, 14%) and at eps = 0.5 the merge scale is on sqrt(nu/s) (0.050 vs 0.062): amplified self-induction does not
carry that seam through the floor. The twist@0.05 keeps rising at the clock in both runs (0.033, 0.040 - 2-3x the
unforced peak): the one clause pointing "through", undecided because the clock expires on it (t = 1.40, 1.60). max|w|
is not accelerating in the last quarter in either run. At 256^3 the runs leave the clock before t = 1 and say nothing.

**C28 (registered 2026-09-08, before the run). The seam flips again below the cut.** Reconnection leaves threads that
are antiparallel to each other at a smaller scale (Hussain's bridging; the Yao-Hussain 2020 reconnection cascade). In
the Lagrangian mode, FLIP = the fraction of tagged particles whose sign of omega . xi_ref has reversed since seeding.
At nu = 2e-3, 320^3, seeded at t = 1.0: the flip fraction rises from ~0 to 0.10-0.30 across the merge (t ~ 1.7), and
the material gap, after re-opening, closes a second time before the clock (2.4) - a second, smaller seam. Refuted by:
a flip fraction below 0.05 at the clock (a clean cut), or no second closing inside the clock. Why it matters: a seam
that flips and is cut at every level is the turbulence cascade as a staircase of reconnections (the safe outcome);
a level that flips and is *not* cut would have to arrive with lambda <= 1/2 (Theorem 7).
**Result (2026-09-08, v10, 320^3, nu = 2e-3): between.** Flip fraction 0.000 through the merge (t = 1.55), then
rising steadily as the material vorticity falls and the jump collapses: 0.01 (1.8), 0.04 (2.0), 0.09 (2.2), 0.127 at
the clock (2.4), ~0.2 per time unit and still climbing - the cut is a sign reversal on material fluid, one particle in
eight by the clock (first clause met). The material gap re-opens monotonically after the merge (0.156 -> 0.188) and
does not close again inside the clock: no second seam among the tagged fluid at this Reynolds number (second clause
fails). C26 passes a third time (lambda 0.99, jump x1.04). The staircase's first step exists; its exponent lambda_stair
= log(r_s)/log(r_t) needs a second step, beyond ~400^3.

**C29 (registered 2026-09-08, before the run). The sheets close under their own induction.** In the Lagrangian mode,
for every nearest A-B pair, the closing rate along the separation from the full velocity, (u_A - u_B) . d_hat, and from
the pair's own Biot-Savart field u_H (the velocity induced by the smoothly-masked high set alone). On the descent
(t = 1.0 to the merge) at nu = 2e-3 and 1e-3, 320^3, the self-induced rate is >= 70% of the measured closing rate: the
gap law g = 0.51 (T* - t) is Helmholtz - the pair pressing itself - not a fit and not external strain. Refuted by: a
self-induced share below 40%. Also reported: where the flipped fluid sits after the merge (distance to the other sheet
over the gap: ~0 bridges, ~1 threads).

**Not claimed.** Anything about the regularity of 3-D Navier-Stokes. `THEORY.md`, Theorem 4, records why the
structure used here cannot decide it.

## Where we are weakest (poke here first)

1. **One seed everywhere.** Random initial conditions (2-D decaying, ABC perturbation, Kolmogorov) use seed 0. No
   error bars.
2. **The upwind baseline is the weakest reasonable comparator.** First-order upwind is what textbooks use to show
   numerical dissipation; a high-order WENO or a standard energy-conserving finite-volume scheme would be the fair
   opponent. The point being made (dissipation hides growth) survives, but the size of the gap in the figures is
   against an easy target.
3. **The learned baseline in the Kolmogorov experiments is small** (a 3-layer CNN, 16-step unroll in v2). Stronger
   learned simulators exist (Kochkov et al. 2021, FNO variants). "Learned drifts" only counts against a strong one.
4. **Resolution.** 128^3 is small by modern standards; the t = 4 Taylor-Green numbers are a lower bound on a curve
   still rising with resolution, as stated. Nothing here approaches the resolutions of published singularity searches.
5. ~~RK4 is not exactly conservative.~~ Closed: `midpoint.py` (implicit midpoint on the skew part) brings the 3-D drift
   from 3.5e-8 to 5.7e-12. The RK4 numbers elsewhere stand as reported; the midpoint integrator is available.
6. **The 2/3 rule and the skew form.** Energy conservation on the dealiased modes is standard; we have not separately
   verified the aliasing error bound at the highest retained modes for the 3-D runs.
7. **Comparison with the literature is qualitative.** "In range" of Brachet's curves is not a digit-for-digit match;
   a digit-for-digit comparison at matched resolution and time step is the obvious next check.
8. **Kida-Pelz is not resolved here.** The flow has octahedral symmetry and the literature runs it on 1/64 of the box
   (effective thousands cubed). Our full-box 128^3 is their 32^3; the analyticity-strip clock in `criteria3d.py`
   says where each run stops being reliable. Symmetry reduction is the next real step, not a bigger box.
9. **The theory note's Proposition 3 (memory) is about a different system** than the fluid experiments; it is in the
   note because the same inequality governs both, not because the fluid runs test it.
