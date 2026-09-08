# The seam, measured: where Kelvin hands over to Stokes

*Varun Sharma — draft note, 2026-09-08. Repository: github.com/varun29ankuS/navier-stokes-solution-search, release
v2026.09.08. Every claim below was registered before its run and scored afterwards (`CLAIMS.md`); numbers are
reported only inside the analyticity-strip clock. A search, not a solution.*

## 1. What was measured

An adversarial search for the fastest enstrophy growth in 3-D Navier-Stokes on the periodic box produces, every
time, a pair of antiparallel vortex sheets pressed together - the configuration of Kerr (1993) and Moffatt-Kimura,
found rather than built. We call the surface between them, across which the vorticity reverses and $|\omega|$
passes through zero, the *seam*. An instrument was built to watch it: a pseudo-spectral solver with zero numerical
dissipation (skew transport, exact viscous integrating factor, energy 1.000000 at $\nu = 0$), a GPU port certified
against it to three digits, and a clock - the analyticity-strip width $\delta(t)$ from the spectrum's tail; a number
is reported only while $\delta > 2\,dx$.

The seam is made visible by one measure. With $\xi = \omega/|\omega|$, the signed angle $\beta = 1 - \xi\cdot\xi'$
between a point and a neighbour at a fixed physical separation is $> 1$ across a reversal; the **strong twist** at
separation $\sigma$ is the enstrophy-weighted mean of $(\beta - 1)_+^2$ over neighbours at that separation. Measured
at $\sigma = 0.56, 0.40, 0.28, 0.20, 0.14, 0.10, 0.07, 0.05$ on the fastest field at $256^3$ and $320^3$, for
$\nu = 2\times10^{-3}, 10^{-3}, 5\times10^{-4}$.

## 2. Five numbers

1. **The reversal descends through the scales as a wave.** The twist at each separation rises and falls in turn,
   largest separation first; the peak times agree across the three viscosities to one sampling step (0.05) at every
   separation resolved. The descent is inviscid above the viscous thickness. (C22, passes.)
2. **The gap relaxes exponentially onto the sheets' thickness, driven by their own induction, which cancels as they
   meet.** In Kelvin's frame $dg/dt = -k\,(g - \delta)$ with $k = 1.7 \pm 0.1$ at $
u = 2	imes10^{-3}, 10^{-3},
   5	imes10^{-4}$, constant to 3% within a run (C32); the sheets supply 84% of the closing and the far field 16% at
   all three viscosities; the gap never reaches zero. A tube pair has its own $k pprox 0.16$ and closes against the
   far field. The straight line $g = 0.51\,(1.78 - t)$ was its tangent over a
   short window, and the shrinking time per octave ($\approx 0.55$ to $\approx 0.25$) is the relaxation on a log axis. In Kelvin's frame the sheets' own Biot-Savart field supplies 84%
   of the closing early and goes to zero at the merge - antiparallel sheets' fields cancel as they approach - with
   the far field (16%, rising to 29%) carrying the end. The curvature law (Da Rios) is refuted: $\kappa$ rises
   $0.75 \to 3.6$ while the induced closing falls. A hand-built tube pair shows a mild late speed-up (3.7, 2.7, 2.4
   then 0.6-0.7 per octave). (C22, C23, C29-C31; an earlier reading of self-accelerating induction is withdrawn.)
3. **The pair merges on the viscous scale.** The gap arm (inviscid, linear) meets the sheets' thickness arm
   ($\delta(t)$, viscous, exponential) at $t = 1.69$ and $1.72$; the meeting scale is $0.052$ vs $\sqrt{\nu/s} =
   0.051$ at $\nu = 2\times10^{-3}$ and $0.031$ vs $0.028$ at $10^{-3}$ ($320^3$). The merge is the cut: enstrophy
   and $\max|\omega|$ peak and decay inside the resolved window. (C20; an earlier reading "above the floor at
   $10^{-3}$" was a biased fit and is retracted.)
4. **In Kelvin's frame the jump is bounded.** Fluid tagged on both sheets at $t = 1$ and followed: the material gap
   closes with $\lambda = 0.99$ at both viscosities, its minimum coinciding with the Eulerian merge to 0.04; the
   velocity jump across the pair rises by at most 9% before the merge, the same at both $\nu$; the material vorticity,
   however, never follows the inviscid budget: the strains overshoot its growth from the first row by about half the
   compression rate, with the sign and viscosity-scaling of diffusion across the sheets (C34) - Kelvin's frame holds
   for circulation, not for vorticity. (C26, passes. An Eulerian estimator had
   suggested the jump grows as $\nu$ falls; the direct measurement says it does not, and that reading is withdrawn.)
5. **A steady smooth force does not carry the seam through the floor.** With $f = \varepsilon\times$(the large scales)
   held on - the smooth force Fefferman's breakdown statements (C)/(D) admit - the twist turns over and the seam is
   cut exactly as unforced; growth continues only by re-supply, the pumped large scales rebuilding seams. A driven
   flow, not a collapse. A force that *tracks* the collapse (the pair's own induction, amplified) speeds the descent without
   changing its shape - gap slope 0.51 → 0.72 → 0.88, the same twist peak arriving earlier, $\max|\omega|$ 4.6× →
   8.7× → 15× - yet the seam at separation 0.1 still turns over on the viscous scale; a thinner reversal at 0.05 is
   still rising when the clock expires, undecided at $320^3$. (C24 refuted as registered; C27 between.)

Controls: Kida-Pelz (coherent, symmetric focusing) shows no wave; a Gross-Pitaevskii pair (a fluid with a fixed
floor, the healing length) reconnects on contact with the literature's gap law and never returns.

## 3. One theorem, one lemma, one open bound

**Theorem 7 (the race).** Under Type I strain, $\sup|\nabla u| \le A/(T - t)$, a length $\ell \le L(T - t)^\lambda$
satisfies $\ell/\sqrt{\nu/\sup|\nabla u|} \le L\sqrt{A/\nu}\,(T - t)^{\lambda - 1/2}$: for $\lambda > \tfrac12$ the
structure is inside its viscous scale before $T$; at $\lambda = \tfrac12$ the ratio is a Reynolds number of the
structure; below $\tfrac12$ nothing follows. *The measured closing is exponential - inside its viscous scale faster than any power law - and did not outrun
viscosity. The dangerous class is $\lambda = \tfrac12$ exactly, which is where the Ganeshram-Duruisseaux-
Anandkumar travelling profile sits by construction.*

**Lemma 8 (Kelvin).** The jump across a vortex sheet is circulation per unit length; circulation on a material loop
is conserved; hence $\tfrac{d}{dt}\log\Delta u = -\hat t\cdot S\,\hat t$ - the jump changes only by compression of
the sheet *along itself*, never by thinning across it.

**Conjecture C25, reformulated.** For the seam class, $\sup_t \Delta u < \infty$ iff the total in-plane compression
$\int_0^T (\hat t\cdot S\,\hat t)_-\,dt$ is finite. This replaces a global hypothesis (bounded velocity: the whole
problem, by Serrin) with a geometric one on one structure - *a sheet cannot be compressed along itself forever* - and
the runs above are two solutions on which it holds (9% in 0.65 time units). What would close the seam class is a
two-sheet Biot-Savart bound on that compression: pen-and-paper mathematics of the 1858-1871 kind. Not proved here.

## 4. What this is not

Not a proof, not a claim about (A)/(B), not a claim about the forced route. One adversarial field (plus a second,
a tube pair, Kida-Pelz), $\mathrm{Re} \approx 500$-$1500$, $\lambda$ measured on the resolved window only. What it is:
the handover from the inviscid nineteenth-century laws (Cauchy 1815, Biot-Savart, Helmholtz 1858, Kelvin 1869, the
sheet instability of 1868-71) to the viscous term of 1822/1845, put under a ruler for the first time on the structure
that every singularity candidate is built from.

For a century and a half the field asked how *big* a solution could get, when the answer was in what *shape* it took.
The norm era was not wrong; it was aimed at the wrong variable for this one problem, and productive enough elsewhere
that nobody noticed the exception. A singularity, if it exists, is a shape. The proof, if it comes, will be about a
shape. The tools for shapes predate the detour.
