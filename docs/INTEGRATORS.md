# Which other REBOUND integrators port well

An assessment of REBOUND's other integrators for the same treatment
IAS15 got here: arithmetic through `cft.h`, binary64 first as a
bit-for-bit gate against REBOUND, then binary128 and binary256, with
the per-lane work as sequencer programs where it fits. Three questions
per integrator, from reading the pinned source
(`third_party/rebound/src/integrator_*.c`, REBOUND 5.1.1 at
bdfda4bd), not from running anything:

1. **Is it accuracy-limited or speed-limited?** A method whose error
   floor is round-off benefits from a wider format; a method whose
   error is its own truncation, or whose purpose is throughput, does
   not.
2. **Does it need transcendentals?** libcft provides them correctly
   rounded (`cft_sin`, `cft_exp`, `cft_pow`, ...), so the question is
   cost and whether they sit inside the per-step loop, not whether
   they exist.
3. **Is the per-step state small enough to stay resident?** On the
   tile the state that must persist between runs is what a lane
   carries in its scratch block (256 format-width slots a lane today -
   cft-fp256's docs/SEQUENCER.md, not this project's docs/HARDWARE.md,
   which counts what IAS15 puts in those slots rather than how many
   there are) or what the host can hold in device-resident buffers.

The IAS15 measurements in docs/VALIDATION.md are the calibration for
the first question: a wider format lowers the round-off floor exactly
as far as the format's precision says, but only when the method's
truncation error is already below it, and IAS15's 15th order reaches
that only by shrinking the step and by iterating its corrector longer.

## WHFast

Wisdom-Holman with symplectic correctors and kernels (Rein & Tamayo
2015; Rein, Tamayo & Brown 2019). **Accuracy-limited by round-off in
its intended regime.** For a nearly-Keplerian system with a fixed step
the truncation error is a bounded oscillation set by the step and the
corrector order, while round-off accumulates as Brouwer's t^(1/2) walk
in energy and secularly in phase; Rein & Tamayo spent a whole section
on keeping that walk unbiased at binary64 because it is what limits a
Gyr integration. A wider format lowers that walk by 2^(53-p) and does
nothing to the bounded part, which is the useful shape: no smaller
step is needed to collect the gain, unlike IAS15. **No
transcendentals.** The Kepler step uses Stumpff series with
quarter-angle recurrences and a Newton/Halley iteration on the
universal anomaly - multiplies, adds, one square root per particle-step
and a handful of divides. The 1/n! table is written as `1./5040.` etc.
in the source, i.e. derived by the compiler from exact integers rather
than typed as decimals, which is the right kind of literal; at binary256
the same integers divide exactly the same way through `cft_div`. **State
is tiny**: Jacobi (or heliocentric, WHDS, barycentric) coordinates of N
particles, 6N values plus masses - well inside one lane's scratch even
per particle. The divides in the Kepler solver are the one obstacle to
a single sequencer program, exactly as for IAS15 (docs/HARDWARE.md);
the FMA-form treatment applied here to IAS15 applies to them
unchanged, with reciprocal constants for the fixed divisors and the
seed-plus-Newton route for the data-dependent ones.

## WHFast512

The same algorithm hand-written in AVX-512 assembly for at most 8
planets, integrating 10^6 steps per call. **Speed-limited by
construction and binary64 by construction**: its reason to exist is
the eight-wide double register. It could not even be built in this
tree with MinGW (ref/whfast512_stub.c). Nothing about it is worth
porting that is not already WHFast; a tile running WHFast over eight
lanes at binary64 is the same idea in a different machine.

## SABA and SABAC

Laskar & Robutel's SABA1-4, the corrected SABAC variants, and the
generalised-order SABA(10,4), SABA(8,6,4), SABA(10,6,4) compositions,
all on the same drift-kick kernel as WHFast. **Accuracy-limited the
same way as WHFast, with a much lower truncation floor at the same
step**, since the higher-order compositions push the leading error to
O(eps dt^10 + eps^2 dt^6 ...) for perturbation eps. That makes them the
natural partner of a wide format: at binary256 one can shrink the step
until a SABA(10,6,4) truncation error meets a 2^-237 round-off floor
without the 15th-order penalty IAS15 pays in corrector iterations.
**No transcendentals** beyond the Kepler step's square roots. **State
tiny.** The one thing to get right is the coefficient table:
`integrator_saba.c` carries the stage coefficients as 40-digit decimal
literals, which is above binary64 and below binary256, and some are
algebraic (SABA2's 1/2 - sqrt(3)/6, SABA3's 1/2 - sqrt(15)/10, SABA4's
Gauss-Legendre nodes) while the generalised-order ones are roots of
the order conditions. They must be re-derived at 130 digits by the
same rule tools/gen_constants.py applies to the Gauss-Radau nodes, and
checked by reproducing the published doubles.

## MERCURIUS

Hybrid WHFast with an IAS15 switch-over during close encounters,
Chambers-style. **Accuracy-limited in the encounter, speed-limited
outside it**, and the switching itself is a modelling error the
format cannot touch. **Transcendentals inside the step**: the default
switching function is polynomial, but `L_infinity` is `exp(-1/x)` and
the critical radii use `cbrt` - both correctly rounded in libcft, both
cheap at the encounter count. **State is moderate** (two particle
backups, the switching radii, the encounter map), and the control flow
is data-dependent: which particles are in an encounter changes per
step, which on a lane machine is a `SETACT` mask and on the host is a
loop. Worth doing only after WHFast and IAS15 both exist at the wide
formats, since it is exactly their composition; the gain over running
IAS15 alone at binary256 is speed, not accuracy.

## TRACE

Hernandez & Dehnen's time-reversible hybrid, WHFast with BS or IAS15
for encounters and pericentre passages, with a reversibility check that
can reject and redo a step. The same profile as MERCURIUS with heavier
host control flow (the step can be repeated with a different
partition), no transcendentals of its own beyond what its pieces use,
and the largest per-step state of the family (three particle backups).
Last on the list: everything it needs is either WHFast or IAS15, and
its own contribution is bookkeeping the host does anyway.

## JANUS

Rein & Tamayo's bit-wise time-reversible integrator: positions and
velocities live on an integer lattice (`scale_pos = scale_vel =
1e-16`, `int64_t` coordinates) and only the force is floating point.
**Its accuracy is set by the lattice spacing, not the format**, so at
binary64 a wider float buys nothing; the interesting move is the other
way round - a finer lattice needs wider integers, and the tile's
integer group (`IADD`, `ISUB`, `ISHL`, `ISHR` at the full format width,
256 bits at binary256) is that wider integer for free. A JANUS whose
lattice is 2^-200 with a binary256 force is a genuinely new object,
exactly reversible bit for bit, which is this project's own kind of
claim. **No transcendentals.** **State tiny** (6N integers, 6N floats).
The lattice rounding (float to integer and back, `to_int`/`to_double`)
is a convert the sequencer does not have as an opcode, but at fixed
scale it is one add of a large constant and one subtract, which it
does have. High on the list for what it says rather than for what it
computes.

## EOS

Rein's embedded operator splitting: an outer splitting (leapfrog,
Yoshida LF4/LF6/LF8, the modified-kick PMLF4/PMLF6, PLF7_6_4) with an
inner one applied n times per outer step. **Accuracy-limited like the
SABA family and designed to reach machine precision cheaply**, so a
wide format is what lets it go further. **No transcendentals**; the
modified kicks need the jerk (the force's time derivative), which is
more arithmetic of the same kind. **State tiny.** Its obstacle is the
same table problem as SABA, worse: `integrator_eos.c` holds its
composition coefficients as decimal literals of 16 to 40 digits
(`pmlf6_z`, `pmlf6_y`, `plf7_6_4_a`, ...), some visibly truncated at
16, and they are numerical solutions of order conditions rather than
closed forms. A binary256 EOS needs those systems re-solved at 130
digits, a real but bounded job, and until then the port would be
capped at the table's own precision - the exact trap the IAS15
constants were regenerated to avoid.

## LEAPFROG

Stormer-Verlet and its 4th, 6th and 8th order compositions. cft-fp256
already carries this on libcft as `cft-orbits` with a 300-digit oracle
(its docs/ORBITS.md), and measured what a wide format is worth here:
the energy error is the method's and identical at every format, the
angular-momentum drift is the arithmetic's and falls as 2^-p. Nothing
left to learn from porting it again; its value is as the calibration
of the others. Trivially feasible, low value.

## SEI

The symplectic epicycle integrator for the shearing sheet (planetary
rings). **Speed-limited**: ring simulations are collisional and
dissipative at 10^4-10^6 particles, and round-off is nowhere near
their error budget. The sines and tangents of OMEGA dt are computed
once per step size and cached, so transcendentals are a set-up cost,
not a loop cost. The drift-kick part would run on a lane machine
beautifully - particles are lanes - but the collision search and the
tree gravity are host work, and nothing in the physics wants
binary128. Feasible, no reason.

## BS

Gragg-Bulirsch-Stoer with Richardson extrapolation, adaptive in step
and order. **Accuracy-limited** and general (it integrates arbitrary
ODEs), but its cost for d digits grows with the extrapolation depth
roughly as d^2, and its step control uses `pow(error, -1/(2k+1))` -
correctly rounded in libcft, and cheap at one per step. The
extrapolation tableau is sequential, division-heavy and data-dependent
in depth, so it is a host loop over `cft_run` calls, not a program.
Moderate value for the non-N-body ODE surface; moderate feasibility.

## Ranking, value times feasibility

| rank | integrator | value at a wide format | feasibility | note |
|---|---|---|---|---|
| 1 | **WHFast** | high: round-off is its long-term limit and the gain needs no smaller step | high: tiny state, no transcendentals, Kepler step in FMA form plus one sqrt | the workhorse; port this first |
| 2 | **SABA / SABAC** | high: same kernel, truncation floor low enough to meet a 2^-237 walk | high, once the coefficient tables are re-derived at 130 digits | the constant tables are the whole job |
| 3 | **JANUS** | high for what it demonstrates: a 256-bit lattice, exactly reversible | high: the tile's integer group is the wider lattice | a new object, not a faster old one |
| 4 | **EOS** | high, same family as SABA | medium: coefficients are numerical solutions of order conditions, some truncated at 16 digits in the source | re-solve, then it is SABA's job again |
| 5 | **MERCURIUS** | medium: accuracy comes from its IAS15 half, which exists | medium: data-dependent switching, exp and cbrt in the loop | after 1 and IAS15 |
| 6 | **BS** | medium: general ODEs | medium: sequential tableau, host-driven | not a program |
| 7 | **TRACE** | medium | low: the heaviest host control flow of the family | last |
| 8 | **LEAPFROG** | low: already measured in cft-fp256 | trivial | calibration only |
| 9 | **SEI** | none: dissipative physics | high for the kernel, host for the rest | skip |
| 10 | **WHFast512** | none beyond WHFast | not buildable here; binary64 by design | port WHFast instead |

**Recommendation.** WHFast next, in exactly the way IAS15 was done
here: a bit-for-bit binary64 gate against REBOUND first, then the
formats, with the Kepler solver's fixed divisors turned into
reciprocal constants for the sequencer and its data-dependent divide
handled by the seed-and-Newton route the tile already runs. Then SABA
on the same kernel, which is mostly a constant-derivation task that
tools/gen_constants.py already knows how to gate. JANUS third, because
a 256-bit lattice is the one thing on this list that no CPU port of
REBOUND can offer.
