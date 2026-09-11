# The horizon: when a binary64 answer stops being worth having

The owner's question, answered with measurements: after how much
simulated time does binary64 break down where a wider format still
holds, what kind of problem shortens that time, and where binary64 is
entirely sufficient - which is most places. Every number below is in
docs/VALIDATION.md with the run that produced it; this file is the
argument built from them. Binary64 numbers come from REBOUND's own
IAS15 (`build/ias15_ref`), which tools/check_equivalence.py shows is
the port at binary64 bit for bit and which runs 900 times faster; the
wide-format numbers come from the port on the software backend.

## The criterion

A long integration exists to say where a body is. The error that
grows fastest and that the user feels is the **phase**: how late or
early the body is on its orbit. tools/oracle.py measures it against
the closed two-body solution as the along-track timing error, the
position error projected on the true velocity, in units of the
period - so that it does not depend on where in the orbit a sample
happens to land (a phase error shows up in `|dx|/a` 200 times larger
at the pericentre of an e = 0.99 orbit than at its apocentre). The
energy is kept beside it: it is the invariant everyone reports, and
it is precisely the quantity that hides the problem, because a
random-walking energy error of 1e-13 is a phase error of 1e-6 orbits
by the time it has walked for a million orbits.

"Stops being worth having" is the user's threshold, so the tables give
the crossing of three: **1e-6 of an orbit** (30 seconds in a year: an
ephemeris-class demand), **1e-3** (hours in a year: the body is where
you said it would be, to a resolution element), and **0.1** (the phase
is being lost). The percentile is the 90th of an ensemble, for the
reason in the next section.

## Two growth laws, and why a single run cannot measure the first

**Round-off is a random walk.** Brouwer (1937): unbiased rounding
makes the energy error grow as t^1/2 and the phase error, which
integrates the period error, as t^3/2. A single run is one draw from
that distribution and does not follow the law: the same binary64
Kepler problem run to a million orbits at three eccentricities gave
apparent phase exponents of 0.75, 1.9 and 1.9 with 0.3 to 0.5 dex of
scatter, and one of them (e = 0.9) sat *below* the e = 1/2 run at the
same orbit count by luck of its realisation. An ensemble is the
measurement: 64 copies of the e = 1/2 problem with the planet's
starting x shifted by 2k ulps, each run for 1e7 steps at REBOUND's
default epsilon (2e5 orbits), give at every sample a distribution
whose 90th percentile is a power law to 0.02 dex:

    p90 phase error = 2.11e-16 * orbits^1.525
    median          = 2.07e-16 * orbits^1.446
    p90 energy error ~ orbits^0.5, 1.3e-13 at 1.8e5 orbits

The one long single run (1e8 steps, 2e6 orbits) landed at 9.9e-7 of
an orbit, where the p90 law says 1.1e-6; a billion-step run (1e7
orbits at the fixed step 1/16) landed at 8.2e-7, a low draw against
that step's median law of 2.6e-6, with its energy error still in the
1e-14 to 1e-13 band. This is the law of the arithmetic and no setting
moves it: a smaller step adds steps and round-off with them - the
fixed-step ensemble at dt = 1/16, 100 steps an orbit against the
adaptive run's 51, has p90 = 3.21e-16 * orbits^1.471, a coefficient
1.5x larger for a walk 1.4x longer.

**Truncation is a ruler.** Every wide-format record of the sweep -
where the method's own error is far above the format's floor - has
energy error proportional to t^1.00 and phase error to t^2.00, both
to 0.00 dex over 200 orbits: a systematic per-orbit error, the same
every orbit because the orbit is the same, integrated twice. Its
coefficient is the step control's:

    binary128, epsilon 1e-9  (REBOUND's default, 51 steps an orbit):  phase = 1.71e-21 * orbits^2
    binary128, epsilon 1e-12 (138 steps an orbit):                    phase = 3.86e-28 * orbits^2
    binary128, epsilon 1e-16 (513 steps an orbit, at its own floor):  phase = 1.9e-36  * orbits^2.07
    binary128, fixed dt 0.05 (126 steps an orbit):                    phase = 7.04e-24 * orbits^2

binary256 gives the same numbers to every digit at every one of these
settings (docs/VALIDATION.md, the sweep): it differs from binary128
only below binary128's floor, which none of these runs reaches in the
phase. So "a wider format" means binary128 here, and binary256 is for
the one regime the sweep found (the outer solar system at 10-20 day
steps) and for ensembles on the tile, which is a binary256 engine.

## The horizon of a regular orbit

Kepler, e = 1/2, the crossing of the 90th-percentile phase error
(binary64: measured to 2e5 orbits by the ensemble and to 2e6 by the
single run, extrapolated on the t^1.5 law beyond; binary128: the t^2
law, measured to 200 orbits, extrapolated beyond):

| phase error reaches | binary64, any epsilon | binary128, epsilon 1e-9 | binary128, epsilon 1e-12 | binary128, epsilon 1e-16 |
|---|---|---|---|---|
| 1e-9 orbits | 2.4e4 orbits (measured 2.7e4) | 7.6e5 | 5.1e7 | 2e13 |
| 1e-6 orbits | 2.2e6 (one run measured 2.0e6) | 2.4e7 | 5.1e10 | 2e14 |
| 1e-3 orbits | 2.0e8 | 7.6e8 | 1.6e12 | 7e15 |
| 0.1 orbits | 4.2e9 | 7.6e9 | 1.6e13 | 7e16 |
| cost per orbit, force evaluations, relative | 1 (51 steps x 3.4 passes) | 2.5 (51 x 8.6) | 4.9 (138 x 6.1) | 12.5 (513 x 4.2) |

Read it across. **At REBOUND's default tolerance the wider format buys
an order of magnitude in horizon at the tightest threshold and almost
nothing at the loosest**: 2.4e7 against 2.2e6 orbits at 1e-6, 7.6e9
against 4.2e9 at 0.1, because the method's t^2 truncation growth
overtakes the arithmetic's t^1.5 walk (the two laws cross at 1.5e10
orbits at this epsilon). The gain is bought with the step control:
at epsilon 1e-12, 2.7 times the steps, binary128 holds 1e-6 of an
orbit for 5e10 orbits, 23,000 times longer than binary64 can at any
setting, and at 1e-16 it is at its own floor and the horizon is
beyond any run anyone will make. That is the shape of the answer:
binary64's horizon is fixed by the arithmetic; binary128's is set by
the method and is for sale, at 5 to 12 times the force evaluations per
orbit plus the format's own cost per operation (1.3x on the software
backend, nothing on the tile, whose rate is format-blind).

For scale: 2e6 orbits is 2 Myr of an Earth-like orbit or 25 Myr of a
Jupiter-like one; 2e8 orbits is 200 Myr and 2.4 Gyr. The solar
system's own age is 4.5e9 Earth orbits.

## What shortens the horizon

**Eccentricity.** The energy of an e = 0.99 orbit is a small
difference of a kinetic and a potential term each 200 times larger
than it at pericentre ((1+e)/(1-e)), and the step control spends 159
steps an orbit there against 51 at e = 1/2. Both laws are worse for
it. The round-off floor, measured as the binary64-minus-binary128
difference on one replayed sequence of 20,000 steps, is 1.7e-13 in
energy and 1.7e-11 orbits of phase after 126 orbits at e = 0.99,
against 3e-15 and a few 1e-13 at e = 1/2: about fifty times, for
three times the steps an orbit. The method's truncation ruler at
REBOUND's default tolerance is 1.57e-20 orbits^2 at e = 0.99 against
1.72e-21 at e = 1/2, nine times. One long binary64 run at the default
settings (a single draw, the caveat above) had a phase error of
1.9e-6 of an orbit after 6.3e5 orbits, where the e = 1/2 ensemble's
p90 law gives 1.3e-7. So a highly eccentric orbit reaches every
threshold ten to fifty times sooner and costs three times more per
orbit to get there, at either format - the binary128 floor at
e = 0.99 is 9.7e-32 in energy against 2.3e-33 at e = 1/2, the same
factor of fifty, and 2^60 under binary64's at either eccentricity.
The eccentricity is the orbit's amplifier, not the format's.

**Chaos.** Past its Lyapunov time a system forgets its initial
condition exponentially, and the round-off floor is the perturbation
it forgets from. Two runs of the same chaotic system that differ only
in their arithmetic diverge from the floor at the system's own rate,
so a wider format buys exactly its extra bits of divergence before
the solution is lost. Measured on Burrau's Pythagorean problem
(masses 3, 4, 5 at rest on a 3-4-5 triangle; close encounters, a step
spanning four orders of magnitude, a binary that ejects the third
body near t = 60-70) with three formats on literally the same 6,000
steps, the binary64 solution differs from the binary256 one by 1e-12
at t = 6.9, 1e-9 at t = 23, 1e-6 at t = 46, 1e-3 at t = 60 and is
unrelated to it from t = 66, before the ejection that is the whole
point of the problem; the binary128 solution differs from binary256
by 1e-30, 1e-27, 1e-24 and 1e-21 at t = 6.9, 23, 46 and 61 - the
same curve, eighteen decades (2^60) lower, crossing each threshold at
the same time to within one sample - and ends the evolution at t = 73
three parts in 1e21 from binary256. The curve is a staircase, flat
between encounters and a jump at each, averaging a decade every four
time units; one close passage costs binary64 three orders of
magnitude of energy (4e-14 to 7e-11 in a step) and costs binary128 an
energy error that is binary256's to every digit, i.e. the method's.

The method's own error does not appear in any of those differences,
so it was measured separately at one time: two step sequences, at
epsilon 1e-9 and at 1e-11, cut to end at exactly t = 60 and both
replayed at binary256, differ by 1.3e-13 there, while the binary64
run of the epsilon 1e-9 steps is 3.4e-4 from the same reference.
The truncation-seeded error is nine decades *below* the round-off-
seeded one at the moment the binary64 solution is lost - a
fifteenth-order method's truncation error is smooth and lies along
the flow, a time shift the encounters do not amplify, whereas
round-off is random and has the transverse components they do. (The
first draft of this document inferred the opposite from the sweep's
energy errors and was wrong; docs/VALIDATION.md keeps both.) So on
a chaotic problem at REBOUND's default tolerance a wider format buys
its first nine decades of trajectory accuracy for nothing but the
format's cost, and the tighter step control that the regular orbit
needs is only for the science that wants the trajectory further
still.

And the ensemble says the same thing from the other side. Eight
copies of the problem with one coordinate displaced by 1, 2, 4, ...
64 binary64 ulps, integrated together on the same steps: at
binary256 their divergences from the unperturbed member stay in the
ratios 2.00, 4.00, 8.00, 16.00, 32.00, 64.00 at every sample through
the whole evolution while the common amplification climbs to 2e8 -
the linear-response regime a Lyapunov measurement assumes, seen
directly; at binary64 the same eight are within a factor of four of
each other in no order from t = 10 on, the seeds drowned in the
arithmetic's noise. At binary64 an ulp-scale ensemble measures the
arithmetic; at binary256 it measures the dynamics.

**Close encounters** are where the amplifying is done, in bursts, and
they are what a wide format is *for* in a chaotic system: not a longer
Lyapunov horizon in itself, but the encounter passed with the energy
still the method's and the ensemble still a measurement. **Long
secular runs** are where the t^1.5 law gets its time.

## Where binary64 is entirely sufficient

Most places, and it is worth being plain about it:

- Any regular system integrated for fewer than about a million orbits,
  if 1e-6 of an orbit of phase is enough, or a hundred million if 1e-3
  is. That is nearly every N-body integration whose purpose is not a
  Gyr-scale phase. The outer solar system's energy error at 40-day
  steps stays between 2e-14 and 5e-14 for 1e7 years (the ledger), and
  its angular momentum at 2e-14: nothing a model of the real system
  needs is limited by it.
- Any chaotic system beyond its Lyapunov horizon, which is where the
  interesting ones are run: the phase is unpredictable there at any
  precision, the science is statistical, an ensemble is the tool, and
  binary64 members are as good as any.
- Any problem whose physics (masses, initial conditions, the model) is
  uncertain at the 1e-10 level, which is all real ones; the arithmetic
  is not the weakest link.
- Speed-limited work of every kind. On this host REBOUND's binary64
  IAS15 takes 5 to 9 us a step; the port at binary128 on the software
  backend takes 25 ms, and the tile's binary256 rate is a few thousand
  system-steps a second per tile however wide the ensemble
  (docs/HARDWARE.md).

## The niche

Where a wider format pays is the intersection: a regular or weakly
chaotic system, integrated for long enough that a t^1.5 walk from
1e-16 reaches the phase precision the science needs, with the step
control tightened enough that the method's own error does not get
there first, and with the patience for 5 to 12 times the force
evaluations. Concretely: **phases of a regular system past a few
million orbits at 1e-6, or past a hundred million at 1e-3, at epsilon
1e-12 or tighter** - a Gyr of a planetary system that is not chaotic,
a long-lived resonance whose libration phase matters, a stability map
where the question is a phase and not a survival. And, separately,
**a certificate**: a result that is exact rather than close, that the
card and the software backend reproduce bit for bit, that an ensemble
reproduces bit for bit against its members, and whose round-off floor
is measured by differencing two formats rather than estimated.

## What exactness permits that a tolerance does not

Three things, all done here. Two integrations can be compared *without
a tolerance*: the difference between a binary64 and a binary256 run of
the same steps is the binary64 round-off, exactly, which is how every
floor in the ledger was measured. An ensemble can be gated: 64 or
1,000 systems in one vector are demonstrably each the solo
integration, not approximately (tools/check_ensemble.py). And a run can
be replayed: the recorded step sequence, fed back, reproduces the record
to the bit - the `--dt-file` entry of the ledger, 2026-09-10, on Kepler
and the Pythagorean problem at 2,000 steps each. That one is a
measurement rather than a gate: nothing in `make check` replays a step
sequence, and `tools/check_ensemble.py` - which an earlier version of
this line credited for it - has no `--dt-file` case.

What exactness does not give IAS15 is an exact return. Integrated
forward N steps and back N steps from the exact recorded state
(tools/reversal.py), it does not come back bit for bit at any format:
its Gauss-Radau nodes include the start of the step and not the end,
so the backward polynomial is a different polynomial and the method
is not time-symmetric even in exact arithmetic. What the test shows
instead is *which error is which*: on Kepler at dt = 1/16, 1,000 steps
out and back, binary128 and binary256 return to the same 7.81e-19 of
the orbit - the method's own asymmetry, the arithmetic invisible under
it - while binary64 returns to 1.24e-13, all of it arithmetic; on the
outer solar system 1.16e-31 and 1.15e-31 against 7.9e-16. A bit-exact
return needs a symmetric or a lattice scheme, which is why JANUS is
ranked where it is in docs/INTEGRATORS.md.
