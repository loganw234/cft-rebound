# results/horizon: the records and tables of the 2026-09-10 horizon campaign

Every number in docs/HORIZON.md and in the 2026-09-10 entries of
docs/VALIDATION.md comes from a file here (or from a gate run whose
output is quoted in the ledger). Records (`.rec`) are the programs'
own output, every value an exact hex float; tables (`.csv`) are
tools/oracle.py's scoring of a record (`rel_phase_err` is the along-
track timing error in orbits) or tools/compare_formats.py's difference
of two records on the same steps or tools/divergence.py's member-
against-member distances. The commands are in the ledger entries;
the shapes are:

    fp64_single/     REBOUND's own IAS15 (build/ias15_ref, bit for bit the port at binary64):
                     Kepler e = 1/2 at the fixed step 1/16 for 1e4, 1e5, 1e6 and 1e7 orbits;
                     e = 1/2, 0.9 and 0.99 at epsilon 1e-9 for 1e8 steps; the outer solar system
                     at 40-day steps for 1e6 and 1e7 years
    kepler_ulps_E64.txt   the 64-member ensemble file (planet x shifted by 2k ulps)
    ens64_adapt/     the 64 members' oracle tables, epsilon 1e-9, 1e7 steps each (2e5 orbits);
                     tools/horizon.py results/horizon/ens64_adapt/m*.csv --fit-from 20000
    ens64_fixed/     the same at the fixed step 1/16 (1e5 orbits)
    pythagorean/     the fp256 adaptive seed (epsilon 1e-9) and its step sequence, the
                     binary64-rounded sequence dt64.txt, the three replays on it and their
                     scores, the format differences div_*.csv, the eight-member ulp ensemble
                     file and its records and divergences at fp64 and fp256, the epsilon 1e-11
                     seed and the three T* = 60 replays (the truncation measurement), and
                     REBOUND's own adaptive binary64 run to t = 1,142
    eccentric/       the fp256 seeds at e = 0.99 and e = 1/2 (20,000 steps), their rounded
                     sequences, the six replays and their scores, and the floors
    reversal/        forward and backward records and the turn-around and return problem
                     files of tools/reversal.py, Kepler and outer, three formats
    dtfile/          the --dt-file self-check: an adaptive run, its recorded sequence, and
                     the replay that reproduces it
    throughput/      the trailers of the three throughput scans (scan1 the first ensemble
                     build, scan2 after the convergence-test and step-control change, scan3
                     the final source) and the script that runs one

The 64-member records themselves (7.7 MB, reproducible in twenty
minutes with ens64_run.sh, whose paths are the scratch directory it ran in) are not
committed; their tables are.
