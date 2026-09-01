# Checked-in global summary

`summary.json` merges `worker_001.json` through `worker_018.json`. The saved
aggregate fields are

```text
verified = 117083
discarded = 15222
bisected = 132289
unverified = 0
veigs_certified = 117083
max_depth = 29
min_certified_margin = 8.0698658297961856e-6
worker_count = 18
initial_retained_all = 16
complete = true
```

The run used source commit `b9944f600ba63bed3ddbbde9890febebb613b3ac`,
MATLAB R2023b Update 5, INTLAB V12, and `veigs` commit
`6556d39a0d9819bb172d232062b698aa76e420f6` on `liulab-hpc2023`.

Each of the 16 retained initial boxes ran directly through
`qn_global_certified_cover`. The 18 output slots are kept for the worker
launcher convention; workers 17 and 18 have no assigned initial box.

These files are historical summaries, not a current finite-cover certificate.
This source predates the explicit lower-`lambda_1` record required by the
current implementation of (44). Its floating-point bisection can also leave a
rounding-size gap between children, and the leaf boxes were not saved, so the
old run cannot be re-audited for (47)--(49). Recompute the global calculation
with the current source before using it as the paper's global certificate.

`RUN_PROVENANCE.txt` records the environment. `SHA256SUMS` covers every saved
result and provenance file except itself.
