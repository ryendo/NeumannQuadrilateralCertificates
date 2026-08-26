# Global result directory

`summary.json` merges `worker_001.json` through `worker_018.json`. The saved
certificate satisfies

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

Seventeen root boxes ran directly through `qn_global_certified_cover`. The
remaining root (worker 12) was split along longest sides into an exact dyadic
cover and evaluated in parallel with the same discard and certificate
predicates. `parallel_parts_012.json` records the 50 disjoint completed
subcovers used in the merge; their dyadic intervals cover the root exactly,
and the 49 fixed bisections are included in `worker_012.json` and the global
total. `parallel_driver_012.m` is the run driver retained for auditability.

`RUN_PROVENANCE.txt` records the environment. `SHA256SUMS` covers every saved
result and provenance file except itself.
