# Global result directory

`summary.json` and `worker_001.json` through `worker_018.json` are the output
of the current MATLAB/INTLAB/veigs run of `qn_global_certified_cover`. A
complete certificate must have:

```text
complete = true
unverified = 0
min_certified_margin > 0
```

The checked-in summary satisfies all three conditions:

```text
verified = 166928
discarded = 25285
bisected = 192195
unverified = 0
veigs_certified = 166928
max_depth = 30
min_certified_margin = 2.6545819961754091e-5
worker_count = 18
initial_retained_all = 18
complete = true
```

This liulab run used source commit
`7a0d1b42e0c9660fd66a54feca2e38555a940e34` and `veigs` commit
`6556d39a0d9819bb172d232062b698aa76e420f6`. The global source is unchanged
in repository commit `2d4d4ce4d40202614f98107f3a284c632ce05c13`.
`RUN_PROVENANCE.txt` records the deployed revisions, and `SHA256SUMS` covers
the summary and all 18 worker outputs.

The worker IDs, completion flags, counts, maximum depth, and minimum margin
were independently recomputed from the JSON files. These current results
supersede the earlier one-/two-vector interval Rayleigh-test counts.
