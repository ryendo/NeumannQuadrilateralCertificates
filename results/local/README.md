# Current local completion record

This directory records the 32-worker local run from clean source commit
`671ddd9f6fb34dbb74b4474a3bdf71d11c53b766` (PBS job
`295.liulab-hpc2023`). `summary.json` reports

```text
top-level direction boxes       = 13824 / 13824
failures                        = 0
maximum subdivision depth       = 1
S lower minimum                 = 0.010109281357063793
lambda_1 lower minimum          = 7.6300046469977758
lambda_2 upper maximum          = 13.267374191855003
lambda_min(M) lower minimum     = 0.052100929203335616
lambda_3 lower minimum          = 17.834653422811055
lambda_3 - 3*pi^2/2 lower min   = 3.0302468211770126
verified                        = true
```

Every `meta_NNN.json` records the common source, `veigs`, MATLAB, and INTLAB
provenance and binds `res_NNN.csv` and `done_NNN.txt` by SHA-256. The
summarizer checked the exact IDs `1:13824`, all completion markers, every
recorded bound, and all 32 artifact bindings. `SHA256SUMS` authenticates every
checked-in result, metadata, completion, and provenance file in this
directory.
