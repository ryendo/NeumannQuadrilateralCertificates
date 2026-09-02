# Current global completion record

This directory records the 16-worker global run from clean source commit
`671ddd9f6fb34dbb74b4474a3bdf71d11c53b766` (PBS job
`294.liulab-hpc2023`). `summary.json` reports

```text
initial retained boxes = 16
accepted               = 114627
discarded              = 15623
bisections              = 130234
unresolved              = 0
maximum depth           = 30
q lower minimum         = 0.8273243052571837
lambda_1 lower minimum  = 3.4247176582671237
Delta lower minimum     = 2.6545819956425021e-5
complete                = true
```

The merger checked the exact initial IDs
`[41,42,50,51,52,53,54,68,69,71,72,77,78,79,80,81]`, with no missing or
duplicate representative, as well as the binary-tree accounting identity,
the required radius, positive certified bounds, and common provenance.

The worker JSON files are aggregate completion records; they do not list every
accepted/discarded leaf. Rerun the pinned source to reconstruct the complete
subdivision. `RUN_PROVENANCE.txt` records the environment and `SHA256SUMS`
authenticates every checked-in result and provenance file in this directory.
