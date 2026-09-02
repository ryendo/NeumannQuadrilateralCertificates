# Source and computation provenance

## Verified numerical dependencies

Both certificate paths use MATLAB R2023b, INTLAB V12, and the MIT-licensed
[`yuuka-math/veigs`](https://github.com/yuuka-math/veigs) package at commit
`6556d39a0d9819bb172d232062b698aa76e420f6`. The dependency is not bundled;
set `VEIGS_ROOT`, or pass its checkout as the second argument to `qn_setup`.

`qn_veigs_indices` symmetrizes each interval pencil and targets every requested
index separately. It accepts a result only when the certified index data
contains that index. If `veigs` cannot separate a target, the wrapper calls
the verified small-pencil `veig` routine from the same pinned package. The
local calculation targets indices 1 and 3 and requires
`inf(lambda_3) > 16`; the global calculation targets index 1.

## Checked-in computations (2026-09-01--02)

The checked-in results were recomputed from the clean source commit
`671ddd9f6fb34dbb74b4474a3bdf71d11c53b766`.
The executable-tree digest for `src/`, `scripts/`, `tests/`, and `data/` is
`a063015a056a186841d0d040f649c9e24b4baa9012e601ddf31c40aeadb75d12`.

| computation | workers | completion | certified minimum |
|---|---:|---|---:|
| local | 32 | 13,824/13,824 top-level boxes, 0 failures | `S >= 0.010109` |
| global | 16 | 16/16 initial representatives, 0 unresolved boxes | `Delta >= 2.6545e-5` |

The local run additionally records
`lambda_1 >= 7.630`,
`lambda_2 <= 13.268`,
`lambda_min(M) >= 0.0521`, and
`lambda_3 >= 17.834`. Its 32 metadata files record the common
source/dependency provenance and the SHA-256 digests of the corresponding CSV
and completion marker. The summarizer checks the exact box IDs `1:13824`, no
duplicates, all 32 completion markers, every bound, and every artifact binding
before reporting `verified=true`.

The global run accepted 114,627 boxes, discarded 15,623, performed 130,234
bisections, reached depth 30, and left no unresolved box. It records
`q >= 0.8273` and
`lambda_1 >= 3.4247` on accepted boxes. The merger verifies the
exact 16 initial-box IDs, the binary-tree accounting identity, the required
radius, all five conditions represented in the current schema, and common
source/dependency provenance before reporting `complete=true`. The worker JSON
files are aggregate completion records, not a saved list of every leaf; the
full leaf subdivision is reconstructed by rerunning the source.

Both jobs used MATLAB R2023b Update 5. Their recorded INTLAB startup digest is
`fd313a5a13bca7153627f9c1f875ecc27a75efb36a48b06edf4343103067d4b5`,
and the INTLAB tree digest (excluding INTLAB's generated startup cache) is
`5f5ab318797270747ff92a3e3ef10f1f021a3a99603ccf2d926e08341a0a132a`.
The per-certificate `RUN_PROVENANCE.txt` files record the run ID, software
environment, start/end times, and clean source state. The
`SHA256SUMS` manifests cover all checked-in result, metadata, completion, and
provenance files; diagnostic worker logs were retrieved separately and are not
part of the repository certificate artifacts.
