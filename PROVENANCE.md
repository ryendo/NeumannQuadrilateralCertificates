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

## Checked-in computations (2026-08-26--27)

The checked-in results were computed on `liulab-hpc2023` from source commit
`b9944f600ba63bed3ddbbde9890febebb613b3ac`.

They are checked-in outputs, not outputs regenerated from the current
working tree. The current local code uses the literal formulas (34) and (39).
The current global code tests all five conditions in (44) and records the
aggregate lower bounds for `q`, `lambda_1`, and `Delta`. A new local and global
run is therefore required to attach current-source provenance to new numerical
summaries. In particular, the old global JSON schema does not record the lower
bound for `lambda_1` on every accepted box.

| computation | workers | recorded completion | recorded minimum |
|---|---:|---|---:|
| local | 60 | 13,824/13,824 top-level boxes, 0 failures | `S >= 0.0044848550155087707` |
| global | 18 | 16/16 retained initial boxes, 0 unverified boxes | margin `>= 8.0698658297961856e-6` |

The local audit found exactly the box IDs `1:13824`, without duplicates,
all 60 completion markers, and the recorded ten-column schema. It independently
recomputed the saved extrema, including
`lambda_min(M) >= 0.052100929203335616` and
`lambda_3 >= 17.834653422811055`.

The global audit found worker IDs `1:18`, recomputed every sum and extremum in
`summary.json`, and checked the recorded `complete=true` and `unverified=0`
fields for every worker. Each of the 16 retained initial boxes was processed
directly by `qn_global_certified_cover`; workers 17 and 18 have no assigned
box. This audit checks the aggregate files only. The old bisection used
floating-point child descriptors that can leave a rounding-size gap, and the
leaf boxes were not saved, so the recorded global run is not independently
auditable as the finite cover in (49). It must be recomputed with the current
outward-covering bisection.

The per-certificate `RUN_PROVENANCE.txt` files record the exact environment,
and the `SHA256SUMS` manifests cover all saved outputs.

## Mathematical source alignment

The local Taylor coefficients originated in
`ryendo/dq2-quadrilateral-certificate` commit
`6d44991af03102a2338e7af69c67767a32e55178`. The implementation in the current
working tree keeps Taylor-remainder accumulation and Gershgorin radii in
interval arithmetic and uses the exact two-root identity described in the
appendix. These source changes postdate the saved run identified above.

The global implementation is a mathematical port of the supplied
`quad_neumann_global` code. For the fixed ambient trigonometric trial space,

```text
grad_ref(psi o Phi) = DPhi' grad_x(psi),
```

so the inverse metric cancels and the stiffness integrand is evaluated as
`(grad_x psi_i . grad_x psi_j) J`. This algebraically regular form removes the
source representation's artificial `1/J` pole without changing the trial
space, Rayleigh--Ritz pencil, box test, or certified conclusion. The interval
assembly uses the common exact kernel, INTLAB automatic gradients, and
Bernstein-ellipse Gauss--Legendre truncation padding.
