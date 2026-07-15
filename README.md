# Certified Neumann-eigenvalue bounds for quadrilaterals

> **Private supplementary-code repository. Do not make this repository public.**

This repository contains the MATLAB + INTLAB computations accompanying

> Braxton Osting and Ryoki Endo, *Maximizing Laplace--Neumann eigenvalues on
> quadrilaterals* (2026).

The paper proves that, among convex planar quadrilaterals of fixed area, the
square uniquely maximizes the first nonzero Laplace--Neumann eigenvalue. The
computer-assisted part uses a five-dimensional trigonometric Rayleigh--Ritz
space on a fixed reference square and has two certified components:

1. a local second-order difference-quotient certificate on a ball about the
   square; and
2. a global adaptive box cover outside a smaller seam ball.

Both components now use MATLAB and INTLAB. Python is not required.

## Certified statements

Let `p=(a,b,c,d)` be the paper's quadrilateral parameter and let
`K(p)x = lambda M(p)x` be the 5-by-5 Rayleigh--Ritz pencil.

### Local step

For

```text
rho_local = rho# = 3232/(27*pi^6) approximately 0.1245,
```

the DQ2 computation certifies

```text
f_2(p) > pi^(-2),       0 < ||p|| <= rho_local.
```

Equivalently, the final interval sign test proves `S(t,e)>0` on a finite cover
of `S^3 x (0,rho_local]`, where `p=t e`.

The saved result in `results/local/` has been checked as:

| quantity | certified saved value |
|---|---:|
| top-level boxes | 13,824 |
| certified boxes | 13,824 |
| failures | 0 |
| box-ID coverage | complete, with no duplicates |
| worker completion markers | 48/48 |
| minimum lower bound for `S` | 0.00493500422156500917 |
| minimum lower bound for `lambda_1` | 7.58452255183086610 |
| maximum upper bound for `lambda_2` | 13.2685132182832994 |
| maximum subdivision depth | 2 |
| wall time | 13.41 h on 48 workers |

### Global step

The global computation certifies

```text
|Q_p| lambda_1(p) < pi^2
for p in P_K \ B(0,rho_seam),

rho_seam = rho#/2 = 1616/(27*pi^6) approximately 0.062256.
```

Because the trial space is contained in the zero-mean Sobolev space,
`lambda_1(p) >= mu_1(Q_p)`, so the same strict inequality follows for the
first nonzero Neumann eigenvalue. The local and global regions overlap on
`rho#/2 <= ||p|| <= rho#`.

The source Python code called the seam `RHO_SHARP`; this repository uses the
paper-aligned names `rho_local` and `rho_seam` to remove that ambiguity.

> **Recalculation status.** No fresh full-cover result is checked in yet. The supplied Python
> extract is not the run-of-record driver named in its own README: it omits
> `eigbound_shrink.py`, the parent `quad_neumann_release/` tree, and
> `GLOBAL_STEP_AUDIT.md`. A faithful run of the supplied `run_cover(n_init=3)`
> reaches `qn:Jacobian` on boxes converging to the degenerate boundary of
> `P_K`. The current MATLAB/INTLAB assembly removes that artificial obstruction
> by using the exactly equivalent physical-gradient formula described below.
> On liulab, the revised smoke test certifies a box crossing the triangle face
> `c_3=0` by the 2-vector route with positive margin (approximately `1.3612`). A new full global run is
> still required before the global certificate can be declared complete.

The reference global run reported 1,420 verified boxes, 337 discarded boxes,
1,683 bisections, no unverified boxes, maximum depth 26, and positive certified
slack. A fresh MATLAB/INTLAB run writes its own record to
`results/global/summary.json`; do not substitute the reference counts for a
new run's output.

## Requirements

- MATLAB R2023b or later (the local computation was tested with R2023b).
- INTLAB 12. INTLAB is not bundled.
- A POSIX shell only for the optional launch scripts.

All non-integer proof constants and decimal bounds are parsed from strings,
for example `intval('0.5')`, `intval('pi')`, and `intval('1e-12')`. Floating
point is used only for choices that are mathematically valid for any choice:
the center test frame and the bisection heuristic. Every accept/reject proof
decision is made from INTLAB interval endpoints.

## Project structure

```text
.
├── QuadrilateralProofRunner.m       unified entry point
├── src/
│   ├── local/                       DQ2 local certificate
│   │   ├── dq2_run_certificate.m
│   │   ├── dq2_algorithm1_box.m
│   │   └── qn_summarize_local_results.m
│   └── global/                      adaptive global certificate
│       ├── qn_run_global_cover.m
│       ├── qn_certify_box.m
│       ├── qn_km_enclosure.m
│       └── qn_gl_pad.m
├── data/local/taylor_coefficients.mat
├── results/
│   ├── local/                       saved 48-worker certified run
│   └── global/                      output of a fresh global run
├── tests/qn_smoke_test.m
├── scripts/
│   ├── run_local_workers.sh
│   ├── run_global.sh
│   ├── run_global_workers.sh
│   └── run_smoke.sh
└── docs/PROVENANCE.md
```

## Setup and quick check

In MATLAB:

```matlab
r = QuadrilateralProofRunner('/path/to/Intlab_V12');
r.setup();
report = r.smokeTest();
```

The smoke test checks that the certified global pencil enclosure contains the
floating-point center value, exercises a representative per-box test, and
verifies the saved local CSV set.

From a shell:

```bash
export INTLAB_ROOT=/path/to/Intlab_V12
./scripts/run_smoke.sh
```

## Check the saved local result

```matlab
r = QuadrilateralProofRunner('/path/to/Intlab_V12');
r.setup();
summary = r.summarizeLocal();
assert(summary.verified)
```

`verified=true` requires every saved row to have `ok=1` and the minimum saved
lower endpoint for `S` to be positive. The proof decision itself was made in
INTLAB before the decimal endpoint was written to CSV.

## Recompute the local certificate

On one machine:

```bash
export INTLAB_ROOT=/path/to/Intlab_V12
./scripts/run_local_workers.sh 16 12 results/local_new
```

On liulab, use separate INTLAB copies for independent MATLAB processes:

```bash
export INTLAB_ROOT_PATTERN='/home/rendo/Code_Endo/Intlab_Group/Intlab_V12_no%d'
./scripts/run_local_workers.sh 48 12 results/local_new
```

After all `done_*.txt` files appear:

```matlab
r.summarizeLocal(fullfile(r.Root,'results','local_new'))
```

## Recompute the global certificate

```bash
export INTLAB_ROOT=/path/to/Intlab_V12
./scripts/run_global.sh
```

or in MATLAB:

```matlab
result = r.runGlobal();
assert(result.complete)
assert(result.unverified == 0)
assert(result.min_certified_margin > 0)
```

On liulab, split the independent root-box forests across INTLAB copies:

```bash
export INTLAB_ROOT_PATTERN='/home/rendo/Code_Endo/Intlab_Group/Intlab_V12_no%d'
./scripts/run_global_workers.sh 10
```

After `worker_001.json` through `worker_010.json` have appeared:

```matlab
result = qn_merge_global_results('results/global',10,'results/global/summary.json');
assert(result.complete)
```

The cover starts from the paper's bounding cube `P_C`, retains one
representative from each `D_4` orbit, discards boxes proved to lie inside the
local seam ball or outside `P_K`, and applies the following certified test:

```text
Q(B) Lambda(V;B) < pi^2,
Q(B) = sup_B |Q_p|.
```

The center spectral gap selects either the bottom eigenvector or the bottom
two-dimensional frame. The 2-vector route uses the guarded closed form for the
smaller generalized eigenvalue and the source algorithm's 90-direction
Rayleigh sweep when the interval discriminant crosses zero. Failed boxes are
bisected using the source slack-driven coordinate rule; `qn_interval_box`
outward-rounds each child without imposing a minimum box width.

### Boundary-regular stiffness assembly

For the ambient trigonometric trial functions used in the paper,
`grad_ref(psi o Phi) = DPhi' * grad_x(psi)`. Hence the inverse-metric factors
in the generic pullback cancel exactly and

```text
K_ij(p) = integral_square
          (grad_x psi_i . grad_x psi_j)(Phi_p(u,v)) J(u,v;p) du dv.
```

The implementation assembles this form directly. It is algebraically
identical to the source Rayleigh--Ritz pencil in the interior, but contains no
division by `J` and remains regular when one corner value `J=c_i` vanishes on a
triangle face. The Bernstein-ellipse GL pad is correspondingly an
entire-integrand bound; it also encloses the quadrature errors in the
box-uniform parameter gradients used by the mean-value form.

## Paper-to-code map

| paper item | implementation |
|---|---|
| local single-box Algorithm 1 | `src/local/dq2_algorithm1_box.m` |
| local cover and subdivision | `src/local/dq2_run_certificate.m` |
| local Taylor remainder | `src/local/dq2_bound_taylor_remainder.m` |
| Proposition `p:box-bound` | `src/global/qn_certify_box.m` |
| interval pencil `K(B),M(B)` | `src/global/qn_km_enclosure.m` |
| GL truncation enclosure | `src/global/qn_gl_pad.m` |
| Theorem `t:box-cover-terminates` | `src/global/qn_run_global_cover.m` |
| `P_K`, `P_C`, seam and `D_4` logic | `src/global/qn_*box*.m`, `qn_global_constants.m` |

## Trust boundary

- INTLAB outward rounding and the explicit Taylor/GL remainder bounds are
  soundness-critical.
- A local box is accepted only when INTLAB proves `inf(S)>0`.
- A global box is accepted only when INTLAB proves positivity of the mass
  matrix by interval LDL and `sup(Q(B)*Lambda(V;B)) < inf(pi^2)`.
- The 20-by-20 GL quadrature is not treated as exact: every stiffness, raw
  mass, and mean entry, and every parameter derivative used in its mean-value
  enclosure, is widened by the Bernstein-ellipse truncation bound.
- Center eigensolves and finite-difference sensitivities are non-certified,
  but they only select a frame and split coordinate; the box theorem is valid
  for any such choices.

See [docs/PROVENANCE.md](docs/PROVENANCE.md) for source commits, hashes, and the
Arb-to-INTLAB porting record.

## Confidentiality

There is intentionally no open-source license. See [PRIVATE.md](PRIVATE.md).
Keep the GitHub repository private and do not publish releases or mirrors.
