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

Both components use MATLAB, INTLAB, and the verified generalized-eigenvalue
solver `veigs`. Python is not required.

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

This saved run predates the 2026-07-16 rounding hardening and the mandatory
`veigs` index-1 check. The hardening keeps
the same mathematical algorithm but accumulates every Taylor-remainder cell
bound and every Gershgorin row radius with outward-rounded INTLAB operations.
It also records the coarse `eta/beta` bound separately from the exact
root-identity refinement.  Until the scheduled full rerun finishes, the table
above is retained as historical reference data rather than a refreshed result
of the hardened, `veigs`-checked code. A fresh complete local run is required
before the table can be promoted to the current run of record.

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

> **Recalculation status.** No `veigs`-based full-cover result is checked in
> yet. A complete pre-`veigs` MATLAB/INTLAB run exists on liulab, but its
> per-box acceptance used the former one-/two-vector Rayleigh tests and is not
> a run of record for the present code. The supplied Python
> extract is not the run-of-record driver named in its own README: it omits
> `eigbound_shrink.py`, the parent `quad_neumann_release/` tree, and
> `GLOBAL_STEP_AUDIT.md`. A faithful run of the supplied `run_cover(n_init=3)`
> reaches `qn:Jacobian` on boxes converging to the degenerate boundary of
> `P_K`. The current MATLAB/INTLAB assembly removes that artificial obstruction
> by using the exactly equivalent physical-gradient formula described below.
> On liulab, the revised smoke test certifies representative interior and
> triangle-boundary boxes using `veigs` with index 1. A new full global run is
> required before the current global certificate can be declared complete.

The reference global run reported 1,420 verified boxes, 337 discarded boxes,
1,683 bisections, no unverified boxes, maximum depth 26, and positive certified
slack. A fresh MATLAB/INTLAB run writes its own record to
`results/global/summary.json`; do not substitute the reference counts for a
new run's output.

## Core Libraries & Dependencies

This project relies on specialized libraries for verified numerical
computation:

1. **MATLAB** R2023b or later (the local computation was tested with R2023b).
2. **INTLAB**: The fundamental toolbox for rigorous interval arithmetic in
   MATLAB. INTLAB is not bundled.
   * **Source:** [http://www.tuhh.de/ti3/intlab/](http://www.tuhh.de/ti3/intlab/)
     [INTLAB V12 was used for the computation.]
3. **veigs**: Used for solving generalized matrix eigenvalue problems with
   rigorous error bounds with the information of indices. `veigs` is not
   bundled.
   * **Source:** [https://github.com/yuuka-math/veigs](https://github.com/yuuka-math/veigs)
     [2025/12/13]
   * **Pinned revision:** `6556d39a0d9819bb172d232062b698aa76e420f6`
4. A POSIX shell, required only for the optional launch scripts.

All non-integer proof constants and decimal bounds are parsed from strings,
for example `intval('0.5')`, `intval('pi')`, and `intval('1e-12')`. Floating
point is used only for the non-certified bisection heuristic. Every
accept/reject proof decision is made from INTLAB interval endpoints.

## Project structure

```text
.
├── QuadrilateralProofRunner.m       unified entry point
├── src/
│   ├── common/
│   │   └── qn_veigs_smallest.m      verified index-1 eigensolver wrapper
│   ├── local/                       DQ2 local certificate
│   │   ├── dq2_run_certificate.m
│   │   ├── dq2_algorithm1_box.m
│   │   ├── dq2_algorithm1_scalars.m
│   │   ├── dq2_evaluate_taylor_coefficients_vectorized.m
│   │   ├── dq2_bound_taylor_remainder_vectorized.m
│   │   └── qn_summarize_local_results.m
│   └── global/                      adaptive global certificate
│       ├── qn_run_global_cover.m
│       ├── qn_certify_box.m
│       ├── qn_km_enclosure.m
│       ├── qn_assemble_interval_center.m
│       ├── qn_assemble_interval_grad.m
│       └── qn_gl_pad.m
├── data/local/taylor_coefficients.mat
├── results/
│   ├── local/                       saved 48-worker certified run
│   └── global/                      output of a fresh global run
├── tests/qn_smoke_test.m
├── tests/qn_global_vectorization_test.m
├── tests/dq2_vectorization_test.m
├── tests/dq2_remainder_vectorization_test.m
├── scripts/
│   ├── matlab_runner.sh
│   ├── run_local_workers.sh
│   ├── run_global.sh
│   ├── run_global_workers.sh
│   └── run_smoke.sh
└── docs/PROVENANCE.md
```

## Setup and quick check

Install the pinned `veigs` revision:

```bash
git clone https://github.com/yuuka-math/veigs.git /path/to/veigs
git -C /path/to/veigs checkout 6556d39a0d9819bb172d232062b698aa76e420f6
```

In MATLAB:

```matlab
r = QuadrilateralProofRunner('/path/to/Intlab_V12','/path/to/veigs');
r.setup();
report = r.smokeTest();
```

The smoke test checks that the certified global pencil enclosure contains the
floating-point center value, exercises a representative per-box test, and
verifies the saved local CSV set.

The focused local proof-path regression test additionally exercises rigorous
Taylor-remainder accumulation, the coarse Algorithm-1 bound, and the exact
root-identity refinement:

```matlab
report = qn_local_rigour_test(r.Root);
global_fast = qn_global_vectorization_test(r.Root);
local_fast = dq2_vectorization_test(r.Root);
remainder_fast = dq2_remainder_vectorization_test(r.Root);
```

The vectorization regressions retain the former scalar INTLAB kernels as
independent references. They compare every stiffness, raw-mass, mean, gradient,
and local Hessian interval, and exercise a global box crossing `c_3=0`.

From a shell:

```bash
export INTLAB_ROOT=/path/to/Intlab_V12
export VEIGS_ROOT=/path/to/veigs
./scripts/run_smoke.sh
```

## Check the saved local result

```matlab
r = QuadrilateralProofRunner('/path/to/Intlab_V12','/path/to/veigs');
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
export VEIGS_ROOT=/path/to/veigs
./scripts/run_local_workers.sh 16 12 results/local_new
```

On liulab, use separate INTLAB copies for independent MATLAB processes:

```bash
export INTLAB_ROOT_PATTERN='/home/rendo/Code_Endo/Intlab_Group/Intlab_V12_no%d'
export VEIGS_ROOT=/path/to/veigs
./scripts/run_local_workers.sh 48 12 results/local_new
```

After all `done_*.txt` files appear:

```matlab
r.summarizeLocal(fullfile(r.Root,'results','local_new'))
```

## Recompute the global certificate

```bash
export INTLAB_ROOT=/path/to/Intlab_V12
export VEIGS_ROOT=/path/to/veigs
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
export VEIGS_ROOT=/path/to/veigs
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
Q(B) sup(lambda_1(B)) < pi^2,
Q(B) = sup_B |Q_p|.
```

For every retained box, `veigs(K(B),M(B),1,'sa')` returns a rigorous interval
for the generalized eigenvalue with index 1. A box is accepted only when the
returned index data contains `1` and the upper endpoint proves
`sup(Q(B))*sup(lambda_1(B)) < inf(pi^2)`. Floating center eigensolves are used
only by the source slack-driven bisection heuristic. `qn_interval_box`
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

### Proof-preserving vectorization

The global assembly evaluates the same 20-by-20 interval GL rule as the scalar
reference, but batches its 400 nodes in INTLAB arrays. Parameter derivatives
are the explicit chain-rule formulas for `X`, `Y`, `J`, the five modes, and
their physical gradients; no derivative is approximated by floating point.
On a representative liulab box, the complete warm per-box test decreased from
about 25.2 seconds to 0.24 seconds. The scalar implementation remains in the
repository solely for interval-overlap regression tests.

The local code similarly evaluates the existing exact monomial coefficient
table by exponent groups and batches the four Taylor-remainder integration
cells. Generic Hessian `0^0` is deliberately avoided: exponent zero is assigned
as the exact constant one and positive powers are built recursively. On the
focused local regression, the certified lower bound increased slightly while
the first-call wall time decreased from about 17.3 seconds to 8.3 seconds.
These timing values are diagnostics, not proof inputs.

## Paper-to-code map

| paper item | implementation |
|---|---|
| local single-box Algorithm 1 | `src/local/dq2_algorithm1_box.m` |
| local cover and subdivision | `src/local/dq2_run_certificate.m` |
| local Taylor remainder | `src/local/dq2_bound_taylor_remainder_vectorized.m` |
| verified index-1 generalized eigenvalue | `src/common/qn_veigs_smallest.m` |
| Proposition `p:box-bound` | `src/global/qn_certify_box.m` |
| interval pencil `K(B),M(B)` | `src/global/qn_km_enclosure.m` |
| GL truncation enclosure | `src/global/qn_gl_pad.m` |
| Theorem `t:box-cover-terminates` | `src/global/qn_run_global_cover.m` |
| `P_K`, `P_C`, seam and `D_4` logic | `src/global/qn_*box*.m`, `qn_global_constants.m` |

## Trust boundary

- INTLAB outward rounding and the explicit Taylor/GL remainder bounds are
  soundness-critical.
- A local box is accepted only when INTLAB proves `inf(S)>0` and `veigs`
  returns a certified enclosure whose index data contains eigenvalue 1. The
  DQ2 and `veigs` enclosures for `lambda_1` are intersected.
- Cellwise Taylor-remainder magnitudes and Gershgorin row radii are summed as
  intervals; floating-point endpoint sums are not used as certified bounds.
- Vectorized kernels use only INTLAB operations for proof quantities. Their
  scalar reference implementations remain covered by componentwise interval
  overlap/subset tests.
- The coarse `eta/beta` lower bound and the exact root-identity improvement are
  both retained in the per-box result structure for audit.
- A global box is accepted only when INTLAB proves positivity of the mass
  matrix by interval LDL, `veigs` certifies index 1, and
  `sup(Q(B))*sup(lambda_1(B)) < inf(pi^2)`.
- The 20-by-20 GL quadrature is not treated as exact: every stiffness, raw
  mass, and mean entry, and every parameter derivative used in its mean-value
  enclosure, is widened by the Bernstein-ellipse truncation bound.
- Center eigensolves and finite-difference sensitivities are non-certified,
  but they only select the split coordinate and do not enter acceptance.

See [docs/PROVENANCE.md](docs/PROVENANCE.md) for source commits, hashes, and the
Arb-to-INTLAB porting record.

## Confidentiality

There is intentionally no open-source license. See [PRIVATE.md](PRIVATE.md).
Keep the GitHub repository private and do not publish releases or mirrors.
