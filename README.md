# Certified Neumann-eigenvalue bounds for quadrilaterals

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

This is equation (20) of the paper. Equivalently, the final interval sign test
proves equation (27), `S(t,e)>0`, on a finite cover of
`S^3 x (0,rho_local]`, where `p=t e` as in equation (21).

The certified results reported in Appendix B, Table 4 of the paper are:

| quantity | certified saved value |
|---|---:|
| initial products / maximum depth | 13,824 × 9 / 3 |
| certified lower bound for `S` | 0.0039623492860414444 |
| certified bounds for (`lambda_1`, `lambda_2`) | `lambda_1 >= 7.584`, `lambda_2 <= 13.269` |
| outcome | verified on every final box |

The numerical realization in Section 5 records the more precise cluster
bounds `lambda_1 >= 7.5846307638015062` and
`lambda_2 <= 13.268331793903352`. The checked-in CSV files and summary provide
the machine-readable output for the single-box `S>0` test and the index-1
`veigs` check.

### Global step

The global computation certifies

```text
|Q_p| lambda_1(p) < pi^2
for p in P_K \ B(0,rho_seam),

rho_seam = rho#/2 = 1616/(27*pi^6) approximately 0.062256.
```

This is equation (16) of the paper. Because the trial space is contained in
the zero-mean Sobolev space,
`lambda_1(p) >= mu_1(Q_p)`, so the same strict inequality follows for the
first nonzero Neumann eigenvalue. The local and global regions overlap on
`rho#/2 <= ||p|| <= rho#`.

The source Python code called the seam `RHO_SHARP`; this repository uses the
paper-aligned names `rho_local` and `rho_seam` to remove that ambiguity.

The certified output reported in Appendix C, equation (58) of the paper is:

| quantity | certified saved value |
|---|---:|
| initial boxes / maximum depth | 18 / 30 |
| accepted / discarded / unresolved boxes | 166,928 / 25,285 / 0 |
| bisections | 192,195 |
| `Delta_*` | 2.6545819961754091e-5 |
| conclusion | `|Q_p| lambda_1(p) < pi^2` on `Omega_II` |

Every accepted box is certified by `veigs`, and no box is left unresolved.
The checked-in worker JSON files provide the machine-readable output for this
finite cover.

## Core Libraries & Dependencies

Both components use MATLAB, INTLAB, and the verified generalized-eigenvalue
solver `veigs`.

## Project structure

```text
.
├── QuadrilateralProofRunner.m       unified entry point
├── src/                             all proof implementation
│   ├── qn_single_box_certificate.m  paper Algorithm 1
│   ├── qn_global_certified_cover.m  paper Appendix C finite cover
│   ├── qn_*.m                       certificate orchestration/shared code
│   └── dq2_*.m                      low-level local DQ2 kernels
├── data/taylor_coefficients.mat
├── results/
│   ├── local/                       current local veigs run
│   └── global/                      current global veigs run
├── tests/                           smoke and regression tests
├── scripts/                         shell launchers
└── PROVENANCE.md
```

## Setup and quick check

Install `veigs`:

```bash
git clone https://github.com/yuuka-math/veigs.git /path/to/veigs
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
Taylor-remainder accumulation, the coarse single-box bound, and the exact
root-identity refinement:

```matlab
report = qn_local_rigour_test(r.Root);
global_fast = qn_global_vectorization_test(r.Root);
local_fast = dq2_vectorization_test(r.Root);
remainder_fast = dq2_remainder_vectorization_test(r.Root);
```

The vectorization regressions check that box evaluations contain their center
values and center Hessians, that Taylor-cell subdivision does not enlarge the
remainder bounds, and that the global assembly certifies a box crossing
`c_3=0`.

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

To use separate INTLAB copies for independent MATLAB processes:

```bash
export INTLAB_ROOT_PATTERN='/path/to/Intlab_V12_no%d'
export VEIGS_ROOT=/path/to/veigs
./scripts/run_local_workers.sh 40 12 results/local_new
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
result = r.runGlobalCertifiedCover();
assert(result.complete)
assert(result.unverified == 0)
assert(result.min_certified_margin > 0)
```

To split the independent root-box forests across INTLAB copies:

```bash
export INTLAB_ROOT_PATTERN='/path/to/Intlab_V12_no%d'
export VEIGS_ROOT=/path/to/veigs
./scripts/run_global_workers.sh 18
```

After `worker_001.json` through `worker_018.json` have appeared:

```matlab
result = qn_merge_global_results('results/global',18,'results/global/summary.json');
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
`sup(Q(B))*sup(lambda_1(B)) < inf(pi^2)`, implementing the conditions in
equation (52) and the conclusion in equation (53). Floating center
eigensolves are used only by the source slack-driven bisection heuristic.
`qn_interval_box` outward-rounds each child without imposing a minimum box
width.

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
On a representative box, the complete warm per-box test decreased from about
25.2 seconds to 0.24 seconds during development. The former scalar reference
implementation has since been removed.

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
| Algorithm 1 (Appendix B.4) | `src/qn_single_box_certificate.m` |
| local cover and subdivision | `src/qn_local_certificate_cover.m` |
| local Taylor remainder | `src/dq2_bound_taylor_remainder_vectorized.m` |
| verified index-1 generalized eigenvalue | `src/qn_veigs_smallest.m` |
| Proposition 6.1 | `src/qn_certify_box.m` |
| interval pencil `K(B),M(B)` | `src/qn_km_enclosure.m` |
| GL truncation enclosure | `src/qn_gl_pad.m` |
| Appendix C finite cover, Theorem 6.2 | `src/qn_global_certified_cover.m` |
| `P_K`, `P_C`, seam and `D_4` logic | `src/qn_*box*.m`, `src/qn_global_constants.m` |

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

See [PROVENANCE.md](PROVENANCE.md) for source provenance and the Arb-to-INTLAB
porting record.

## Confidentiality

There is intentionally no open-source license. See [PRIVATE.md](PRIVATE.md).
Keep the GitHub repository private and do not publish releases or mirrors.
