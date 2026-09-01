# Computer-Assisted Proof for the Neumann Eigenvalue Problem on Quadrilaterals

This repository contains the MATLAB source code and saved computation summaries
for the paper

> Ryoki Endo and Braxton Osting, *Maximizing the fundamental Laplace--Neumann
> eigenvalue on quadrilaterals* (2026).

The equation, theorem, algorithm, and table numbers below refer to the current
manuscript. The paper-to-code table also gives the corresponding LaTeX labels
so that each reference remains identifiable when the manuscript is revised.

## Background

For the quadrilateral parameter

$$
\mathbf p=(a,b,c,d), \qquad Q_{\mathbf 0}=\square,
$$

the area is

$$
q(\mathbf p)=|Q_{\mathbf p}|=1-a^2-d^2
$$

by equation (2), `e:Area`. Sections 3 and 4.3 define the five-dimensional
Rayleigh--Ritz pencil

$$
K(\mathbf p)x=\lambda M(\mathbf p)x,
\qquad
\lambda_1(\mathbf p)\le\cdots\le\lambda_5(\mathbf p),
$$

using integrals over the fixed reference square. The Rayleigh--Ritz comparison
in equation (12), `e:KeyInequality`, reduces the proof of Theorem 1.1,
`t:Main`, to two certified computations.

1. **Step (I), equation (21), `eq:local-statement`.** For

   <div align="center">

   $`\displaystyle
   \rho^\sharp=\frac{3232}{27\pi^6},
   `$

   </div>

   certify

   <div align="center">

   $`\displaystyle
   f_2(\mathbf p)>\pi^{-2}
   \qquad
   (0<\|\mathbf p\|\le\rho^\sharp).
   `$

   </div>

   With $\mathbf p=t\mathbf e$, this is reduced by equations (22)--(26) to
   the inequality $S(t,\mathbf e)>0$. Appendix A, `app:diff-quotient`, gives
   the interval conditions (34)--(35), the enclosures (36), the refinement formulas
   (38)--(39), and the finite cover (40)--(42). Algorithm 1 in Appendix A.4 is
   implemented by `src/qn_single_box_certificate.m`.

2. **Step (II), equation (20), `e:global-target`.** On

   <div align="center">

   $`\displaystyle
   \Omega_{\mathrm{II}}
   =\mathcal P_{\mathrm K}
   \setminus B\!\left(\mathbf 0,\frac{\rho^\sharp}{2}\right),
   `$

   </div>

   certify

   <div align="center">

   $`\displaystyle
   |Q_{\mathbf p}|\lambda_1(\mathbf p)<\pi^2.
   `$

   </div>

   Proposition 6.1, `p:box-bound`, gives the boxwise implication. Appendix B,
   `app:implementation`,
   defines the interval box in (43), the accepted-box conditions and conclusion
   in (44)--(45), the discard and bisection rules in (46)--(47), and the finite
   completion condition in (48). Theorem 6.2, `t:box-cover-terminates`, then
   yields the required inequality on $\Omega_{\mathrm{II}}$.

The two computations overlap on

$$
\frac{\rho^\sharp}{2}\le\|\mathbf p\|\le\rho^\sharp.
$$

The code uses the paper notation through the fields `rho_sharp` and
`rho_sharp_over_2`.

## Core Libraries & Dependencies

This project relies on specialized libraries for verified numerical computation:

1. **INTLAB**: The fundamental toolbox for rigorous interval arithmetic in MATLAB.
   - **Source:** [http://www.tuhh.de/ti3/intlab/](http://www.tuhh.de/ti3/intlab/) [INTLAB_V12, INTLAB_V14 were used for the computation.]
2. Revised version of **VFEM2D**: Used for rigorous finite element matrix assembly and high-precision eigenvalue bounds (Lehmann–Goerisch method).
   - **Source:** [https://github.com/xfliu/VFEM2D](https://github.com/xfliu/VFEM2D) [2025/12/13]
3. **veigs**: Used for solving generalized matrix eigenvalue problems with rigorous error bounds with the information of indices.
   - **Source:** [https://github.com/yuuka-math/veigs](https://github.com/yuuka-math/veigs) [2025/12/13]

The source code in this repository directly calls INTLAB and `veigs`. It does
not call VFEM2D: the matrices $K(\mathbf p)$ and $M(\mathbf p)$ are assembled
from the fixed-reference-square integrals of Section 4.3. VFEM2D is listed as a
related verified finite-element library, not as an executable dependency of
this repository.

The saved computations recorded in `PROVENANCE.md` used MATLAB R2023b,
INTLAB V12, and `veigs` at commit
`6556d39a0d9819bb172d232062b698aa76e420f6`.

## Project Structure

```text
.
├── src/
│   ├── qn_single_box_certificate.m   # Algorithm 1, Appendix A.4
│   ├── qn_single_box_blocks.m        # K, M, and C_t blocks in (27)--(28)
│   ├── qn_single_box_quantities.m    # F_t and coefficients in (28), (31)--(32)
│   ├── qn_exact_root_refinement_corrections.m # corrections in (39)
│   ├── qn_local_certificate_cover.m  # finite cover in (40)--(42)
│   ├── qn_certify_box.m              # certified per-box implementation of (43)--(45)
│   ├── qn_global_code_test.m          # all five accepted-box conditions in (44)
│   ├── qn_global_certified_cover.m   # finite cover in (46)--(49)
│   ├── qn_km_enclosure.m             # interval matrices in Appendix B
│   ├── qn_km_float.m                 # coordinate choice before (47), not a proof bound
│   ├── qn_veigs_indices.m            # verified eigenvalue indices
│   ├── qn_quadrilateral_kernels.m    # integrands of Section 4.3
│   └── dq2_*.m                       # Taylor coefficients used in Appendix A
├── data/
│   └── taylor_coefficients.mat       # coefficients for equation (17)
├── results/
│   ├── local/                        # per-box output from the recorded Algorithm 1 run
│   └── global/                       # aggregate summaries from the recorded global run
├── tests/                            # representative and saved-output checks
├── scripts/                          # MATLAB launch scripts
├── PROVENANCE.md                     # software versions and run information
└── CITATION.cff
```

## Installation & Setup

### 1. Clone the repository

```bash
git clone https://github.com/ryendo/NeumannQuadrilateralCertificates.git
cd NeumannQuadrilateralCertificates
```

### 2. Install INTLAB and `veigs`

Install INTLAB from the source listed above. Clone the `veigs` version used by
the saved computation:

```bash
git clone https://github.com/yuuka-math/veigs.git /path/to/veigs
git -C /path/to/veigs checkout 6556d39a0d9819bb172d232062b698aa76e420f6
```

### 3. Initialize MATLAB

```matlab
addpath('src');
repo_root = qn_setup('/path/to/Intlab_V12', '/path/to/veigs');
```

`qn_setup` starts INTLAB, adds `src/` and `tests/` to the MATLAB path, and
checks that both `veigs` and `veig` are available.

## Usage

### 1. Representative verification

Run the representative local and global checks and verify the saved local
summary:

```matlab
report = qn_smoke_test(repo_root);
assert(report.ok)
```

From a shell:

```bash
export INTLAB_ROOT=/path/to/Intlab_V12
export VEIGS_ROOT=/path/to/veigs
./scripts/run_smoke.sh
```

This test does not recompute either finite cover.

### 2. Step (I): equations (21)--(42) and Algorithm 1

For each product box defined by (40)--(41), Algorithm 1 verifies the
conditions (34)--(35), computes (36), optionally intersects it with
(38)--(39), and substitutes the resulting intervals into (23). Acceptance
requires the inequality (26). The implementation also obtains verified
enclosures of generalized-eigenvalue indices 1 and 3 and requires
$\inf[\lambda_3]>16$, as reported in Table 3.

#### Run the complete finite cover

```bash
export INTLAB_ROOT=/path/to/Intlab_V12
export VEIGS_ROOT=/path/to/veigs
./scripts/run_local_workers.sh 16 12 results/local_new
```

The arguments are the number of MATLAB processes, the number of subintervals
per coordinate in (40), and the output directory. To use separate INTLAB
copies for simultaneous MATLAB processes:

```bash
export INTLAB_ROOT_PATTERN='/path/to/Intlab_V12_no%d'
export VEIGS_ROOT=/path/to/veigs
./scripts/run_local_workers.sh 40 12 results/local_new
```

After all processes finish, verify the output:

```matlab
summary = qn_summarize_local_results( ...
    fullfile(repo_root, 'results', 'local_new'));
assert(summary.verified)
```

#### Check the saved local results

```matlab
summary = qn_summarize_local_results( ...
    fullfile(repo_root, 'results', 'local'));
assert(summary.verified)
```

The saved output gives

```text
direction boxes obtained from (40) / intervals in (41):  13,824 / 9
maximum subdivision depth:             1
inf S:                                 0.0044848550155087707
inf lambda_1:                          7.6300046469977758
sup lambda_2:                          13.267374191855003
inf lambda_min(M):                     0.052100929203335616
inf lambda_3:                          17.834653422811055
inf (lambda_3 - 3*pi^2/2):             3.0302468211770126
uncertified products:                  0
```

These are the unrounded values underlying the rounded entries in Table 3.
They were generated from source commit
`b9944f600ba63bed3ddbbde9890febebb613b3ac`, as recorded in
`results/local/RUN_PROVENANCE.txt`. The assertion above audits the saved CSV
rows and completion markers; it does not claim that they were regenerated by
the current working tree. The current implementation uses the literal
factor $t^2$ in $-\partial_\nu F_t$ from (34) and the literal
$\widehat\nu_i$ terms in (39), so a new `results/local_new` run is required
to produce results with the current source code. This does not invalidate the saved local
proof: its derivative test implies the condition in (34) for
$0\le t<1$, and its earlier implementation of (39) used a wider interval.

### 3. Step (II): Proposition 6.1, Theorem 6.2, and equations (43)--(49)

For a box $\mathcal B$ in (43), the code encloses
$`[K]_{\mathcal B}`$, $`[M]_{\mathcal B}`$, and $`[q]_{\mathcal B}`$,
requires $\underline q_{\mathcal B}>0$, certifies
$`[M]_{\mathcal B}\succ0`$, and obtains a verified enclosure
$`[\lambda]_{\mathcal B}`$ whose index data contains 1. It also requires
$\underline\lambda_{\mathcal B}>0$ and evaluates

$$
\Delta_{\mathcal B}
=\underline{\pi^2}
-\overline q_{\mathcal B}\,\overline\lambda_{\mathcal B},
$$

as in equation (44). These five conditions give both inequalities in (45).
The remaining boxes are treated by (46)--(47), and completion is the condition
(48).

The coordinate used in (47) is selected from floating-point changes in the
center pencil. If that coordinate is too short, the code selects a longest
side, ensuring $h_r\ge\tfrac12\max_s h_s$. This choice does not enter the
interval proof; every child box is tested independently.

Because the stored centers and half-widths are binary floating-point numbers,
`qn_bisect_box.m` rounds the children outward. The implemented relation is the
conservative covering
$\mathcal B\subseteq\mathcal B^-\cup\mathcal B^+$; the children can extend
by a few rounding units beyond the exact halves displayed in (47).

#### Run the complete finite cover

The complete computation is intended for multiple MATLAB processes. With one
INTLAB copy per process, write to a new directory:

```bash
export INTLAB_ROOT_PATTERN='/path/to/Intlab_V12_no%d'
export VEIGS_ROOT=/path/to/veigs
./scripts/run_global_workers.sh 16 results/global_new
```

After all workers finish, merge the summaries and test (48):

```matlab
result = qn_merge_global_results( ...
    fullfile(repo_root, 'results', 'global_new'), 16, ...
    fullfile(repo_root, 'results', 'global_new', 'summary.json'));
assert(result.complete)
```

For a single MATLAB process and a new output directory:

```matlab
output_dir = fullfile(repo_root, 'results', 'global_new');
if ~isfolder(output_dir), mkdir(output_dir); end
result = qn_global_certified_cover( ...
    3, 60, true, fullfile(output_dir, 'summary.json'));
assert(result.complete)
assert(result.unverified == 0)
assert(result.q_lower_min > 0)
assert(result.lambda_lower_min > 0)
assert(result.delta_star_lower > 0)
```

Here `3` is the subdivision count in each coordinate of the initial family
$\mathcal G_0$ in Appendix B, and `60` is the maximum subdivision depth.
The JSON files are aggregate run summaries. They record counts and certified
minima, but not every accepted and discarded leaf box; reproducing the full
finite cover therefore requires rerunning the code. The current schema also
records each worker's assigned initial-box IDs; the merge rejects missing or
duplicate assignments before reporting `complete=true`.

#### Verify the saved output

```matlab
result = jsondecode(fileread( ...
    fullfile(repo_root, 'results', 'global', 'summary.json')));
assert(result.unverified == 0)
assert(result.min_certified_margin > 0)
```

These two assertions check only the aggregate fields recorded by the saved
run. They do not upgrade the old JSON schema to the current five-condition
test in (44), or establish the finite cover in (49).

The saved output gives

```text
initial boxes retained after (46) and D_4 reduction:  16
accepted boxes:                                       117,083
discarded boxes:                                      15,222
unresolved boxes:                                     0
bisections:                                           132,289
maximum subdivision depth:                            29
Delta_* in (48):                                      8.0698658297961856e-6
```

The saved summary records $\mathcal U=\varnothing$ and $\Delta_*>0$, the
two conditions in (48). Its box counts come from the run recorded in
`results/global/RUN_PROVENANCE.txt`. That run predates the explicit
$\underline\lambda_{\mathcal B}>0$ branch now implemented for (44), and its
JSON files do not store the lower eigenvalue bound for every accepted box.
It also used a binary-floating-point bisection that can leave a rounding-size
gap between the two children in (47). Because the leaf boxes were not saved,
the old subdivision cannot be re-audited. Consequently, these files are
historical run summaries, not a current finite-cover certificate; the global
calculation must be recomputed with the present code.

The current manuscript presently reports the following values from an earlier
run:

```text
initial active boxes after (46) and D_4 reduction:  18
accepted boxes:                                     166,928
discarded boxes:                                    25,285
unresolved boxes:                                   0
bisections:                                         192,195
maximum subdivision depth:                          30
Delta_* in (48):                                    at least 2.654e-5
```

That earlier run used an obsolete $D_4$-representative rule. On the
$3^4$ initial grid, the rule selected 18 boxes but omitted one of the 16
distinct orbits and selected three other orbits twice. The omitted orbit meets
$\mathcal P_{\mathrm K}\setminus B(0,\rho^\sharp/2)$, so those 18-box values
must not be used as a finite-cover certificate. The current rule selects
exactly one representative from every orbit; `qn_smoke_test.m` checks this and
the resulting count of 16. The checked-in 16-box run uses this corrected rule,
but a new run is still required to record every condition in the current (44)
schema.

### 4. Verify the saved files

The SHA-256 manifests cover the saved numerical output and run information:

```bash
(cd results/local && sha256sum -c SHA256SUMS)
(cd results/global && sha256sum -c SHA256SUMS)
```

## Paper-to-Code Correspondence

| Paper item | Formula or condition | Implementation | Saved output or check |
|---|---|---|---|
| Section 4.3, `app:integral-rep` | Fixed-reference-square integrals for $K(\mathbf p)$ and $M(\mathbf p)$ | `src/qn_quadrilateral_kernels.m`, `src/qn_assemble_interval.m` | `tests/qn_smoke_test.m` |
| Equation (17), `e:KM-expansion` | Taylor coefficients at $\mathbf p=\mathbf 0$ | `data/taylor_coefficients.mat`, `src/dq2_load_taylor_coefficients.m` | `tests/dq2_vectorization_test.m` |
| Equations (21)--(26), `eq:local-statement`--`eq:S-positive-target` | $f_2(\mathbf p)>\pi^{-2}$ reduced to $S(t,\mathbf e)>0$ | `src/qn_single_box_certificate.m`, `src/qn_local_certificate_cover.m` | `results/local/summary.json` |
| Equation (27), `eq:Ct-definition` | $C_t(\nu)$ and its Taylor enclosure | `src/qn_single_box_blocks.m`, `src/dq2_evaluate_taylor_coefficients_vectorized.m`, `src/dq2_bound_taylor_remainder_vectorized.m` | `tests/qn_local_rigour_test.m` |
| Equation (28), `eq:Ft-definition`; equations (31)--(36), `eq:d1-continuous-extension`--`eq:symmetric-single-box-enclosures` | $F_t$, $\widehat F_t$, and the enclosures of $L$ and $\nu_1\nu_2$ | `src/qn_single_box_quantities.m`, `src/qn_single_box_certificate.m` | `tests/qn_local_rigour_test.m` |
| Equations (38)--(39), `eq:exact-root-refinement`, `eq:refined-symmetric-quantities` | Refinement when $`[\nu_1]\cap[\nu_2]=\varnothing`$ | `src/qn_single_box_certificate.m`, `src/qn_exact_root_refinement_corrections.m` | `tests/qn_local_rigour_test.m`, `tests/qn_smoke_test.m` |
| Equations (40)--(42), `eq:sphere-chart`--`eq:subdivision-inclusion` | Cover of $S^3\times[0,\rho^\sharp]$ and subdivision | `src/dq2_face_direction.m`, `src/qn_local_certificate_cover.m` | `results/local/` |
| Algorithm 1, `alg:single-box-certificate` | Conditions (34)--(35), equations (36), (38)--(39), and $S>0$ | `src/qn_single_box_certificate.m` | Table 3, `tab:local-second-order-certificate`; saved run at source commit `b9944f6` |
| Proposition 6.1, `p:box-bound`; equations (43)--(45), `eq:global-box-notation`--`eq:global-accepted-box-conclusion` | The five conditions in (44) and the conclusion in (45) | `src/qn_km_enclosure.m`, `src/qn_interval_ldl_pd.m`, `src/qn_certify_box.m`, `src/qn_global_code_test.m`, `src/qn_veigs_indices.m` | `tests/qn_smoke_test.m`; current global recomputation required |
| Equations (46)--(47), `eq:global-discard-test`, `eq:global-bisection` | Discard and bisection rules | `src/qn_box_inside_ball.m`, `src/qn_box_outside_pk.m`, `src/qn_bisect_box.m` | `tests/qn_smoke_test.m`; current global recomputation required |
| Equation (48), `eq:global-completion-test` | $\mathcal U=\varnothing$ and $\Delta_*>0$ | `src/qn_global_certified_cover.m`, `src/qn_merge_global_results.m` | `tests/qn_smoke_test.m`; saved global summary uses the earlier schema |
| Equation (49), `eq:global-final-cover`; Theorem 6.2, `t:box-cover-terminates` | Finite cover of $\Omega_{\mathrm{II}}$ | `src/qn_global_certified_cover.m` | current recomputation required; saved JSON files are aggregate summaries |

The current `qn_certify_box.m` tests all five conditions displayed in (44).
In particular, `qn_km_enclosure.m` checks
$\underline q_{\mathcal B}>0$ before forming $`[M]_{\mathcal B}`$, and
`qn_global_code_test.m` separately checks
$\underline\lambda_{\mathcal B}>0$ before testing
$\Delta_{\mathcal B}>0$.

## Citation and Provenance

Use `CITATION.cff` to cite this repository and the accompanying paper. Exact
software versions, source commits, worker counts, hosts, and output hashes are
recorded in `PROVENANCE.md` and in the `RUN_PROVENANCE.txt` files under
`results/local/` and `results/global/`.
