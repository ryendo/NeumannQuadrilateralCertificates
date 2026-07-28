# Source provenance

The local certificate was imported from `ryendo/dq2-quadrilateral-certificate`
at commit `6d44991af03102a2338e7af69c67767a32e55178` (2026-06-18). The saved
48-worker result was copied byte-for-byte before the paths were reorganized.

On 2026-07-16 the integrated copy was hardened in two proof-critical places:
Taylor-remainder bounds from the spatial cells are now accumulated as INTLAB
intervals before extracting their upper endpoints, and the Gershgorin row
radii used for `M(te)>0` are summed in interval arithmetic.  The exact 2-by-2
root-identity refinement already used by the imported code is now exposed in
the result diagnostics and documented as a proved refinement of Algorithm 1
in the paper.  The historical saved result predates this hardening and a fresh
48-worker local run is required before replacing it.

The global MATLAB/INTLAB implementation is a line-by-line mathematical port of
the local `quad_neumann_global` Python/Arb extract supplied with the paper on
2026-07-15. Important source SHA-256 values are:

| source | SHA-256 |
|---|---|
| `geometry.py` | `7faf190ed03dc732e7e36b2767d1b559bdb2fc209ef0399eb26e87b3b005649b` |
| `global_step.py` | `bce0b526324e64db0bd07b4805428d9d24f361003600b6aeac88880e7215eeca` |
| `pencil.py` | `7d672d6daae302a4ef163047f9817c3a5dd199da35cf856b9b2143f78bbffcb9` |
| `gl_pad_v2.py` | `9b3d39f751c1db2c9b0455f4f7303f22a14498c3f5585e47689c313ef180f8df` |
| `quad_truncation.py` | `e658c59332007136e10c702404f95a5a15b9e155bf2e9e1c09ab9cfca7e0e0f4` |

The reference layout was `ryendo/DirichletSimplicityClustered` at commit
`9331995405ac7fdbcb97e2a9a3ab26c0c9f6bf2b` (2025-12-14).

## Verified generalized eigensolver

Both certificate paths use the MIT-licensed
[`yuuka-math/veigs`](https://github.com/yuuka-math/veigs) package at commit
`6556d39a0d9819bb172d232062b698aa76e420f6` (2025-12-15). The repository does
not bundle the dependency; `VEIGS_ROOT` or the second
`QuadrilateralProofRunner` constructor argument must point to that exact
checkout.

The wrapper `src/common/qn_veigs_smallest.m` symmetrizes the interval pencil
by taking the componentwise hull with its transpose, calls
`veigs(K,M,1,'sa')`, and rejects results whose returned index information does
not contain index 1. In the local DQ2 path this enclosure is intersected with
the Schur-complement enclosure for `lambda_1`; the DQ2 argument remains
responsible for the double cluster and the sign of `S(t,e)`. In the global
path the upper endpoint returned by `veigs` replaces the former one-/two-vector
Rayleigh acceptance test.

The historical local CSV files and the completed pre-`veigs` liulab global
run were generated before this dependency became mandatory. They remain
useful regression data but are not run-of-record outputs for the present
source.

## Missing global run-of-record source

The supplied `quad_neumann_global/README.md` attributes the paper's global
counts to a development driver named `eigbound_shrink.py` in a parent
`quad_neumann_release/` tree. It also refers to `GLOBAL_STEP_AUDIT.md` and a
`python_flint_audit/` directory. None of those files is present in the supplied
extract, elsewhere in the paper workspace, or in the accessible GitHub source.

This distinction is observable in the cover itself: the supplied
`run_cover(n_init=3)` retains 18 root boxes, whereas the paper reports 74 roots
through `1420 + 337 - 1683 = 74`. On liulab, a faithful MATLAB/INTLAB run of the
18-root extract reached depth 60 with `qn:Jacobian` on boxes approaching the
degenerate boundary of `P_K`. The run was stopped and no global summary was
accepted. `qn_run_global_cover` records the center, half-widths, depth, and
failure reason for every such box so that a supplied run-of-record algorithm
can be compared directly.

## Algebraically regular stiffness representation

The Python source evaluates the generic inverse-metric pullback of the
stiffness form and therefore introduces a factor `1/J`. For the paper's actual
trial space this factor is removable: the functions are restrictions of fixed
ambient trigonometric functions, so

```text
grad_ref(psi o Phi) = DPhi' grad_x(psi)
```

and the inverse metric cancels to give `(grad_x psi_i . grad_x psi_j) J`.
The MATLAB/INTLAB code evaluates this algebraically identical expression
directly. This is not a change to the trial space, Rayleigh--Ritz pencil, box
test, or subdivision algorithm. It removes the source representation's
artificial pole at `c_i=J(corner)=0` and permits an entire-integrand GL
remainder bound on triangle faces.

The boundary-regular smoke test encloses the analytic square pencil, an
interior pencil, and a parameter box crossing the genuine triangle face
`c_3=0`. In the current source the interior and boundary boxes must also be
certified by `veigs` with returned index information containing `1`.

Porting rules:

- The global cover geometry, D4 representative rule, and directional
  subdivision fallback are kept. The former 1-vector and 2-vector acceptance
  tests are replaced by the verified index-1 upper bound from `veigs`. The paper's stated
  longest-side fallback is applied when the maximum-slack axis is already less
  than half the longest half-width; this restores the diameter-to-zero premise
  of the termination proof that the Python extract's bare `argmax(slack)` did
  not enforce.
- Arb intervals are replaced by INTLAB `intval` intervals.
- Decimal proof constants are constructed with `intval('...')`.
- The pencil uses a centered mean-value enclosure. A lightweight explicit
  forward-mode jet evaluates the INTLAB gradient over the complete parameter
  box, so `f(c)+Df(B)(B-c)` includes all parameter dependence directly and
  needs no Taylor remainder; no exposed
  floating derivative is used. Bernstein-ellipse Gauss--Legendre truncation
  padding is added to both the center values and the parameter gradients in
  the mean-value form. The spatial rule uses order 20 because the
  source kernel itself notes that global/far boxes require order at least about
  20; its default call left the order at 12, whose non-vanishing truncation pad
  cannot certify skewed interior points. Floating-point values choose the
  split coordinate only.
