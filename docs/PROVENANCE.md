# Source provenance

The local certificate was imported from `ryendo/dq2-quadrilateral-certificate`
at commit `6d44991af03102a2338e7af69c67767a32e55178` (2026-06-18). The saved
48-worker result was copied byte-for-byte before the paths were reorganized.

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

Porting rules:

- The global cover geometry, D4 representative rule, gap threshold, 1-vector
  and 2-vector tests, and directional fallback are kept. The paper's stated
  longest-side fallback is applied when the maximum-slack axis is already less
  than half the longest half-width; this restores the diameter-to-zero premise
  of the termination proof that the Python extract's bare `argmax(slack)` did
  not enforce.
- Arb intervals are replaced by INTLAB `intval` intervals.
- Decimal proof constants are constructed with `intval('...')`.
- The pencil uses a centered mean-value enclosure. A lightweight explicit
  forward-mode jet evaluates the INTLAB gradient over the complete parameter
  box, so `f(c)+Df(B)(B-c)` includes all
  parameter dependence directly and needs no Taylor remainder; no exposed
  floating derivative is used. The same Bernstein-ellipse Gauss--Legendre
  truncation padding is added. The spatial rule uses order 20 because the source kernel itself notes
  that global/far boxes require order at least about 20; its default call left
  the order at 12, whose non-vanishing truncation pad cannot certify skewed
  interior points. Floating-point values choose frames and split coordinates
  only.
