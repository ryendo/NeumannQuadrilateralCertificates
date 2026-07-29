# Global result directory

`summary.json` is created only by a fresh MATLAB/INTLAB/veigs run of
`qn_global_certified_cover`. A complete certificate must have:

```text
complete = true
unverified = 0
min_certified_margin > 0
```

No `veigs`-based `summary.json` is checked in at present. A completed
pre-`veigs` liulab run exists, but its former one-/two-vector acceptance test
does not make it a result for the current source. The supplied Python extract omits
the paper's run-of-record driver. Its unsimplified `1/J` stiffness assembly
reaches degenerate-boundary boxes with `reason = qn:Jacobian`; the current
MATLAB/INTLAB code instead uses the algebraically equivalent, boundary-regular
physical-gradient assembly. A fresh full run of this revised enclosure is
still required. See `PROVENANCE.md` before attempting a paper-count
reproduction.

The historical Python/Arb counts quoted in the main README are reference data,
not a substitute for this MATLAB/INTLAB output.
