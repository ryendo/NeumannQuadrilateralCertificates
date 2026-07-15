# Global result directory

`summary.json` is created only by a fresh MATLAB/INTLAB run of
`qn_run_global_cover`. A complete certificate must have:

```text
complete = true
unverified = 0
min_certified_margin > 0
```

No `summary.json` is checked in at present. The supplied Python extract omits
the paper's run-of-record driver. Its unsimplified `1/J` stiffness assembly
reaches degenerate-boundary boxes with `reason = qn:Jacobian`; the current
MATLAB/INTLAB code instead uses the algebraically equivalent, boundary-regular
physical-gradient assembly. A fresh full run of this revised enclosure is
still required. See `docs/PROVENANCE.md` before attempting a paper-count
reproduction.

The historical Python/Arb counts quoted in the main README are reference data,
not a substitute for this MATLAB/INTLAB output.
