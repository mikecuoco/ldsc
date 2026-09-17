# ldsc-rs Python API

This package exposes the computation-only Rust APIs for LD-score calculation,
heritability, and genetic correlation. It does not invoke the `ldsc` CLI or
parse its output.

The first release intentionally accepts in-memory, already harmonised columns.
File-oriented `munge_sumstats`, `estimate_h2`, and `estimate_rg` workflows will
be added after their CLI orchestration is separated from formatting and file
writing in the Rust crate.

```python
from ldsc_rs import fit_h2

result = fit_h2(
    chi2,
    ref_ld,
    weight_ld,
    sample_size,
    m_snps=1_190_321,
)
print(result.heritability.value, result.heritability.standard_error)
```

`fit_h2` follows the scalar CLI default and uses the two-step estimator at a
cutoff of 30. Its default `two_step="auto"` disables two-step estimation when
an intercept is constrained; pass a number or `None` to override that choice.
`fit_rg` follows the `rg` CLI default and does not enable two-step estimation
automatically. A constrained intercept is incompatible with two-step
estimation.

Inputs to `fit_h2` and `fit_rg` are equal-length finite numeric sequences.
For `fit_rg`, SNP joining and allele harmonisation must happen before the call;
reference and weight LD scores must be non-negative. A constrained intercept's
standard error is represented as `None` in Python.

This initial API is intended for small and medium in-memory integrations.
Python sequences and BED bytes are copied into Rust-owned buffers. For
genome-scale jobs, use the existing CLI until a path-based or zero-copy array
API is available.
