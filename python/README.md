# ldsc-rs Python API

This package exposes the Rust `ldsc` crate's LD-score, heritability, and
genetic correlation computations to Python, in two layers:

- **In-memory functions** (`fit_h2`, `fit_h2_partitioned`, `fit_rg`,
  `fit_rg_partitioned`, `compute_ld_scores_from_bytes`) take already-loaded,
  already-harmonised columns and do no file I/O.
- **File-oriented functions** (`estimate_h2`, `estimate_rg`,
  `estimate_ldscore`, `munge_sumstats`) take file paths, mirroring the
  `ldsc h2`/`rg`/`l2`/`munge-sumstats` CLI subcommands, and do the same
  loading, merging, and filtering — but return structured Python objects
  instead of printing. They cover the CLI's core options, including
  stratified LD Score regression (S-LDSC): partitioned LD score computation
  (`estimate_ldscore(..., annot=...)`) and overlap-corrected enrichment
  (`estimate_h2(..., overlap_annot=True)`). `--h2-cts` and jackknife-
  diagnostics printing remain out of scope for the Python API.

```python
from ldsc_rs import fit_h2, estimate_h2

# In-memory: caller already has aligned chi2/ref_ld/weight_ld/N columns.
result = fit_h2(chi2, ref_ld, weight_ld, sample_size, m_snps=1_190_321)
print(result.heritability.value, result.heritability.standard_error)

# File-oriented: caller has sumstats + LD score files on disk, like the CLI.
result = estimate_h2(
    "trait.sumstats.gz",
    ref_ld_chr="eur_w_ld_chr/",
    w_ld_chr="eur_w_ld_chr/",
)
print(result.heritability.value, result.n_snps)
```

## Scalar vs. partitioned (K=1 vs. K>1)

`fit_h2`/`fit_rg` are scalar (a single LD-score column). For partitioned or
stratified LD scores (multiple annotation columns), use `fit_h2_partitioned`/
`fit_rg_partitioned`, which take a 2-D `ref_ld` of shape `(n_snps, k)` and a
per-annotation `m_vec`. `fit_h2_partitioned` reports a per-annotation h2
breakdown (`per_annotation`) alongside the total; `fit_rg_partitioned`
reports totals only, matching `fit_rg` and the `rg` CLI's own multi-annotation
behavior. Partitioned h2 does not support the two-step estimator, matching
the CLI's own K>1 restriction.

`estimate_h2` dispatches between the scalar and partitioned path
automatically, based on how many annotation columns the loaded reference LD
scores carry — the returned `H2FileResult`'s `mean_chi2`/`lambda_gc`/`ratio`
fields are populated for K==1, and `l2_cols`/`per_annotation`/`m_values` for
K>1 (never both).

## Two-step and constrained intercepts

`fit_h2` follows the scalar CLI default and uses the two-step estimator at a
cutoff of 30. Its default `two_step="auto"` disables two-step estimation when
an intercept is constrained; pass a number or `None` to override that choice.
`fit_rg` follows the `rg` CLI default and does not enable two-step estimation
automatically. A constrained intercept is incompatible with two-step
estimation, and its standard error is represented as `None` in Python.

## Input validation

Inputs to the in-memory functions are equal-length finite numeric sequences;
reference and weight LD scores must be non-negative. For `fit_rg`/
`fit_rg_partitioned`, SNP joining and allele harmonisation must happen before
the call. The file-oriented functions do this joining/filtering themselves
(matching the CLI), and additionally clamp negative LD scores to zero rather
than rejecting them, since real per-annotation LD scores can dip slightly
below zero from estimation noise — the same defensive clamp the `rg` CLI
subcommand already applies.

## NumPy support

Every numeric sequence parameter (1-D or 2-D) accepts either a plain Python
sequence or a NumPy `float64` array. Passing a NumPy array avoids PyO3's
per-element Python-sequence extraction in favor of a direct buffer copy —
meaningfully faster for large arrays (whole-genome LD-score columns, up to
~1.6M SNPs). This does not make the eventual Rust-side allocation
zero-copy — the underlying linear-algebra type always owns its buffer — it
only skips the slow generic-sequence path on the way in.

## LD score computation from files

`estimate_ldscore(bfile, ...)` reads `{bfile}.bed`/`.bim`/`.fam` (a PLINK
`--bfile` prefix) from disk and is otherwise identical to
`compute_ld_scores_from_bytes` for the scalar (K==1) case — same parameters
(`window`, `chunk_size`, `dtype`, `sketch`, `sketch_maf_aware`,
`snp_level_masking`, `pq_exp`), no `--extract`/`--keep` filtering.

Pass `annot=<path or prefix>` for partitioned (S-LDSC) LD scores — one
column per annotation category, rows aligned 1:1 with the BIM (`thin_annot`
selects the thin `.annot` format with no CHR/SNP/BP/CM metadata columns).
The result's `l2_by_annot`/`annot_names`/`m_values`/`m_values_5_50` fields
are then populated; `ld_score` stays the first annotation's column for
backward compatibility with the scalar path. `annot` mirrors the CLI's
`--annot` for a single fileset (an explicit path, or a prefix auto-resolved
the same way): call `estimate_ldscore` once per chromosome for real S-LDSC
workflows, just like `ldsc l2 --annot`. Mutually exclusive with `pq_exp`.

Pass `out=<prefix>` to additionally write `{out}.l2.ldscore.gz`/`.l2.M`/
`.l2.M_5_50`, in exactly the format the `l2` CLI writes for a single
(non-chromosome-looped) `--bfile`/`--out` pair. This is the only way to
produce files `estimate_h2`/`estimate_rg` (or the CLI itself) can load
directly — omitting `out` returns the same `LdScoreResult` but does no
file I/O:

```python
from ldsc_rs import estimate_ldscore, estimate_h2

estimate_ldscore("1000G.chr1", out="1000G.chr1")
result = estimate_h2("trait.sumstats.gz", ref_ld="1000G.chr1", w_ld="1000G.chr1")
```

## Stratified LD Score regression (S-LDSC)

`estimate_h2(..., overlap_annot=True)` computes overlap-corrected
enrichment (Finucane et al. 2015), matching the `h2` CLI's
`--overlap-annot`. It requires a partitioned (K>1) fit and reads the
`.annot[.gz|.bz2]` files at the same location as `ref_ld`/`ref_ld_chr`
(one annotation file per chromosome for `ref_ld_chr`). Unless
`not_m_5_50=True`, it also requires `frqfile` (with `ref_ld`) or
`frqfile_chr` (with `ref_ld_chr`) to restrict M-counts to
0.05 < MAF < 0.95, exactly like the CLI. Results land in
`H2FileResult.overlap_enrichment`, a tuple of `OverlapEnrichmentCategory`
(one per category, in reference-LD-score-column order): `prop_snps`,
`prop_h2`, `enrichment` (each an `Estimate` where applicable),
`enrichment_diff_p`, and `coefficient`.

Pass `out=<prefix>` (only meaningful with `overlap_annot=True`) to
additionally write `{out}.results`, in exactly the format the CLI's
`h2 --overlap-annot --out {out}` writes (`print_coefficients=True` adds
coefficient/SE/z-score columns, matching `--print-coefficients`). Omitting
`out` does no file I/O; `overlap_enrichment` is populated either way:

```python
estimate_h2(
    "trait.sumstats.gz",
    ref_ld_chr="baselineLD.",
    w_ld="weights.",
    overlap_annot=True,
    frqfile_chr="1000G.frq.",
    out="trait",
)
```

## Scope and limitations

- `--h2-cts` and `--print-cov`/`--print-delete-vals` (jackknife
  diagnostics) are CLI-only; not exposed here.
- `estimate_ldscore`/`compute_ld_scores_from_bytes` don't support
  `--extract`/`--keep` filtering.
- `estimate_rg` computes `sumstats[0]` against every other trait in
  `sumstats`, like `ldsc rg --rg a,b,c`, returning one result per pair in
  the same order as `sumstats[1:]`.
- `munge_sumstats` writes `{out}.sumstats.gz`, exactly like the CLI, so its
  output is directly consumable by `estimate_h2`/`estimate_rg`.
- Liability-scale conversion is available on `estimate_h2` via
  `samp_prev`/`pop_prev` (populates `liability_heritability`); there is no
  rg equivalent, matching the CLI (which validates but never applies
  liability-scale conversion for `rg`).
