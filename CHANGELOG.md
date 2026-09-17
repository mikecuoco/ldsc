# Changelog

All notable changes to this project will be documented here.
Format: [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).
Versioning: [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- **Python API (`ldsc_rs`): partitioned h2/rg, file-oriented workflows, NumPy
  input.** Extends the computation-only PyO3 bindings introduced in 0.5.0:
  - `fit_h2_partitioned`/`fit_rg_partitioned` — K>1 (partitioned/stratified
    LD score) heritability and genetic correlation on in-memory columns.
    `run_hsq_ldsc`/`run_rg_ldsc` (core crate) generalized to accept any K;
    `H2Result` gained optional per-annotation fields; a new
    `run_h2_ldsc_partitioned`/`PartitionedH2Result` reuses the `h2`
    subcommand's own regression exactly for CLI parity.
  - `estimate_h2`/`estimate_rg`/`estimate_ldscore`/`munge_sumstats` —
    file-oriented counterparts to the `h2`/`rg`/`l2`/`munge-sumstats` CLI
    subcommands: same file loading, merging, and filtering, structured
    Python results instead of stdout. Required splitting
    `run_h2`/`run_rg`/`munge::run`'s file I/O from their pure computation
    in the core crate (`estimate_h2_from_files`, `estimate_rg_from_files`,
    `munge_sumstats_df`/`munge_sumstats_from_files`,
    `compute_l2_from_bfile`). `estimate_h2` auto-dispatches between
    scalar (K==1) and partitioned (K>1) based on the loaded LD scores;
    `estimate_ldscore` is scalar-only, matching
    `compute_ld_scores_from_bytes`.
  - Every numeric input (1-D and 2-D) now accepts a NumPy `float64` array,
    via `python/src/pyarray.rs`'s `FloatColumn`/`FloatMatrix`, avoiding
    PyO3's per-element Python-sequence extraction for large arrays.
  - See `python/README.md` for the full API and its scope relative to the
    CLI (`--h2-cts`/jackknife-diagnostics printing are out of scope).
- **Python API: full S-LDSC support (`--annot` LD scores, `--overlap-annot`
  enrichment).** Closes the two remaining gaps for running Stratified LDSC
  from Python:
  - `estimate_ldscore` gained `annot`/`thin_annot` for partitioned LD score
    computation, matching `ldsc l2 --annot`. `L2Config`/`L2Output` (core
    crate) gained `annot`/`annot_names` and `l2_by_annot`/`annot_names`/
    `m_vec`/`m_vec_5_50` respectively; `compute_l2_from_bfile` gained
    `annot`/`thin_annot` parameters. `LdScoreResult` gained
    `l2_by_annot`/`annot_names`/`m_values`/`m_values_5_50`.
  - `estimate_h2` gained `overlap_annot`/`frqfile`/`frqfile_chr` for
    overlap-corrected enrichment (Finucane et al. 2015), matching
    `ldsc h2 --overlap-annot`. The CLI's `write_overlap_results` (core
    crate) was split into a pure `compute_overlap_enrichment` plus a thin
    file-writing wrapper, so `estimate_h2`/`H2FileOptions` reuse the exact
    same math via the new public `OverlapEnrichmentResult`.
    `H2FileResult.overlap_enrichment` exposes one `OverlapEnrichmentCategory`
    per annotation (`prop_snps`/`prop_h2`/`enrichment`/`enrichment_diff_p`/
    `coefficient`).

## [0.5.0] — 2026-05-12

### Added
- **`--sketch d --snp-level-masking` is now the recommended high-accuracy
  fast mode.** Combining sketch with per-SNP exact windows (the LDSC
  paper's mathematical definition of `ℓ_j = Σ_k r²_jk`) reaches the
  exact-mode h² at sketch speed:
  - Biobank N=50K: `--sketch 1000 --snp-level-masking` matches
    `--snp-level-masking --fast-f32` h² within 0.001 across BMI/Height
    GWAS, at ~17× the speed of exact (19s vs 361s).
  - 1000G EUR N=503: `--sketch 480 --snp-level-masking` matches GCTA-
    quality LD scores (Pearson r = 0.994 vs GCTA 0.993) with the lowest
    max-error of any method tested (39 vs 115-150 for chunked methods).
  - The combo was always supported but never tested at biobank scale —
    sweep results in `docs/perf-log.md` 2026-05-12 entry.
- **`--sketch-maf-aware`** — opt-in kurtosis-corrected quadratic bias
  formula. Mathematically tightens the asymptotic correction
  `(1−r²)(1−2r²)/d` with a per-SNP `−2r²(κ_j+κ_k)/d` term. Empirically
  zero impact on h² across all tested datasets/scales (the correction
  is O(1/N), below the regression measurement floor at any realistic
  N), but free at runtime. Kept for partitioned-LD use cases where
  rare-variant-stratified annotations might surface the term. See
  `docs/countsketch-math-analysis.md` §14.

### Changed
- **CLAUDE.md and `docs/countsketch-math-analysis.md` updated** with the
  sketch+masking guidance as the recommended high-accuracy path. Old
  recommendation was "use `--snp-level-masking --fast-f32` if you care
  about exactness" — at biobank scale that's now 17× slower than
  necessary.

### Removed
- **`--sketch-hybrid` and `--sketch-hybrid-fused`.** The hybrid variants
  ran the within-chunk B×B GEMM exact while keeping cross-window A×B
  sketched. Empirically gave a ~30% accuracy lift over sketch-only at
  ~1.7× the cost on 1000G — but `--sketch d --snp-level-masking` (which
  was always supported) gives a *larger* accuracy lift at *lower* cost.
  Hybrid was strictly dominated on both axes. Code complexity (separate
  B×B GEMM paths in F32+F64+GPU branches, side-buffer write in fused
  projector, mmap-vs-non-fused interaction handling) wasn't justified.
  Preserved on the `experiment/hybrid-and-maf-aware` branch for
  reference.

## [0.4.0] — 2026-03-14

### Added
- **Fused CountSketch** (`--sketch <d> --sketch-method countsketch`) — hash-based dimensionality
  reduction with fused BED-decode-normalize-scatter-add kernel. O(N) regardless of sketch
  dimension d. 20× speedup at biobank scale (N=50K), 3.6× faster than Rademacher sketch.
- **Rademacher sketch** (`--sketch <d>`) with ratio estimator for bias correction.
  Dense ±1/√d projection, 101× vs Python at 1000G scale (d=50).
- **Stochastic trace estimation** (`--stochastic <T>`) — Hutchinson's method with T random
  probes. ~43× vs Python at T=50. Best at small N; avoid T>50 (cache thrashing).
- **Random subsampling** (`--subsample <N>`) — subsample N' individuals before LD computation.
  Combinable with sketch modes for compound speedups.
- **mmap BED reader** (`--mmap`) — memory-mapped I/O with `MADV_SEQUENTIAL`/`MADV_WILLNEED`
  for GPFS and networked HPC filesystems. Zero-copy for fused CountSketch path.
- GPU native f64 matmul (`--gpu-f64`).
- AVX2+FMA SIMD intrinsics for normalization sum/sum-of-squares reduction.
- Branchless byte-level LUT for fused CountSketch scatter-add (eliminates NaN branching).
- PGO build script (`scripts/build_pgo.sh`).
- Comprehensive AWS Batch benchmark infrastructure (`scripts/aws-bench.sh`, `scripts/biobank-bench.sh`).
- Parity sweep script (`scripts/verify_parity_all.sh`) — 8-mode exact/sketch × f64/f32 vs Python.
- LDSC paper and supplementary note in `docs/`.

### Changed
- l2 module refactored from single `l2.rs` (1765 lines) into `l2/` directory with 6 submodules
  (compute, window, normalize, snp_stats, io, mod).
- `--sketch` now errors on out-of-bounds dimension (was a silent warning that fell back to exact).
- `--sketch-method` validates input; rejects unknown methods (was silently defaulting to Rademacher).
- Hot-path performance: bounds-check elimination, pre-computed r²_unbiased constants,
  column-major accumulation, iterator-based loops.
- h2/rg internals: `Par::Seq` for small matmuls, O(n) median selection.
- Removed unused `tokio` dependency.

### Performance (AWS EPYC 7R13, 1.66M SNPs, `--ld-wind-kb 1000`)
- **1000G (N=2,490)**: exact 37.7×, stochastic-50 43×, sketch-50 101× vs Python (1548s)
- **Biobank (N=50,000)**: countsketch-50 20.1×, countsketch-200 19.7× vs exact-f64 (666s)

## [0.3.1] — 2026-03-08

### Added
- GPU flex32 precision (`--gpu-flex32`) — half-precision compute with f32 accumulation.

### Removed
- Dropped `half` crate dependency (flex32 uses CubeCL native support).

## [0.3.0] — 2026-03-07

### Added
- Runtime f32 dispatch (`--fast-f32`) — runtime flag instead of compile-time feature.
- Optional GPU acceleration (`--gpu`) via CubeCL CUDA backend.

## [0.2.1] — 2026-03-04

### Added
- GitHub release builds no longer require BLAS feature flags; faer-only build across platforms.

### Fixed
- Crates.io packaging excludes local repo copies (`ldsc.wiki`, `ldsc_py`, `bed-reader`) and special-character filenames.
- Minor L2 trace timing/format and clippy cleanup.

## [0.2.0] — 2026-03-04

### Added
- Internal SNP-major BED reader to remove the external bed-reader dependency.
- `faer`-based linear algebra helpers (`src/la.rs`) to standardize matrix operations.
- L2 parity/perf tooling scripts for quick validation and benchmarking.
- Experimental `fast-f32` feature for L2 matmuls (f32 compute, f64 accumulation; not parity-safe).

### Changed
- Replaced the ndarray/OpenBLAS backend with `faer`; removed BLAS build dependencies and `--blas-threads`.
- L2 defaults and hot-path optimizations (ring alignment, faster MAF prefilter, tuned chunk size).

### Removed
- OpenBLAS/vcpkg build plumbing (`build.rs`, vcpkg manifests).

## [0.1.3] — 2026-02-23

### Added
- `cts-annot` continuous-annotation binning (Python `--cts-bin`) and `.annot` output.
- Cell-type-specific `h2-cts` support with `.ldcts` inputs and `--print-all-cts`.
- Overlap-aware partitioned h2 via `--overlap-annot` with frequency inputs.
- Windows (MSVC) system-BLAS support via vcpkg (OpenBLAS + Clapack).
- GitHub Release artifacts for Linux/macOS/Windows with checksums.

### Changed
- Release packaging uses system OpenBLAS on Linux/macOS.
- README reorganized for a faster install/usage path.

## [0.1.0] — 2026-02-19

### Added
- `munge-sumstats` — Polars LazyFrame streaming pipeline; full Python API parity
  (column-name overrides, `--signed-sumstats`, `--info-list`, `--nstudy`, `--a1-inc`, etc.)
- `l2` — ring-buffer DGEMM loop; global sequential pass matching Python's
  cross-chromosome LD window behaviour; 7× speedup over Python on 1000G data
  - `--annot` (partitioned LD scores, per-chromosome annot files)
  - `--extract`, `--print-snps`, `--maf`, `--keep`, `--per-allele`
- `h2` — scalar and partitioned heritability; `--two-step`, `--intercept-h2`,
  `--no-intercept`, `--print-coefficients`, liability-scale output
- `rg` — genetic correlation for all trait pairs; `--two-step`, `--intercept-gencov`,
  `--no-intercept`, `--samp-prev`/`--pop-prev`
- `make-annot` — BED-file and gene-set annotation generators
- Statically linked OpenBLAS (no runtime library needed)
- 40 unit tests; integration smoke test for l2

[Unreleased]: https://github.com/sharifhsn/ldsc/compare/v0.5.0...HEAD
[0.5.0]: https://github.com/sharifhsn/ldsc/compare/v0.4.0...v0.5.0
[0.4.0]: https://github.com/sharifhsn/ldsc/compare/v0.3.1...v0.4.0
[0.3.1]: https://github.com/sharifhsn/ldsc/compare/v0.3.0...v0.3.1
[0.3.0]: https://github.com/sharifhsn/ldsc/compare/v0.2.1...v0.3.0
[0.2.1]: https://github.com/sharifhsn/ldsc/compare/v0.2.0...v0.2.1
[0.2.0]: https://github.com/sharifhsn/ldsc/compare/v0.1.3...v0.2.0
[0.1.3]: https://github.com/sharifhsn/ldsc/releases/tag/v0.1.3
[0.1.0]: https://github.com/sharifhsn/ldsc/releases/tag/v0.1.0
