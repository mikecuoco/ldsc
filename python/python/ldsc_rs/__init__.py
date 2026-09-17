"""Typed Python interface to the ldsc-rs numerical core."""

from __future__ import annotations

from dataclasses import dataclass
from math import isnan
from typing import Any, Literal, Optional, Sequence, Tuple, Union

from . import _native

__version__ = "0.5.0"
_AUTO_TWO_STEP = "auto"

#: A 2-D `(n_snps, k)` input: a NumPy array or a list of `k`-length rows.
FloatMatrix = Union[Sequence[Sequence[float]], Any]


@dataclass(frozen=True)
class Estimate:
    value: float
    standard_error: Optional[float]


@dataclass(frozen=True)
class HeritabilityResult:
    heritability: Estimate
    intercept: Estimate
    mean_chi2: float
    lambda_gc: float
    ratio: Optional[Estimate]


@dataclass(frozen=True)
class GeneticCovarianceResult:
    covariance: Estimate
    intercept: Estimate
    mean_z1z2: float


@dataclass(frozen=True)
class GeneticCorrelationResult:
    h2_1: HeritabilityResult
    h2_2: HeritabilityResult
    genetic_covariance: GeneticCovarianceResult
    correlation: Estimate
    z: float
    p: float
    n_snps: int


@dataclass(frozen=True)
class LdScoreResult:
    """Result of :func:`estimate_ldscore`/:func:`compute_ld_scores_from_bytes`.
    `ld_score` is always the first annotation column (identical to the
    scalar result when there was no `annot`). `l2_by_annot`/`annot_names`/
    `m_values`/`m_values_5_50` are populated only when `annot` was given."""

    snp: tuple[str, ...]
    chromosome: tuple[int, ...]
    base_pair: tuple[int, ...]
    centimorgan: tuple[float, ...]
    ld_score: tuple[float, ...]
    maf: tuple[float, ...]
    wall_seconds: float
    l2_by_annot: Optional[tuple[tuple[float, ...], ...]] = None
    annot_names: Optional[tuple[str, ...]] = None
    m_values: Optional[tuple[float, ...]] = None
    m_values_5_50: Optional[tuple[float, ...]] = None


@dataclass(frozen=True)
class PartitionedHeritabilityResult:
    """Partitioned (K>=1) heritability: a total plus a per-annotation breakdown."""

    heritability: Estimate
    intercept: Estimate
    per_annotation: tuple[Estimate, ...]
    m_values: tuple[float, ...]
    n_snps: int


@dataclass(frozen=True)
class OverlapEnrichmentCategory:
    """One category's overlap-corrected enrichment (Finucane et al. 2015
    `--overlap-annot`), as returned within :attr:`H2FileResult.overlap_enrichment`."""

    name: str
    prop_snps: float
    prop_h2: Estimate
    enrichment: Estimate
    enrichment_diff_p: Optional[float]
    coefficient: Estimate


@dataclass(frozen=True)
class H2FileResult:
    """Result of :func:`estimate_h2`. `per_annotation`/`l2_cols`/`m_values` are
    populated only when the loaded reference LD scores carried more than one
    annotation column; otherwise `mean_chi2`/`lambda_gc`/`ratio` are.
    `overlap_enrichment` is populated only when `overlap_annot=True` was
    passed and the estimate was partitioned (K>1)."""

    heritability: Estimate
    intercept: Estimate
    n_snps: int
    mean_chi2: Optional[float] = None
    lambda_gc: Optional[float] = None
    ratio: Optional[Estimate] = None
    l2_cols: Optional[tuple[str, ...]] = None
    per_annotation: Optional[tuple[Estimate, ...]] = None
    m_values: Optional[tuple[float, ...]] = None
    liability_heritability: Optional[float] = None
    overlap_enrichment: Optional[tuple[OverlapEnrichmentCategory, ...]] = None


@dataclass(frozen=True)
class MungeSummary:
    rows_in: int
    rows_out: int
    duplicates_removed: int
    out_path: str


def _h2_result(result: object) -> HeritabilityResult:
    ratio = None
    if result.ratio is not None:
        ratio = Estimate(result.ratio[0], result.ratio[1])
    return HeritabilityResult(
        heritability=Estimate(result.h2, result.h2_se),
        intercept=Estimate(
            result.intercept,
            None if isnan(result.intercept_se) else result.intercept_se,
        ),
        mean_chi2=result.mean_chi2,
        lambda_gc=result.lambda_gc,
        ratio=ratio,
    )


def _rg_h2_result(result: object, trait: int) -> HeritabilityResult:
    prefix = f"h2_{trait}"
    ratio_value = getattr(result, f"{prefix}_ratio")
    ratio = None if ratio_value is None else Estimate(ratio_value[0], ratio_value[1])
    return HeritabilityResult(
        heritability=Estimate(
            getattr(result, prefix), getattr(result, f"{prefix}_se")
        ),
        intercept=Estimate(
            getattr(result, f"{prefix}_intercept"),
            None
            if isnan(getattr(result, f"{prefix}_intercept_se"))
            else getattr(result, f"{prefix}_intercept_se"),
        ),
        mean_chi2=getattr(result, f"{prefix}_mean_chi2"),
        lambda_gc=getattr(result, f"{prefix}_lambda_gc"),
        ratio=ratio,
    )


def fit_h2(
    chi2: Sequence[float],
    ref_ld: Sequence[float],
    weight_ld: Sequence[float],
    sample_size: Sequence[float],
    *,
    m_snps: float,
    n_blocks: int = 200,
    two_step: Union[float, None, Literal["auto"]] = _AUTO_TWO_STEP,
    intercept: Optional[float] = None,
) -> HeritabilityResult:
    """Fit scalar LDSC h2 to aligned, in-memory columns."""
    if two_step == _AUTO_TWO_STEP:
        two_step = None if intercept is not None else 30.0
    return _h2_result(
        _native.fit_h2(
            chi2,
            ref_ld,
            weight_ld,
            sample_size,
            m_snps,
            n_blocks,
            two_step,
            intercept,
        )
    )


def _rg_result(result: object) -> GeneticCorrelationResult:
    return GeneticCorrelationResult(
        h2_1=_rg_h2_result(result, 1),
        h2_2=_rg_h2_result(result, 2),
        genetic_covariance=GeneticCovarianceResult(
            covariance=Estimate(result.gencov, result.gencov_se),
            intercept=Estimate(
                result.gencov_intercept,
                None
                if isnan(result.gencov_intercept_se)
                else result.gencov_intercept_se,
            ),
            mean_z1z2=result.mean_z1z2,
        ),
        correlation=Estimate(result.rg, result.rg_se),
        z=result.z,
        p=result.p,
        n_snps=result.n_snps,
    )


def fit_rg(
    z1: Sequence[float],
    z2: Sequence[float],
    ref_ld: Sequence[float],
    weight_ld: Sequence[float],
    sample_size1: Sequence[float],
    sample_size2: Sequence[float],
    *,
    m_snps: float,
    n_blocks: int = 200,
    two_step: Optional[float] = None,
    intercept_h2: Tuple[Optional[float], Optional[float]] = (None, None),
    intercept_gencov: Optional[float] = None,
) -> GeneticCorrelationResult:
    """Fit bivariate LDSC to SNP-aligned and allele-harmonised columns."""
    if len(intercept_h2) != 2:
        raise ValueError("intercept_h2 must contain exactly two values")
    result = _native.fit_rg(
        z1,
        z2,
        ref_ld,
        weight_ld,
        sample_size1,
        sample_size2,
        m_snps,
        n_blocks,
        two_step,
        intercept_h2[0],
        intercept_h2[1],
        intercept_gencov,
    )
    return _rg_result(result)


def fit_h2_partitioned(
    chi2: Sequence[float],
    ref_ld: FloatMatrix,
    weight_ld: Sequence[float],
    sample_size: Sequence[float],
    *,
    m_vec: Sequence[float],
    n_blocks: int = 200,
    intercept: Optional[float] = None,
) -> PartitionedHeritabilityResult:
    """Fit partitioned (K>=1) LDSC h2 to aligned, in-memory columns.

    `ref_ld` is 2-D: shape ``(n_snps, k)`` — a NumPy array or a list of
    `k`-length rows, one per SNP. `m_vec` carries one M value per
    annotation (length `k`). Does not support the two-step estimator
    (matching the `h2` CLI's own K>1 restriction).
    """
    result = _native.fit_h2_partitioned(
        chi2, ref_ld, weight_ld, sample_size, m_vec, n_blocks, intercept
    )
    return PartitionedHeritabilityResult(
        heritability=Estimate(result.h2_total, result.h2_total_se),
        intercept=Estimate(result.intercept, result.intercept_se),
        per_annotation=tuple(
            Estimate(h, se)
            for h, se in zip(result.h2_per_annot, result.h2_per_annot_se)
        ),
        m_values=tuple(result.m_vec),
        n_snps=result.n_snps,
    )


def fit_rg_partitioned(
    z1: Sequence[float],
    z2: Sequence[float],
    ref_ld: FloatMatrix,
    weight_ld: Sequence[float],
    sample_size1: Sequence[float],
    sample_size2: Sequence[float],
    *,
    m_vec: Sequence[float],
    n_blocks: int = 200,
    two_step: Optional[float] = None,
    intercept_h2: Tuple[Optional[float], Optional[float]] = (None, None),
    intercept_gencov: Optional[float] = None,
) -> GeneticCorrelationResult:
    """Fit bivariate LDSC with partitioned (K>=1) LD scores.

    `ref_ld` is 2-D like :func:`fit_h2_partitioned`'s. Reports totals only
    (matching :func:`fit_rg` and the `rg` CLI's own multi-annotation
    behavior) — no per-annotation breakdown.
    """
    if len(intercept_h2) != 2:
        raise ValueError("intercept_h2 must contain exactly two values")
    result = _native.fit_rg_partitioned(
        z1,
        z2,
        ref_ld,
        weight_ld,
        sample_size1,
        sample_size2,
        m_vec,
        n_blocks,
        two_step,
        intercept_h2[0],
        intercept_h2[1],
        intercept_gencov,
    )
    return _rg_result(result)


def _ldscore_result(result: object) -> LdScoreResult:
    l2_by_annot = None
    if result.l2_by_annot is not None:
        l2_by_annot = tuple(tuple(col) for col in result.l2_by_annot)
    return LdScoreResult(
        snp=tuple(result.snp),
        chromosome=tuple(result.chromosome),
        base_pair=tuple(result.base_pair),
        centimorgan=tuple(result.centimorgan),
        ld_score=tuple(result.ld_score),
        maf=tuple(result.maf),
        wall_seconds=result.wall_seconds,
        l2_by_annot=l2_by_annot,
        annot_names=tuple(result.annot_names) if result.annot_names is not None else None,
        m_values=tuple(result.m_vec) if result.m_vec is not None else None,
        m_values_5_50=tuple(result.m_vec_5_50) if result.m_vec_5_50 is not None else None,
    )


def compute_ld_scores_from_bytes(
    bed: bytes,
    bim: str,
    fam: str,
    *,
    window: Tuple[str, float] = ("kb", 1_000.0),
    chunk_size: int = 200,
    dtype: str = "float64",
    sketch: Optional[int] = None,
    sketch_maf_aware: bool = False,
    snp_level_masking: bool = False,
    pq_exp: Optional[float] = None,
) -> LdScoreResult:
    """Compute scalar LD scores from in-memory PLINK BED/BIM/FAM contents."""
    if len(window) != 2:
        raise ValueError("window must contain exactly (unit, value)")
    result = _native.compute_ld_scores_from_bytes(
        bed,
        bim,
        fam,
        window[0],
        window[1],
        chunk_size,
        dtype,
        sketch,
        sketch_maf_aware,
        snp_level_masking,
        pq_exp,
    )
    return _ldscore_result(result)


def estimate_ldscore(
    bfile: str,
    *,
    annot: Optional[str] = None,
    thin_annot: bool = False,
    out: Optional[str] = None,
    window: Tuple[str, float] = ("kb", 1_000.0),
    chunk_size: int = 200,
    dtype: str = "float64",
    sketch: Optional[int] = None,
    sketch_maf_aware: bool = False,
    snp_level_masking: bool = False,
    pq_exp: Optional[float] = None,
) -> LdScoreResult:
    """File-oriented, computation-only counterpart to the `l2` CLI subcommand.

    Reads `{bfile}.bed`/`.bim`/`.fam` (a PLINK `--bfile` prefix) from disk.
    No `--extract`/`--keep` filtering.

    Pass `annot` (a `.annot`/`.annot.gz` path, or a prefix auto-resolved the
    same way the CLI resolves a single-fileset `--annot`) for partitioned
    (S-LDSC) LD scores — one column per annotation category, rows aligned
    1:1 with the BIM. The result's `l2_by_annot`/`annot_names`/`m_values`/
    `m_values_5_50` are then populated; `ld_score` stays the first
    annotation's column for backward compatibility. Set `thin_annot=True`
    if the annotation file has no CHR/SNP/BP/CM metadata columns.

    Mutually exclusive with `pq_exp` (matching the CLI's per-chromosome
    `--annot` usage, call this once per chromosome bfile/annot pair).

    Pass `out` to additionally write `{out}.l2.ldscore.gz`/`.l2.M`/
    `.l2.M_5_50`, in exactly the format the `l2` CLI writes for a single
    (non-chromosome-looped) `--bfile`/`--out` pair — so the files are
    directly loadable by `estimate_h2`/`estimate_rg` (or the CLI itself)
    without writing them yourself. Omitting `out` does no file I/O; the
    returned `LdScoreResult` is unaffected either way.
    """
    if len(window) != 2:
        raise ValueError("window must contain exactly (unit, value)")
    result = _native.estimate_ldscore(
        bfile,
        annot=annot,
        thin_annot=thin_annot,
        out=out,
        window_unit=window[0],
        window_value=window[1],
        chunk_size=chunk_size,
        dtype=dtype,
        sketch=sketch,
        sketch_maf_aware=sketch_maf_aware,
        snp_level_masking=snp_level_masking,
        pq_exp=pq_exp,
    )
    return _ldscore_result(result)


def estimate_h2(
    sumstats: str,
    *,
    ref_ld: Optional[str] = None,
    ref_ld_chr: Optional[str] = None,
    w_ld: Optional[str] = None,
    w_ld_chr: Optional[str] = None,
    m_snps: Optional[float] = None,
    not_m_5_50: bool = False,
    n_blocks: int = 200,
    two_step: Optional[float] = None,
    intercept_h2: Optional[float] = None,
    no_intercept: bool = False,
    chisq_max: Optional[float] = None,
    samp_prev: Optional[float] = None,
    pop_prev: Optional[float] = None,
    overlap_annot: bool = False,
    frqfile: Optional[str] = None,
    frqfile_chr: Optional[str] = None,
    out: Optional[str] = None,
    print_coefficients: bool = False,
) -> H2FileResult:
    """File-oriented, computation-only counterpart to the `h2` CLI subcommand.

    Loads `sumstats` plus reference/weight LD scores from disk, then fits
    scalar (K==1) or partitioned (K>1) h2 depending on how many annotation
    columns the loaded reference LD scores carry. Exactly one of
    `ref_ld`/`ref_ld_chr` and one of `w_ld`/`w_ld_chr` must be given.

    Pass `overlap_annot=True` for overlap-corrected enrichment (Finucane et
    al. 2015 `--overlap-annot`), read from the `.annot[.gz|.bz2]` files at
    the same location as `ref_ld`/`ref_ld_chr`. K>1 only. Unless
    `not_m_5_50=True`, this also requires `frqfile` (with `ref_ld`) or
    `frqfile_chr` (with `ref_ld_chr`) to restrict M-counts to
    0.05 < MAF < 0.95, matching the CLI.

    `out`/`print_coefficients` are only meaningful together with
    `overlap_annot=True`: they additionally write `{out}.results`, in
    exactly the format the CLI's `h2 --overlap-annot --out {out}` writes
    (`--print-coefficients` adds coefficient/SE/z-score columns).
    Omitting `out` does no file I/O; `H2FileResult.overlap_enrichment` is
    populated either way — like `estimate_ldscore`, the returned result
    doesn't depend on whether you also asked for a file.
    """
    result = _native.estimate_h2(
        sumstats,
        ref_ld=ref_ld,
        ref_ld_chr=ref_ld_chr,
        w_ld=w_ld,
        w_ld_chr=w_ld_chr,
        m_snps=m_snps,
        not_m_5_50=not_m_5_50,
        n_blocks=n_blocks,
        two_step=two_step,
        intercept_h2=intercept_h2,
        no_intercept=no_intercept,
        chisq_max=chisq_max,
        samp_prev=samp_prev,
        pop_prev=pop_prev,
        overlap_annot=overlap_annot,
        frqfile=frqfile,
        frqfile_chr=frqfile_chr,
        out=out,
        print_coefficients=print_coefficients,
    )
    ratio = None if result.ratio is None else Estimate(result.ratio[0], result.ratio[1])
    per_annotation = None
    if result.h2_per_annot is not None:
        per_annotation = tuple(
            Estimate(h, se)
            for h, se in zip(result.h2_per_annot, result.h2_per_annot_se)
        )
    overlap_enrichment = None
    if result.overlap_enrichment is not None:
        oe = result.overlap_enrichment
        overlap_enrichment = tuple(
            OverlapEnrichmentCategory(
                name=name,
                prop_snps=prop_snps,
                prop_h2=Estimate(prop_h2, prop_h2_se),
                enrichment=Estimate(enrichment, enrichment_se),
                enrichment_diff_p=diff_p,
                coefficient=Estimate(coef, coef_se),
            )
            for name, prop_snps, prop_h2, prop_h2_se, enrichment, enrichment_se, diff_p, coef, coef_se in zip(
                oe.category_names,
                oe.prop_m_overlap,
                oe.prop_h2_overlap,
                oe.prop_h2_overlap_se,
                oe.enrichment,
                oe.enrichment_se,
                oe.enrichment_diff_p,
                oe.coefficient,
                oe.coefficient_se,
            )
        )
    return H2FileResult(
        heritability=Estimate(result.h2, result.h2_se),
        intercept=Estimate(result.intercept, result.intercept_se),
        n_snps=result.n_snps,
        mean_chi2=result.mean_chi2,
        lambda_gc=result.lambda_gc,
        ratio=ratio,
        l2_cols=None if result.l2_cols is None else tuple(result.l2_cols),
        per_annotation=per_annotation,
        m_values=None if result.m_vec is None else tuple(result.m_vec),
        liability_heritability=result.liability_h2,
        overlap_enrichment=overlap_enrichment,
    )


def estimate_rg(
    sumstats: Sequence[str],
    *,
    ref_ld: Optional[str] = None,
    ref_ld_chr: Optional[str] = None,
    w_ld: Optional[str] = None,
    w_ld_chr: Optional[str] = None,
    m_snps: Optional[float] = None,
    not_m_5_50: bool = False,
    n_blocks: int = 200,
    two_step: Optional[float] = None,
    chisq_max: Optional[float] = None,
    no_check_alleles: bool = False,
    no_intercept: bool = False,
    intercept_h2: Sequence[float] = (),
    intercept_gencov: Sequence[float] = (),
) -> Tuple[GeneticCorrelationResult, ...]:
    """File-oriented, computation-only counterpart to the `rg` CLI subcommand.

    Computes `sumstats[0]` vs every other trait in `sumstats`, in the same
    order as `sumstats[1:]`. Exactly one of `ref_ld`/`ref_ld_chr` and one of
    `w_ld`/`w_ld_chr` must be given.
    """
    results = _native.estimate_rg(
        list(sumstats),
        ref_ld=ref_ld,
        ref_ld_chr=ref_ld_chr,
        w_ld=w_ld,
        w_ld_chr=w_ld_chr,
        m_snps=m_snps,
        not_m_5_50=not_m_5_50,
        n_blocks=n_blocks,
        two_step=two_step,
        chisq_max=chisq_max,
        no_check_alleles=no_check_alleles,
        no_intercept=no_intercept,
        intercept_h2=list(intercept_h2),
        intercept_gencov=list(intercept_gencov),
    )
    return tuple(_rg_result(result) for result in results)


def munge_sumstats(
    sumstats: str,
    out: str,
    *,
    merge_alleles: Optional[str] = None,
    daner: bool = False,
    daner_n: bool = False,
    n_min: float = 0.0,
    maf: float = 0.01,
    info_min: float = 0.9,
    n: Optional[float] = None,
    n_cas: Optional[float] = None,
    n_con: Optional[float] = None,
    snp_col: Optional[str] = None,
    n_col: Optional[str] = None,
    n_cas_col: Optional[str] = None,
    n_con_col: Optional[str] = None,
    a1_col: Optional[str] = None,
    a2_col: Optional[str] = None,
    p_col: Optional[str] = None,
    frq_col: Optional[str] = None,
    info_col: Optional[str] = None,
    signed_sumstats: Optional[str] = None,
    ignore: Optional[str] = None,
    keep_maf: bool = False,
    a1_inc: bool = False,
    no_alleles: bool = False,
    info_list: Optional[str] = None,
    nstudy: Optional[str] = None,
    nstudy_min: Optional[int] = None,
) -> MungeSummary:
    """File-oriented, computation-only counterpart to `munge_sumstats.py`.

    Writes `{out}.sumstats.gz`.
    """
    result = _native.munge_sumstats(
        sumstats,
        out,
        merge_alleles=merge_alleles,
        daner=daner,
        daner_n=daner_n,
        n_min=n_min,
        maf=maf,
        info_min=info_min,
        n=n,
        n_cas=n_cas,
        n_con=n_con,
        snp_col=snp_col,
        n_col=n_col,
        n_cas_col=n_cas_col,
        n_con_col=n_con_col,
        a1_col=a1_col,
        a2_col=a2_col,
        p_col=p_col,
        frq_col=frq_col,
        info_col=info_col,
        signed_sumstats=signed_sumstats,
        ignore=ignore,
        keep_maf=keep_maf,
        a1_inc=a1_inc,
        no_alleles=no_alleles,
        info_list=info_list,
        nstudy=nstudy,
        nstudy_min=nstudy_min,
    )
    return MungeSummary(
        rows_in=result.rows_in,
        rows_out=result.rows_out,
        duplicates_removed=result.duplicates_removed,
        out_path=f"{out}.sumstats.gz",
    )


__all__ = [
    "Estimate",
    "GeneticCorrelationResult",
    "GeneticCovarianceResult",
    "H2FileResult",
    "HeritabilityResult",
    "LdScoreResult",
    "MungeSummary",
    "OverlapEnrichmentCategory",
    "PartitionedHeritabilityResult",
    "compute_ld_scores_from_bytes",
    "estimate_h2",
    "estimate_ldscore",
    "estimate_rg",
    "fit_h2",
    "fit_h2_partitioned",
    "fit_rg",
    "fit_rg_partitioned",
    "munge_sumstats",
]
