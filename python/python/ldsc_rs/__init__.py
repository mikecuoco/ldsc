"""Typed Python interface to the ldsc-rs numerical core."""

from __future__ import annotations

from dataclasses import dataclass
from math import isnan
from typing import Literal, Optional, Sequence, Tuple, Union

from . import _native

__version__ = "0.5.0"
_AUTO_TWO_STEP = "auto"


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
    snp: tuple[str, ...]
    chromosome: tuple[int, ...]
    base_pair: tuple[int, ...]
    centimorgan: tuple[float, ...]
    ld_score: tuple[float, ...]
    maf: tuple[float, ...]
    wall_seconds: float


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
    return LdScoreResult(
        snp=tuple(result.snp),
        chromosome=tuple(result.chromosome),
        base_pair=tuple(result.base_pair),
        centimorgan=tuple(result.centimorgan),
        ld_score=tuple(result.ld_score),
        maf=tuple(result.maf),
        wall_seconds=result.wall_seconds,
    )


__all__ = [
    "Estimate",
    "GeneticCorrelationResult",
    "GeneticCovarianceResult",
    "HeritabilityResult",
    "LdScoreResult",
    "compute_ld_scores_from_bytes",
    "fit_h2",
    "fit_rg",
]
