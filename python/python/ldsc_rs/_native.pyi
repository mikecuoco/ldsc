from typing import List, Optional, Sequence, Tuple

class NativeH2Result:
    h2: float
    h2_se: float
    intercept: float
    intercept_se: float
    mean_chi2: float
    lambda_gc: float
    ratio: Optional[Tuple[float, float]]

class NativeRgResult:
    h2_1: float
    h2_1_se: float
    h2_1_intercept: float
    h2_1_intercept_se: float
    h2_1_mean_chi2: float
    h2_1_lambda_gc: float
    h2_1_ratio: Optional[Tuple[float, float]]
    h2_2: float
    h2_2_se: float
    h2_2_intercept: float
    h2_2_intercept_se: float
    h2_2_mean_chi2: float
    h2_2_lambda_gc: float
    h2_2_ratio: Optional[Tuple[float, float]]
    gencov: float
    gencov_se: float
    gencov_intercept: float
    gencov_intercept_se: float
    mean_z1z2: float
    rg: float
    rg_se: float
    z: float
    p: float
    n_snps: int

class NativeLdScoreResult:
    snp: List[str]
    chromosome: List[int]
    base_pair: List[int]
    centimorgan: List[float]
    ld_score: List[float]
    maf: List[float]
    wall_seconds: float

def fit_h2(
    chi2: Sequence[float],
    ref_ld: Sequence[float],
    weight_ld: Sequence[float],
    sample_size: Sequence[float],
    m_snps: float,
    n_blocks: int = ...,
    two_step: Optional[float] = ...,
    intercept: Optional[float] = ...,
) -> NativeH2Result: ...

def fit_rg(
    z1: Sequence[float],
    z2: Sequence[float],
    ref_ld: Sequence[float],
    weight_ld: Sequence[float],
    sample_size1: Sequence[float],
    sample_size2: Sequence[float],
    m_snps: float,
    n_blocks: int = ...,
    two_step: Optional[float] = ...,
    intercept_h2_1: Optional[float] = ...,
    intercept_h2_2: Optional[float] = ...,
    intercept_gencov: Optional[float] = ...,
) -> NativeRgResult: ...

def compute_ld_scores_from_bytes(
    bed: bytes,
    bim: str,
    fam: str,
    window_unit: str = ...,
    window_value: float = ...,
    chunk_size: int = ...,
    dtype: str = ...,
    sketch: Optional[int] = ...,
    sketch_maf_aware: bool = ...,
    snp_level_masking: bool = ...,
    pq_exp: Optional[float] = ...,
) -> NativeLdScoreResult: ...
