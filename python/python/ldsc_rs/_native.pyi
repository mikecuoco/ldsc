from typing import Any, List, Optional, Sequence, Tuple, Union

FloatArray = Union[Sequence[float], Any]  # Any covers numpy.ndarray
FloatMatrix = Union[Sequence[Sequence[float]], Any]  # Any covers numpy.ndarray (2-D)

class NativeH2Result:
    h2: float
    h2_se: float
    intercept: float
    intercept_se: float
    mean_chi2: float
    lambda_gc: float
    ratio: Optional[Tuple[float, float]]

class NativePartitionedH2Result:
    h2_per_annot: List[float]
    h2_per_annot_se: List[float]
    h2_total: float
    h2_total_se: float
    intercept: float
    intercept_se: Optional[float]
    m_vec: List[float]
    n_snps: int

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
    l2_by_annot: Optional[List[List[float]]]
    annot_names: Optional[List[str]]
    m_vec: Optional[List[float]]
    m_vec_5_50: Optional[List[float]]

class NativeOverlapEnrichmentResult:
    category_names: List[str]
    prop_m_overlap: List[float]
    prop_h2_overlap: List[float]
    prop_h2_overlap_se: List[float]
    enrichment: List[float]
    enrichment_se: List[float]
    enrichment_diff_p: List[Optional[float]]
    coefficient: List[float]
    coefficient_se: List[float]

class NativeH2FileResult:
    h2: float
    h2_se: float
    intercept: float
    intercept_se: Optional[float]
    mean_chi2: Optional[float]
    lambda_gc: Optional[float]
    ratio: Optional[Tuple[float, float]]
    l2_cols: Optional[List[str]]
    h2_per_annot: Optional[List[float]]
    h2_per_annot_se: Optional[List[float]]
    m_vec: Optional[List[float]]
    overlap_enrichment: Optional[NativeOverlapEnrichmentResult]
    n_snps: int
    liability_h2: Optional[float]

class NativeMungeSummary:
    rows_in: int
    rows_out: int
    duplicates_removed: int

def fit_h2(
    chi2: FloatArray,
    ref_ld: FloatArray,
    weight_ld: FloatArray,
    sample_size: FloatArray,
    m_snps: float,
    n_blocks: int = ...,
    two_step: Optional[float] = ...,
    intercept: Optional[float] = ...,
) -> NativeH2Result: ...

def fit_h2_partitioned(
    chi2: FloatArray,
    ref_ld: FloatMatrix,
    weight_ld: FloatArray,
    sample_size: FloatArray,
    m_vec: FloatArray,
    n_blocks: int = ...,
    intercept: Optional[float] = ...,
) -> NativePartitionedH2Result: ...

def fit_rg(
    z1: FloatArray,
    z2: FloatArray,
    ref_ld: FloatArray,
    weight_ld: FloatArray,
    sample_size1: FloatArray,
    sample_size2: FloatArray,
    m_snps: float,
    n_blocks: int = ...,
    two_step: Optional[float] = ...,
    intercept_h2_1: Optional[float] = ...,
    intercept_h2_2: Optional[float] = ...,
    intercept_gencov: Optional[float] = ...,
) -> NativeRgResult: ...

def fit_rg_partitioned(
    z1: FloatArray,
    z2: FloatArray,
    ref_ld: FloatMatrix,
    weight_ld: FloatArray,
    sample_size1: FloatArray,
    sample_size2: FloatArray,
    m_vec: FloatArray,
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

def estimate_ldscore(
    bfile: str,
    *,
    annot: Optional[str] = ...,
    thin_annot: bool = ...,
    window_unit: str = ...,
    window_value: float = ...,
    chunk_size: int = ...,
    dtype: str = ...,
    sketch: Optional[int] = ...,
    sketch_maf_aware: bool = ...,
    snp_level_masking: bool = ...,
    pq_exp: Optional[float] = ...,
) -> NativeLdScoreResult: ...

def estimate_h2(
    sumstats: str,
    *,
    ref_ld: Optional[str] = ...,
    ref_ld_chr: Optional[str] = ...,
    w_ld: Optional[str] = ...,
    w_ld_chr: Optional[str] = ...,
    m_snps: Optional[float] = ...,
    not_m_5_50: bool = ...,
    n_blocks: int = ...,
    two_step: Optional[float] = ...,
    intercept_h2: Optional[float] = ...,
    no_intercept: bool = ...,
    chisq_max: Optional[float] = ...,
    samp_prev: Optional[float] = ...,
    pop_prev: Optional[float] = ...,
    overlap_annot: bool = ...,
    frqfile: Optional[str] = ...,
    frqfile_chr: Optional[str] = ...,
) -> NativeH2FileResult: ...

def estimate_rg(
    sumstats: Sequence[str],
    *,
    ref_ld: Optional[str] = ...,
    ref_ld_chr: Optional[str] = ...,
    w_ld: Optional[str] = ...,
    w_ld_chr: Optional[str] = ...,
    m_snps: Optional[float] = ...,
    not_m_5_50: bool = ...,
    n_blocks: int = ...,
    two_step: Optional[float] = ...,
    chisq_max: Optional[float] = ...,
    no_check_alleles: bool = ...,
    no_intercept: bool = ...,
    intercept_h2: Sequence[float] = ...,
    intercept_gencov: Sequence[float] = ...,
) -> List[NativeRgResult]: ...

def munge_sumstats(
    sumstats: str,
    out: str,
    *,
    merge_alleles: Optional[str] = ...,
    daner: bool = ...,
    daner_n: bool = ...,
    n_min: float = ...,
    maf: float = ...,
    info_min: float = ...,
    n: Optional[float] = ...,
    n_cas: Optional[float] = ...,
    n_con: Optional[float] = ...,
    snp_col: Optional[str] = ...,
    n_col: Optional[str] = ...,
    n_cas_col: Optional[str] = ...,
    n_con_col: Optional[str] = ...,
    a1_col: Optional[str] = ...,
    a2_col: Optional[str] = ...,
    p_col: Optional[str] = ...,
    frq_col: Optional[str] = ...,
    info_col: Optional[str] = ...,
    signed_sumstats: Optional[str] = ...,
    ignore: Optional[str] = ...,
    keep_maf: bool = ...,
    a1_inc: bool = ...,
    no_alleles: bool = ...,
    info_list: Optional[str] = ...,
    nstudy: Optional[str] = ...,
    nstudy_min: Optional[int] = ...,
) -> NativeMungeSummary: ...
