mod pyarray;

use ldsc::h2::{H2Result, run_h2_ldsc};
use ldsc::l2::{L2Config, WindowMode, compute_l2_from_bytes};
use ldsc::la::col_from_vec;
use ldsc::munge::{MungeOptions, MungeSummary, munge_sumstats_from_files};
use ldsc::regressions::{
    GeneticCorrelationResult, H2Estimate, H2FileOptions, H2FileResult, PartitionedH2Result,
    RgFileOptions, estimate_h2_from_files, estimate_rg_from_files, run_h2_ldsc_partitioned,
    run_rg_ldsc,
};
use pyarray::{FloatColumn, FloatMatrix};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyBytes;

fn value_error(error: impl std::fmt::Display) -> PyErr {
    PyValueError::new_err(format!("{error:#}"))
}

#[pyclass(frozen, get_all, module = "ldsc_rs._native", skip_from_py_object)]
#[derive(Clone)]
struct NativeH2Result {
    h2: f64,
    h2_se: f64,
    intercept: f64,
    intercept_se: f64,
    mean_chi2: f64,
    lambda_gc: f64,
    ratio: Option<(f64, f64)>,
}

impl From<H2Result> for NativeH2Result {
    fn from(result: H2Result) -> Self {
        Self {
            h2: result.h2,
            h2_se: result.h2_se,
            intercept: result.intercept,
            intercept_se: result.intercept_se,
            mean_chi2: result.mean_chi2,
            lambda_gc: result.lambda_gc,
            ratio: result.ratio,
        }
    }
}

#[pyclass(frozen, get_all, module = "ldsc_rs._native")]
struct NativePartitionedH2Result {
    h2_per_annot: Vec<f64>,
    h2_per_annot_se: Vec<f64>,
    h2_total: f64,
    h2_total_se: f64,
    intercept: f64,
    intercept_se: Option<f64>,
    m_vec: Vec<f64>,
    n_snps: usize,
}

impl From<PartitionedH2Result> for NativePartitionedH2Result {
    fn from(result: PartitionedH2Result) -> Self {
        Self {
            h2_per_annot: result.h2_per_annot,
            h2_per_annot_se: result.h2_per_annot_se,
            h2_total: result.h2_total,
            h2_total_se: result.h2_total_se,
            intercept: result.intercept,
            intercept_se: result.intercept_se,
            m_vec: result.m_vec,
            n_snps: result.n_snps,
        }
    }
}

#[pyclass(frozen, get_all, module = "ldsc_rs._native")]
struct NativeRgResult {
    h2_1: f64,
    h2_1_se: f64,
    h2_1_intercept: f64,
    h2_1_intercept_se: f64,
    h2_1_mean_chi2: f64,
    h2_1_lambda_gc: f64,
    h2_1_ratio: Option<(f64, f64)>,
    h2_2: f64,
    h2_2_se: f64,
    h2_2_intercept: f64,
    h2_2_intercept_se: f64,
    h2_2_mean_chi2: f64,
    h2_2_lambda_gc: f64,
    h2_2_ratio: Option<(f64, f64)>,
    gencov: f64,
    gencov_se: f64,
    gencov_intercept: f64,
    gencov_intercept_se: f64,
    mean_z1z2: f64,
    rg: f64,
    rg_se: f64,
    z: f64,
    p: f64,
    n_snps: usize,
}

impl From<GeneticCorrelationResult> for NativeRgResult {
    fn from(result: GeneticCorrelationResult) -> Self {
        Self {
            h2_1: result.h2_1.h2,
            h2_1_se: result.h2_1.h2_se,
            h2_1_intercept: result.h2_1.intercept,
            h2_1_intercept_se: result.h2_1.intercept_se,
            h2_1_mean_chi2: result.h2_1.mean_chi2,
            h2_1_lambda_gc: result.h2_1.lambda_gc,
            h2_1_ratio: result.h2_1.ratio,
            h2_2: result.h2_2.h2,
            h2_2_se: result.h2_2.h2_se,
            h2_2_intercept: result.h2_2.intercept,
            h2_2_intercept_se: result.h2_2.intercept_se,
            h2_2_mean_chi2: result.h2_2.mean_chi2,
            h2_2_lambda_gc: result.h2_2.lambda_gc,
            h2_2_ratio: result.h2_2.ratio,
            gencov: result.genetic_covariance.covariance,
            gencov_se: result.genetic_covariance.covariance_se,
            gencov_intercept: result.genetic_covariance.intercept,
            gencov_intercept_se: result.genetic_covariance.intercept_se,
            mean_z1z2: result.genetic_covariance.mean_z1z2,
            rg: result.rg,
            rg_se: result.rg_se,
            z: result.z,
            p: result.p,
            n_snps: result.n_snps,
        }
    }
}

#[pyclass(frozen, get_all, module = "ldsc_rs._native")]
struct NativeLdScoreResult {
    snp: Vec<String>,
    chromosome: Vec<u8>,
    base_pair: Vec<u32>,
    centimorgan: Vec<f64>,
    ld_score: Vec<f64>,
    maf: Vec<f64>,
    wall_seconds: f64,
}

#[pyfunction]
#[pyo3(signature = (
    chi2,
    ref_ld,
    weight_ld,
    sample_size,
    m_snps,
    n_blocks=200,
    two_step=Some(30.0),
    intercept=None
))]
#[allow(clippy::too_many_arguments)]
fn fit_h2<'py>(
    py: Python<'py>,
    chi2: FloatColumn<'py>,
    ref_ld: FloatColumn<'py>,
    weight_ld: FloatColumn<'py>,
    sample_size: FloatColumn<'py>,
    m_snps: f64,
    n_blocks: usize,
    two_step: Option<f64>,
    intercept: Option<f64>,
) -> PyResult<NativeH2Result> {
    let chi2 = chi2.to_vec()?;
    let ref_ld = ref_ld.to_vec()?;
    let weight_ld = weight_ld.to_vec()?;
    let sample_size = sample_size.to_vec()?;
    let n = chi2.len();
    if ref_ld.len() != n || weight_ld.len() != n || sample_size.len() != n {
        return Err(value_error("run_h2_ldsc: input length mismatch"));
    }
    if !m_snps.is_finite() || m_snps <= 0.0 {
        return Err(value_error("m_snps must be finite and > 0"));
    }
    if n_blocks <= 1 {
        return Err(value_error("n_blocks must be > 1"));
    }
    let n_blocks = n_blocks.min(n);
    if n_blocks <= 1 {
        return Err(value_error("at least two observations are required"));
    }
    for i in 0..n {
        if !chi2[i].is_finite()
            || !ref_ld[i].is_finite()
            || !weight_ld[i].is_finite()
            || !sample_size[i].is_finite()
        {
            return Err(value_error("all inputs must be finite"));
        }
        if ref_ld[i] < 0.0 || weight_ld[i] < 0.0 {
            return Err(value_error("LD scores must be non-negative"));
        }
        if sample_size[i] <= 0.0 {
            return Err(value_error("sample sizes must be > 0"));
        }
    }
    if two_step.is_some_and(|value| !value.is_finite()) {
        return Err(value_error("two_step must be finite"));
    }
    if intercept.is_some_and(|value| !value.is_finite()) {
        return Err(value_error("intercept must be finite"));
    }
    let result = py
        .detach(move || {
            run_h2_ldsc(
                &col_from_vec(chi2),
                &col_from_vec(ref_ld),
                &col_from_vec(weight_ld),
                &col_from_vec(sample_size),
                m_snps,
                n_blocks,
                two_step,
                intercept,
            )
        })
        .map_err(value_error)?;
    Ok(result.into())
}

#[pyfunction]
#[pyo3(signature = (
    z1,
    z2,
    ref_ld,
    weight_ld,
    sample_size1,
    sample_size2,
    m_snps,
    n_blocks=200,
    two_step=None,
    intercept_h2_1=None,
    intercept_h2_2=None,
    intercept_gencov=None
))]
#[allow(clippy::too_many_arguments)]
fn fit_rg<'py>(
    py: Python<'py>,
    z1: FloatColumn<'py>,
    z2: FloatColumn<'py>,
    ref_ld: FloatColumn<'py>,
    weight_ld: FloatColumn<'py>,
    sample_size1: FloatColumn<'py>,
    sample_size2: FloatColumn<'py>,
    m_snps: f64,
    n_blocks: usize,
    two_step: Option<f64>,
    intercept_h2_1: Option<f64>,
    intercept_h2_2: Option<f64>,
    intercept_gencov: Option<f64>,
) -> PyResult<NativeRgResult> {
    let z1 = z1.to_vec()?;
    let z2 = z2.to_vec()?;
    let ref_ld = ref_ld.to_vec()?;
    let weight_ld = weight_ld.to_vec()?;
    let sample_size1 = sample_size1.to_vec()?;
    let sample_size2 = sample_size2.to_vec()?;
    let result = py
        .detach(move || {
            run_rg_ldsc(
                &col_from_vec(z1),
                &col_from_vec(z2),
                &[col_from_vec(ref_ld)],
                &col_from_vec(weight_ld),
                &col_from_vec(sample_size1),
                &col_from_vec(sample_size2),
                &[m_snps],
                n_blocks,
                two_step,
                intercept_h2_1,
                intercept_h2_2,
                intercept_gencov,
            )
        })
        .map_err(value_error)?;
    Ok(result.into())
}

/// Partitioned (K>=1) LD Score regression h2 on aligned in-memory columns.
/// `ref_ld` is 2-D: shape `(n_snps, k)` — a NumPy array or a list of `k`-length
/// rows. `m_vec` carries one M value per annotation (length `k`).
#[pyfunction]
#[pyo3(signature = (
    chi2,
    ref_ld,
    weight_ld,
    sample_size,
    m_vec,
    n_blocks=200,
    intercept=None
))]
#[allow(clippy::too_many_arguments)]
fn fit_h2_partitioned<'py>(
    py: Python<'py>,
    chi2: FloatColumn<'py>,
    ref_ld: FloatMatrix<'py>,
    weight_ld: FloatColumn<'py>,
    sample_size: FloatColumn<'py>,
    m_vec: FloatColumn<'py>,
    n_blocks: usize,
    intercept: Option<f64>,
) -> PyResult<NativePartitionedH2Result> {
    let chi2 = chi2.to_vec()?;
    let ref_l2_k = ref_ld.to_columns()?;
    let weight_ld = weight_ld.to_vec()?;
    let sample_size = sample_size.to_vec()?;
    let m_vec = m_vec.to_vec()?;
    let result = py
        .detach(move || {
            let ref_l2_k: Vec<_> = ref_l2_k.into_iter().map(col_from_vec).collect();
            run_h2_ldsc_partitioned(
                &col_from_vec(chi2),
                &ref_l2_k,
                &col_from_vec(weight_ld),
                &col_from_vec(sample_size),
                &m_vec,
                n_blocks,
                intercept,
            )
        })
        .map_err(value_error)?;
    Ok(result.into())
}

/// Partitioned (K>=1) bivariate LD Score regression on aligned in-memory
/// columns. `ref_ld` is 2-D like [`fit_h2_partitioned`]'s. Reports totals
/// only (matching [`fit_rg`] and the `rg` CLI's own multi-annotation
/// behavior) — no per-annotation breakdown.
#[pyfunction]
#[pyo3(signature = (
    z1,
    z2,
    ref_ld,
    weight_ld,
    sample_size1,
    sample_size2,
    m_vec,
    n_blocks=200,
    two_step=None,
    intercept_h2_1=None,
    intercept_h2_2=None,
    intercept_gencov=None
))]
#[allow(clippy::too_many_arguments)]
fn fit_rg_partitioned<'py>(
    py: Python<'py>,
    z1: FloatColumn<'py>,
    z2: FloatColumn<'py>,
    ref_ld: FloatMatrix<'py>,
    weight_ld: FloatColumn<'py>,
    sample_size1: FloatColumn<'py>,
    sample_size2: FloatColumn<'py>,
    m_vec: FloatColumn<'py>,
    n_blocks: usize,
    two_step: Option<f64>,
    intercept_h2_1: Option<f64>,
    intercept_h2_2: Option<f64>,
    intercept_gencov: Option<f64>,
) -> PyResult<NativeRgResult> {
    let z1 = z1.to_vec()?;
    let z2 = z2.to_vec()?;
    let ref_l2_k = ref_ld.to_columns()?;
    let weight_ld = weight_ld.to_vec()?;
    let sample_size1 = sample_size1.to_vec()?;
    let sample_size2 = sample_size2.to_vec()?;
    let m_vec = m_vec.to_vec()?;
    let result = py
        .detach(move || {
            let ref_l2_k: Vec<_> = ref_l2_k.into_iter().map(col_from_vec).collect();
            run_rg_ldsc(
                &col_from_vec(z1),
                &col_from_vec(z2),
                &ref_l2_k,
                &col_from_vec(weight_ld),
                &col_from_vec(sample_size1),
                &col_from_vec(sample_size2),
                &m_vec,
                n_blocks,
                two_step,
                intercept_h2_1,
                intercept_h2_2,
                intercept_gencov,
            )
        })
        .map_err(value_error)?;
    Ok(result.into())
}

#[pyfunction]
#[pyo3(signature = (
    bed,
    bim,
    fam,
    window_unit="kb",
    window_value=1000.0,
    chunk_size=200,
    dtype="float64",
    sketch=None,
    sketch_maf_aware=false,
    snp_level_masking=false,
    pq_exp=None
))]
#[allow(clippy::too_many_arguments)]
fn compute_ld_scores_from_bytes(
    py: Python<'_>,
    bed: &Bound<'_, PyBytes>,
    bim: String,
    fam: String,
    window_unit: &str,
    window_value: f64,
    chunk_size: usize,
    dtype: &str,
    sketch: Option<usize>,
    sketch_maf_aware: bool,
    snp_level_masking: bool,
    pq_exp: Option<f64>,
) -> PyResult<NativeLdScoreResult> {
    if !window_value.is_finite() || window_value <= 0.0 {
        return Err(value_error("window value must be finite and > 0"));
    }
    if chunk_size == 0 {
        return Err(value_error("chunk_size must be > 0"));
    }
    let mode = match window_unit {
        "cm" => WindowMode::Cm(window_value),
        "kb" => WindowMode::Kb(window_value),
        "snps" if window_value.fract() == 0.0 => WindowMode::Snp(window_value as usize),
        "snps" => return Err(value_error("SNP window value must be an integer")),
        _ => return Err(value_error("window unit must be 'cm', 'kb', or 'snps'")),
    };
    let use_f32 = match dtype {
        "float64" => false,
        "float32" => true,
        _ => return Err(value_error("dtype must be 'float64' or 'float32'")),
    };
    if sketch_maf_aware && sketch.is_none() {
        return Err(value_error("sketch_maf_aware requires sketch"));
    }
    if sketch == Some(0) {
        return Err(value_error("sketch dimension must be > 0"));
    }
    let bed = bed.as_bytes().to_vec();
    let result = py
        .detach(move || {
            compute_l2_from_bytes(
                bed,
                &bim,
                &fam,
                L2Config {
                    mode,
                    chunk_size,
                    use_f32,
                    sketch,
                    sketch_maf_aware,
                    snp_level_masking,
                    yes_really: true,
                    pq_exp,
                    verbose_timing: false,
                },
            )
        })
        .map_err(value_error)?;

    Ok(NativeLdScoreResult {
        snp: result.snps.iter().map(|snp| snp.snp.clone()).collect(),
        chromosome: result.snps.iter().map(|snp| snp.chr).collect(),
        base_pair: result.snps.iter().map(|snp| snp.bp).collect(),
        centimorgan: result.snps.iter().map(|snp| snp.cm).collect(),
        ld_score: result.l2,
        maf: result.maf,
        wall_seconds: result.wall_seconds,
    })
}

/// File-oriented h2 result: fields are populated for the scalar (K==1)
/// case, or the partitioned (K>1) case, depending on which annotation
/// count the loaded reference LD scores carried — never both.
#[pyclass(frozen, get_all, module = "ldsc_rs._native")]
struct NativeH2FileResult {
    h2: f64,
    h2_se: f64,
    intercept: f64,
    intercept_se: Option<f64>,
    // Scalar (K==1) only.
    mean_chi2: Option<f64>,
    lambda_gc: Option<f64>,
    ratio: Option<(f64, f64)>,
    // Partitioned (K>1) only.
    l2_cols: Option<Vec<String>>,
    h2_per_annot: Option<Vec<f64>>,
    h2_per_annot_se: Option<Vec<f64>>,
    m_vec: Option<Vec<f64>>,
    // Always present.
    n_snps: usize,
    liability_h2: Option<f64>,
}

impl From<H2FileResult> for NativeH2FileResult {
    fn from(result: H2FileResult) -> Self {
        let liability_h2 = result.liability_h2;
        let n_snps = result.n_snps;
        match result.estimate {
            H2Estimate::Scalar(r) => Self {
                h2: r.h2,
                h2_se: r.h2_se,
                intercept: r.intercept,
                intercept_se: if r.intercept_se.is_nan() {
                    None
                } else {
                    Some(r.intercept_se)
                },
                mean_chi2: Some(r.mean_chi2),
                lambda_gc: Some(r.lambda_gc),
                ratio: r.ratio,
                l2_cols: None,
                h2_per_annot: None,
                h2_per_annot_se: None,
                m_vec: None,
                n_snps,
                liability_h2,
            },
            H2Estimate::Partitioned(p, l2_cols) => Self {
                h2: p.h2_total,
                h2_se: p.h2_total_se,
                intercept: p.intercept,
                intercept_se: p.intercept_se,
                mean_chi2: None,
                lambda_gc: None,
                ratio: None,
                l2_cols: Some(l2_cols),
                h2_per_annot: Some(p.h2_per_annot),
                h2_per_annot_se: Some(p.h2_per_annot_se),
                m_vec: Some(p.m_vec),
                n_snps,
                liability_h2,
            },
        }
    }
}

#[pyclass(frozen, get_all, module = "ldsc_rs._native")]
struct NativeMungeSummary {
    rows_in: usize,
    rows_out: usize,
    duplicates_removed: usize,
}

impl From<MungeSummary> for NativeMungeSummary {
    fn from(s: MungeSummary) -> Self {
        Self {
            rows_in: s.rows_in,
            rows_out: s.rows_out,
            duplicates_removed: s.duplicates_removed,
        }
    }
}

/// File-oriented, computation-only counterpart to the `h2` CLI subcommand.
/// See [`ldsc::regressions::estimate_h2_from_files`].
#[pyfunction]
#[pyo3(signature = (
    sumstats,
    *,
    ref_ld=None,
    ref_ld_chr=None,
    w_ld=None,
    w_ld_chr=None,
    m_snps=None,
    not_m_5_50=false,
    n_blocks=200,
    two_step=None,
    intercept_h2=None,
    no_intercept=false,
    chisq_max=None,
    samp_prev=None,
    pop_prev=None
))]
#[allow(clippy::too_many_arguments)]
fn estimate_h2(
    py: Python<'_>,
    sumstats: String,
    ref_ld: Option<String>,
    ref_ld_chr: Option<String>,
    w_ld: Option<String>,
    w_ld_chr: Option<String>,
    m_snps: Option<f64>,
    not_m_5_50: bool,
    n_blocks: usize,
    two_step: Option<f64>,
    intercept_h2: Option<f64>,
    no_intercept: bool,
    chisq_max: Option<f64>,
    samp_prev: Option<f64>,
    pop_prev: Option<f64>,
) -> PyResult<NativeH2FileResult> {
    let opts = H2FileOptions {
        m_snps,
        not_m_5_50,
        n_blocks,
        two_step,
        intercept_h2,
        no_intercept,
        chisq_max,
        samp_prev,
        pop_prev,
    };
    let result = py
        .detach(move || {
            estimate_h2_from_files(
                &sumstats,
                ref_ld.as_deref(),
                ref_ld_chr.as_deref(),
                w_ld.as_deref(),
                w_ld_chr.as_deref(),
                &opts,
            )
        })
        .map_err(value_error)?;
    Ok(result.into())
}

/// File-oriented, computation-only counterpart to the `rg` CLI subcommand.
/// One [`NativeRgResult`] per non-reference trait, in the same order as
/// `sumstats[1..]`. See [`ldsc::regressions::estimate_rg_from_files`].
#[pyfunction]
#[pyo3(signature = (
    sumstats,
    *,
    ref_ld=None,
    ref_ld_chr=None,
    w_ld=None,
    w_ld_chr=None,
    m_snps=None,
    not_m_5_50=false,
    n_blocks=200,
    two_step=None,
    chisq_max=None,
    no_check_alleles=false,
    no_intercept=false,
    intercept_h2=Vec::new(),
    intercept_gencov=Vec::new()
))]
#[allow(clippy::too_many_arguments)]
fn estimate_rg(
    py: Python<'_>,
    sumstats: Vec<String>,
    ref_ld: Option<String>,
    ref_ld_chr: Option<String>,
    w_ld: Option<String>,
    w_ld_chr: Option<String>,
    m_snps: Option<f64>,
    not_m_5_50: bool,
    n_blocks: usize,
    two_step: Option<f64>,
    chisq_max: Option<f64>,
    no_check_alleles: bool,
    no_intercept: bool,
    intercept_h2: Vec<f64>,
    intercept_gencov: Vec<f64>,
) -> PyResult<Vec<NativeRgResult>> {
    let opts = RgFileOptions {
        m_snps,
        not_m_5_50,
        n_blocks,
        two_step,
        chisq_max,
        no_check_alleles,
        no_intercept,
        intercept_h2,
        intercept_gencov,
    };
    let results = py
        .detach(move || {
            estimate_rg_from_files(
                &sumstats,
                ref_ld.as_deref(),
                ref_ld_chr.as_deref(),
                w_ld.as_deref(),
                w_ld_chr.as_deref(),
                &opts,
            )
        })
        .map_err(value_error)?;
    Ok(results.into_iter().map(NativeRgResult::from).collect())
}

/// File-oriented, computation-only counterpart to `munge_sumstats.py`.
/// Writes `{out}.sumstats.gz`. See [`ldsc::munge::munge_sumstats_from_files`].
#[pyfunction]
#[pyo3(signature = (
    sumstats,
    out,
    *,
    merge_alleles=None,
    daner=false,
    daner_n=false,
    n_min=0.0,
    maf=0.01,
    info_min=0.9,
    n=None,
    n_cas=None,
    n_con=None,
    snp_col=None,
    n_col=None,
    n_cas_col=None,
    n_con_col=None,
    a1_col=None,
    a2_col=None,
    p_col=None,
    frq_col=None,
    info_col=None,
    signed_sumstats=None,
    ignore=None,
    keep_maf=false,
    a1_inc=false,
    no_alleles=false,
    info_list=None,
    nstudy=None,
    nstudy_min=None
))]
#[allow(clippy::too_many_arguments)]
fn munge_sumstats(
    py: Python<'_>,
    sumstats: String,
    out: String,
    merge_alleles: Option<String>,
    daner: bool,
    daner_n: bool,
    n_min: f64,
    maf: f64,
    info_min: f64,
    n: Option<f64>,
    n_cas: Option<f64>,
    n_con: Option<f64>,
    snp_col: Option<String>,
    n_col: Option<String>,
    n_cas_col: Option<String>,
    n_con_col: Option<String>,
    a1_col: Option<String>,
    a2_col: Option<String>,
    p_col: Option<String>,
    frq_col: Option<String>,
    info_col: Option<String>,
    signed_sumstats: Option<String>,
    ignore: Option<String>,
    keep_maf: bool,
    a1_inc: bool,
    no_alleles: bool,
    info_list: Option<String>,
    nstudy: Option<String>,
    nstudy_min: Option<u64>,
) -> PyResult<NativeMungeSummary> {
    let opts = MungeOptions {
        daner,
        daner_n,
        n_min,
        maf,
        info_min,
        n,
        n_cas,
        n_con,
        snp_col,
        n_col,
        n_cas_col,
        n_con_col,
        a1_col,
        a2_col,
        p_col,
        frq_col,
        info_col,
        signed_sumstats,
        ignore,
        keep_maf,
        a1_inc,
        no_alleles,
        info_list,
        nstudy,
        nstudy_min,
    };
    let out_path = format!("{}.sumstats.gz", out);
    let summary = py
        .detach(move || {
            munge_sumstats_from_files(&sumstats, opts, merge_alleles.as_deref(), &out_path)
        })
        .map_err(value_error)?;
    Ok(summary.into())
}

#[pymodule]
fn _native(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<NativeH2Result>()?;
    module.add_class::<NativePartitionedH2Result>()?;
    module.add_class::<NativeRgResult>()?;
    module.add_class::<NativeLdScoreResult>()?;
    module.add_class::<NativeH2FileResult>()?;
    module.add_class::<NativeMungeSummary>()?;
    module.add_function(wrap_pyfunction!(fit_h2, module)?)?;
    module.add_function(wrap_pyfunction!(fit_h2_partitioned, module)?)?;
    module.add_function(wrap_pyfunction!(fit_rg, module)?)?;
    module.add_function(wrap_pyfunction!(fit_rg_partitioned, module)?)?;
    module.add_function(wrap_pyfunction!(compute_ld_scores_from_bytes, module)?)?;
    module.add_function(wrap_pyfunction!(estimate_h2, module)?)?;
    module.add_function(wrap_pyfunction!(estimate_rg, module)?)?;
    module.add_function(wrap_pyfunction!(munge_sumstats, module)?)?;
    Ok(())
}
