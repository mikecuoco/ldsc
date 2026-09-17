use ldsc::h2::{H2Result, run_h2_ldsc};
use ldsc::l2::{L2Config, WindowMode, compute_l2_from_bytes};
use ldsc::la::col_from_vec;
use ldsc::regressions::{GeneticCorrelationResult, run_rg_ldsc};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyBytes;

fn value_error(error: impl std::fmt::Display) -> PyErr {
    PyValueError::new_err(format!("{error:#}"))
}

#[pyclass(frozen, get_all, module = "ldsc_rs._native")]
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
fn fit_h2(
    py: Python<'_>,
    chi2: Vec<f64>,
    ref_ld: Vec<f64>,
    weight_ld: Vec<f64>,
    sample_size: Vec<f64>,
    m_snps: f64,
    n_blocks: usize,
    two_step: Option<f64>,
    intercept: Option<f64>,
) -> PyResult<NativeH2Result> {
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
fn fit_rg(
    py: Python<'_>,
    z1: Vec<f64>,
    z2: Vec<f64>,
    ref_ld: Vec<f64>,
    weight_ld: Vec<f64>,
    sample_size1: Vec<f64>,
    sample_size2: Vec<f64>,
    m_snps: f64,
    n_blocks: usize,
    two_step: Option<f64>,
    intercept_h2_1: Option<f64>,
    intercept_h2_2: Option<f64>,
    intercept_gencov: Option<f64>,
) -> PyResult<NativeRgResult> {
    let result = py
        .detach(move || {
            run_rg_ldsc(
                &col_from_vec(z1),
                &col_from_vec(z2),
                &col_from_vec(ref_ld),
                &col_from_vec(weight_ld),
                &col_from_vec(sample_size1),
                &col_from_vec(sample_size2),
                m_snps,
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

#[pymodule]
fn _native(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<NativeH2Result>()?;
    module.add_class::<NativeRgResult>()?;
    module.add_class::<NativeLdScoreResult>()?;
    module.add_function(wrap_pyfunction!(fit_h2, module)?)?;
    module.add_function(wrap_pyfunction!(fit_rg, module)?)?;
    module.add_function(wrap_pyfunction!(compute_ld_scores_from_bytes, module)?)?;
    Ok(())
}
