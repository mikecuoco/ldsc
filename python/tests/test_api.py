from __future__ import annotations

from dataclasses import FrozenInstanceError
from math import isfinite, sqrt

import pytest

from ldsc_rs import compute_ld_scores_from_bytes, fit_h2, fit_rg


def synthetic_h2_columns(n: int = 100, h2: float = 0.25):
    ref_ld = [1.0 + 0.02 * i for i in range(n)]
    sample_size = [10_000.0] * n
    m_snps = 1_000.0
    chi2 = [1.0 + sample_size[i] * h2 * ref_ld[i] / m_snps for i in range(n)]
    weight_ld = [max(1.0, value) for value in ref_ld]
    return chi2, ref_ld, weight_ld, sample_size, m_snps


def test_fit_h2_recovers_exact_linear_signal():
    chi2, ref_ld, weight_ld, sample_size, m_snps = synthetic_h2_columns()
    result = fit_h2(
        chi2,
        ref_ld,
        weight_ld,
        sample_size,
        m_snps=m_snps,
        n_blocks=10,
        intercept=1.0,
    )

    assert result.heritability.value == pytest.approx(0.25, abs=1e-8)
    assert result.intercept.value == 1.0
    assert result.intercept.standard_error is None
    assert result.ratio is None
    with pytest.raises(FrozenInstanceError):
        result.mean_chi2 = 1.0


def test_fit_h2_rejects_mismatched_columns():
    with pytest.raises(ValueError, match="length mismatch"):
        fit_h2([1.0, 2.0], [1.0], [1.0, 1.0], [100.0, 100.0], m_snps=10)

    with pytest.raises(ValueError, match="non-negative"):
        fit_h2([1.0, 2.0], [-1.0, 1.0], [1.0, 1.0], [100.0, 100.0], m_snps=10)


def test_fit_rg_identical_traits_is_one():
    chi2, ref_ld, weight_ld, sample_size, m_snps = synthetic_h2_columns()
    z = [sqrt(value) for value in chi2]
    result = fit_rg(
        z,
        z,
        ref_ld,
        weight_ld,
        sample_size,
        sample_size,
        m_snps=m_snps,
        n_blocks=10,
        two_step=None,
        intercept_h2=(1.0, 1.0),
    )

    assert result.correlation.value == pytest.approx(1.0, abs=1e-8)
    assert result.h2_1.heritability.value == pytest.approx(0.25, abs=1e-8)
    assert result.h2_2.heritability.value == pytest.approx(0.25, abs=1e-8)
    assert result.n_snps == len(z)
    assert isfinite(result.genetic_covariance.covariance.value)


def test_fit_rg_rejects_invalid_inputs():
    with pytest.raises(ValueError, match="length mismatch"):
        fit_rg([1.0, 2.0], [1.0], [1.0, 1.0], [1.0, 1.0], [100.0] * 2, [100.0] * 2, m_snps=10)
    with pytest.raises(ValueError, match="at least two observations"):
        fit_rg([1.0], [1.0], [1.0], [1.0], [100.0], [100.0], m_snps=10)
    with pytest.raises(ValueError, match="non-negative"):
        fit_rg([1.0, 1.0], [1.0, 1.0], [-1.0, 1.0], [1.0, 1.0], [100.0] * 2, [100.0] * 2, m_snps=10)
    with pytest.raises(ValueError, match="exactly two"):
        fit_rg([1.0, 1.0], [1.0, 1.0], [1.0, 1.0], [1.0, 1.0], [100.0] * 2, [100.0] * 2, m_snps=10, intercept_h2=(1.0,))


def test_compute_ld_scores_microbed():
    bed = bytes([0x6C, 0x1B, 0x01, 0xF0, 0xF0])
    bim = "1\trs1\t0\t100\tA\tG\n1\trs2\t0\t200\tA\tG\n"
    fam = (
        "F1\tI1\t0\t0\t1\t-9\n"
        "F2\tI2\t0\t0\t1\t-9\n"
        "F3\tI3\t0\t0\t1\t-9\n"
        "F4\tI4\t0\t0\t1\t-9\n"
    )

    result = compute_ld_scores_from_bytes(
        bed,
        bim,
        fam,
        window=("kb", 100.0),
        chunk_size=2,
    )

    assert result.snp == ("rs1", "rs2")
    assert result.maf == pytest.approx((0.5, 0.5))
    assert result.ld_score == pytest.approx((2.0, 2.0))
    assert result.wall_seconds >= 0.0


def test_compute_ld_scores_rejects_truncated_bed():
    bed = bytes([0x6C, 0x1B, 0x01, 0xF0])
    bim = "1\trs1\t0\t100\tA\tG\n1\trs2\t0\t200\tA\tG\n"
    fam = "F1\tI1\t0\t0\t1\t-9\n"

    with pytest.raises(ValueError, match="BED"):
        compute_ld_scores_from_bytes(bed, bim, fam)


def test_compute_ld_scores_validates_configuration():
    bed = bytes([0x6C, 0x1B, 0x01, 0xF0])
    bim = "1\trs1\t0\t100\tA\tG\n"
    fam = "F1\tI1\t0\t0\t1\t-9\n"

    with pytest.raises(ValueError, match="exactly"):
        compute_ld_scores_from_bytes(bed, bim, fam, window=("kb",))
    with pytest.raises(ValueError, match="integer"):
        compute_ld_scores_from_bytes(bed, bim, fam, window=("snps", 1.5))
    with pytest.raises(ValueError, match="dtype"):
        compute_ld_scores_from_bytes(bed, bim, fam, dtype="float16")
    with pytest.raises(ValueError, match="dimension"):
        compute_ld_scores_from_bytes(bed, bim, fam, sketch=0)
