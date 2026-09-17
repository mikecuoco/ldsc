from __future__ import annotations

import gzip
from dataclasses import FrozenInstanceError
from math import isfinite, sqrt

import numpy as np
import pytest

from ldsc_rs import (
    compute_ld_scores_from_bytes,
    estimate_h2,
    estimate_rg,
    fit_h2,
    fit_h2_partitioned,
    fit_rg,
    fit_rg_partitioned,
    munge_sumstats,
)


def synthetic_h2_columns(n: int = 100, h2: float = 0.25):
    ref_ld = [1.0 + 0.02 * i for i in range(n)]
    sample_size = [10_000.0] * n
    m_snps = 1_000.0
    chi2 = [1.0 + sample_size[i] * h2 * ref_ld[i] / m_snps for i in range(n)]
    weight_ld = [max(1.0, value) for value in ref_ld]
    return chi2, ref_ld, weight_ld, sample_size, m_snps


def synthetic_partitioned_h2_columns(n: int = 100, h2_1: float = 0.15, h2_2: float = 0.10):
    # ref_ld1 (linear in i) and ref_ld2 (an unrelated modular pattern) are
    # deliberately non-collinear so the K=2 design matrix is well-conditioned.
    ref_ld1 = [1.0 + 0.02 * i for i in range(n)]
    ref_ld2 = [1.0 + 0.03 * ((i * 7) % 13) for i in range(n)]
    sample_size = [10_000.0] * n
    m1, m2 = 500.0, 500.0
    chi2 = [
        1.0
        + sample_size[i] * (h2_1 * ref_ld1[i] / m1 + h2_2 * ref_ld2[i] / m2)
        for i in range(n)
    ]
    weight_ld = [max(1.0, ref_ld1[i] + ref_ld2[i]) for i in range(n)]
    ref_ld = [[a, b] for a, b in zip(ref_ld1, ref_ld2)]
    return chi2, ref_ld, weight_ld, sample_size, [m1, m2]


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


def test_fit_h2_numpy_input_matches_list_input():
    chi2, ref_ld, weight_ld, sample_size, m_snps = synthetic_h2_columns()
    kwargs = dict(m_snps=m_snps, n_blocks=10, intercept=1.0)

    from_lists = fit_h2(chi2, ref_ld, weight_ld, sample_size, **kwargs)
    from_arrays = fit_h2(
        np.array(chi2),
        np.array(ref_ld),
        np.array(weight_ld),
        np.array(sample_size),
        **kwargs,
    )

    assert from_arrays.heritability.value == pytest.approx(from_lists.heritability.value)
    assert from_arrays.intercept.value == from_lists.intercept.value


def test_fit_h2_partitioned_recovers_total_and_sums_correctly():
    chi2, ref_ld, weight_ld, sample_size, m_vec = synthetic_partitioned_h2_columns()

    result = fit_h2_partitioned(
        chi2,
        ref_ld,
        weight_ld,
        sample_size,
        m_vec=m_vec,
        n_blocks=10,
        intercept=1.0,
    )

    assert len(result.per_annotation) == 2
    assert result.heritability.value == pytest.approx(0.25, abs=1e-6)
    total_from_parts = sum(e.value for e in result.per_annotation)
    assert total_from_parts == pytest.approx(result.heritability.value, abs=1e-9)
    assert result.per_annotation[0].value == pytest.approx(0.15, abs=0.05)
    assert result.per_annotation[1].value == pytest.approx(0.10, abs=0.05)
    assert result.intercept.value == 1.0
    assert result.intercept.standard_error is None
    assert result.m_values == tuple(m_vec)
    assert result.n_snps == len(chi2)


def test_fit_h2_partitioned_accepts_numpy_2d_ref_ld():
    chi2, ref_ld, weight_ld, sample_size, m_vec = synthetic_partitioned_h2_columns()

    from_lists = fit_h2_partitioned(
        chi2, ref_ld, weight_ld, sample_size, m_vec=m_vec, n_blocks=10, intercept=1.0
    )
    from_arrays = fit_h2_partitioned(
        np.array(chi2),
        np.array(ref_ld),
        np.array(weight_ld),
        np.array(sample_size),
        m_vec=np.array(m_vec),
        n_blocks=10,
        intercept=1.0,
    )

    assert from_arrays.heritability.value == pytest.approx(from_lists.heritability.value)
    assert from_arrays.per_annotation[0].value == pytest.approx(
        from_lists.per_annotation[0].value
    )
    assert from_arrays.per_annotation[1].value == pytest.approx(
        from_lists.per_annotation[1].value
    )


def test_fit_rg_partitioned_identical_traits_is_one():
    chi2, ref_ld, weight_ld, sample_size, m_vec = synthetic_partitioned_h2_columns()
    z = [sqrt(value) for value in chi2]

    result = fit_rg_partitioned(
        z,
        z,
        ref_ld,
        weight_ld,
        sample_size,
        sample_size,
        m_vec=m_vec,
        n_blocks=10,
        two_step=None,
        intercept_h2=(1.0, 1.0),
    )

    assert result.correlation.value == pytest.approx(1.0, abs=1e-6)
    assert result.h2_1.heritability.value == pytest.approx(0.25, abs=1e-6)
    assert result.h2_2.heritability.value == pytest.approx(0.25, abs=1e-6)
    assert result.n_snps == len(z)


def _write_lines(path, lines):
    path.write_text("\n".join(lines) + "\n")


def test_munge_sumstats_writes_output_file(tmp_path):
    sumstats_path = tmp_path / "raw.sumstats"
    _write_lines(
        sumstats_path,
        [
            "SNP\tA1\tA2\tN\tZ",
            "rs1\tA\tG\t10000\t1.5",
            "rs2\tA\tG\t10000\t-2.0",
            "rs3\tA\tG\t10000\t0.5",
        ],
    )
    out_prefix = str(tmp_path / "munged")

    summary = munge_sumstats(str(sumstats_path), out_prefix)

    assert summary.rows_in == 3
    assert summary.rows_out == 3
    assert summary.duplicates_removed == 0
    assert summary.out_path == f"{out_prefix}.sumstats.gz"
    with gzip.open(summary.out_path, "rt") as f:
        header = f.readline().split()
    assert header == ["SNP", "A1", "A2", "Z", "N"]


def test_estimate_h2_scalar_matches_fit_h2(tmp_path):
    chi2, ref_ld, weight_ld, sample_size, m_snps = synthetic_h2_columns(n=60)
    snps = [f"rs{i}" for i in range(len(chi2))]
    z = [sqrt(v) for v in chi2]

    sumstats_path = tmp_path / "trait.sumstats"
    _write_lines(
        sumstats_path,
        ["SNP\tN\tZ"] + [f"{s}\t{n}\t{zv}" for s, n, zv in zip(snps, sample_size, z)],
    )
    ref_ld_path = tmp_path / "trait.l2.ldscore"
    _write_lines(
        ref_ld_path,
        ["SNP\tL2"] + [f"{s}\t{v}" for s, v in zip(snps, ref_ld)],
    )
    w_ld_path = tmp_path / "trait.w_ld.ldscore"
    _write_lines(
        w_ld_path,
        ["SNP\tL2"] + [f"{s}\t{v}" for s, v in zip(snps, weight_ld)],
    )

    result = estimate_h2(
        str(sumstats_path),
        ref_ld=str(ref_ld_path),
        w_ld=str(w_ld_path),
        m_snps=m_snps,
        n_blocks=10,
        intercept_h2=1.0,
    )
    expected = fit_h2(
        chi2, ref_ld, weight_ld, sample_size, m_snps=m_snps, n_blocks=10, intercept=1.0
    )

    assert result.heritability.value == pytest.approx(expected.heritability.value)
    assert result.mean_chi2 is not None
    assert result.l2_cols is None
    assert result.n_snps == len(chi2)
    assert result.liability_heritability is None


def test_estimate_rg_matches_fit_rg(tmp_path):
    chi2, ref_ld, weight_ld, sample_size, m_snps = synthetic_h2_columns(n=60)
    snps = [f"rs{i}" for i in range(len(chi2))]
    z = [sqrt(v) for v in chi2]

    ref_ld_path = tmp_path / "trait.l2.ldscore"
    _write_lines(
        ref_ld_path,
        ["SNP\tL2"] + [f"{s}\t{v}" for s, v in zip(snps, ref_ld)],
    )
    w_ld_path = tmp_path / "trait.w_ld.ldscore"
    _write_lines(
        w_ld_path,
        ["SNP\tL2"] + [f"{s}\t{v}" for s, v in zip(snps, weight_ld)],
    )

    def write_trait(name):
        path = tmp_path / f"{name}.sumstats"
        _write_lines(
            path,
            ["SNP\tN\tZ"] + [f"{s}\t{n}\t{zv}" for s, n, zv in zip(snps, sample_size, z)],
        )
        return path

    trait1 = write_trait("trait1")
    trait2 = write_trait("trait2")

    results = estimate_rg(
        [str(trait1), str(trait2)],
        ref_ld=str(ref_ld_path),
        w_ld=str(w_ld_path),
        m_snps=m_snps,
        n_blocks=10,
        no_check_alleles=True,
        no_intercept=True,
    )
    expected = fit_rg(
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
        intercept_gencov=0.0,
    )

    assert len(results) == 1
    assert results[0].correlation.value == pytest.approx(expected.correlation.value, abs=1e-6)
