//! NumPy-aware parameter types for the PyO3 boundary.
//!
//! PyO3's default `Vec<f64>: FromPyObject` extracts a Python sequence by
//! iterating and converting element-by-element — fine for small inputs, but
//! a real, avoidable cost for genome-scale arrays (millions of SNPs).
//! [`FloatColumn`]/[`FloatMatrix`] instead borrow a NumPy `float64` array's
//! buffer directly when one is passed, falling back to the generic
//! sequence path (a plain Python list, unchanged behavior) otherwise.
//!
//! This does not make the eventual `ColF`/`MatF` (`faer::Mat`) allocation
//! zero-copy — `faer::Mat` always owns its buffer — it only skips PyO3's
//! slow per-element extraction on the way in.

use numpy::{PyReadonlyArray1, PyReadonlyArray2};
use pyo3::Borrowed;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

/// A 1-D numeric column: either a borrowed NumPy `float64` array, or a
/// `Vec<f64>` copied from any other Python sequence.
pub enum FloatColumn<'py> {
    Array(PyReadonlyArray1<'py, f64>),
    Owned(Vec<f64>),
}

impl FloatColumn<'_> {
    pub fn to_vec(&self) -> PyResult<Vec<f64>> {
        match self {
            FloatColumn::Array(arr) => arr
                .as_slice()
                .map(<[f64]>::to_vec)
                .or_else(|_| Ok(arr.as_array().iter().copied().collect::<Vec<f64>>())),
            FloatColumn::Owned(v) => Ok(v.clone()),
        }
    }
}

impl<'py> FromPyObject<'_, 'py> for FloatColumn<'py> {
    type Error = PyErr;

    fn extract(obj: Borrowed<'_, 'py, PyAny>) -> PyResult<Self> {
        if let Ok(arr) = obj.extract::<PyReadonlyArray1<'py, f64>>() {
            return Ok(FloatColumn::Array(arr));
        }
        Ok(FloatColumn::Owned(obj.extract()?))
    }
}

/// A 2-D numeric matrix, row-major (shape `(n_rows, k_cols)`): either a
/// borrowed NumPy `float64` array, or a `Vec<Vec<f64>>` copied from any
/// other Python sequence of sequences.
pub enum FloatMatrix<'py> {
    Array(PyReadonlyArray2<'py, f64>),
    Owned(Vec<Vec<f64>>),
}

impl FloatMatrix<'_> {
    /// Transpose into `k` column vectors, one per annotation, each of
    /// length `n_rows`. Errors if the owned (list-of-rows) form has
    /// inconsistent row lengths.
    pub fn to_columns(&self) -> PyResult<Vec<Vec<f64>>> {
        match self {
            FloatMatrix::Array(arr) => {
                let view = arr.as_array();
                let (n, k) = (view.shape()[0], view.shape()[1]);
                let mut cols: Vec<Vec<f64>> = (0..k).map(|_| Vec::with_capacity(n)).collect();
                for i in 0..n {
                    for (j, col) in cols.iter_mut().enumerate() {
                        col.push(view[[i, j]]);
                    }
                }
                Ok(cols)
            }
            FloatMatrix::Owned(rows) => {
                if rows.is_empty() {
                    return Ok(Vec::new());
                }
                let k = rows[0].len();
                for row in rows {
                    if row.len() != k {
                        return Err(PyValueError::new_err(
                            "all rows of a 2-D input must have the same length",
                        ));
                    }
                }
                let mut cols: Vec<Vec<f64>> =
                    (0..k).map(|_| Vec::with_capacity(rows.len())).collect();
                for row in rows {
                    for (j, col) in cols.iter_mut().enumerate() {
                        col.push(row[j]);
                    }
                }
                Ok(cols)
            }
        }
    }
}

impl<'py> FromPyObject<'_, 'py> for FloatMatrix<'py> {
    type Error = PyErr;

    fn extract(obj: Borrowed<'_, 'py, PyAny>) -> PyResult<Self> {
        if let Ok(arr) = obj.extract::<PyReadonlyArray2<'py, f64>>() {
            return Ok(FloatMatrix::Array(arr));
        }
        Ok(FloatMatrix::Owned(obj.extract()?))
    }
}
