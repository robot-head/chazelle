//! C/C++ bindings for Chazelle polygon triangulation library with algorithm selection.

use std::ffi::c_char;
use std::slice;
use zerocopy::{FromBytes, IntoBytes, KnownLayout, Immutable};
use crate::geometry::Point;
use crate::monotone::TriangulationError;
use crate::{Algorithm, triangulate_points_with_algorithm};

/// 2D point representation matching C layout and google/zerocopy traits.
#[repr(C)]
#[derive(Debug, Clone, Copy, PartialEq, FromBytes, IntoBytes, KnownLayout, Immutable)]
pub struct ChazellePoint {
    pub x: f64,
    pub y: f64,
}

impl ChazellePoint {
    pub const fn new(x: f64, y: f64) -> Self {
        Self { x, y }
    }
}

/// Triangle representation matching C layout and google/zerocopy traits.
#[repr(C)]
#[derive(Debug, Clone, Copy, PartialEq, Eq, FromBytes, IntoBytes, KnownLayout, Immutable)]
pub struct ChazelleTriangle {
    pub a: usize,
    pub b: usize,
    pub c: usize,
}

impl ChazelleTriangle {
    pub const fn new(a: usize, b: usize, c: usize) -> Self {
        Self { a, b, c }
    }

    pub fn slice_as_bytes(tris: &[ChazelleTriangle]) -> &[u8] {
        tris.as_bytes()
    }

    pub fn slice_from_bytes(bytes: &[u8]) -> Option<&[ChazelleTriangle]> {
        <[ChazelleTriangle]>::ref_from_bytes(bytes).ok()
    }
}

#[repr(C)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ChazelleStatus {
    Success = 0,
    ErrorInvalidArgument = 1,
    ErrorPolygonTooSmall = 2,
    ErrorDegenerate = 3,
    ErrorInternal = 4,
    ErrorBufferTooSmall = 5,
}

/// Choice of polygon triangulation algorithm for FFI callers.
#[repr(C)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ChazelleAlgorithm {
    /// Bernard Chazelle's deterministic linear-time O(n) algorithm.
    Chazelle = 0,
    /// Classic plane-sweep monotone polygon decomposition in O(n log n).
    MonotoneSweep = 1,
    /// Raimund Seidel's randomized incremental algorithm in O(n log* n).
    Seidel = 2,
}

impl From<ChazelleAlgorithm> for Algorithm {
    fn from(a: ChazelleAlgorithm) -> Self {
        match a {
            ChazelleAlgorithm::Chazelle => Algorithm::Chazelle,
            ChazelleAlgorithm::MonotoneSweep => Algorithm::MonotoneSweep,
            ChazelleAlgorithm::Seidel => Algorithm::Seidel,
        }
    }
}

/// In-place zero-allocation triangulation with algorithm selection.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn chazelle_triangulate_into_with_algorithm(
    points: *const ChazellePoint,
    num_points: usize,
    out_triangles: *mut ChazelleTriangle,
    out_capacity: usize,
    out_num_triangles: *mut usize,
    algorithm: ChazelleAlgorithm,
) -> ChazelleStatus {
    if points.is_null() || out_triangles.is_null() || out_num_triangles.is_null() {
        return ChazelleStatus::ErrorInvalidArgument;
    }

    if num_points < 3 {
        return ChazelleStatus::ErrorPolygonTooSmall;
    }

    let required_triangles = num_points - 2;
    if out_capacity < required_triangles {
        return ChazelleStatus::ErrorBufferTooSmall;
    }

    let res = std::panic::catch_unwind(|| {
        let pts_slice = unsafe { slice::from_raw_parts(points as *const Point, num_points) };

        match triangulate_points_with_algorithm(pts_slice, algorithm.into()) {
            Ok(tris) => {
                let count = tris.len();
                if count > out_capacity {
                    return ChazelleStatus::ErrorBufferTooSmall;
                }
                unsafe {
                    let out_slice = slice::from_raw_parts_mut(out_triangles, count);
                    for (i, t) in tris.iter().enumerate() {
                        out_slice[i] = ChazelleTriangle {
                            a: t[0],
                            b: t[1],
                            c: t[2],
                        };
                    }
                    *out_num_triangles = count;
                }
                ChazelleStatus::Success
            }
            Err(TriangulationError::PolygonTooSmall) => ChazelleStatus::ErrorPolygonTooSmall,
            Err(TriangulationError::DegeneratePolygon) => ChazelleStatus::ErrorDegenerate,
            Err(TriangulationError::InternalError(_)) => ChazelleStatus::ErrorInternal,
        }
    });

    match res {
        Ok(status) => status,
        Err(_) => ChazelleStatus::ErrorInternal,
    }
}

/// In-place zero-allocation triangulation using default Chazelle algorithm.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn chazelle_triangulate_into(
    points: *const ChazellePoint,
    num_points: usize,
    out_triangles: *mut ChazelleTriangle,
    out_capacity: usize,
    out_num_triangles: *mut usize,
) -> ChazelleStatus {
    unsafe {
        chazelle_triangulate_into_with_algorithm(
            points,
            num_points,
            out_triangles,
            out_capacity,
            out_num_triangles,
            ChazelleAlgorithm::Chazelle,
        )
    }
}

/// Triangulate with algorithm selection and dynamic memory allocation.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn chazelle_triangulate_with_algorithm(
    points: *const ChazellePoint,
    num_points: usize,
    out_triangles: *mut *mut ChazelleTriangle,
    out_num_triangles: *mut usize,
    algorithm: ChazelleAlgorithm,
) -> ChazelleStatus {
    if points.is_null() || out_triangles.is_null() || out_num_triangles.is_null() {
        return ChazelleStatus::ErrorInvalidArgument;
    }

    if num_points < 3 {
        return ChazelleStatus::ErrorPolygonTooSmall;
    }

    let res = std::panic::catch_unwind(|| {
        let pts_slice = unsafe { slice::from_raw_parts(points as *const Point, num_points) };

        match triangulate_points_with_algorithm(pts_slice, algorithm.into()) {
            Ok(mut tris) => {
                let count = tris.len();
                let c_tris: Vec<ChazelleTriangle> = unsafe {
                    let ptr = tris.as_mut_ptr() as *mut ChazelleTriangle;
                    let len = tris.len();
                    let cap = tris.capacity();
                    std::mem::forget(tris);
                    Vec::from_raw_parts(ptr, len, cap)
                };

                let mut boxed = c_tris.into_boxed_slice();
                let ptr = boxed.as_mut_ptr();
                std::mem::forget(boxed);

                unsafe {
                    *out_triangles = ptr;
                    *out_num_triangles = count;
                }
                ChazelleStatus::Success
            }
            Err(TriangulationError::PolygonTooSmall) => ChazelleStatus::ErrorPolygonTooSmall,
            Err(TriangulationError::DegeneratePolygon) => ChazelleStatus::ErrorDegenerate,
            Err(TriangulationError::InternalError(_)) => ChazelleStatus::ErrorInternal,
        }
    });

    match res {
        Ok(status) => status,
        Err(_) => ChazelleStatus::ErrorInternal,
    }
}

/// Standard allocating C API using default Chazelle algorithm.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn chazelle_triangulate(
    points: *const ChazellePoint,
    num_points: usize,
    out_triangles: *mut *mut ChazelleTriangle,
    out_num_triangles: *mut usize,
) -> ChazelleStatus {
    unsafe {
        chazelle_triangulate_with_algorithm(
            points,
            num_points,
            out_triangles,
            out_num_triangles,
            ChazelleAlgorithm::Chazelle,
        )
    }
}

#[unsafe(no_mangle)]
pub unsafe extern "C" fn chazelle_free_triangles(
    triangles: *mut ChazelleTriangle,
    num_triangles: usize,
) {
    if !triangles.is_null() && num_triangles > 0 {
        unsafe {
            let _ = Vec::from_raw_parts(triangles, num_triangles, num_triangles);
        }
    }
}

#[unsafe(no_mangle)]
pub extern "C" fn chazelle_version() -> *const c_char {
    c"0.1.0".as_ptr()
}
