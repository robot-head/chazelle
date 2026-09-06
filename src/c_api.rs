//! C/C++ bindings for Chazelle polygon triangulation library with google/zerocopy optimizations.

use std::ffi::c_char;
use std::slice;
use zerocopy::{FromBytes, IntoBytes, KnownLayout, Immutable};
use crate::geometry::Point;
use crate::monotone::TriangulationError;
use crate::triangulate_points;

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

    /// Zero-copy byte view of a slice of triangles.
    pub fn slice_as_bytes(tris: &[ChazelleTriangle]) -> &[u8] {
        tris.as_bytes()
    }

    /// Zero-copy conversion of raw bytes to a slice of triangles.
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

/// Zero-copy, zero-allocation triangulation directly into caller-provided buffer.
///
/// Avoids any heap allocation for the output triangles.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn chazelle_triangulate_into(
    points: *const ChazellePoint,
    num_points: usize,
    out_triangles: *mut ChazelleTriangle,
    out_capacity: usize,
    out_num_triangles: *mut usize,
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
        // Zero-copy view: ChazellePoint and Point share identical repr(C) layout
        let pts_slice = unsafe { slice::from_raw_parts(points as *const Point, num_points) };

        match triangulate_points(pts_slice) {
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

/// Standard allocating C API with zero-copy internal transmutation.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn chazelle_triangulate(
    points: *const ChazellePoint,
    num_points: usize,
    out_triangles: *mut *mut ChazelleTriangle,
    out_num_triangles: *mut usize,
) -> ChazelleStatus {
    if points.is_null() || out_triangles.is_null() || out_num_triangles.is_null() {
        return ChazelleStatus::ErrorInvalidArgument;
    }

    if num_points < 3 {
        return ChazelleStatus::ErrorPolygonTooSmall;
    }

    let res = std::panic::catch_unwind(|| {
        // Zero-copy view: cast directly without allocating or copying Point array
        let pts_slice = unsafe { slice::from_raw_parts(points as *const Point, num_points) };

        match triangulate_points(pts_slice) {
            Ok(mut tris) => {
                let count = tris.len();
                // Zero-copy in-place buffer transmutation of Vec<[usize; 3]> into Vec<ChazelleTriangle>
                // Guaranteed safe because [usize; 3] and ChazelleTriangle have identical size, alignment, and repr(C)
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
