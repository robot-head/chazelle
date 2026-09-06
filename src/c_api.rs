//! C/C++ bindings for Chazelle polygon triangulation library.

use std::ffi::c_char;
use std::slice;
use crate::geometry::Point;
use crate::monotone::TriangulationError;
use crate::triangulate_points;

#[repr(C)]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ChazellePoint {
    pub x: f64,
    pub y: f64,
}

#[repr(C)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ChazelleTriangle {
    pub a: usize,
    pub b: usize,
    pub c: usize,
}

#[repr(C)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ChazelleStatus {
    Success = 0,
    ErrorInvalidArgument = 1,
    ErrorPolygonTooSmall = 2,
    ErrorDegenerate = 3,
    ErrorInternal = 4,
}

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
        let pts_slice = unsafe { slice::from_raw_parts(points, num_points) };
        let poly: Vec<Point> = pts_slice.iter().map(|p| Point::new(p.x, p.y)).collect();

        match triangulate_points(&poly) {
            Ok(tris) => {
                let count = tris.len();
                let mut c_tris: Vec<ChazelleTriangle> = Vec::with_capacity(count);
                for t in tris {
                    c_tris.push(ChazelleTriangle {
                        a: t[0],
                        b: t[1],
                        c: t[2],
                    });
                }
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
