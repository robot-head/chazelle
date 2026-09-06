//! # Chazelle's Linear Time Polygon Triangulation Algorithm
//!
//! Implementation of Bernard Chazelle's deterministic $O(n)$ algorithm for triangulating
//! a simple polygon, as described in:
//!
//! > Bernard Chazelle, "Triangulating a Simple Polygon in Linear Time",
//! > *Discrete & Computational Geometry* 6:485-524 (1991).
//! > <https://www.cs.princeton.edu/~chazelle/pubs/polygon-triang.pdf>
//!
//! ## Optimizations with `google/zerocopy`
//!
//! This library integrates `zerocopy` to achieve zero-copy data flow across:
//! - **Memory Casts:** Safe, zero-copy casting between `ChazellePoint`, `Point`, and raw byte buffers.
//! - **FFI Boundary:** Zero-copy passing of vertex slices from C/C++ without heap allocation or copying.
//! - **Zero-Allocation Output:** In-place triangulation via [`chazelle_triangulate_into`] directly writing
//!   into caller-provided buffers.
//! - **Buffer Transmutation:** In-place transmutation of triangle vectors into C FFI types without element copies.

pub mod geometry;
pub mod double_boundary;
pub mod submap;
pub mod fusion;
pub mod conformality;
pub mod oracles;
pub mod up_phase;
pub mod down_phase;
pub mod monotone;
pub mod c_api;
pub mod datasets;

pub use geometry::Point;
pub use monotone::TriangulationError;
pub use c_api::{ChazellePoint, ChazelleTriangle, ChazelleStatus, chazelle_triangulate, chazelle_triangulate_into, chazelle_free_triangles};

use geometry::signed_polygon_area;
use up_phase::UpPhaseHierarchy;
use down_phase::run_down_phase;
use monotone::triangulate_from_visibility_map;

/// Triangulate a simple polygon given as a slice of (x, y) coordinate pairs.
///
/// # Example
/// ```
/// use chazelle::triangulate;
///
/// let poly = vec![(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)];
/// let triangles = triangulate(&poly).unwrap();
/// assert_eq!(triangles.len(), 2);
/// ```
pub fn triangulate(polygon: &[(f64, f64)]) -> Result<Vec<[usize; 3]>, TriangulationError> {
    let pts: Vec<Point> = polygon.iter().map(|&(x, y)| Point::new(x, y)).collect();
    triangulate_points(&pts)
}

/// Triangulate a simple polygon given as a slice of [`Point`]s with zero unnecessary copies.
pub fn triangulate_points(polygon: &[Point]) -> Result<Vec<[usize; 3]>, TriangulationError> {
    let n = polygon.len();
    if n < 3 {
        return Err(TriangulationError::PolygonTooSmall);
    }
    if n == 3 {
        return Ok(vec![[0, 1, 2]]);
    }

    let raw_area = signed_polygon_area(polygon);
    if raw_area.abs() <= Point::EPSILON {
        return Err(TriangulationError::DegeneratePolygon);
    }

    let is_ccw = raw_area > 0.0;
    let (ccw_polygon, index_map): (Vec<Point>, Vec<usize>) = if is_ccw {
        (polygon.to_vec(), (0..n).collect())
    } else {
        let mut pts = Vec::with_capacity(n);
        let mut map = Vec::with_capacity(n);
        for i in (0..n).rev() {
            pts.push(polygon[i]);
            map.push(i);
        }
        (pts, map)
    };

    // Step 1: Execute bottom-up Up-Phase
    let hierarchy = UpPhaseHierarchy::build(&ccw_polygon);

    // Step 2: Execute top-down Down-Phase to obtain full visibility map V(P)
    let visibility_map = run_down_phase(&hierarchy, &ccw_polygon);

    // Step 3: Triangulate from the visibility map
    let ccw_triangles = triangulate_from_visibility_map(&ccw_polygon, &visibility_map)?;

    // Remap indices back to the original polygon ordering
    let mapped_triangles = ccw_triangles
        .into_iter()
        .map(|[a, b, c]| [index_map[a], index_map[b], index_map[c]])
        .collect();

    Ok(mapped_triangles)
}

/// Zero-copy parsing and triangulation from a raw byte buffer containing `Point` structs.
pub fn triangulate_from_bytes(bytes: &[u8]) -> Result<Vec<ChazelleTriangle>, TriangulationError> {
    let pts = Point::slice_from_bytes(bytes).ok_or_else(|| {
        TriangulationError::InternalError("Invalid byte alignment or length for Point slice".into())
    })?;

    let mut tris = triangulate_points(pts)?;
    let c_tris = unsafe {
        let ptr = tris.as_mut_ptr() as *mut ChazelleTriangle;
        let len = tris.len();
        let cap = tris.capacity();
        std::mem::forget(tris);
        Vec::from_raw_parts(ptr, len, cap)
    };
    Ok(c_tris)
}
