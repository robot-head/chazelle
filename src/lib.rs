//! # Chazelle's Linear Time Polygon Triangulation Algorithm
//!
//! Implementation of Bernard Chazelle's deterministic $O(n)$ algorithm for triangulating
//! a simple polygon, as described in:
//!
//! > Bernard Chazelle, "Triangulating a Simple Polygon in Linear Time",
//! > *Discrete & Computational Geometry* 6:485-524 (1991).
//! > <https://www.cs.princeton.edu/~chazelle/pubs/polygon-triang.pdf>
//!
//! ## Alternative Triangulation Solvers Available:
//! - **Classic Plane Sweep Monotone Decomposition ($O(n \log n)$):**
//!   Textbook algorithm (de Berg et al., Preparata & Shamos) using plane-sweep vertex classification
//!   and linear-time stack-based $y$-monotone polygon triangulation.
//! - **Seidel's Randomized Algorithm ($O(n \log^* n)$):**
//!   Raimund Seidel's incremental randomized trapezoidal decomposition and monotone mountain triangulation.
//!
//! ## Optimizations with `google/zerocopy`
//! - Zero-copy memory casting between `ChazellePoint` and `Point` without heap copies.
//! - In-place triangulation via [`chazelle_triangulate_into`] directly writing into caller buffers.
//! - In-place buffer transmutation between `[usize; 3]` and `ChazelleTriangle`.
//! - Zero-copy binary serialization and deserialization via `Point::slice_from_bytes`.

pub mod geometry;
pub mod double_boundary;
pub mod submap;
pub mod fusion;
pub mod conformality;
pub mod oracles;
pub mod up_phase;
pub mod down_phase;
pub mod monotone;
pub mod monotone_sweep;
pub mod seidel;
pub mod datasets;
pub mod c_api;

pub use geometry::Point;
pub use monotone::TriangulationError;
pub use monotone_sweep::triangulate_monotone_sweep;
pub use seidel::triangulate_seidel;
pub use c_api::{
    ChazellePoint, ChazelleTriangle, ChazelleStatus, ChazelleAlgorithm,
    chazelle_triangulate, chazelle_triangulate_into, chazelle_free_triangles,
    chazelle_triangulate_with_algorithm, chazelle_triangulate_into_with_algorithm,
};

use geometry::signed_polygon_area;
use up_phase::UpPhaseHierarchy;
use monotone::triangulate_from_visibility_map;

/// Available polygon triangulation algorithms in the library.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Algorithm {
    /// Bernard Chazelle's deterministic linear-time $O(n)$ algorithm (default).
    #[default]
    Chazelle,
    /// Classic textbook plane-sweep monotone decomposition in $O(n \log n)$.
    MonotoneSweep,
    /// Raimund Seidel's randomized incremental algorithm in $O(n \log^* n)$.
    Seidel,
}

/// Triangulate a simple polygon given as a slice of (x, y) coordinate pairs using the default Chazelle algorithm.
pub fn triangulate(polygon: &[(f64, f64)]) -> Result<Vec<[usize; 3]>, TriangulationError> {
    let pts: Vec<Point> = polygon.iter().map(|&(x, y)| Point::new(x, y)).collect();
    triangulate_points(&pts)
}

/// Triangulate a simple polygon using a specified algorithm.
pub fn triangulate_with_algorithm(
    polygon: &[(f64, f64)],
    algorithm: Algorithm,
) -> Result<Vec<[usize; 3]>, TriangulationError> {
    let pts: Vec<Point> = polygon.iter().map(|&(x, y)| Point::new(x, y)).collect();
    triangulate_points_with_algorithm(&pts, algorithm)
}

/// Triangulate a simple polygon given as a slice of [`Point`]s using Chazelle's linear-time algorithm.
pub fn triangulate_points(polygon: &[Point]) -> Result<Vec<[usize; 3]>, TriangulationError> {
    triangulate_points_with_algorithm(polygon, Algorithm::Chazelle)
}

/// Triangulate a simple polygon given as a slice of [`Point`]s using the selected algorithm.
pub fn triangulate_points_with_algorithm(
    polygon: &[Point],
    algorithm: Algorithm,
) -> Result<Vec<[usize; 3]>, TriangulationError> {
    let n = polygon.len();
    if n < 3 {
        return Err(TriangulationError::PolygonTooSmall);
    }
    if n == 3 {
        return Ok(vec![[0, 1, 2]]);
    }

    match algorithm {
        Algorithm::Chazelle => {
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

            let oracle = oracles::RayShootingOracle::new(&ccw_polygon);
            let hierarchy = UpPhaseHierarchy::build_with_oracle(&ccw_polygon, &oracle);
            let visibility_map = down_phase::run_down_phase_with_oracle(&hierarchy, &ccw_polygon, &oracle);
            let ccw_triangles = triangulate_from_visibility_map(&ccw_polygon, &visibility_map)?;

            let mapped_triangles = ccw_triangles
                .into_iter()
                .map(|[a, b, c]| [index_map[a], index_map[b], index_map[c]])
                .collect();

            Ok(mapped_triangles)
        }
        Algorithm::MonotoneSweep => triangulate_monotone_sweep(polygon),
        Algorithm::Seidel => triangulate_seidel(polygon),
    }
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
