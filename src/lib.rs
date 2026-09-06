//! # Chazelle's Linear Time Polygon Triangulation Algorithm
//!
//! Implementation of Bernard Chazelle's deterministic $O(n)$ algorithm for triangulating
//! a simple polygon, as described in:
//!
//! > Bernard Chazelle, "Triangulating a Simple Polygon in Linear Time",
//! > *Discrete & Computational Geometry* 6:485-524 (1991).
//! > <https://www.cs.princeton.edu/~chazelle/pubs/polygon-triang.pdf>
//!
//! ## Overview of the Algorithm
//!
//! The algorithm computes a horizontal visibility map (trapezoidal map) of the simple polygon
//! through two balanced phases:
//!
//! 1. **Up-Phase (Bottom-Up, Section 4.1):**
//!    The polygon boundary is partitioned into dyadic chains across grades $\lambda = 0, \dots, p$.
//!    For each chain, canonical conformal submaps with bounded granularity are computed and merged
//!    along balanced binary trees using local ray shooting and the polygon-cutting theorem.
//!
//! 2. **Merging Two Submaps (Section 3):**
//!    - **Fusion (Section 3.1):** Boundary-walking procedure using local shooting in conformal regions.
//!    - **Restoring Conformality (Section 3.2):** Regions with $> 4$ arcs are split using chords discovered
//!      via centroid tree search.
//!    - **Maintaining Granularity (Section 3.3):** Superfluous chords are removed to keep submap size bounded.
//!
//! 3. **Down-Phase (Top-Down, Section 4.2):**
//!    The canonical submap for the entire polygon is refined top-down, restoring missing visibility chords
//!    until the full horizontal visibility map $V(P)$ emerges at granularity 1.
//!
//! 4. **Triangulation (Section 4.3):**
//!    The horizontal chords from $V(P)$ decompose the polygon into monotone mountains, which are
//!    triangulated in linear time using a stack.
//!
//! Rust Edition 2024 is used throughout, along with complete C/C++ FFI bindings.

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

pub use geometry::Point;
pub use monotone::TriangulationError;
pub use c_api::{ChazellePoint, ChazelleTriangle, ChazelleStatus};

use geometry::signed_polygon_area;
use up_phase::UpPhaseHierarchy;
use down_phase::run_down_phase;
use monotone::triangulate_from_visibility_map;

/// Triangulate a simple polygon given as a slice of (x, y) coordinate pairs.
///
/// The vertices can be in clockwise or counter-clockwise order.
///
/// # Returns
/// A vector of triangles, where each triangle is an array of 3 indices referencing
/// the vertices of the input polygon.
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

/// Triangulate a simple polygon given as a slice of [`Point`]s.
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

    // Normalize polygon to CCW orientation while keeping track of original indices
    let is_ccw = raw_area > 0.0;
    let (ccw_polygon, index_map): (Vec<Point>, Vec<usize>) = if is_ccw {
        let pts = polygon.to_vec();
        let map = (0..n).collect();
        (pts, map)
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
