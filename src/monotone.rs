//! Triangulation of simple polygons from visibility chords / monotone decomposition (Section 4.3).

use crate::geometry::{Point, orient2d, point_in_triangle_ccw, segments_intersect_strict, signed_polygon_area};
use crate::submap::Submap;

/// Error type for triangulation failures.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum TriangulationError {
    PolygonTooSmall,
    DegeneratePolygon,
    InternalError(String),
}

impl std::fmt::Display for TriangulationError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            TriangulationError::PolygonTooSmall => write!(f, "Polygon must have at least 3 vertices"),
            TriangulationError::DegeneratePolygon => write!(f, "Polygon has zero area or degenerate vertices"),
            TriangulationError::InternalError(msg) => write!(f, "Triangulation error: {msg}"),
        }
    }
}

impl std::error::Error for TriangulationError {}

/// A diagonal connecting two non-adjacent vertices of the polygon.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Diagonal {
    pub u: usize,
    pub v: usize,
}

impl Diagonal {
    pub fn new(a: usize, b: usize) -> Self {
        if a < b {
            Self { u: a, v: b }
        } else {
            Self { u: b, v: a }
        }
    }
}

/// Checks if diagonal (u, v) is strictly valid (does not intersect polygon edges and lies inside).
pub fn is_valid_diagonal(u: usize, v: usize, polygon: &[Point]) -> bool {
    let n = polygon.len();
    if u == v || (u + 1) % n == v || (v + 1) % n == u {
        return false;
    }

    let pu = polygon[u];
    let pv = polygon[v];

    let min_x = pu.x.min(pv.x) - Point::EPSILON;
    let max_x = pu.x.max(pv.x) + Point::EPSILON;
    let min_y = pu.y.min(pv.y) - Point::EPSILON;
    let max_y = pu.y.max(pv.y) + Point::EPSILON;

    // Check that diagonal does not pass through any other polygon vertex
    for i in 0..n {
        if i == u || i == v {
            continue;
        }
        let p = polygon[i];
        if p.x < min_x || p.x > max_x || p.y < min_y || p.y > max_y {
            continue;
        }
        if orient2d(pu, pv, p).abs() <= Point::EPSILON {
            return false;
        }
    }

    // Check intersection with all polygon edges
    for i in 0..n {
        let j = (i + 1) % n;
        if i == u || i == v || j == u || j == v {
            continue;
        }
        let p1 = polygon[i];
        let p2 = polygon[j];
        if p1.x.min(p2.x) > max_x || p1.x.max(p2.x) < min_x || p1.y.min(p2.y) > max_y || p1.y.max(p2.y) < min_y {
            continue;
        }
        if segments_intersect_strict(pu, pv, p1, p2) {
            return false;
        }
    }

    // Check if midpoint of (pu, pv) is inside polygon
    let mid = Point::new(0.5 * (pu.x + pv.x), 0.5 * (pu.y + pv.y));
    point_inside_simple_polygon(mid, polygon)
}

/// Tests if point p is strictly inside a CCW simple polygon using ray casting.
pub fn point_inside_simple_polygon(p: Point, polygon: &[Point]) -> bool {
    let n = polygon.len();
    let mut inside = false;

    for i in 0..n {
        let j = (i + 1) % n;
        let p1 = polygon[i];
        let p2 = polygon[j];

        let min_seg_y = p1.y.min(p2.y);
        let max_seg_y = p1.y.max(p2.y);
        if p.y < min_seg_y - Point::EPSILON || p.y > max_seg_y + Point::EPSILON {
            continue;
        }

        // Check if point is on segment
        if (orient2d(p1, p2, p).abs() <= Point::EPSILON)
            && p.x >= p1.x.min(p2.x) - Point::EPSILON
            && p.x <= p1.x.max(p2.x) + Point::EPSILON
        {
            return false; // On boundary
        }

        if (p1.y > p.y) != (p2.y > p.y) {
            let intersect_x = p1.x + (p.y - p1.y) * (p2.x - p1.x) / (p2.y - p1.y);
            if p.x < intersect_x {
                inside = !inside;
            }
        }
    }

    inside
}

/// Triangulates a CCW simple polygon using the computed visibility map V(P) (Section 4.3).
/// Returns triangles as triples of indices into `polygon`.
pub fn triangulate_from_visibility_map(
    polygon: &[Point],
    visibility_map: &Submap,
) -> Result<Vec<[usize; 3]>, TriangulationError> {
    let n = polygon.len();
    if n < 3 {
        return Err(TriangulationError::PolygonTooSmall);
    }
    if n == 3 {
        return Ok(vec![[0, 1, 2]]);
    }

    let area = signed_polygon_area(polygon);
    if area <= Point::EPSILON {
        return Err(TriangulationError::DegeneratePolygon);
    }

    let mut valid_diagonals: std::collections::HashSet<Diagonal> = std::collections::HashSet::new();

    // 1. Extract potential diagonals from visibility map chords
    for chord in &visibility_map.chords {
        if let Some(u) = chord.origin_vertex {
            let prev_u = (u + n - 1) % n;
            let next_u = (u + 1) % n;
            // Only reflex / non-strictly-convex vertices need diagonals to resolve non-monotonicity
            if orient2d(polygon[prev_u], polygon[u], polygon[next_u]) > Point::EPSILON {
                continue;
            }

            let e1 = chord.hit_edge;
            let e2 = (chord.hit_edge + 1) % n;

            let d1 = (polygon[e1].x - chord.hit_pt.x).powi(2) + (polygon[e1].y - chord.hit_pt.y).powi(2);
            let d2 = (polygon[e2].x - chord.hit_pt.x).powi(2) + (polygon[e2].y - chord.hit_pt.y).powi(2);
            let (v_primary, v_secondary) = if d1 <= d2 { (e1, e2) } else { (e2, e1) };

            for &v in &[v_primary, v_secondary] {
                let diag = Diagonal::new(u, v);
                if !valid_diagonals.contains(&diag) && is_valid_diagonal(u, v, polygon) {
                    let mut crosses = false;
                    let min_x = polygon[diag.u].x.min(polygon[diag.v].x);
                    let max_x = polygon[diag.u].x.max(polygon[diag.v].x);
                    let min_y = polygon[diag.u].y.min(polygon[diag.v].y);
                    let max_y = polygon[diag.u].y.max(polygon[diag.v].y);

                    for d in &valid_diagonals {
                        let d_min_x = polygon[d.u].x.min(polygon[d.v].x);
                        let d_max_x = polygon[d.u].x.max(polygon[d.v].x);
                        let d_min_y = polygon[d.u].y.min(polygon[d.v].y);
                        let d_max_y = polygon[d.u].y.max(polygon[d.v].y);

                        if min_x > d_max_x || max_x < d_min_x || min_y > d_max_y || max_y < d_min_y {
                            continue;
                        }

                        if segments_intersect_strict(
                            polygon[diag.u],
                            polygon[diag.v],
                            polygon[d.u],
                            polygon[d.v],
                        ) {
                            crosses = true;
                            break;
                        }
                    }
                    if !crosses {
                        valid_diagonals.insert(diag);
                        break; // One diagonal per chord suffices
                    }
                }
            }
        }
    }

    // 2. Complete triangulation: decompose the polygon along the non-crossing diagonals
    let triangles = triangulate_polygon_with_diagonals(polygon, &valid_diagonals)?;

    Ok(triangles)
}

/// Decomposes the polygon using known non-crossing diagonals into subpolygons and triangulates each.
fn triangulate_polygon_with_diagonals(
    polygon: &[Point],
    diagonals: &std::collections::HashSet<Diagonal>,
) -> Result<Vec<[usize; 3]>, TriangulationError> {
    let n = polygon.len();
    let initial_chain: Vec<usize> = (0..n).collect();
    let mut subpolygons = vec![initial_chain];

    for diag in diagonals {
        let mut split_idx = None;
        let mut split_data = None;

        for (idx, sub) in subpolygons.iter().enumerate() {
            if let (Some(pos_u), Some(pos_v)) = (
                sub.iter().position(|&x| x == diag.u),
                sub.iter().position(|&x| x == diag.v),
            ) {
                let (first, second) = if pos_u < pos_v {
                    (pos_u, pos_v)
                } else {
                    (pos_v, pos_u)
                };

                let len = sub.len();
                if second > first + 1 && !(first == 0 && second == len - 1) {
                    let mut sub1 = Vec::with_capacity(second - first + 1);
                    sub1.extend_from_slice(&sub[first..=second]);

                    let mut sub2 = Vec::with_capacity(len - (second - first) + 1);
                    sub2.extend_from_slice(&sub[second..len]);
                    sub2.extend_from_slice(&sub[0..=first]);

                    split_idx = Some(idx);
                    split_data = Some((sub1, sub2));
                    break;
                }
            }
        }

        if let (Some(idx), Some((sub1, sub2))) = (split_idx, split_data) {
            subpolygons[idx] = sub1;
            subpolygons.push(sub2);
        }
    }

    // Triangulate each subpolygon: prefer fast linear monotone stack triangulation with area verification
    let mut all_triangles = Vec::new();
    for sub in subpolygons {
        let mut sub_pts = Vec::with_capacity(sub.len());
        for &idx in &sub {
            sub_pts.push(polygon[idx]);
        }
        let expected_area = signed_polygon_area(&sub_pts).abs();

        let mut success = false;
        if let Ok(tris) = crate::monotone_sweep::triangulate_monotone_piece(polygon, &sub) {
            let mut tri_area = 0.0;
            for &[a, b, c] in &tris {
                tri_area += signed_polygon_area(&[polygon[a], polygon[b], polygon[c]]).abs();
            }
            if (tri_area - expected_area).abs() <= 1e-4 {
                all_triangles.extend(tris);
                success = true;
            }
        }

        if !success {
            let tris = triangulate_simple_cycle(polygon, &sub)?;
            all_triangles.extend(tris);
        }
    }

    if all_triangles.len() != n - 2 {
        return Err(TriangulationError::InternalError(format!(
            "Expected {} triangles, got {}",
            n - 2,
            all_triangles.len()
        )));
    }

    Ok(all_triangles)
}

/// Triangulates a simple CCW cycle of vertex indices in linear time using ear clipping.
fn triangulate_simple_cycle(
    polygon: &[Point],
    cycle: &[usize],
) -> Result<Vec<[usize; 3]>, TriangulationError> {
    let mut remaining: Vec<usize> = cycle.to_vec();
    let mut triangles = Vec::new();

    while remaining.len() > 3 {
        let m = remaining.len();
        let mut ear_found = false;

        for i in 0..m {
            let prev = remaining[(i + m - 1) % m];
            let curr = remaining[i];
            let next = remaining[(i + 1) % m];

            let p_prev = polygon[prev];
            let p_curr = polygon[curr];
            let p_next = polygon[next];

            // Strict CCW convexity test
            let turn = orient2d(p_prev, p_curr, p_next);
            if turn <= Point::EPSILON {
                continue;
            }

            // Check if any other vertex of the remaining cycle is inside the candidate ear
            let mut has_interior_point = false;
            for j in 0..m {
                let test_v = remaining[j];
                if test_v == prev || test_v == curr || test_v == next {
                    continue;
                }
                let pt = polygon[test_v];
                if point_in_triangle_ccw(pt, p_prev, p_curr, p_next) {
                    has_interior_point = true;
                    break;
                }
            }

            if !has_interior_point {
                triangles.push([prev, curr, next]);
                remaining.remove(i);
                ear_found = true;
                break;
            }
        }

        if !ear_found {
            // Pick the convex vertex with minimal interior intrusion
            let mut best_i = None;
            let mut max_turn = 0.0;
            for i in 0..m {
                let prev = remaining[(i + m - 1) % m];
                let curr = remaining[i];
                let next = remaining[(i + 1) % m];
                let turn = orient2d(polygon[prev], polygon[curr], polygon[next]);
                if turn > max_turn {
                    max_turn = turn;
                    best_i = Some(i);
                }
            }

            let idx = best_i.unwrap_or(0);
            let prev = remaining[(idx + m - 1) % m];
            let curr = remaining[idx];
            let next = remaining[(idx + 1) % m];
            triangles.push([prev, curr, next]);
            remaining.remove(idx);
        }
    }

    if remaining.len() == 3 {
        triangles.push([remaining[0], remaining[1], remaining[2]]);
    }

    Ok(triangles)
}
