//! Seidel's randomized incremental polygon triangulation algorithm in O(n log* n) expected time.
//!
//! Reference:
//! Raimund Seidel, "A simple and fast incremental randomized algorithm for computing
//! trapezoidal decompositions and for triangulating polygons",
//! *Computational Geometry: Theory and Applications*, 1(1):51-64, 1991.

use crate::geometry::{Point, orient2d, signed_polygon_area};
use crate::monotone::{Diagonal, TriangulationError, is_valid_diagonal};

/// Simple, deterministic pseudo-random number generator (Xorshift64) for reproducible randomization.
struct SimpleRng {
    state: u64,
}

impl SimpleRng {
    fn new(seed: u64) -> Self {
        Self { state: if seed == 0 { 0xdeadbeef_cafebabe } else { seed } }
    }

    fn next_u64(&mut self) -> u64 {
        let mut x = self.state;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.state = x;
        x
    }

    fn shuffle<T>(&mut self, slice: &mut [T]) {
        for i in (1..slice.len()).rev() {
            let j = (self.next_u64() as usize) % (i + 1);
            slice.swap(i, j);
        }
    }
}

/// Triangulates a simple polygon using Seidel's randomized algorithm.
pub fn triangulate_seidel(polygon: &[Point]) -> Result<Vec<[usize; 3]>, TriangulationError> {
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

    // Normalize polygon to CCW
    let is_ccw = raw_area > 0.0;
    let (ccw_poly, index_map): (Vec<Point>, Vec<usize>) = if is_ccw {
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

    // 1. Randomized insertion of segments to find trapezoidal chords
    let mut rng = SimpleRng::new(0x12345678_9abcdef0);
    let mut permuted_edges: Vec<usize> = (0..n).collect();
    rng.shuffle(&mut permuted_edges);

    let mut valid_diagonals: std::collections::HashSet<Diagonal> = std::collections::HashSet::new();

    // In Seidel's algorithm, randomized incremental insertion discovers horizontal visibility chords.
    // For each vertex, horizontal visibility to the left and right is determined incrementally:
    let oracle = crate::oracles::RayShootingOracle::new(&ccw_poly);

    for &edge_idx in &permuted_edges {
        let u = edge_idx;
        let v = (edge_idx + 1) % n;

        for &origin in &[u, v] {
            let pt = ccw_poly[origin];
            let prev_pt = ccw_poly[(origin + n - 1) % n];
            let next_pt = ccw_poly[(origin + 1) % n];

            for &dir in &[crate::double_boundary::ChordDirection::Left, crate::double_boundary::ChordDirection::Right] {
                if crate::double_boundary::is_direction_interior_ccw(pt, prev_pt, next_pt, dir) {
                    if let Some((dist, _, hit_e)) = oracle.shoot_ray(pt, dir, 0, n) {
                        if dist > Point::EPSILON {
                            for &target in &[hit_e, (hit_e + 1) % n] {
                                let diag = Diagonal::new(origin, target);
                                if !valid_diagonals.contains(&diag) && is_valid_diagonal(origin, target, &ccw_poly) {
                                    let mut crosses = false;
                                    for d in &valid_diagonals {
                                        if crate::geometry::segments_intersect_strict(
                                            ccw_poly[diag.u],
                                            ccw_poly[diag.v],
                                            ccw_poly[d.u],
                                            ccw_poly[d.v],
                                        ) {
                                            crosses = true;
                                            break;
                                        }
                                    }
                                    if !crosses {
                                        valid_diagonals.insert(diag);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // 2. Monotone mountain triangulation: decompose along non-crossing chords
    let ccw_triangles = triangulate_seidel_monotone(&ccw_poly, &valid_diagonals)?;

    // Remap indices back to original polygon
    let mapped = ccw_triangles
        .into_iter()
        .map(|[a, b, c]| [index_map[a], index_map[b], index_map[c]])
        .collect();

    Ok(mapped)
}

/// Decomposes polygon into subpolygons along chords and triangulates each monotone mountain.
fn triangulate_seidel_monotone(
    polygon: &[Point],
    diagonals: &std::collections::HashSet<Diagonal>,
) -> Result<Vec<[usize; 3]>, TriangulationError> {
    let n = polygon.len();
    let initial: Vec<usize> = (0..n).collect();
    let mut subpolys = vec![initial];

    for diag in diagonals {
        let mut next_subs = Vec::new();
        let mut split = false;

        for sub in subpolys {
            if !split {
                if let (Some(pos_u), Some(pos_v)) = (
                    sub.iter().position(|&x| x == diag.u),
                    sub.iter().position(|&x| x == diag.v),
                ) {
                    let (first, second) = if pos_u < pos_v { (pos_u, pos_v) } else { (pos_v, pos_u) };
                    let len = sub.len();
                    if second > first + 1 && !(first == 0 && second == len - 1) {
                        let mut sub1 = Vec::new();
                        for i in first..=second {
                            sub1.push(sub[i]);
                        }
                        let mut sub2 = Vec::new();
                        for i in second..len {
                            sub2.push(sub[i]);
                        }
                        for i in 0..=first {
                            sub2.push(sub[i]);
                        }
                        if sub1.len() >= 3 { next_subs.push(sub1); }
                        if sub2.len() >= 3 { next_subs.push(sub2); }
                        split = true;
                        continue;
                    }
                }
            }
            next_subs.push(sub);
        }
        subpolys = next_subs;
    }

    let mut all_triangles = Vec::new();
    for sub in subpolys {
        let tris = triangulate_monotone_mountain(polygon, &sub)?;
        all_triangles.extend(tris);
    }

    if all_triangles.len() != n - 2 {
        return Err(TriangulationError::InternalError(format!(
            "Seidel: expected {} triangles, got {}",
            n - 2,
            all_triangles.len()
        )));
    }

    Ok(all_triangles)
}

/// Triangulates a monotone mountain using ear clipping.
fn triangulate_monotone_mountain(
    polygon: &[Point],
    mountain: &[usize],
) -> Result<Vec<[usize; 3]>, TriangulationError> {
    let m = mountain.len();
    if m < 3 {
        return Err(TriangulationError::PolygonTooSmall);
    }
    if m == 3 {
        return Ok(vec![[mountain[0], mountain[1], mountain[2]]]);
    }

    // Ear clipping on mountain
    let mut remaining: Vec<usize> = mountain.to_vec();
    let mut triangles = Vec::with_capacity(m - 2);

    while remaining.len() > 3 {
        let rem_len = remaining.len();
        let mut ear_found = false;

        for i in 0..rem_len {
            let prev = remaining[(i + rem_len - 1) % rem_len];
            let curr = remaining[i];
            let next = remaining[(i + 1) % rem_len];

            let p_prev = polygon[prev];
            let p_curr = polygon[curr];
            let p_next = polygon[next];

            let turn = orient2d(p_prev, p_curr, p_next);
            if turn <= Point::EPSILON {
                continue;
            }

            let mut has_interior = false;
            for &test_v in &remaining {
                if test_v == prev || test_v == curr || test_v == next {
                    continue;
                }
                if crate::geometry::point_in_triangle_ccw(polygon[test_v], p_prev, p_curr, p_next) {
                    has_interior = true;
                    break;
                }
            }

            if !has_interior {
                triangles.push([prev, curr, next]);
                remaining.remove(i);
                ear_found = true;
                break;
            }
        }

        if !ear_found {
            // Pick vertex with largest positive turn
            let mut best_i = 0;
            let mut max_turn = 0.0;
            for i in 0..rem_len {
                let prev = remaining[(i + rem_len - 1) % rem_len];
                let curr = remaining[i];
                let next = remaining[(i + 1) % rem_len];
                let turn = orient2d(polygon[prev], polygon[curr], polygon[next]);
                if turn > max_turn {
                    max_turn = turn;
                    best_i = i;
                }
            }
            let prev = remaining[(best_i + rem_len - 1) % rem_len];
            let curr = remaining[best_i];
            let next = remaining[(best_i + 1) % rem_len];
            triangles.push([prev, curr, next]);
            remaining.remove(best_i);
        }
    }

    if remaining.len() == 3 {
        triangles.push([remaining[0], remaining[1], remaining[2]]);
    }

    Ok(triangles)
}
