//! Classic textbook O(n log n) polygon triangulation via plane-sweep monotone decomposition.
//!
//! Reference:
//! Mark de Berg, Otfried Cheong, Marc van Kreveld, Mark Overmars,
//! "Computational Geometry: Algorithms and Applications", Chapter 3 (Polygon Triangulation).

use crate::geometry::{Point, orient2d, signed_polygon_area};
use crate::monotone::{Diagonal, TriangulationError, is_valid_diagonal};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum VertexType {
    Start,
    End,
    Split,
    Merge,
    RegularLeft,
    RegularRight,
}

#[inline]
fn is_above(p1: Point, idx1: usize, p2: Point, idx2: usize) -> bool {
    if (p1.y - p2.y).abs() > Point::EPSILON {
        p1.y > p2.y
    } else if (p1.x - p2.x).abs() > Point::EPSILON {
        p1.x < p2.x
    } else {
        idx1 < idx2
    }
}

/// Classifies a vertex of a CCW simple polygon into Start, End, Split, Merge, or Regular.
fn classify_vertex(i: usize, polygon: &[Point]) -> VertexType {
    let n = polygon.len();
    let prev_idx = (i + n - 1) % n;
    let next_idx = (i + 1) % n;

    let v = polygon[i];
    let prev = polygon[prev_idx];
    let next = polygon[next_idx];

    let prev_is_below = is_above(v, i, prev, prev_idx);
    let next_is_below = is_above(v, i, next, next_idx);

    let turn = orient2d(prev, v, next);

    if prev_is_below && next_is_below {
        // Both neighbors are below
        if turn > Point::EPSILON {
            VertexType::Start // Interior angle < 180 deg
        } else {
            VertexType::Split // Reflex (interior angle > 180 deg)
        }
    } else if !prev_is_below && !next_is_below {
        // Both neighbors are above
        if turn > Point::EPSILON {
            VertexType::End
        } else {
            VertexType::Merge // Reflex
        }
    } else {
        // One above, one below
        if prev_is_below {
            // Edge from prev to v goes up, interior is to the left
            VertexType::RegularLeft
        } else {
            // Edge from prev to v goes down, interior is to the right
            VertexType::RegularRight
        }
    }
}

/// Triangulates a simple polygon using the classic O(n log n) plane sweep monotone decomposition.
pub fn triangulate_monotone_sweep(polygon: &[Point]) -> Result<Vec<[usize; 3]>, TriangulationError> {
    let n = polygon.len();
    if n < 3 {
        return Err(TriangulationError::PolygonTooSmall);
    }
    if n == 3 {
        return Ok(vec![[0, 1, 2]]);
    }

    let area = signed_polygon_area(polygon);
    if area.abs() <= Point::EPSILON {
        return Err(TriangulationError::DegeneratePolygon);
    }

    // Ensure CCW polygon
    let is_ccw = area > 0.0;
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

    // 1. Classify all vertices
    let mut vertex_types = Vec::with_capacity(n);
    for i in 0..n {
        vertex_types.push(classify_vertex(i, &ccw_poly));
    }

    // 2. Sort vertices for plane sweep (decreasing y, tie-break by increasing x)
    let mut event_order: Vec<usize> = (0..n).collect();
    event_order.sort_by(|&a, &b| {
        if is_above(ccw_poly[a], a, ccw_poly[b], b) {
            std::cmp::Ordering::Less
        } else {
            std::cmp::Ordering::Greater
        }
    });

    let mut diagonals: std::collections::HashSet<Diagonal> = std::collections::HashSet::new();

    // Active edges in sweep line: represented by edge index e (segment ccw_poly[e] -> ccw_poly[(e+1)%n])
    // along with the current helper vertex.
    #[derive(Clone)]
    struct ActiveEdge {
        edge_idx: usize,
        helper: usize,
    }

    let mut active_edges: Vec<ActiveEdge> = Vec::new();

    // Helper: compute x-coordinate of edge at current sweep-line y
    let edge_x_at_y = |edge_idx: usize, y: f64| -> f64 {
        let p1 = ccw_poly[edge_idx];
        let p2 = ccw_poly[(edge_idx + 1) % n];
        if (p1.y - p2.y).abs() <= Point::EPSILON {
            p1.x.min(p2.x)
        } else {
            let t = (y - p1.y) / (p2.y - p1.y);
            p1.x + t * (p2.x - p1.x)
        }
    };

    // Find active edge immediately to the left of point v
    let find_left_edge = |active: &[ActiveEdge], v: Point| -> Option<usize> {
        let mut best_idx = None;
        let mut max_x = f64::NEG_INFINITY;

        for (idx, ae) in active.iter().enumerate() {
            let x = edge_x_at_y(ae.edge_idx, v.y);
            if x <= v.x + Point::EPSILON && x > max_x {
                max_x = x;
                best_idx = Some(idx);
            }
        }
        best_idx
    };

    let add_diagonal = |u: usize, v: usize, diags: &mut std::collections::HashSet<Diagonal>, poly: &[Point]| {
        if u != v && (u + 1) % n != v && (v + 1) % n != u {
            let d = Diagonal::new(u, v);
            if !diags.contains(&d) && is_valid_diagonal(u, v, poly) {
                diags.insert(d);
            }
        }
    };

    // Process vertices in sweep order
    for &v_idx in &event_order {
        let v = ccw_poly[v_idx];
        let v_type = vertex_types[v_idx];
        let prev_e = (v_idx + n - 1) % n;

        match v_type {
            VertexType::Start => {
                // Insert edge e_i (v_idx -> v_idx + 1)
                active_edges.push(ActiveEdge {
                    edge_idx: v_idx,
                    helper: v_idx,
                });
            }
            VertexType::End => {
                // If helper(e_{i-1}) is a merge vertex, add diagonal
                if let Some(pos) = active_edges.iter().position(|e| e.edge_idx == prev_e) {
                    let helper = active_edges[pos].helper;
                    if vertex_types[helper] == VertexType::Merge {
                        add_diagonal(v_idx, helper, &mut diagonals, &ccw_poly);
                    }
                    active_edges.remove(pos);
                }
            }
            VertexType::Split => {
                // Find edge immediately to the left of v
                if let Some(left_idx) = find_left_edge(&active_edges, v) {
                    let helper = active_edges[left_idx].helper;
                    add_diagonal(v_idx, helper, &mut diagonals, &ccw_poly);
                    active_edges[left_idx].helper = v_idx;
                }
                active_edges.push(ActiveEdge {
                    edge_idx: v_idx,
                    helper: v_idx,
                });
            }
            VertexType::Merge => {
                if let Some(pos) = active_edges.iter().position(|e| e.edge_idx == prev_e) {
                    let helper = active_edges[pos].helper;
                    if vertex_types[helper] == VertexType::Merge {
                        add_diagonal(v_idx, helper, &mut diagonals, &ccw_poly);
                    }
                    active_edges.remove(pos);
                }
                if let Some(left_idx) = find_left_edge(&active_edges, v) {
                    let helper = active_edges[left_idx].helper;
                    if vertex_types[helper] == VertexType::Merge {
                        add_diagonal(v_idx, helper, &mut diagonals, &ccw_poly);
                    }
                    active_edges[left_idx].helper = v_idx;
                }
            }
            VertexType::RegularRight => {
                // Interior to the right: edge prev_e is replaced by v_idx
                if let Some(pos) = active_edges.iter().position(|e| e.edge_idx == prev_e) {
                    let helper = active_edges[pos].helper;
                    if vertex_types[helper] == VertexType::Merge {
                        add_diagonal(v_idx, helper, &mut diagonals, &ccw_poly);
                    }
                    active_edges[pos] = ActiveEdge {
                        edge_idx: v_idx,
                        helper: v_idx,
                    };
                } else {
                    active_edges.push(ActiveEdge {
                        edge_idx: v_idx,
                        helper: v_idx,
                    });
                }
            }
            VertexType::RegularLeft => {
                // Interior to the left: check edge to left
                if let Some(left_idx) = find_left_edge(&active_edges, v) {
                    let helper = active_edges[left_idx].helper;
                    if vertex_types[helper] == VertexType::Merge {
                        add_diagonal(v_idx, helper, &mut diagonals, &ccw_poly);
                    }
                    active_edges[left_idx].helper = v_idx;
                }
            }
        }
    }

    // 3. Decompose polygon along diagonals into monotone subpolygons and triangulate each
    let ccw_triangles = triangulate_subpolygons_with_diagonals(&ccw_poly, &diagonals)?;

    // Remap indices back to original polygon
    let result = ccw_triangles
        .into_iter()
        .map(|[a, b, c]| [index_map[a], index_map[b], index_map[c]])
        .collect();

    Ok(result)
}

/// Decomposes polygon into subpolygons using diagonals and triangulates each $y$-monotone piece.
fn triangulate_subpolygons_with_diagonals(
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
        let tris = triangulate_monotone_piece(polygon, &sub)?;
        all_triangles.extend(tris);
    }

    if all_triangles.len() != n - 2 {
        return Err(TriangulationError::InternalError(format!(
            "Monotone sweep: expected {} triangles, got {}",
            n - 2,
            all_triangles.len()
        )));
    }

    Ok(all_triangles)
}

/// Triangulates a y-monotone piece using the classic linear stack algorithm (Garey et al. 1978).
pub fn triangulate_monotone_piece(
    polygon: &[Point],
    piece: &[usize],
) -> Result<Vec<[usize; 3]>, TriangulationError> {
    let m = piece.len();
    if m < 3 {
        return Err(TriangulationError::PolygonTooSmall);
    }
    if m == 3 {
        return Ok(vec![[piece[0], piece[1], piece[2]]]);
    }

    // Find top and bottom vertices in y
    let mut top_pos = 0;
    let mut bot_pos = 0;
    for i in 1..m {
        if is_above(polygon[piece[i]], piece[i], polygon[piece[top_pos]], piece[top_pos]) {
            top_pos = i;
        }
        if is_above(polygon[piece[bot_pos]], piece[bot_pos], polygon[piece[i]], piece[i]) {
            bot_pos = i;
        }
    }

    // Classify vertices as left chain or right chain
    // Traversal from top_pos to bot_pos in cycle order is one chain; remaining is other chain
    let mut is_left_chain = vec![false; m];
    let mut curr = top_pos;
    while curr != bot_pos {
        is_left_chain[curr] = true;
        curr = (curr + 1) % m;
    }
    is_left_chain[bot_pos] = true;

    // Sort vertices of this piece by decreasing y (merge of two sorted chains)
    let mut sorted_indices: Vec<usize> = (0..m).collect();
    sorted_indices.sort_by(|&a, &b| {
        if is_above(polygon[piece[a]], piece[a], polygon[piece[b]], piece[b]) {
            std::cmp::Ordering::Less
        } else {
            std::cmp::Ordering::Greater
        }
    });

    let mut stack: Vec<usize> = Vec::new();
    let mut triangles = Vec::new();

    stack.push(sorted_indices[0]);
    stack.push(sorted_indices[1]);

    for j in 2..m {
        let v_local = sorted_indices[j];
        let v_orig = piece[v_local];
        let top_local = *stack.last().unwrap();

        if is_left_chain[v_local] != is_left_chain[top_local] {
            // Different chain: pop all and connect diagonals
            let last = *stack.last().unwrap();
            while stack.len() > 1 {
                let u = stack.pop().unwrap();
                let next_u = *stack.last().unwrap();
                triangles.push([piece[u], piece[next_u], v_orig]);
            }
            stack.clear();
            stack.push(last);
            stack.push(v_local);
        } else {
            // Same chain: pop while internal turn
            let mut last_popped = stack.pop().unwrap();
            while let Some(&next_top) = stack.last() {
                let p1 = polygon[piece[next_top]];
                let p2 = polygon[piece[last_popped]];
                let p3 = polygon[v_orig];

                let turn = orient2d(p1, p2, p3);
                let valid_turn = if is_left_chain[v_local] {
                    turn > Point::EPSILON // Left chain requires left turn
                } else {
                    turn < -Point::EPSILON // Right chain requires right turn
                };

                if valid_turn {
                    triangles.push([piece[next_top], piece[last_popped], v_orig]);
                    last_popped = stack.pop().unwrap();
                } else {
                    break;
                }
            }
            stack.push(last_popped);
            stack.push(v_local);
        }
    }

    if triangles.len() == m - 2 {
        Ok(triangles)
    } else {
        Err(TriangulationError::InternalError("Piece not strictly monotone".into()))
    }
}
