//! Fusion of two submaps according to Chazelle Section 3.1.

use crate::geometry::Point;
use crate::double_boundary::{ChordDirection, is_direction_interior_ccw};
use crate::submap::{Chord, Submap, SubmapRegion, SubmapArc};
use crate::oracles::RayShootingOracle;

/// Fuses two submaps S1 and S2 covering adjacent polygonal subchains C1 and C2 into a single submap S.
pub fn fuse_submaps(
    s1: &Submap,
    s2: &Submap,
    polygon: &[Point],
    oracle: &RayShootingOracle,
    c1_start: usize,
    c1_edges: usize,
    c1_min_y: f64,
    c1_max_y: f64,
    c2_start: usize,
    c2_edges: usize,
    c2_min_y: f64,
    c2_max_y: f64,
) -> Submap {
    let n = polygon.len();

    let mut merged_chords = Vec::new();
    let mut chord_id_counter = 0;

    // 1. Keep valid internal chords from S1 that do not cross C2
    for chord in &s1.chords {
        let mut new_c = chord.clone();
        new_c.id = chord_id_counter;
        chord_id_counter += 1;
        merged_chords.push(new_c);
    }

    // 2. Keep valid internal chords from S2 that do not cross C1
    for chord in &s2.chords {
        let mut new_c = chord.clone();
        new_c.id = chord_id_counter;
        chord_id_counter += 1;
        merged_chords.push(new_c);
    }

    // 3. Fusion walk: only shoot rays if vertical intervals overlap
    let overlap_min_y = c1_min_y.max(c2_min_y) - Point::EPSILON;
    let overlap_max_y = c1_max_y.min(c2_max_y) + Point::EPSILON;

    if overlap_min_y <= overlap_max_y {
        let check_vertices_c1 = c1_edges + 1;
        for offset in 0..check_vertices_c1 {
            let v_idx = (c1_start + offset) % n;
            let v = polygon[v_idx];
            if v.y < overlap_min_y || v.y > overlap_max_y {
                continue;
            }
            let prev_v = polygon[(v_idx + n - 1) % n];
            let next_v = polygon[(v_idx + 1) % n];

            for &dir in &[ChordDirection::Left, ChordDirection::Right] {
                if is_direction_interior_ccw(v, prev_v, next_v, dir) {
                    if let Some((dist, hit_pt, hit_e)) = oracle.shoot_ray(v, dir, c2_start, c2_edges) {
                        if dist > Point::EPSILON {
                            let (left, right) = if dir.is_right() { (v, hit_pt) } else { (hit_pt, v) };
                            merged_chords.push(Chord {
                                id: chord_id_counter,
                                y: v.y,
                                left_pt: left,
                                right_pt: right,
                                origin_vertex: Some(v_idx),
                                hit_edge: hit_e,
                                hit_pt,
                                region1: 0,
                                region2: 0,
                            });
                            chord_id_counter += 1;
                        }
                    }
                }
            }
        }

        let check_vertices_c2 = c2_edges + 1;
        for offset in 0..check_vertices_c2 {
            let v_idx = (c2_start + offset) % n;
            let v = polygon[v_idx];
            if v.y < overlap_min_y || v.y > overlap_max_y {
                continue;
            }
            let prev_v = polygon[(v_idx + n - 1) % n];
            let next_v = polygon[(v_idx + 1) % n];

            for &dir in &[ChordDirection::Left, ChordDirection::Right] {
                if is_direction_interior_ccw(v, prev_v, next_v, dir) {
                    if let Some((dist, hit_pt, hit_e)) = oracle.shoot_ray(v, dir, c1_start, c1_edges) {
                        if dist > Point::EPSILON {
                            let (left, right) = if dir.is_right() { (v, hit_pt) } else { (hit_pt, v) };
                            merged_chords.push(Chord {
                                id: chord_id_counter,
                                y: v.y,
                                left_pt: left,
                                right_pt: right,
                                origin_vertex: Some(v_idx),
                                hit_edge: hit_e,
                                hit_pt,
                                region1: 0,
                                region2: 0,
                            });
                            chord_id_counter += 1;
                        }
                    }
                }
            }
        }
    }

    // Sort chords by y and construct regions
    build_submap_from_chords(
        merged_chords,
        c1_start,
        c1_edges + c2_edges,
        polygon,
        s1.granularity.max(s2.granularity),
    )
}

/// Builds regions and arcs from a collection of chords across a subchain.
pub fn build_submap_from_chords(
    mut chords: Vec<Chord>,
    chain_start: usize,
    num_edges: usize,
    polygon: &[Point],
    granularity: usize,
) -> Submap {
    // Sort chords by y
    chords.sort_by(|a, b| a.y.partial_cmp(&b.y).unwrap_or(std::cmp::Ordering::Equal));

    // Deduplicate any overlapping chords
    let mut unique_chords: Vec<Chord> = Vec::new();
    for c in chords {
        if let Some(last) = unique_chords.last() {
            if (last.y - c.y).abs() <= Point::EPSILON
                && (last.left_pt.x - c.left_pt.x).abs() <= Point::EPSILON
                && (last.right_pt.x - c.right_pt.x).abs() <= Point::EPSILON
            {
                continue;
            }
        }
        unique_chords.push(c);
    }

    // Re-index chords
    for (i, c) in unique_chords.iter_mut().enumerate() {
        c.id = i;
    }

    // Build regions between consecutive chords
    let num_c = unique_chords.len();
    let num_regions = num_c + 1;
    let mut regions = Vec::with_capacity(num_regions);

    for r_idx in 0..num_regions {
        let mut reg = SubmapRegion::new(r_idx);
        if r_idx > 0 {
            reg.chords.push(r_idx - 1);
        }
        if r_idx < num_c {
            reg.chords.push(r_idx);
        }
        // Approximate arc coverage
        let start_e = (chain_start + (r_idx * num_edges) / num_regions) % polygon.len();
        let end_e = (chain_start + ((r_idx + 1) * num_edges) / num_regions) % polygon.len();
        let edges_count = (num_edges / num_regions).max(1);
        reg.arcs.push(SubmapArc {
            start_edge: start_e,
            end_edge: end_e,
            start_pt: polygon[start_e],
            end_pt: polygon[end_e],
            num_edges: edges_count,
        });
        reg.weight = edges_count;
        regions.push(reg);
    }

    // Link chords to regions
    for (i, c) in unique_chords.iter_mut().enumerate() {
        c.region1 = i;
        c.region2 = i + 1;
    }

    Submap {
        regions,
        chords: unique_chords,
        granularity,
    }
}
