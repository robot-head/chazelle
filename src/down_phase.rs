//! Top-down refinement phase (Section 4.2): refines submaps to full visibility map V(P).

use crate::geometry::Point;
use crate::double_boundary::{ChordDirection, is_direction_interior_ccw};
use crate::submap::{Chord, Submap};
use crate::up_phase::UpPhaseHierarchy;
use crate::oracles::RayShootingOracle;
use crate::fusion::build_submap_from_chords;

/// Executes the down-phase, refining the coarse submaps from `hierarchy` into the full visibility map V(P).
pub fn run_down_phase(hierarchy: &UpPhaseHierarchy, polygon: &[Point]) -> Submap {
    let n = polygon.len();
    let oracle = RayShootingOracle::new(polygon);

    // Collect all chords from the top submap
    let mut all_chords = Vec::new();
    if let Some(top) = hierarchy.top_submap() {
        all_chords.extend(top.chords.clone());
    }

    // Identify which vertices already have their interior horizontal chords
    let mut vertex_has_chord = vec![false; n];
    for c in &all_chords {
        if let Some(v) = c.origin_vertex {
            if v < n {
                vertex_has_chord[v] = true;
            }
        }
    }

    // For any remaining vertices, shoot horizontal rays to complete V(P)
    let mut chord_id_counter = all_chords.len();
    for v_idx in 0..n {
        let v = polygon[v_idx];
        let prev_v = polygon[(v_idx + n - 1) % n];
        let next_v = polygon[(v_idx + 1) % n];

        for &dir in &[ChordDirection::Left, ChordDirection::Right] {
            if is_direction_interior_ccw(v, prev_v, next_v, dir) {
                // Shoot ray across polygon edges
                if let Some((dist, hit_pt, hit_e)) = oracle.shoot_ray(v, dir, 0, n) {
                    if dist > Point::EPSILON {
                        let (left, right) = if dir.is_right() { (v, hit_pt) } else { (hit_pt, v) };
                        all_chords.push(Chord {
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

    // Build the fully refined granularity-1 visibility map
    build_submap_from_chords(all_chords, 0, n, polygon, 1)
}
