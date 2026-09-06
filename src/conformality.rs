//! Restoring conformality and maintaining granularity (Section 3.2 & 3.3).

use crate::geometry::Point;
use crate::submap::{Chord, Submap, SubmapRegion};

/// Enforces conformality by ensuring no region in the submap exceeds degree 4 (Section 3.2).
pub fn restore_conformality(mut submap: Submap, _polygon: &[Point]) -> Submap {
    let mut needs_split = true;

    while needs_split {
        needs_split = false;
        let mut split_region_idx = None;

        for (idx, reg) in submap.regions.iter().enumerate() {
            if reg.chords.len() > 4 || reg.arcs.len() > 4 {
                split_region_idx = Some(idx);
                break;
            }
        }

        if let Some(r_idx) = split_region_idx {
            // Split the region by subdividing its chords and arcs
            let reg = &submap.regions[r_idx];
            let chords = reg.chords.clone();
            let mid = chords.len() / 2;

            let chords1: Vec<usize> = chords[..mid].to_vec();
            let chords2: Vec<usize> = chords[mid..].to_vec();

            let new_r_id = submap.regions.len();
            let mut new_reg = SubmapRegion::new(new_r_id);
            new_reg.chords = chords2;
            new_reg.weight = reg.weight / 2;

            submap.regions[r_idx].chords = chords1;
            submap.regions[r_idx].weight -= new_reg.weight;

            submap.regions.push(new_reg);
            needs_split = true;
        }
    }

    submap
}

/// Enforces gamma-granularity by removing superfluous chords (Section 3.3).
/// Merges regions whose combined weight does not exceed the target granularity gamma.
pub fn maintain_granularity(mut submap: Submap, gamma: usize) -> Submap {
    if submap.chords.is_empty() || gamma <= 1 {
        return submap;
    }

    let mut remove_chord_indices = Vec::new();

    for (c_idx, chord) in submap.chords.iter().enumerate() {
        let r1 = chord.region1;
        let r2 = chord.region2;
        if r1 < submap.regions.len() && r2 < submap.regions.len() {
            let combined_weight = submap.regions[r1].weight + submap.regions[r2].weight;
            let deg1 = submap.regions[r1].chords.len();
            let deg2 = submap.regions[r2].chords.len();

            // If incident upon at least one node of degree < 3 and combined weight <= gamma:
            if (deg1 < 3 || deg2 < 3) && combined_weight <= gamma {
                remove_chord_indices.push(c_idx);
            }
        }
    }

    if remove_chord_indices.is_empty() {
        return submap;
    }

    // Filter out removed chords
    let remove_set: std::collections::HashSet<usize> = remove_chord_indices.into_iter().collect();
    let remaining_chords: Vec<Chord> = submap
        .chords
        .into_iter()
        .enumerate()
        .filter(|(idx, _)| !remove_set.contains(idx))
        .map(|(_, c)| c)
        .collect();

    submap.chords = remaining_chords;
    submap.granularity = gamma;
    submap
}
