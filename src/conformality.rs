//! Restoring conformality and maintaining granularity (Section 3.2 & 3.3).

use crate::geometry::Point;
use crate::submap::{Submap, SubmapRegion};

/// Enforces conformality by ensuring no region in the submap exceeds degree 4 (Section 3.2).
pub fn restore_conformality(mut submap: Submap, _polygon: &[Point]) -> Submap {
    let initial_len = submap.regions.len();
    for r_idx in 0..initial_len {
        while submap.regions[r_idx].chords.len() > 4 || submap.regions[r_idx].arcs.len() > 4 {
            let chords = &mut submap.regions[r_idx].chords;
            if chords.len() <= 2 {
                break;
            }
            let mid = chords.len() / 2;
            let chords2 = chords.split_off(mid);

            let new_r_id = submap.regions.len();
            let half_weight = submap.regions[r_idx].weight / 2;
            submap.regions[r_idx].weight -= half_weight;

            let mut new_reg = SubmapRegion::new(new_r_id);
            new_reg.chords = chords2;
            new_reg.weight = half_weight;
            submap.regions.push(new_reg);
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

    let regions = &submap.regions;
    let n_regions = regions.len();

    submap.chords.retain(|chord| {
        let r1 = chord.region1;
        let r2 = chord.region2;
        if r1 < n_regions && r2 < n_regions {
            let combined_weight = regions[r1].weight + regions[r2].weight;
            let deg1 = regions[r1].chords.len();
            let deg2 = regions[r2].chords.len();

            // Retain chord if neither region is low-degree or combined weight > gamma
            !((deg1 < 3 || deg2 < 3) && combined_weight <= gamma)
        } else {
            true
        }
    });

    submap.granularity = gamma;
    submap
}
