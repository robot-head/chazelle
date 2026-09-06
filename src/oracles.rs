//! Ray-shooting and arc-cutting oracles (Section 3.4).

use crate::geometry::{Point, ray_h_intersect_segment};
use crate::double_boundary::ChordDirection;
use crate::submap::Submap;

/// Ray-shooting oracle for a polygonal curve subchain.
pub struct RayShootingOracle<'a> {
    pub polygon: &'a [Point],
}

impl<'a> RayShootingOracle<'a> {
    pub fn new(polygon: &'a [Point]) -> Self {
        Self { polygon }
    }

    /// Shoot a horizontal ray from `origin` in `dir` against a subchain `[start_edge..=end_edge]`.
    pub fn shoot_ray(
        &self,
        origin: Point,
        dir: ChordDirection,
        start_edge: usize,
        num_edges: usize,
    ) -> Option<(f64, Point, usize)> {
        let n = self.polygon.len();
        let mut closest_dist = f64::INFINITY;
        let mut best_hit = None;

        for step in 0..num_edges {
            let e = (start_edge + step) % n;
            let next_e = (e + 1) % n;
            let p1 = self.polygon[e];
            let p2 = self.polygon[next_e];

            if let Some((dist, hit_pt)) = ray_h_intersect_segment(origin, dir.is_right(), p1, p2) {
                if dist < closest_dist {
                    closest_dist = dist;
                    best_hit = Some((dist, hit_pt, e));
                }
            }
        }

        best_hit
    }

    /// Ray shooting using a precomputed submap (Lemma 3.6).
    pub fn shoot_ray_with_submap(
        &self,
        submap: &Submap,
        region_id: usize,
        origin: Point,
        dir: ChordDirection,
    ) -> Option<(f64, Point, usize)> {
        submap.local_ray_shoot(region_id, origin, dir, self.polygon)
    }
}

/// Arc-cutting oracle: partitions an arc of the submap into canonical subchains.
pub struct ArcCuttingOracle;

impl ArcCuttingOracle {
    /// Cuts a subchain `[start_edge, start_edge + num_edges)` into chunks of size at most `max_chunk_size`.
    pub fn cut_arc(
        start_edge: usize,
        num_edges: usize,
        max_chunk_size: usize,
    ) -> Vec<(usize, usize)> {
        let mut chunks = Vec::new();
        let chunk_size = max_chunk_size.max(1);
        let mut rem = num_edges;
        let mut curr = start_edge;

        while rem > 0 {
            let take = rem.min(chunk_size);
            chunks.push((curr, take));
            curr += take;
            rem -= take;
        }

        chunks
    }
}
