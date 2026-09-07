//! Ray-shooting and arc-cutting oracles (Section 3.4).

use crate::geometry::{Point, ray_h_intersect_segment};
use crate::double_boundary::ChordDirection;
use crate::submap::Submap;

/// Ray-shooting oracle for a polygonal curve subchain with spatial slab binning.
pub struct RayShootingOracle<'a> {
    pub polygon: &'a [Point],
    min_y: f64,
    max_y: f64,
    inv_bin_h: f64,
    num_bins: usize,
    bin_offsets: Vec<usize>,
    bin_edges: Vec<usize>,
}

impl<'a> RayShootingOracle<'a> {
    pub fn new(polygon: &'a [Point]) -> Self {
        let n = polygon.len();
        if n == 0 {
            return Self {
                polygon,
                min_y: 0.0,
                max_y: 0.0,
                inv_bin_h: 1.0,
                num_bins: 0,
                bin_offsets: vec![0],
                bin_edges: Vec::new(),
            };
        }

        let mut min_y = polygon[0].y;
        let mut max_y = polygon[0].y;
        for p in polygon {
            if p.y < min_y { min_y = p.y; }
            if p.y > max_y { max_y = p.y; }
        }

        let height = max_y - min_y;
        let num_bins = if n < 16 || height <= Point::EPSILON {
            1
        } else {
            (n / 4).clamp(32, 2048)
        };

        let inv_bin_h = if height <= Point::EPSILON || num_bins <= 1 {
            1.0
        } else {
            (num_bins as f64) / height
        };

        let mut counts = vec![0usize; num_bins];

        // Pass 1: count edges per bin
        for e in 0..n {
            let p1 = polygon[e];
            let p2 = polygon[(e + 1) % n];
            let y_low = (p1.y.min(p2.y) - Point::EPSILON).max(min_y);
            let y_high = (p1.y.max(p2.y) + Point::EPSILON).min(max_y);

            let b_low = if num_bins <= 1 {
                0
            } else {
                (((y_low - min_y) * inv_bin_h).floor() as isize).clamp(0, num_bins as isize - 1) as usize
            };
            let b_high = if num_bins <= 1 {
                0
            } else {
                (((y_high - min_y) * inv_bin_h).floor() as isize).clamp(0, num_bins as isize - 1) as usize
            };

            for b in b_low..=b_high {
                counts[b] += 1;
            }
        }

        // Prefix sum for offsets
        let mut bin_offsets = Vec::with_capacity(num_bins + 1);
        bin_offsets.push(0);
        let mut total = 0;
        for &c in &counts {
            total += c;
            bin_offsets.push(total);
        }

        let mut bin_edges = vec![0usize; total];
        let mut cursor = bin_offsets.clone();

        // Pass 2: fill bin_edges
        for e in 0..n {
            let p1 = polygon[e];
            let p2 = polygon[(e + 1) % n];
            let y_low = (p1.y.min(p2.y) - Point::EPSILON).max(min_y);
            let y_high = (p1.y.max(p2.y) + Point::EPSILON).min(max_y);

            let b_low = if num_bins <= 1 {
                0
            } else {
                (((y_low - min_y) * inv_bin_h).floor() as isize).clamp(0, num_bins as isize - 1) as usize
            };
            let b_high = if num_bins <= 1 {
                0
            } else {
                (((y_high - min_y) * inv_bin_h).floor() as isize).clamp(0, num_bins as isize - 1) as usize
            };

            for b in b_low..=b_high {
                bin_edges[cursor[b]] = e;
                cursor[b] += 1;
            }
        }

        // Pass 3: sort edges in each bin by min_x for directional early-exit pruning
        for b in 0..num_bins {
            let start = bin_offsets[b];
            let end = bin_offsets[b + 1];
            if end > start + 1 {
                bin_edges[start..end].sort_unstable_by(|&e1, &e2| {
                    let min_x1 = polygon[e1].x.min(polygon[(e1 + 1) % n].x);
                    let min_x2 = polygon[e2].x.min(polygon[(e2 + 1) % n].x);
                    min_x1.partial_cmp(&min_x2).unwrap_or(std::cmp::Ordering::Equal)
                });
            }
        }

        Self {
            polygon,
            min_y,
            max_y,
            inv_bin_h,
            num_bins,
            bin_offsets,
            bin_edges,
        }
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
        if num_edges == 0 || n == 0 {
            return None;
        }

        // For small subchains, direct linear scan is fastest
        if num_edges <= 16 || self.num_bins <= 1 {
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

            return best_hit;
        }

        let y = origin.y;
        if y < self.min_y - Point::EPSILON || y > self.max_y + Point::EPSILON {
            return None;
        }

        let b = (((y - self.min_y) * self.inv_bin_h).floor() as isize)
            .clamp(0, self.num_bins as isize - 1) as usize;

        let mut closest_dist = f64::INFINITY;
        let mut best_hit = None;

        let b_start = if b > 0 && (y - (self.min_y + b as f64 / self.inv_bin_h)).abs() <= Point::EPSILON {
            b - 1
        } else {
            b
        };
        let b_end = if b + 1 < self.num_bins
            && ((self.min_y + (b + 1) as f64 / self.inv_bin_h) - y).abs() <= Point::EPSILON
        {
            b + 1
        } else {
            b
        };

        if dir.is_right() {
            for curr_b in b_start..=b_end {
                let start_idx = self.bin_offsets[curr_b];
                let end_idx = self.bin_offsets[curr_b + 1];

                for idx in start_idx..end_idx {
                    let e = self.bin_edges[idx];
                    let next_e = (e + 1) % n;
                    let p1 = self.polygon[e];
                    let p2 = self.polygon[next_e];

                    let min_x = p1.x.min(p2.x);
                    if min_x > origin.x + closest_dist {
                        break; // Sorted by min_x: all subsequent edges are farther than closest hit
                    }

                    let max_x = p1.x.max(p2.x);
                    if max_x < origin.x - Point::EPSILON {
                        continue; // Completely behind ray
                    }

                    if num_edges < n {
                        let diff = (e + n - start_edge) % n;
                        if diff >= num_edges {
                            continue;
                        }
                    }

                    if let Some((dist, hit_pt)) = ray_h_intersect_segment(origin, true, p1, p2) {
                        if dist < closest_dist {
                            closest_dist = dist;
                            best_hit = Some((dist, hit_pt, e));
                        }
                    }
                }
            }
        } else {
            // Left ray: scan from right to left (reverse order)
            for curr_b in b_start..=b_end {
                let start_idx = self.bin_offsets[curr_b];
                let end_idx = self.bin_offsets[curr_b + 1];

                for idx in (start_idx..end_idx).rev() {
                    let e = self.bin_edges[idx];
                    let next_e = (e + 1) % n;
                    let p1 = self.polygon[e];
                    let p2 = self.polygon[next_e];

                    let max_x = p1.x.max(p2.x);
                    if max_x < origin.x - closest_dist {
                        break; // All remaining edges are too far to the left
                    }

                    let min_x = p1.x.min(p2.x);
                    if min_x > origin.x + Point::EPSILON {
                        continue; // Completely behind ray
                    }

                    if num_edges < n {
                        let diff = (e + n - start_edge) % n;
                        if diff >= num_edges {
                            continue;
                        }
                    }

                    if let Some((dist, hit_pt)) = ray_h_intersect_segment(origin, false, p1, p2) {
                        if dist < closest_dist {
                            closest_dist = dist;
                            best_hit = Some((dist, hit_pt, e));
                        }
                    }
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
