//! Normal-form representation of visibility submaps and centroid tree decomposition (Section 2.2-2.4).

use crate::geometry::{Point, ray_h_intersect_segment};
use crate::double_boundary::ChordDirection;

/// Represents a horizontal chord in a visibility submap.
#[derive(Debug, Clone, PartialEq)]
pub struct Chord {
    pub id: usize,
    pub y: f64,
    pub left_pt: Point,
    pub right_pt: Point,
    /// Origin vertex index in the input polygon, if chord emanates from a vertex.
    pub origin_vertex: Option<usize>,
    /// Polygon edge (i, (i+1)%n) that was hit by the horizontal ray.
    pub hit_edge: usize,
    pub hit_pt: Point,
    /// Adjacent regions in the submap tree.
    pub region1: usize,
    pub region2: usize,
}

impl Chord {
    pub fn length(&self) -> f64 {
        (self.right_pt.x - self.left_pt.x).abs()
    }
}

/// An arc is a contiguous sequence of polygon boundary edges between two chord endpoints.
#[derive(Debug, Clone, PartialEq)]
pub struct SubmapArc {
    pub start_edge: usize,
    pub end_edge: usize,
    pub start_pt: Point,
    pub end_pt: Point,
    pub num_edges: usize,
}

/// A region in a visibility submap (dual to a node in the submap tree).
/// In a conformal submap, each region has at most 4 chords and 4 arcs.
#[derive(Debug, Clone)]
pub struct SubmapRegion {
    pub id: usize,
    pub arcs: Vec<SubmapArc>,
    pub chords: Vec<usize>,
    pub weight: usize,
}

impl SubmapRegion {
    pub fn new(id: usize) -> Self {
        Self {
            id,
            arcs: Vec::new(),
            chords: Vec::new(),
            weight: 0,
        }
    }

    pub fn is_conformal(&self) -> bool {
        self.chords.len() <= 4 && self.arcs.len() <= 4
    }
}

/// Centroid tree decomposition of a conformal submap tree (Section 2.3, Fig 2.7).
/// Allows O(log m) navigation and binary search over the submap.
#[derive(Debug, Clone)]
pub enum CentroidTreeNode {
    Internal {
        chord_id: usize,
        left_child: Box<CentroidTreeNode>,
        right_child: Box<CentroidTreeNode>,
    },
    Leaf {
        region_id: usize,
    },
}

/// A visibility submap given in normal form (Section 2.4).
#[derive(Debug, Clone)]
pub struct Submap {
    pub regions: Vec<SubmapRegion>,
    pub chords: Vec<Chord>,
    pub granularity: usize,
}

impl Submap {
    pub fn new(granularity: usize) -> Self {
        Self {
            regions: Vec::new(),
            chords: Vec::new(),
            granularity,
        }
    }

    pub fn num_regions(&self) -> usize {
        self.regions.len()
    }

    pub fn num_chords(&self) -> usize {
        self.chords.len()
    }

    /// Performs local ray-shooting within a region against its bounding arcs.
    /// Since the submap is conformal, at most 4 arcs (and their edges) are checked.
    pub fn local_ray_shoot(
        &self,
        region_id: usize,
        origin: Point,
        dir: ChordDirection,
        polygon: &[Point],
    ) -> Option<(f64, Point, usize)> {
        if region_id >= self.regions.len() {
            return None;
        }
        let region = &self.regions[region_id];
        let n = polygon.len();

        let mut closest_dist = f64::INFINITY;
        let mut best_hit = None;

        for arc in &region.arcs {
            let edge_idx = arc.start_edge;
            let count = arc.num_edges;
            for step in 0..count {
                let e = (edge_idx + step) % n;
                let next_e = (e + 1) % n;
                let p1 = polygon[e];
                let p2 = polygon[next_e];

                if let Some((dist, hit_pt)) = ray_h_intersect_segment(origin, dir.is_right(), p1, p2) {
                    if dist < closest_dist {
                        closest_dist = dist;
                        best_hit = Some((dist, hit_pt, e));
                    }
                }
            }
        }

        best_hit
    }

    /// Builds a hierarchical centroid tree decomposition of the submap (Section 2.3).
    pub fn build_centroid_tree(&self) -> Option<CentroidTreeNode> {
        if self.regions.is_empty() {
            return None;
        }
        if self.regions.len() == 1 {
            return Some(CentroidTreeNode::Leaf { region_id: 0 });
        }

        // Tree of regions: adjacency via chords
        let num_r = self.regions.len();
        let mut adj: Vec<Vec<(usize, usize)>> = vec![Vec::new(); num_r];
        for chord in &self.chords {
            if chord.region1 < num_r && chord.region2 < num_r {
                adj[chord.region1].push((chord.region2, chord.id));
                adj[chord.region2].push((chord.region1, chord.id));
            }
        }

        let active_regions: Vec<usize> = (0..num_r).collect();
        Some(Self::decompose_subtree(&adj, &active_regions, &self.chords))
    }

    fn decompose_subtree(
        adj: &[Vec<(usize, usize)>],
        active: &[usize],
        chords: &[Chord],
    ) -> CentroidTreeNode {
        if active.len() <= 1 {
            let r = active.first().copied().unwrap_or(0);
            return CentroidTreeNode::Leaf { region_id: r };
        }

        // Subtree size from active nodes
        let active_set: std::collections::HashSet<usize> = active.iter().copied().collect();
        let total_size = active.len();

        // Find centroid edge whose removal splits the active subtree into components <= 3/4 total_size
        let mut best_edge = None;
        let mut best_balance = usize::MAX;

        // BFS / DFS to find edge sizes
        for &u in active {
            for &(v, chord_id) in &adj[u] {
                if u < v && active_set.contains(&v) {
                    // Count size of component containing v if edge (u, v) is cut
                    let comp_size = Self::count_component(v, u, adj, &active_set);
                    let other_size = total_size - comp_size;
                    let diff = comp_size.abs_diff(other_size);
                    if diff < best_balance {
                        best_balance = diff;
                        best_edge = Some((u, v, chord_id));
                    }
                }
            }
        }

        if let Some((u, v, chord_id)) = best_edge {
            let comp_u = Self::get_component(u, v, adj, &active_set);
            let comp_v = Self::get_component(v, u, adj, &active_set);

            let left = Self::decompose_subtree(adj, &comp_u, chords);
            let right = Self::decompose_subtree(adj, &comp_v, chords);

            CentroidTreeNode::Internal {
                chord_id,
                left_child: Box::new(left),
                right_child: Box::new(right),
            }
        } else {
            CentroidTreeNode::Leaf { region_id: active[0] }
        }
    }

    fn count_component(
        start: usize,
        blocked: usize,
        adj: &[Vec<(usize, usize)>],
        active: &std::collections::HashSet<usize>,
    ) -> usize {
        let mut visited = std::collections::HashSet::new();
        let mut stack = vec![start];
        visited.insert(start);

        while let Some(curr) = stack.pop() {
            for &(nbr, _) in &adj[curr] {
                if nbr != blocked && active.contains(&nbr) && visited.insert(nbr) {
                    stack.push(nbr);
                }
            }
        }
        visited.len()
    }

    fn get_component(
        start: usize,
        blocked: usize,
        adj: &[Vec<(usize, usize)>],
        active: &std::collections::HashSet<usize>,
    ) -> Vec<usize> {
        let mut visited = std::collections::HashSet::new();
        let mut stack = vec![start];
        visited.insert(start);

        let mut res = Vec::new();
        while let Some(curr) = stack.pop() {
            res.push(curr);
            for &(nbr, _) in &adj[curr] {
                if nbr != blocked && active.contains(&nbr) && visited.insert(nbr) {
                    stack.push(nbr);
                }
            }
        }
        res
    }
}
