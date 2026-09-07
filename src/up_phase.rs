//! Bottom-up phase (Section 4.1): computes canonical submaps along dyadic chain hierarchy.

use crate::geometry::Point;
use crate::submap::Submap;
use crate::fusion::{fuse_submaps, build_submap_from_chords};
use crate::conformality::{restore_conformality, maintain_granularity};

/// A dyadic chain node in the grade hierarchy.
#[derive(Debug, Clone)]
pub struct ChainGradeNode {
    pub grade: usize,
    pub start_vertex: usize,
    pub num_edges: usize,
    pub min_y: f64,
    pub max_y: f64,
    pub submap: Submap,
}

/// Up-phase hierarchy: stores computed canonical submaps across all grades lambda = 0..=p.
pub struct UpPhaseHierarchy {
    pub grades: Vec<Vec<ChainGradeNode>>,
}

impl UpPhaseHierarchy {
    /// Executes the up-phase on polygon P (Section 4.1).
    pub fn build(polygon: &[Point]) -> Self {
        let oracle = crate::oracles::RayShootingOracle::new(polygon);
        Self::build_with_oracle(polygon, &oracle)
    }

    /// Executes the up-phase with a prebuilt spatial oracle.
    pub fn build_with_oracle(polygon: &[Point], oracle: &crate::oracles::RayShootingOracle) -> Self {
        let n = polygon.len();
        if n < 3 {
            return Self { grades: Vec::new() };
        }

        let mut grades: Vec<Vec<ChainGradeNode>> = Vec::new();

        // Grade 0: base chains of small length (16 edges)
        let base_len = 16;
        let mut grade_0 = Vec::new();
        let mut curr = 0;
        while curr < n {
            let edges = (n - curr).min(base_len);
            let mut min_y = polygon[curr].y;
            let mut max_y = polygon[curr].y;
            for offset in 0..=edges {
                let y = polygon[(curr + offset) % n].y;
                if y < min_y { min_y = y; }
                if y > max_y { max_y = y; }
            }

            let submap = build_submap_from_chords(Vec::new(), curr, edges, polygon, 1);
            grade_0.push(ChainGradeNode {
                grade: 0,
                start_vertex: curr,
                num_edges: edges,
                min_y,
                max_y,
                submap,
            });
            curr += edges;
        }
        grades.push(grade_0);

        // Successive grades lambda = 1, 2, ...
        let mut current_grade = 0;
        while grades[current_grade].len() > 1 {
            let prev_nodes = &grades[current_grade];
            let mut next_nodes = Vec::new();
            let gamma = (1usize << ((current_grade + 1) / 4)).max(1);

            let mut i = 0;
            while i < prev_nodes.len() {
                if i + 1 < prev_nodes.len() {
                    let c1 = &prev_nodes[i];
                    let c2 = &prev_nodes[i + 1];

                    let fused = fuse_submaps(
                        &c1.submap,
                        &c2.submap,
                        polygon,
                        oracle,
                        c1.start_vertex,
                        c1.num_edges,
                        c1.min_y,
                        c1.max_y,
                        c2.start_vertex,
                        c2.num_edges,
                        c2.min_y,
                        c2.max_y,
                    );

                    let conformal = restore_conformality(fused, polygon);
                    let granular = maintain_granularity(conformal, gamma);

                    next_nodes.push(ChainGradeNode {
                        grade: current_grade + 1,
                        start_vertex: c1.start_vertex,
                        num_edges: c1.num_edges + c2.num_edges,
                        min_y: c1.min_y.min(c2.min_y),
                        max_y: c1.max_y.max(c2.max_y),
                        submap: granular,
                    });
                    i += 2;
                } else {
                    // Carry forward odd node
                    let mut node = prev_nodes[i].clone();
                    node.grade = current_grade + 1;
                    next_nodes.push(node);
                    i += 1;
                }
            }

            grades.push(next_nodes);
            current_grade += 1;
        }

        Self { grades }
    }

    /// Gets the top-level submap covering the entire polygon.
    pub fn top_submap(&self) -> Option<&Submap> {
        self.grades.last().and_then(|g| g.first()).map(|node| &node.submap)
    }
}
