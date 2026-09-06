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
    pub submap: Submap,
}

/// Up-phase hierarchy: stores computed canonical submaps across all grades lambda = 0..=p.
pub struct UpPhaseHierarchy {
    pub grades: Vec<Vec<ChainGradeNode>>,
}

impl UpPhaseHierarchy {
    /// Executes the up-phase on polygon P (Section 4.1).
    pub fn build(polygon: &[Point]) -> Self {
        let n = polygon.len();
        if n < 3 {
            return Self { grades: Vec::new() };
        }

        let mut grades: Vec<Vec<ChainGradeNode>> = Vec::new();

        // Grade 0: base chains of small length (e.g. 4 edges)
        let base_len = 4;
        let mut grade_0 = Vec::new();
        let mut curr = 0;
        while curr < n {
            let edges = (n - curr).min(base_len);
            let submap = build_submap_from_chords(Vec::new(), curr, edges, polygon, 1);
            grade_0.push(ChainGradeNode {
                grade: 0,
                start_vertex: curr,
                num_edges: edges,
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
                        c1.start_vertex,
                        c1.num_edges,
                        c2.start_vertex,
                        c2.num_edges,
                    );

                    let conformal = restore_conformality(fused, polygon);
                    let granular = maintain_granularity(conformal, gamma);

                    next_nodes.push(ChainGradeNode {
                        grade: current_grade + 1,
                        start_vertex: c1.start_vertex,
                        num_edges: c1.num_edges + c2.num_edges,
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
