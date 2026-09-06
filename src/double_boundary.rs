//! Double boundary representation and vertex chord classification according to Chazelle Section 2.1.

use crate::geometry::{Point, cross2d};

/// Direction of a horizontal chord from a vertex.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ChordDirection {
    Left,
    Right,
}

impl ChordDirection {
    pub fn is_right(&self) -> bool {
        matches!(self, ChordDirection::Right)
    }

    pub fn dx(&self) -> f64 {
        match self {
            ChordDirection::Left => -1.0,
            ChordDirection::Right => 1.0,
        }
    }
}

/// A chord endpoint or chord specification from a polygon vertex.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct VertexChordSpec {
    pub vertex_idx: usize,
    pub direction: ChordDirection,
}

/// Checks whether a given direction vector `d` points into the interior of a CCW polygon at vertex `i`.
///
/// In a CCW polygon, the interior is to the left of the edge `(prev -> curr)` and `(curr -> next)`.
/// Let `u = prev - curr` (pointing backwards along incoming edge)
/// Let `w = next - curr` (pointing forwards along outgoing edge)
/// The interior angular sector is from `w` to `u` counterclockwise.
pub fn is_direction_interior_ccw(curr: Point, prev: Point, next: Point, dir: ChordDirection) -> bool {
    let d = Point::new(dir.dx(), 0.0);
    let u = Point::new(prev.x - curr.x, prev.y - curr.y);
    let w = Point::new(next.x - curr.x, next.y - curr.y);

    // Cross products with origin (0,0)
    let zero = Point::new(0.0, 0.0);
    let w_cross_u = cross2d(zero, w, u);
    let w_cross_d = cross2d(zero, w, d);
    let d_cross_u = cross2d(zero, d, u);

    const EPS: f64 = Point::EPSILON;

    if w_cross_u > EPS {
        // Strictly convex vertex (< 180 deg interior angle)
        // Ray d must be strictly inside the cone from w to u
        w_cross_d > EPS && d_cross_u > EPS
    } else if w_cross_u < -EPS {
        // Reflex vertex (> 180 deg interior angle)
        // Ray d is interior if it is NOT in the exterior cone (which is from u to w)
        !(w_cross_d <= -EPS && d_cross_u <= -EPS)
    } else {
        // Degenerate / collinear incoming and outgoing edges
        // If u and w point in opposite directions (straight line):
        // interior is to the left of (prev -> next)
        let fwd = Point::new(next.x - prev.x, next.y - prev.y);
        cross2d(zero, fwd, d) > EPS
    }
}

/// Represents the double boundary of a subchain C of polygon P (Section 2.1 & 2.4).
#[derive(Debug, Clone)]
pub struct DoubleBoundary<'a> {
    pub vertices: &'a [Point],
    pub chain_start: usize,
    pub chain_end: usize,
}

impl<'a> DoubleBoundary<'a> {
    pub fn new(vertices: &'a [Point], chain_start: usize, chain_end: usize) -> Self {
        Self {
            vertices,
            chain_start,
            chain_end,
        }
    }

    /// Number of edges in the subchain.
    pub fn num_edges(&self) -> usize {
        if self.chain_end >= self.chain_start {
            self.chain_end - self.chain_start
        } else {
            (self.vertices.len() - self.chain_start) + self.chain_end
        }
    }

    /// Gets the polygon vertex at chain offset `i`.
    pub fn vertex(&self, i: usize) -> Point {
        let n = self.vertices.len();
        self.vertices[(self.chain_start + i) % n]
    }
}
