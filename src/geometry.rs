//! Core geometric primitives for 2D polygon triangulation.

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Point {
    pub x: f64,
    pub y: f64,
}

impl Point {
    pub const EPSILON: f64 = 1e-11;

    pub fn new(x: f64, y: f64) -> Self {
        Self { x, y }
    }

    #[inline]
    pub fn approx_eq(&self, other: &Self) -> bool {
        (self.x - other.x).abs() <= Self::EPSILON && (self.y - other.y).abs() <= Self::EPSILON
    }

    #[inline]
    pub fn dist_sq(&self, other: &Self) -> f64 {
        let dx = self.x - other.x;
        let dy = self.y - other.y;
        dx * dx + dy * dy
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Segment {
    pub p1: Point,
    pub p2: Point,
}

impl Segment {
    pub fn new(p1: Point, p2: Point) -> Self {
        Self { p1, p2 }
    }
}

/// Computes the signed area of a 2D simple polygon.
/// Positive if counter-clockwise (CCW), negative if clockwise (CW).
pub fn signed_polygon_area(pts: &[Point]) -> f64 {
    let n = pts.len();
    if n < 3 {
        return 0.0;
    }
    let mut area = 0.0;
    for i in 0..n {
        let j = (i + 1) % n;
        area += pts[i].x * pts[j].y - pts[j].x * pts[i].y;
    }
    0.5 * area
}

/// Returns true if the vertices are in counter-clockwise order.
pub fn is_ccw(pts: &[Point]) -> bool {
    signed_polygon_area(pts) > 0.0
}

/// 2D cross product of vector OA and OB: (A.x - O.x)*(B.y - O.y) - (A.y - O.y)*(B.x - O.x)
#[inline]
pub fn cross2d(o: Point, a: Point, b: Point) -> f64 {
    (a.x - o.x) * (b.y - o.y) - (a.y - o.y) * (b.x - o.x)
}

/// Orientation test: >0 if CCW (left turn), <0 if CW (right turn), 0 if collinear.
#[inline]
pub fn orient2d(o: Point, a: Point, b: Point) -> f64 {
    cross2d(o, a, b)
}

/// Lexicographic comparison by y then x, with index tie-breaking to guarantee total ordering.
#[inline]
pub fn cmp_point_perturbed(p1: Point, idx1: usize, p2: Point, idx2: usize) -> std::cmp::Ordering {
    if (p1.y - p2.y).abs() > Point::EPSILON {
        p1.y.partial_cmp(&p2.y).unwrap_or(std::cmp::Ordering::Equal)
    } else if (p1.x - p2.x).abs() > Point::EPSILON {
        p1.x.partial_cmp(&p2.x).unwrap_or(std::cmp::Ordering::Equal)
    } else {
        idx1.cmp(&idx2)
    }
}

/// Determines if horizontal ray from `origin` in `dir_is_right` intersects segment `seg`.
/// Returns the intersection point and parameter t (distance along ray) if an intersection exists.
pub fn ray_h_intersect_segment(
    origin: Point,
    dir_is_right: bool,
    seg_p1: Point,
    seg_p2: Point,
) -> Option<(f64, Point)> {
    let y = origin.y;
    let y1 = seg_p1.y;
    let y2 = seg_p2.y;

    let min_y = y1.min(y2);
    let max_y = y1.max(y2);

    // If segment is strictly above or below ray y
    if y < min_y - Point::EPSILON || y > max_y + Point::EPSILON {
        return None;
    }

    // If segment is strictly horizontal
    if (y1 - y2).abs() <= Point::EPSILON {
        let min_x = seg_p1.x.min(seg_p2.x);
        let max_x = seg_p1.x.max(seg_p2.x);
        if dir_is_right {
            if max_x >= origin.x - Point::EPSILON {
                let hit_x = if origin.x > min_x { origin.x } else { min_x };
                let dist = hit_x - origin.x;
                return Some((dist.max(0.0), Point::new(hit_x, y)));
            }
        } else {
            if min_x <= origin.x + Point::EPSILON {
                let hit_x = if origin.x < max_x { origin.x } else { max_x };
                let dist = origin.x - hit_x;
                return Some((dist.max(0.0), Point::new(hit_x, y)));
            }
        }
        return None;
    }

    // Linear interpolation for x at y:
    let frac = (y - y1) / (y2 - y1);
    if frac < -Point::EPSILON || frac > 1.0 + Point::EPSILON {
        return None;
    }
    let hit_x = seg_p1.x + frac * (seg_p2.x - seg_p1.x);

    if dir_is_right {
        if hit_x >= origin.x - Point::EPSILON {
            let dist = hit_x - origin.x;
            Some((dist.max(0.0), Point::new(hit_x, y)))
        } else {
            None
        }
    } else {
        if hit_x <= origin.x + Point::EPSILON {
            let dist = origin.x - hit_x;
            Some((dist.max(0.0), Point::new(hit_x, y)))
        } else {
            None
        }
    }
}

/// Checks if segment AB strictly intersects segment CD (excluding endpoints).
pub fn segments_intersect_strict(a: Point, b: Point, c: Point, d: Point) -> bool {
    let o1 = orient2d(a, b, c);
    let o2 = orient2d(a, b, d);
    let o3 = orient2d(c, d, a);
    let o4 = orient2d(c, d, b);

    // Different sides
    if ((o1 > Point::EPSILON && o2 < -Point::EPSILON) || (o1 < -Point::EPSILON && o2 > Point::EPSILON))
        && ((o3 > Point::EPSILON && o4 < -Point::EPSILON) || (o3 < -Point::EPSILON && o4 > Point::EPSILON))
    {
        return true;
    }
    false
}

/// Checks if point p is inside CCW triangle (a, b, c).
#[inline]
pub fn point_in_triangle_ccw(p: Point, a: Point, b: Point, c: Point) -> bool {
    orient2d(a, b, p) >= -Point::EPSILON
        && orient2d(b, c, p) >= -Point::EPSILON
        && orient2d(c, a, p) >= -Point::EPSILON
}
