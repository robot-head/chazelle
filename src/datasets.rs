//! Large dataset generators for benchmarking and stress testing polygon triangulation.
//!
//! Generates orders-of-magnitude larger simple polygons (10k, 50k, 100k, 500k+ vertices)
//! across various geometric topologies:
//! - Harmonic circular waves (star-shaped with multiple frequencies)
//! - Archimedean spiral ribbons (winding corridors with high aspect ratios)
//! - Comb / sawtooth polygons (dense reflex vertices and horizontal visibility)
//! - Star polygons (alternating spikes)
//!
//! Includes zero-copy binary serialization and deserialization via `google/zerocopy`.

use std::fs::File;
use std::io::{Read, Write};
use std::path::Path;
use crate::geometry::Point;

/// Generates a simple harmonic circular polygon with `n` vertices.
///
/// Star-shaped around origin with radius $r(\theta) = R_0 + \sum A_k \sin(f_k \theta)$.
/// Guaranteed to have zero self-intersections as long as radius remains positive.
pub fn generate_harmonic_circle(
    n: usize,
    base_r: f64,
    harmonics: &[(f64, f64)], // (amplitude, frequency)
) -> Vec<Point> {
    let mut pts = Vec::with_capacity(n);
    let two_pi = 2.0 * std::f64::consts::PI;

    for i in 0..n {
        let theta = two_pi * (i as f64) / (n as f64);
        let mut r = base_r;
        for &(amp, freq) in harmonics {
            r += amp * (freq * theta).sin();
        }
        let r_safe = r.max(base_r * 0.1);
        pts.push(Point::new(r_safe * theta.cos(), r_safe * theta.sin()));
    }

    pts
}

/// Generates an Archimedean spiral ribbon polygon with `n` vertices.
///
/// Winding ribbon of `turns` complete rotations and given corridor `width`.
pub fn generate_spiral_ribbon(n: usize, turns: f64, width: f64) -> Vec<Point> {
    let half_n = n / 2;
    let two_pi = 2.0 * std::f64::consts::PI;
    let max_theta = two_pi * turns;
    let r0 = 100.0;
    let c = (width * 3.0) / two_pi; // Safe clearance between turns

    let mut outer = Vec::with_capacity(half_n);
    let mut inner = Vec::with_capacity(half_n);

    for i in 0..half_n {
        let frac = (i as f64) / ((half_n - 1) as f64);
        let theta = max_theta * frac;
        let r_inner = r0 + c * theta;
        let r_outer = r_inner + width;

        outer.push(Point::new(r_outer * theta.cos(), r_outer * theta.sin()));
        inner.push(Point::new(r_inner * theta.cos(), r_inner * theta.sin()));
    }

    inner.reverse();
    outer.extend(inner);
    outer
}

/// Generates a comb / sawtooth polygon with `num_teeth` teeth (approx 3 * num_teeth + 2 vertices).
pub fn generate_comb(num_teeth: usize) -> Vec<Point> {
    let mut pts = Vec::with_capacity(num_teeth * 3 + 2);
    pts.push(Point::new(0.0, 0.0));

    let w = 10.0;
    let tooth_height = 100.0;
    let valley_height = 20.0;

    for i in 0..num_teeth {
        let x0 = (i as f64) * 2.0 * w;
        let x1 = x0 + w;
        pts.push(Point::new(x0, tooth_height));
        pts.push(Point::new(x1, tooth_height));
        pts.push(Point::new(x1, valley_height));
    }

    pts.push(Point::new((num_teeth as f64) * 2.0 * w, 0.0));
    pts
}

/// Generates an alternating star polygon with `n` vertices (must be even).
pub fn generate_star(n: usize, r_inner: f64, r_outer: f64) -> Vec<Point> {
    let mut pts = Vec::with_capacity(n);
    let two_pi = 2.0 * std::f64::consts::PI;

    for i in 0..n {
        let theta = two_pi * (i as f64) / (n as f64);
        let r = if i % 2 == 0 { r_outer } else { r_inner };
        pts.push(Point::new(r * theta.cos(), r * theta.sin()));
    }

    pts
}

/// Zero-copy binary export of a polygon to disk.
pub fn save_polygon_binary<P: AsRef<Path>>(pts: &[Point], path: P) -> std::io::Result<()> {
    let bytes = Point::slice_as_bytes(pts);
    let mut file = File::create(path)?;
    file.write_all(bytes)?;
    Ok(())
}

/// Zero-copy binary load of a polygon from disk into memory.
pub fn load_polygon_binary<P: AsRef<Path>>(path: P) -> std::io::Result<Vec<Point>> {
    let mut file = File::open(path)?;
    let mut buf = Vec::new();
    file.read_to_end(&mut buf)?;

    let pts = Point::slice_from_bytes(&buf).ok_or_else(|| {
        std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            "Buffer size not a multiple of Point struct size",
        )
    })?;

    Ok(pts.to_vec())
}
