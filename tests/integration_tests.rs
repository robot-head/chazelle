use chazelle::{triangulate, Point, geometry::signed_polygon_area};

fn triangle_area(a: (f64, f64), b: (f64, f64), c: (f64, f64)) -> f64 {
    0.5 * ((b.0 - a.0) * (c.1 - a.1) - (b.1 - a.1) * (c.0 - a.0))
}

fn verify_triangulation(poly: &[(f64, f64)], triangles: &[[usize; 3]]) {
    let n = poly.len();
    assert_eq!(triangles.len(), n - 2, "Must produce exactly n - 2 triangles");

    let pts: Vec<Point> = poly.iter().map(|&(x, y)| Point::new(x, y)).collect();
    let poly_area = signed_polygon_area(&pts).abs();

    let mut sum_tri_area = 0.0;
    for &[a, b, c] in triangles {
        assert!(a < n && b < n && c < n, "Indices must be within bounds");
        assert!(a != b && b != c && a != c, "Triangle vertices must be distinct");

        let area = triangle_area(poly[a], poly[b], poly[c]).abs();
        assert!(area > 1e-12, "Triangle area must be positive, got {area}");
        sum_tri_area += area;
    }

    let diff = (sum_tri_area - poly_area).abs();
    assert!(
        diff < 1e-6,
        "Sum of triangle areas ({sum_tri_area}) must equal polygon area ({poly_area}), diff={diff}"
    );
}

#[test]
fn test_triangle() {
    let poly = vec![(0.0, 0.0), (2.0, 0.0), (1.0, 2.0)];
    let tris = triangulate(&poly).expect("Should triangulate triangle");
    verify_triangulation(&poly, &tris);
}

#[test]
fn test_square_ccw() {
    let poly = vec![(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)];
    let tris = triangulate(&poly).expect("Should triangulate square");
    verify_triangulation(&poly, &tris);
}

#[test]
fn test_square_cw() {
    let poly = vec![(0.0, 0.0), (0.0, 1.0), (1.0, 1.0), (1.0, 0.0)];
    let tris = triangulate(&poly).expect("Should triangulate CW square");
    verify_triangulation(&poly, &tris);
}

#[test]
fn test_regular_octagon() {
    let n = 8;
    let mut poly = Vec::new();
    for i in 0..n {
        let angle = 2.0 * std::f64::consts::PI * (i as f64) / (n as f64);
        poly.push((angle.cos(), angle.sin()));
    }
    let tris = triangulate(&poly).expect("Should triangulate octagon");
    verify_triangulation(&poly, &tris);
}

#[test]
fn test_l_shaped_polygon() {
    let poly = vec![
        (0.0, 0.0),
        (3.0, 0.0),
        (3.0, 1.0),
        (1.0, 1.0),
        (1.0, 3.0),
        (0.0, 3.0),
    ];
    let tris = triangulate(&poly).expect("Should triangulate L-shape");
    verify_triangulation(&poly, &tris);
}

#[test]
fn test_five_point_star() {
    // 5-pointed star (concave, 10 vertices)
    let mut poly = Vec::new();
    for i in 0..10 {
        let r = if i % 2 == 0 { 2.0 } else { 0.8 };
        let angle = i as f64 * std::f64::consts::PI / 5.0;
        poly.push((r * angle.cos(), r * angle.sin()));
    }
    let tris = triangulate(&poly).expect("Should triangulate star");
    verify_triangulation(&poly, &tris);
}

#[test]
fn test_comb_polygon() {
    // Comb polygon with teeth
    let mut poly = vec![(0.0, 0.0)];
    for i in 0..5 {
        let x0 = i as f64 * 2.0;
        let x1 = x0 + 1.0;
        poly.push((x0, 5.0));
        poly.push((x1, 5.0));
        poly.push((x1, 1.0));
    }
    poly.push((10.0, 0.0));
    let tris = triangulate(&poly).expect("Should triangulate comb");
    verify_triangulation(&poly, &tris);
}

#[test]
fn test_spiral_polygon() {
    // Spiral ribbon (concave, non-self-intersecting)
    let turns = 2.0;
    let num_pts = 20;
    let width = 0.5;

    let mut outer = Vec::new();
    let mut inner = Vec::new();
    for i in 0..num_pts {
        let theta = 2.0 * std::f64::consts::PI * turns * (i as f64) / ((num_pts - 1) as f64);
        let r1 = 2.0 + theta;
        let r2 = r1 + width;
        outer.push((r2 * theta.cos(), r2 * theta.sin()));
        inner.push((r1 * theta.cos(), r1 * theta.sin()));
    }
    inner.reverse();
    let mut poly = outer;
    poly.extend(inner);

    let tris = triangulate(&poly).expect("Should triangulate spiral ribbon");
    verify_triangulation(&poly, &tris);
}

#[test]
fn test_collinear_and_flat_edges() {
    let poly = vec![
        (0.0, 0.0),
        (2.0, 0.0),
        (4.0, 0.0), // Collinear on bottom edge
        (4.0, 3.0),
        (0.0, 3.0),
    ];
    let tris = triangulate(&poly).expect("Should triangulate polygon with collinear edge");
    verify_triangulation(&poly, &tris);
}

#[test]
fn test_serpentine_polygon() {
    let mut poly = vec![(0.0, 0.0)];
    let steps = 6;
    for i in 0..steps {
        let x = (i * 2) as f64;
        let y_top = 8.0;
        let y_bot = 0.0;
        poly.push((x, y_top));
        poly.push((x + 1.0, y_top));
        poly.push((x + 1.0, y_bot + 1.0));
        poly.push((x + 2.0, y_bot + 1.0));
    }
    poly.push((steps as f64 * 2.0, 0.0));
    let tris = triangulate(&poly).expect("Should triangulate serpentine polygon");
    verify_triangulation(&poly, &tris);
}

#[test]
fn test_error_too_few_vertices() {
    let poly = vec![(0.0, 0.0), (1.0, 1.0)];
    let res = triangulate(&poly);
    assert_eq!(res, Err(chazelle::TriangulationError::PolygonTooSmall));
}

#[test]
fn test_large_polygon_scaling() {
    let n = 200;
    let mut poly = Vec::with_capacity(n);
    for i in 0..n {
        let angle = 2.0 * std::f64::consts::PI * (i as f64) / (n as f64);
        let r = 100.0 + 10.0 * (5.0 * angle).sin();
        poly.push((r * angle.cos(), r * angle.sin()));
    }
    let tris = triangulate(&poly).expect("Should triangulate 200-gon");
    verify_triangulation(&poly, &tris);
}

#[test]
fn test_thousand_vertex_linear_scaling() {
    let n = 1000;
    let mut poly = Vec::with_capacity(n);
    for i in 0..n {
        let angle = 2.0 * std::f64::consts::PI * (i as f64) / (n as f64);
        let r = 500.0 + 50.0 * (10.0 * angle).cos();
        poly.push((r * angle.cos(), r * angle.sin()));
    }
    let start = std::time::Instant::now();
    let tris = triangulate(&poly).expect("Should triangulate 1000-gon");
    let elapsed = start.elapsed();
    verify_triangulation(&poly, &tris);
    println!("1,000-vertex polygon triangulated in {:?}", elapsed);
}

#[test]
fn test_zerocopy_raw_bytes() {
    let pts = vec![
        Point::new(0.0, 0.0),
        Point::new(4.0, 0.0),
        Point::new(4.0, 4.0),
        Point::new(0.0, 4.0),
    ];

    // Zero-copy serialization of Point slice to bytes
    let bytes = Point::slice_as_bytes(&pts);
    assert_eq!(bytes.len(), 4 * std::mem::size_of::<Point>());

    // Zero-copy deserialization and triangulation
    let tris = chazelle::triangulate_from_bytes(bytes).expect("Should triangulate from bytes");
    assert_eq!(tris.len(), 2);

    // Zero-copy serialization of triangles to bytes
    let tri_bytes = chazelle::ChazelleTriangle::slice_as_bytes(&tris);
    assert_eq!(tri_bytes.len(), 2 * std::mem::size_of::<chazelle::ChazelleTriangle>());

    // Zero-copy reading of triangles from bytes
    let tris_recovered = chazelle::ChazelleTriangle::slice_from_bytes(tri_bytes).unwrap();
    assert_eq!(tris_recovered.len(), 2);
    assert_eq!(tris_recovered[0], tris[0]);
    assert_eq!(tris_recovered[1], tris[1]);
}

#[test]
fn test_large_dataset_harmonic_10k() {
    let poly = chazelle::datasets::generate_harmonic_circle(
        10_000,
        1000.0,
        &[(150.0, 16.0), (50.0, 32.0)],
    );
    assert_eq!(poly.len(), 10_000);

    let start = std::time::Instant::now();
    let tris = chazelle::triangulate_points(&poly).expect("Should triangulate 10k harmonic polygon");
    let elapsed = start.elapsed();
    assert_eq!(tris.len(), 9_998);

    println!("10,000-vertex harmonic polygon triangulated in {:?}", elapsed);
}

#[test]
fn test_large_dataset_comb_10k() {
    let poly = chazelle::datasets::generate_comb(3333);
    assert_eq!(poly.len(), 10_001);

    let start = std::time::Instant::now();
    let tris = chazelle::triangulate_points(&poly).expect("Should triangulate 10k comb polygon");
    let elapsed = start.elapsed();
    assert_eq!(tris.len(), 9_999);

    println!("10,000-vertex comb polygon triangulated in {:?}", elapsed);
}

#[test]
fn test_large_dataset_binary_io_10k() {
    let poly = chazelle::datasets::generate_star(10_000, 500.0, 1000.0);
    let tmp_path = std::env::temp_dir().join("chazelle_test_poly_10k.bin");

    // Zero-copy binary export
    chazelle::datasets::save_polygon_binary(&poly, &tmp_path).expect("Save binary");

    // Zero-copy binary reload
    let loaded = chazelle::datasets::load_polygon_binary(&tmp_path).expect("Load binary");
    assert_eq!(loaded.len(), poly.len());
    assert_eq!(loaded, poly);

    let tris = chazelle::triangulate_points(&loaded).expect("Triangulate loaded");
    assert_eq!(tris.len(), 9_998);

    let _ = std::fs::remove_file(tmp_path);
}

#[test]
fn test_all_algorithms_on_concave_shapes() {
    let poly_l = vec![
        (0.0, 0.0),
        (3.0, 0.0),
        (3.0, 1.0),
        (1.0, 1.0),
        (1.0, 3.0),
        (0.0, 3.0),
    ];

    for &algo in &[
        chazelle::Algorithm::Chazelle,
        chazelle::Algorithm::MonotoneSweep,
        chazelle::Algorithm::Seidel,
    ] {
        let tris = chazelle::triangulate_with_algorithm(&poly_l, algo)
            .expect("Should triangulate L-shape");
        verify_triangulation(&poly_l, &tris);
    }

    // 5-point star
    let mut star = Vec::new();
    for i in 0..10 {
        let r = if i % 2 == 0 { 2.0 } else { 0.8 };
        let angle = i as f64 * std::f64::consts::PI / 5.0;
        star.push((r * angle.cos(), r * angle.sin()));
    }

    for &algo in &[
        chazelle::Algorithm::Chazelle,
        chazelle::Algorithm::MonotoneSweep,
        chazelle::Algorithm::Seidel,
    ] {
        let tris = chazelle::triangulate_with_algorithm(&star, algo)
            .expect("Should triangulate star");
        verify_triangulation(&star, &tris);
    }
}
