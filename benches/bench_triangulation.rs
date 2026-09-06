use std::time::Instant;
use chazelle::datasets::{generate_comb, generate_harmonic_circle, generate_spiral_ribbon, generate_star};
use chazelle::geometry::Point;
use chazelle::{triangulate_points, triangulate_from_bytes, ChazelleTriangle, chazelle_triangulate_into, ChazellePoint};

struct BenchmarkResult {
    name: String,
    mean_ms: f64,
    min_ms: f64,
    max_ms: f64,
    throughput_vps: f64,
    triangles: usize,
}

fn bench_case<F>(name: &str, num_vertices: usize, iterations: usize, mut f: F) -> BenchmarkResult
where
    F: FnMut() -> usize,
{
    let mut times_ms = Vec::with_capacity(iterations);
    let mut tris = 0;

    for _ in 0..iterations {
        let t0 = Instant::now();
        tris = f();
        let elapsed = t0.elapsed().as_secs_f64() * 1000.0;
        times_ms.push(elapsed);
    }

    let min_ms = times_ms.iter().copied().fold(f64::INFINITY, f64::min);
    let max_ms = times_ms.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let sum_ms: f64 = times_ms.iter().sum();
    let mean_ms = sum_ms / (iterations as f64);
    let throughput_vps = (num_vertices as f64) / (mean_ms / 1000.0);

    BenchmarkResult {
        name: name.to_string(),
        mean_ms,
        min_ms,
        max_ms,
        throughput_vps,
        triangles: tris,
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let human_readable = args.iter().any(|a| a == "--human" || a == "-h");

    if human_readable {
        eprintln!("Running Chazelle linear-time triangulation benchmarks (orders of magnitude larger datasets)...");
    }

    let mut results = Vec::new();

    // 1. Harmonic Circular Waves (1k, 5k, 10k, 25k, 50k vertices)
    for &n in &[1_000, 5_000, 10_000, 25_000, 50_000] {
        let harmonics = [(150.0, 16.0), (50.0, 32.0)];
        let poly = generate_harmonic_circle(n, 1000.0, &harmonics);
        let iters = if n >= 25_000 { 1 } else { 2 };

        results.push(bench_case(
            &format!("chazelle::triangulate::harmonic_{}k", n / 1000),
            n,
            iters,
            || triangulate_points(&poly).unwrap().len(),
        ));
    }

    // 2. Archimedean Spiral Ribbon (1k, 5k, 10k vertices)
    for &n in &[1_000, 5_000, 10_000] {
        let poly = generate_spiral_ribbon(n, 3.0, 20.0);
        results.push(bench_case(
            &format!("chazelle::triangulate::spiral_{}k", n / 1000),
            n,
            2,
            || triangulate_points(&poly).unwrap().len(),
        ));
    }

    // 3. Comb / Sawtooth Polygons (1k, 3k, 6k, 10k vertices)
    for &(teeth, label) in &[(333, "1k"), (1000, "3k"), (2000, "6k"), (3333, "10k")] {
        let poly = generate_comb(teeth);
        let n = poly.len();
        results.push(bench_case(
            &format!("chazelle::triangulate::comb_{label}"),
            n,
            2,
            || triangulate_points(&poly).unwrap().len(),
        ));
    }

    // 4. Alternating Stars (1k, 5k, 10k vertices)
    for &n in &[1_000, 5_000, 10_000] {
        let poly = generate_star(n, 500.0, 1000.0);
        results.push(bench_case(
            &format!("chazelle::triangulate::star_{}k", n / 1000),
            n,
            2,
            || triangulate_points(&poly).unwrap().len(),
        ));
    }

    // 5. Google/zerocopy In-Place Preallocated Triangulation (10k vertices)
    {
        let n = 10_000;
        let poly = generate_harmonic_circle(n, 1000.0, &[(100.0, 8.0)]);
        let mut out_buffer = vec![ChazelleTriangle::new(0, 0, 0); n - 2];
        let mut num_written = 0;

        results.push(bench_case(
            "chazelle::zerocopy::in_place_10k",
            n,
            2,
            || {
                unsafe {
                    let status = chazelle_triangulate_into(
                        poly.as_ptr() as *const ChazellePoint,
                        n,
                        out_buffer.as_mut_ptr(),
                        out_buffer.len(),
                        &mut num_written,
                    );
                    assert_eq!(status, chazelle::ChazelleStatus::Success);
                }
                num_written
            },
        ));
    }

    // 6. Google/zerocopy Raw Byte Parsing & Triangulation (10k vertices)
    {
        let n = 10_000;
        let poly = generate_harmonic_circle(n, 1000.0, &[(100.0, 8.0)]);
        let bytes = Point::slice_as_bytes(&poly);

        results.push(bench_case(
            "chazelle::zerocopy::raw_bytes_10k",
            n,
            2,
            || triangulate_from_bytes(bytes).unwrap().len(),
        ));
    }

    if human_readable {
        println!("\n{:-<80}", "");
        println!(
            "{:<40} {:>10} {:>12} {:>15}",
            "Benchmark", "Latency", "Throughput", "Triangles"
        );
        println!("{:-<80}", "");
        for r in &results {
            println!(
                "{:<40} {:>7.2} ms {:>9.0} v/s {:>15}",
                r.name, r.mean_ms, r.throughput_vps, r.triangles
            );
        }
        println!("{:-<80}", "");
    } else {
        // Output Bencher Metric Format (BMF) JSON
        println!("{{");
        for (i, r) in results.iter().enumerate() {
            let is_last = i == results.len() - 1;
            println!("  \"{}\": {{", r.name);
            println!("    \"latency\": {{");
            println!("      \"value\": {:.4},", r.mean_ms);
            println!("      \"lower_value\": {:.4},", r.min_ms);
            println!("      \"upper_value\": {:.4}", r.max_ms);
            println!("    }},");
            println!("    \"throughput\": {{");
            println!("      \"value\": {:.1}", r.throughput_vps);
            println!("    }}");
            if is_last {
                println!("  }}");
            } else {
                println!("  }},");
            }
        }
        println!("}}");
    }
}
