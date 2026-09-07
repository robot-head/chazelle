use std::time::Instant;
use chazelle::datasets::{generate_comb, generate_harmonic_circle, generate_spiral_ribbon, generate_star};
use chazelle::geometry::Point;
use chazelle::{
    Algorithm, triangulate_points_with_algorithm, triangulate_from_bytes,
    ChazelleTriangle, chazelle_triangulate_into, ChazellePoint,
};

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
        eprintln!("Benchmarking Chazelle O(n) vs Monotone Sweep O(n log n) vs Seidel O(n log* n)...");
    }

    let mut results = Vec::new();

    // 1. Comparison across Algorithms on Harmonic Waves (1k, 5k, 10k)
    for &n in &[1_000, 5_000, 10_000] {
        let harmonics = [(150.0, 16.0), (50.0, 32.0)];
        let poly = generate_harmonic_circle(n, 1000.0, &harmonics);

        // Chazelle O(n)
        results.push(bench_case(
            &format!("chazelle::harmonic_{}k", n / 1000),
            n,
            2,
            || triangulate_points_with_algorithm(&poly, Algorithm::Chazelle).unwrap().len(),
        ));

        // Monotone Sweep O(n log n)
        results.push(bench_case(
            &format!("monotone_sweep::harmonic_{}k", n / 1000),
            n,
            2,
            || triangulate_points_with_algorithm(&poly, Algorithm::MonotoneSweep).unwrap().len(),
        ));

        // Seidel O(n log* n)
        results.push(bench_case(
            &format!("seidel::harmonic_{}k", n / 1000),
            n,
            2,
            || triangulate_points_with_algorithm(&poly, Algorithm::Seidel).unwrap().len(),
        ));
    }

    // 2. Comparison across Algorithms on Spiral Ribbon (1k, 5k)
    for &n in &[1_000, 5_000] {
        let poly = generate_spiral_ribbon(n, 3.0, 20.0);

        results.push(bench_case(
            &format!("chazelle::spiral_{}k", n / 1000),
            n,
            2,
            || triangulate_points_with_algorithm(&poly, Algorithm::Chazelle).unwrap().len(),
        ));

        results.push(bench_case(
            &format!("monotone_sweep::spiral_{}k", n / 1000),
            n,
            2,
            || triangulate_points_with_algorithm(&poly, Algorithm::MonotoneSweep).unwrap().len(),
        ));

        results.push(bench_case(
            &format!("seidel::spiral_{}k", n / 1000),
            n,
            2,
            || triangulate_points_with_algorithm(&poly, Algorithm::Seidel).unwrap().len(),
        ));
    }

    // 3. Comparison across Algorithms on Comb Polygon (1k, 3k)
    for &(teeth, label) in &[(333, "1k"), (1000, "3k")] {
        let poly = generate_comb(teeth);
        let n = poly.len();

        results.push(bench_case(
            &format!("chazelle::comb_{label}"),
            n,
            2,
            || triangulate_points_with_algorithm(&poly, Algorithm::Chazelle).unwrap().len(),
        ));

        results.push(bench_case(
            &format!("monotone_sweep::comb_{label}"),
            n,
            2,
            || triangulate_points_with_algorithm(&poly, Algorithm::MonotoneSweep).unwrap().len(),
        ));

        results.push(bench_case(
            &format!("seidel::comb_{label}"),
            n,
            2,
            || triangulate_points_with_algorithm(&poly, Algorithm::Seidel).unwrap().len(),
        ));
    }

    // 4. Comparison across Algorithms on Star Polygon (1k, 5k)
    for &n in &[1_000, 5_000] {
        let poly = generate_star(n, 500.0, 1000.0);

        results.push(bench_case(
            &format!("chazelle::star_{}k", n / 1000),
            n,
            2,
            || triangulate_points_with_algorithm(&poly, Algorithm::Chazelle).unwrap().len(),
        ));

        results.push(bench_case(
            &format!("monotone_sweep::star_{}k", n / 1000),
            n,
            2,
            || triangulate_points_with_algorithm(&poly, Algorithm::MonotoneSweep).unwrap().len(),
        ));

        results.push(bench_case(
            &format!("seidel::star_{}k", n / 1000),
            n,
            2,
            || triangulate_points_with_algorithm(&poly, Algorithm::Seidel).unwrap().len(),
        ));
    }

    // 5. Large-scale Chazelle O(n) scaling up to 50k vertices
    for &n in &[25_000, 50_000] {
        let harmonics = [(150.0, 16.0), (50.0, 32.0)];
        let poly = generate_harmonic_circle(n, 1000.0, &harmonics);

        results.push(bench_case(
            &format!("chazelle::scaling_harmonic_{}k", n / 1000),
            n,
            1,
            || triangulate_points_with_algorithm(&poly, Algorithm::Chazelle).unwrap().len(),
        ));
    }

    // 5. Zero-Copy In-Place vs Raw Bytes (10k vertices)
    {
        let n = 10_000;
        let poly = generate_harmonic_circle(n, 1000.0, &[(100.0, 8.0)]);
        let mut out_buffer = vec![ChazelleTriangle::new(0, 0, 0); n - 2];
        let mut num_written = 0;

        results.push(bench_case(
            "zerocopy::in_place_10k",
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

        let bytes = Point::slice_as_bytes(&poly);
        results.push(bench_case(
            "zerocopy::raw_bytes_10k",
            n,
            2,
            || triangulate_from_bytes(bytes).unwrap().len(),
        ));
    }

    // 6. Multithreading / Thread Scaling on Large Polygons (10k Harmonic)
    {
        let n = 10_000;
        let harmonics = [(150.0, 16.0), (50.0, 32.0)];
        let poly = generate_harmonic_circle(n, 1000.0, &harmonics);

        let pool_1 = rayon::ThreadPoolBuilder::new().num_threads(1).build().unwrap();
        results.push(bench_case(
            "chazelle::scaling::1_thread_10k",
            n,
            2,
            || pool_1.install(|| triangulate_points_with_algorithm(&poly, Algorithm::Chazelle).unwrap().len()),
        ));

        let pool_4 = rayon::ThreadPoolBuilder::new().num_threads(4).build().unwrap();
        results.push(bench_case(
            "chazelle::scaling::4_threads_10k",
            n,
            2,
            || pool_4.install(|| triangulate_points_with_algorithm(&poly, Algorithm::Chazelle).unwrap().len()),
        ));

        results.push(bench_case(
            "chazelle::scaling::all_cores_10k",
            n,
            2,
            || triangulate_points_with_algorithm(&poly, Algorithm::Chazelle).unwrap().len(),
        ));
    }

    if human_readable {
        println!("\n{:-<85}", "");
        println!(
            "{:<38} {:>10} {:>14} {:>15}",
            "Benchmark", "Latency", "Throughput", "Triangles"
        );
        println!("{:-<85}", "");
        for r in &results {
            println!(
                "{:<38} {:>7.2} ms {:>11.0} v/s {:>15}",
                r.name, r.mean_ms, r.throughput_vps, r.triangles
            );
        }
        println!("{:-<85}", "");
    } else {
        // Bencher Metric Format (BMF) JSON output
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
