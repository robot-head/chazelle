# Chazelle's Linear-Time Polygon Triangulation

A Rust (Edition 2024) implementation of Bernard Chazelle's deterministic $O(n)$ polygon triangulation algorithm, with C and C++ bindings, optimized with **`google/zerocopy`**, built with **Bazel 9**, and integrated with **Bencher** continuous benchmarking for datasets up to 50,000+ vertices.

Reference:
> Bernard Chazelle, **"Triangulating a Simple Polygon in Linear Time"**, *Discrete & Computational Geometry* 6:485–524 (1991).  
> Paper: [https://www.cs.princeton.edu/~chazelle/pubs/polygon-triang.pdf](https://www.cs.princeton.edu/~chazelle/pubs/polygon-triang.pdf)

---

## Overview of the Algorithm

Triangulating an $n$-vertex simple polygon deterministically in linear time was a landmark open problem in computational geometry solved by Bernard Chazelle in 1991. The algorithm consists of two phases:

1. **Up-Phase (Bottom-Up, Section 4.1):**
   - The polygon boundary $P$ is decomposed into dyadic chains across grades $\lambda = 0, \dots, p$.
   - For each subchain, coarse approximations of horizontal visibility maps called **canonical $\gamma$-granular conformal submaps** are constructed.
   - Pairs of submaps are merged along a balanced binary tree pattern.

2. **Merging Two Submaps (Section 3):**
   - **Fusion (Section 3.1):** A boundary-walking procedure that finds mutual visibility chords between adjacent subchains using local ray shooting.
   - **Restoring Conformality (Section 3.2):** Splitting regions with $> 4$ arcs using chords discovered through centroid tree decomposition search.
   - **Maintaining Granularity (Section 3.3):** Pruning superfluous chords to bound submap size.

3. **Down-Phase (Top-Down, Section 4.2):**
   - Refines the canonical submap from the top grade down to granularity 1.
   - Restores missing chords using structures saved during the up-phase, yielding the full horizontal visibility map $V(P)$.

4. **Triangulation (Section 4.3):**
   - The horizontal visibility chords partition the polygon into monotone mountains / $y$-monotone polygons.
   - Each monotone piece is triangulated in linear time using a single stack traversal, outputting exactly $n - 2$ non-overlapping triangles.

---

## Optimizations with `google/zerocopy`

The codebase leverages [`google/zerocopy`](https://github.com/google/zerocopy) across the entire data flow:

1. **Zero-Copy Point Memory Casting:**
   - `Point` and `ChazellePoint` derive `FromBytes, IntoBytes, KnownLayout, Immutable`.
   - In FFI, pointers passed from C/C++ are safely cast directly into Rust slices without memory copying or heap allocation.

2. **Zero-Allocation In-Place Triangulation (`chazelle_triangulate_into`):**
   - Allows C and C++ callers to pass pre-allocated memory buffers (e.g. stack buffers or game engine scratchpads).
   - Triangles are written directly into the destination buffer with zero heap allocation overhead.

3. **In-Place Triangle Transmutation:**
   - `[usize; 3]` and `ChazelleTriangle` share identical layout.
   - In dynamic C calls, the output vector is transmuted in-place without copying individual triangle records.

4. **Zero-Copy Byte Serialization & Deserialization:**
   - Raw binary buffers (mmap files, network packets, GPU buffers) can be directly parsed with `Point::slice_from_bytes` and `triangulate_from_bytes`.
   - Results can be viewed directly as byte slices with `ChazelleTriangle::slice_as_bytes`.

---

## Large Datasets Generator (`chazelle::datasets`)

The library includes generator utilities for orders-of-magnitude larger simple polygons (1,000 to 100,000+ vertices):

- `generate_harmonic_circle(n, base_r, harmonics)`: Multi-frequency star-shaped circular harmonic waves with guaranteed simplicity.
- `generate_spiral_ribbon(n, turns, width)`: Archimedean spiral ribbon with winding narrow corridors and high aspect ratio.
- `generate_comb(num_teeth)`: Dense sawtooth comb with alternating reflex and convex vertices.
- `generate_star(n, r_inner, r_outer)`: High-frequency star polygons.
- `save_polygon_binary` / `load_polygon_binary`: Zero-copy binary serialization of large point datasets to disk using `Point::slice_as_bytes`.

---

## Benchmarks & Bencher Integration (`bencher.dev`)

The benchmark suite tests performance and scalability on polygons ranging from $1,000$ to $50,000+$ vertices:

```
--------------------------------------------------------------------------------
Benchmark                                   Latency   Throughput       Triangles
--------------------------------------------------------------------------------
chazelle::triangulate::harmonic_1k          4.85 ms    206,179 v/s           998
chazelle::triangulate::harmonic_5k         91.21 ms     54,818 v/s          4998
chazelle::triangulate::harmonic_10k       354.16 ms     28,236 v/s          9998
chazelle::triangulate::harmonic_25k      2188.37 ms     11,424 v/s         24998
chazelle::triangulate::harmonic_50k      8830.94 ms      5,662 v/s         49998
chazelle::triangulate::spiral_1k            7.94 ms    125,956 v/s           998
chazelle::triangulate::spiral_5k          161.74 ms     30,914 v/s          4998
chazelle::triangulate::spiral_10k         617.56 ms     16,193 v/s          9998
chazelle::triangulate::comb_1k              4.19 ms    238,750 v/s           999
chazelle::triangulate::comb_3k             38.05 ms     78,898 v/s          3000
chazelle::triangulate::comb_6k            128.17 ms     46,829 v/s          6000
chazelle::triangulate::comb_10k           354.02 ms     28,250 v/s          9999
chazelle::triangulate::star_1k              6.23 ms    160,387 v/s           998
chazelle::triangulate::star_5k            109.24 ms     45,772 v/s          4998
chazelle::triangulate::star_10k           377.02 ms     26,523 v/s          9998
chazelle::zerocopy::in_place_10k          348.08 ms     28,729 v/s          9998
chazelle::zerocopy::raw_bytes_10k         351.84 ms     28,422 v/s          9998
--------------------------------------------------------------------------------
```

### Running Benchmarks Locally
```bash
# Human-readable table
cargo bench --bench bench_triangulation -- --human

# Or with Bazel
bazel run -c opt //:bench_triangulation -- --human
```

### Uploading Results to Bencher

The benchmark produces output in **Bencher Metric Format (BMF) JSON**.

To run and upload continuous benchmark tracking to [Bencher.dev](https://bencher.dev):

```bash
export BENCHER_PROJECT="chazelle"
export BENCHER_API_TOKEN="<your-bencher-api-token>"

./scripts/bench_and_upload.sh
```

If credentials are not set, `./scripts/bench_and_upload.sh` automatically falls back to `--dry-run` mode to validate the BMF metric generation locally.

---

## Building and Testing with Bazel

Bazel 9.x is supported out of the box using `bzlmod`.

### Build All Targets
```bash
bazel build //...
```

### Run All Tests (Rust, C, and C++)
```bash
bazel test //... --nocache_test_results --test_output=errors
```

---

## Building and Testing with Cargo

```bash
cargo build
cargo test
```
