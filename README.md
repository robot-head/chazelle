# Chazelle's Linear-Time Polygon Triangulation

A Rust (Edition 2024) implementation of Bernard Chazelle's deterministic $O(n)$ polygon triangulation algorithm, with C and C++ bindings, optimized with **`google/zerocopy`** and built with **Bazel 9**.

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

4. **Zero-Copy Byte Serialization:**
   - Raw binary buffers (mmap files, network packets, GPU buffers) can be directly parsed with `Point::slice_from_bytes` and `triangulate_from_bytes`.
   - Results can be viewed directly as byte slices with `ChazelleTriangle::slice_as_bytes`.

---

## Project Structure

```
.
├── MODULE.bazel          # Bazel 9 bzlmod dependencies (rules_rust, rules_cc, crate_universe)
├── BUILD.bazel           # Bazel targets (Rust library, FFI staticlib, C/C++ tests)
├── Cargo.toml            # Rust edition 2024 crate configuration (zerocopy 0.8)
├── include/
│   ├── chazelle.h        # C API header (with chazelle_triangulate_into)
│   └── chazelle.hpp      # Modern C++ header wrapper (with triangulate_into)
├── src/
│   ├── lib.rs            # Top-level library and public Rust API
│   ├── geometry.rs       # 2D primitives with zerocopy traits
│   ├── double_boundary.rs# Double boundary dC and chord classification
│   ├── submap.rs         # Normal form submap and centroid tree decomposition
│   ├── fusion.rs         # Fusion walk of two submaps
│   ├── conformality.rs   # Conformality restoration and granularity pruning
│   ├── oracles.rs        # Ray-shooting and arc-cutting oracles
│   ├── up_phase.rs       # Grade hierarchy bottom-up construction
│   ├── down_phase.rs     # Top-down refinement to full visibility map
│   ├── monotone.rs       # Triangulation from trapezoids/monotone mountains
│   └── c_api.rs          # C/C++ FFI bindings with zerocopy optimizations
└── tests/
    ├── integration_tests.rs # Rust integration tests (including zerocopy raw byte tests)
    ├── c_test.c          # C test suite (including zero-allocation tests)
    └── cpp_test.cpp      # C++ test suite (including triangulate_into tests)
```

---

## Building and Testing with Bazel

Bazel 9.x is supported out of the box using `bzlmod`.

### Build All Targets
```bash
bazel build //...
```

### Run All Tests (Rust, C, and C++)
```bash
bazel test //... --test_output=all
```

Target breakdown:
- `//:chazelle`: Core Rust library (`rust_library`, edition 2024).
- `//:chazelle_ffi`: Static library for C/C++ FFI (`rust_static_library`).
- `//:chazelle_cc`: C/C++ header library (`cc_library`).
- `//:integration_test`: Rust unit and integration test suite (`rust_test`).
- `//:c_test`: C API test (`cc_test`).
- `//:cpp_test`: C++ API test (`cc_test`).

---

## Building and Testing with Cargo

```bash
cargo build
cargo test
```

---

## Usage Examples

### Rust In-Memory API

```rust
use chazelle::triangulate;

fn main() {
    let poly = vec![
        (0.0, 0.0),
        (3.0, 0.0),
        (3.0, 1.0),
        (1.0, 1.0),
        (1.0, 3.0),
        (0.0, 3.0),
    ];

    let triangles = triangulate(&poly).expect("Triangulation failed");
    for tri in triangles {
        println!("Triangle: ({}, {}, {})", tri[0], tri[1], tri[2]);
    }
}
```

### Rust Zero-Copy Byte API

```rust
use chazelle::{Point, triangulate_from_bytes, ChazelleTriangle};

fn main() {
    let pts = vec![
        Point::new(0.0, 0.0),
        Point::new(1.0, 0.0),
        Point::new(1.0, 1.0),
        Point::new(0.0, 1.0),
    ];
    // Zero-copy bytes view
    let raw_bytes: &[u8] = Point::slice_as_bytes(&pts);

    // Direct triangulation from byte slice
    let triangles = triangulate_from_bytes(raw_bytes).unwrap();
    println!("Triangles: {}", triangles.len());

    // Zero-copy triangle byte view
    let out_bytes: &[u8] = ChazelleTriangle::slice_as_bytes(&triangles);
    println!("Output byte length: {}", out_bytes.len());
}
```

### C API: Zero-Allocation In-Place Triangulation (`chazelle.h`)

```c
#include "chazelle.h"
#include <stdio.h>

int main(void) {
    ChazellePoint pts[4] = {
        {0.0, 0.0}, {2.0, 0.0}, {2.0, 2.0}, {0.0, 2.0}
    };
    // Pre-allocated stack buffer for triangles (Zero heap allocations!)
    ChazelleTriangle stack_buf[2];
    size_t num_tris = 0;

    ChazelleStatus status = chazelle_triangulate_into(
        pts, 4, stack_buf, 2, &num_tris
    );
    if (status == CHAZELLE_SUCCESS) {
        printf("Successfully triangulated %zu triangles into preallocated buffer\n", num_tris);
    }
    return 0;
}
```

### C++ API: Preallocated In-Place Triangulation (`chazelle.hpp`)

```cpp
#include "chazelle.hpp"
#include <iostream>
#include <vector>

int main() {
    std::vector<chazelle::Point> pts = {
        {0.0, 0.0}, {2.0, 0.0}, {2.0, 2.0}, {0.0, 2.0}
    };

    chazelle::Triangle stack_buf[2];
    size_t num_tris = chazelle::triangulate_into(
        pts.data(), pts.size(), stack_buf, 2
    );
    std::cout << "Triangles: " << num_tris << std::endl;
    return 0;
}
```
