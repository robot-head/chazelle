# Chazelle's Linear-Time Polygon Triangulation

A Rust (Edition 2024) implementation of Bernard Chazelle's deterministic $O(n)$ polygon triangulation algorithm, with C and C++ bindings and Bazel build integration.

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

## Project Structure

```
.
├── MODULE.bazel          # Bazel 9 bzlmod dependencies (rules_rust, rules_cc)
├── BUILD.bazel           # Bazel targets (Rust library, FFI staticlib, C/C++ tests)
├── Cargo.toml            # Rust edition 2024 crate configuration
├── include/
│   ├── chazelle.h        # C API header
│   └── chazelle.hpp      # Modern C++ header wrapper
├── src/
│   ├── lib.rs            # Top-level library and public Rust API
│   ├── geometry.rs       # 2D primitives, orientation, intersection tests
│   ├── double_boundary.rs# Double boundary dC and chord classification
│   ├── submap.rs         # Normal form submap and centroid tree decomposition
│   ├── fusion.rs         # Fusion walk of two submaps
│   ├── conformality.rs   # Conformality restoration and granularity pruning
│   ├── oracles.rs        # Ray-shooting and arc-cutting oracles
│   ├── up_phase.rs       # Grade hierarchy bottom-up construction
│   ├── down_phase.rs     # Top-down refinement to full visibility map
│   ├── monotone.rs       # Triangulation from trapezoids/monotone mountains
│   └── c_api.rs          # C/C++ FFI bindings (Rust 2024 #[unsafe(no_mangle)])
└── tests/
    ├── integration_tests.rs # Rust integration tests (various topologies & 1k vertices)
    ├── c_test.c          # C binding test suite
    └── cpp_test.cpp      # C++ binding test suite
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

### Rust API

```rust
use chazelle::triangulate;

fn main() {
    // Define simple polygon vertices (CW or CCW)
    let poly = vec![
        (0.0, 0.0),
        (3.0, 0.0),
        (3.0, 1.0),
        (1.0, 1.0),
        (1.0, 3.0),
        (0.0, 3.0),
    ];

    // Returns Vec<[usize; 3]> containing indices of triangle vertices
    let triangles = triangulate(&poly).expect("Triangulation failed");
    for tri in triangles {
        println!("Triangle: ({}, {}, {})", tri[0], tri[1], tri[2]);
    }
}
```

### C API (`chazelle.h`)

```c
#include "chazelle.h"
#include <stdio.h>

int main(void) {
    ChazellePoint pts[] = {
        {0.0, 0.0}, {2.0, 0.0}, {2.0, 2.0}, {0.0, 2.0}
    };
    ChazelleTriangle* tris = NULL;
    size_t num_tris = 0;

    ChazelleStatus status = chazelle_triangulate(pts, 4, &tris, &num_tris);
    if (status == CHAZELLE_SUCCESS) {
        printf("Produced %zu triangles\n", num_tris);
        chazelle_free_triangles(tris, num_tris);
    }
    return 0;
}
```

### C++ API (`chazelle.hpp`)

```cpp
#include "chazelle.hpp"
#include <iostream>
#include <vector>

int main() {
    std::vector<chazelle::Point> pts = {
        {0.0, 0.0}, {2.0, 0.0}, {2.0, 2.0}, {0.0, 2.0}
    };

    try {
        std::vector<chazelle::Triangle> tris = chazelle::triangulate(pts);
        std::cout << "Triangles: " << tris.size() << std::endl;
    } catch (const chazelle::TriangulationError& err) {
        std::cerr << "Error: " << err.what() << std::endl;
    }
    return 0;
}
```
