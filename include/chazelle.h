/**
 * @file chazelle.h
 * @brief C bindings for Chazelle's linear-time polygon triangulation algorithm,
 * with support for Monotone Plane-Sweep and Seidel's randomized triangulation algorithms.
 */

#ifndef CHAZELLE_H
#define CHAZELLE_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    double x;
    double y;
} ChazellePoint;

typedef struct {
    size_t a;
    size_t b;
    size_t c;
} ChazelleTriangle;

typedef enum {
    CHAZELLE_SUCCESS = 0,
    CHAZELLE_ERROR_INVALID_ARGUMENT = 1,
    CHAZELLE_ERROR_POLYGON_TOO_SMALL = 2,
    CHAZELLE_ERROR_DEGENERATE = 3,
    CHAZELLE_ERROR_INTERNAL = 4,
    CHAZELLE_ERROR_BUFFER_TOO_SMALL = 5
} ChazelleStatus;

/**
 * @brief Choice of triangulation algorithm.
 */
typedef enum {
    CHAZELLE_ALGORITHM_CHAZELLE = 0,       /**< Chazelle deterministic linear-time O(n) */
    CHAZELLE_ALGORITHM_MONOTONE_SWEEP = 1,  /**< Classic plane sweep monotone decomposition O(n log n) */
    CHAZELLE_ALGORITHM_SEIDEL = 2          /**< Seidel's randomized incremental O(n log* n) */
} ChazelleAlgorithm;

/**
 * @brief Zero-allocation triangulation writing directly into caller-provided buffer with algorithm selection.
 */
ChazelleStatus chazelle_triangulate_into_with_algorithm(
    const ChazellePoint* points,
    size_t num_points,
    ChazelleTriangle* out_triangles,
    size_t out_capacity,
    size_t* out_num_triangles,
    ChazelleAlgorithm algorithm
);

/**
 * @brief In-place zero-allocation triangulation using default Chazelle algorithm.
 */
ChazelleStatus chazelle_triangulate_into(
    const ChazellePoint* points,
    size_t num_points,
    ChazelleTriangle* out_triangles,
    size_t out_capacity,
    size_t* out_num_triangles
);

/**
 * @brief Triangulate a simple polygon using dynamic memory allocation and algorithm selection.
 */
ChazelleStatus chazelle_triangulate_with_algorithm(
    const ChazellePoint* points,
    size_t num_points,
    ChazelleTriangle** out_triangles,
    size_t* out_num_triangles,
    ChazelleAlgorithm algorithm
);

/**
 * @brief Triangulate a simple polygon using default Chazelle algorithm.
 */
ChazelleStatus chazelle_triangulate(
    const ChazellePoint* points,
    size_t num_points,
    ChazelleTriangle** out_triangles,
    size_t* out_num_triangles
);

void chazelle_free_triangles(ChazelleTriangle* triangles, size_t num_triangles);

const char* chazelle_version(void);

#ifdef __cplusplus
}
#endif

#endif /* CHAZELLE_H */
