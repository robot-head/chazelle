/**
 * @file chazelle.h
 * @brief C bindings for Chazelle's linear-time polygon triangulation algorithm.
 */

#ifndef CHAZELLE_H
#define CHAZELLE_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief 2D point with double-precision floating-point coordinates.
 */
typedef struct {
    double x;
    double y;
} ChazellePoint;

/**
 * @brief Triangle represented by three vertex indices referencing the input polygon.
 */
typedef struct {
    size_t a;
    size_t b;
    size_t c;
} ChazelleTriangle;

/**
 * @brief Status and error codes returned by the Chazelle triangulation library.
 */
typedef enum {
    CHAZELLE_SUCCESS = 0,
    CHAZELLE_ERROR_INVALID_ARGUMENT = 1,
    CHAZELLE_ERROR_POLYGON_TOO_SMALL = 2,
    CHAZELLE_ERROR_DEGENERATE = 3,
    CHAZELLE_ERROR_INTERNAL = 4
} ChazelleStatus;

/**
 * @brief Triangulate a simple polygon in linear time using Chazelle's algorithm.
 *
 * @param points Pointer to an array of vertices forming the simple polygon (in CW or CCW order).
 * @param num_points The number of vertices (must be >= 3).
 * @param out_triangles Pointer to receive the allocated array of output triangles.
 *                      The caller must free this array using chazelle_free_triangles.
 * @param out_num_triangles Pointer to receive the number of triangles (will be num_points - 2 on success).
 * @return ChazelleStatus CHAZELLE_SUCCESS on success, or an error code.
 */
ChazelleStatus chazelle_triangulate(
    const ChazellePoint* points,
    size_t num_points,
    ChazelleTriangle** out_triangles,
    size_t* out_num_triangles
);

/**
 * @brief Free the triangle array allocated by chazelle_triangulate.
 *
 * @param triangles Pointer to the triangle array, or NULL.
 * @param num_triangles The number of triangles in the array.
 */
void chazelle_free_triangles(ChazelleTriangle* triangles, size_t num_triangles);

/**
 * @brief Returns the version string of the Chazelle triangulation library.
 */
const char* chazelle_version(void);

#ifdef __cplusplus
}
#endif

#endif /* CHAZELLE_H */
