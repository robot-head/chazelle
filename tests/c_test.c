#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "include/chazelle.h"

static double triangle_area(ChazellePoint a, ChazellePoint b, ChazellePoint c) {
    return 0.5 * fabs((b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x));
}

int main(void) {
    printf("Chazelle version: %s\n", chazelle_version());

    // Test 1: NULL pointer handling
    {
        ChazelleTriangle* tris = NULL;
        size_t num_tris = 0;
        ChazelleStatus status = chazelle_triangulate(NULL, 0, &tris, &num_tris);
        assert(status == CHAZELLE_ERROR_INVALID_ARGUMENT);
    }

    // Test 2: Polygon too small
    {
        ChazellePoint pts[2] = {{0.0, 0.0}, {1.0, 1.0}};
        ChazelleTriangle* tris = NULL;
        size_t num_tris = 0;
        ChazelleStatus status = chazelle_triangulate(pts, 2, &tris, &num_tris);
        assert(status == CHAZELLE_ERROR_POLYGON_TOO_SMALL);
    }

    // Test 3: Zero-Allocation in-place triangulation (chazelle_triangulate_into)
    {
        ChazellePoint pts[4] = {
            {0.0, 0.0},
            {2.0, 0.0},
            {2.0, 2.0},
            {0.0, 2.0}
        };
        // Pre-allocated stack buffer for triangles
        ChazelleTriangle stack_buffer[2];
        size_t num_tris = 0;

        // Test buffer too small error handling
        ChazelleStatus err_status = chazelle_triangulate_into(pts, 4, stack_buffer, 1, &num_tris);
        assert(err_status == CHAZELLE_ERROR_BUFFER_TOO_SMALL);

        // Test successful in-place zero-allocation triangulation
        ChazelleStatus status = chazelle_triangulate_into(pts, 4, stack_buffer, 2, &num_tris);
        assert(status == CHAZELLE_SUCCESS);
        assert(num_tris == 2);

        double total_area = 0.0;
        for (size_t i = 0; i < num_tris; ++i) {
            assert(stack_buffer[i].a < 4 && stack_buffer[i].b < 4 && stack_buffer[i].c < 4);
            total_area += triangle_area(pts[stack_buffer[i].a], pts[stack_buffer[i].b], pts[stack_buffer[i].c]);
        }
        assert(fabs(total_area - 4.0) < 1e-6);
    }

    // Test 4: Concave L-shape triangulation
    {
        ChazellePoint pts[6] = {
            {0.0, 0.0},
            {3.0, 0.0},
            {3.0, 1.0},
            {1.0, 1.0},
            {1.0, 3.0},
            {0.0, 3.0}
        };
        ChazelleTriangle* tris = NULL;
        size_t num_tris = 0;
        ChazelleStatus status = chazelle_triangulate(pts, 6, &tris, &num_tris);
        assert(status == CHAZELLE_SUCCESS);
        assert(num_tris == 4);

        double total_area = 0.0;
        for (size_t i = 0; i < num_tris; ++i) {
            assert(tris[i].a < 6 && tris[i].b < 6 && tris[i].c < 6);
            total_area += triangle_area(pts[tris[i].a], pts[tris[i].b], pts[tris[i].c]);
        }
        assert(fabs(total_area - 5.0) < 1e-6);

        chazelle_free_triangles(tris, num_tris);
    }

    // Test 5: Triangulation algorithm selection in C
    {
        ChazellePoint square[4] = {
            {0.0, 0.0}, {2.0, 0.0}, {2.0, 2.0}, {0.0, 2.0}
        };
        ChazelleTriangle* tris_sweep = NULL;
        size_t n_sweep = 0;
        ChazelleStatus s1 = chazelle_triangulate_with_algorithm(
            square, 4, &tris_sweep, &n_sweep, CHAZELLE_ALGORITHM_MONOTONE_SWEEP
        );
        assert(s1 == CHAZELLE_SUCCESS);
        assert(n_sweep == 2);
        chazelle_free_triangles(tris_sweep, n_sweep);

        ChazelleTriangle buf_seidel[2];
        size_t n_seidel = 0;
        ChazelleStatus s2 = chazelle_triangulate_into_with_algorithm(
            square, 4, buf_seidel, 2, &n_seidel, CHAZELLE_ALGORITHM_SEIDEL
        );
        assert(s2 == CHAZELLE_SUCCESS);
        assert(n_seidel == 2);
    }

    printf("All C binding tests passed successfully!\n");
    return 0;
}
