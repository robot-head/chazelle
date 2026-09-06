/**
 * @file chazelle.hpp
 * @brief C++ modern wrapper for Chazelle's linear-time polygon triangulation algorithm.
 * Optimized with google/zerocopy for zero-copy memory manipulation and zero-allocation calls.
 */

#ifndef CHAZELLE_HPP
#define CHAZELLE_HPP

#include "chazelle.h"
#include <vector>
#include <array>
#include <stdexcept>
#include <string>

namespace chazelle {

struct Point {
    double x;
    double y;

    Point() : x(0.0), y(0.0) {}
    Point(double x_, double y_) : x(x_), y(y_) {}
};

struct Triangle {
    size_t a;
    size_t b;
    size_t c;

    Triangle() : a(0), b(0), c(0) {}
    Triangle(size_t a_, size_t b_, size_t c_) : a(a_), b(b_), c(c_) {}

    size_t operator[](size_t idx) const {
        if (idx == 0) return a;
        if (idx == 1) return b;
        return c;
    }
};

class TriangulationError : public std::runtime_error {
public:
    ChazelleStatus code;
    explicit TriangulationError(ChazelleStatus c, const std::string& msg)
        : std::runtime_error(msg), code(c) {}
};

/**
 * @brief Zero-allocation triangulation writing directly into caller's pre-allocated buffer.
 *
 * @param points Pointer to points array.
 * @param num_points Number of points.
 * @param out_triangles Destination array with capacity >= num_points - 2.
 * @param out_capacity Capacity of out_triangles.
 * @return size_t Number of triangles written.
 */
inline size_t triangulate_into(
    const Point* points,
    size_t num_points,
    Triangle* out_triangles,
    size_t out_capacity
) {
    if (num_points < 3) {
        throw TriangulationError(CHAZELLE_ERROR_POLYGON_TOO_SMALL, "Polygon must have at least 3 vertices");
    }

    size_t num_written = 0;
    static_assert(sizeof(Point) == sizeof(ChazellePoint), "Point layout mismatch");
    static_assert(sizeof(Triangle) == sizeof(ChazelleTriangle), "Triangle layout mismatch");

    ChazelleStatus status = chazelle_triangulate_into(
        reinterpret_cast<const ChazellePoint*>(points),
        num_points,
        reinterpret_cast<ChazelleTriangle*>(out_triangles),
        out_capacity,
        &num_written
    );

    if (status != CHAZELLE_SUCCESS) {
        throw TriangulationError(status, "Triangulation failed with code " + std::to_string(static_cast<int>(status)));
    }

    return num_written;
}

/**
 * @brief Triangulate a simple polygon in linear time using Chazelle's algorithm.
 *
 * @param vertices A list of 2D points representing the simple polygon.
 * @return std::vector<Triangle> List of (N - 2) triangles.
 */
inline std::vector<Triangle> triangulate(const std::vector<Point>& vertices) {
    if (vertices.size() < 3) {
        throw TriangulationError(CHAZELLE_ERROR_POLYGON_TOO_SMALL, "Polygon must have at least 3 vertices");
    }

    size_t expected_triangles = vertices.size() - 2;
    std::vector<Triangle> result(expected_triangles);

    size_t actual_triangles = triangulate_into(
        vertices.data(),
        vertices.size(),
        result.data(),
        result.size()
    );

    result.resize(actual_triangles);
    return result;
}

inline const char* version() {
    return chazelle_version();
}

} // namespace chazelle

#endif /* CHAZELLE_HPP */
