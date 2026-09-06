/**
 * @file chazelle.hpp
 * @brief C++ modern wrapper for Chazelle's linear-time polygon triangulation algorithm.
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
 * @brief Triangulate a simple polygon in linear time using Chazelle's algorithm.
 *
 * @param vertices A list of 2D points representing the simple polygon.
 * @return std::vector<Triangle> List of (N - 2) triangles.
 * @throws TriangulationError on error.
 */
inline std::vector<Triangle> triangulate(const std::vector<Point>& vertices) {
    if (vertices.size() < 3) {
        throw TriangulationError(CHAZELLE_ERROR_POLYGON_TOO_SMALL, "Polygon must have at least 3 vertices");
    }

    ChazelleTriangle* raw_triangles = nullptr;
    size_t num_triangles = 0;

    static_assert(sizeof(Point) == sizeof(ChazellePoint), "Point layout mismatch");
    ChazelleStatus status = chazelle_triangulate(
        reinterpret_cast<const ChazellePoint*>(vertices.data()),
        vertices.size(),
        &raw_triangles,
        &num_triangles
    );

    if (status != CHAZELLE_SUCCESS) {
        throw TriangulationError(status, "Triangulation failed with code " + std::to_string(static_cast<int>(status)));
    }

    std::vector<Triangle> result;
    result.reserve(num_triangles);
    for (size_t i = 0; i < num_triangles; ++i) {
        result.emplace_back(raw_triangles[i].a, raw_triangles[i].b, raw_triangles[i].c);
    }

    chazelle_free_triangles(raw_triangles, num_triangles);
    return result;
}

inline const char* version() {
    return chazelle_version();
}

} // namespace chazelle

#endif /* CHAZELLE_HPP */
