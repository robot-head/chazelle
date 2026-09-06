#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>
#include "include/chazelle.hpp"

static double triangle_area(const chazelle::Point& a, const chazelle::Point& b, const chazelle::Point& c) {
    return 0.5 * std::abs((b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x));
}

int main() {
    std::cout << "Chazelle C++ wrapper version: " << chazelle::version() << std::endl;

    // Test 1: Exception on small polygon
    try {
        std::vector<chazelle::Point> tiny = {{0.0, 0.0}, {1.0, 1.0}};
        chazelle::triangulate(tiny);
        assert(false && "Should have thrown TriangulationError");
    } catch (const chazelle::TriangulationError& e) {
        assert(e.code == CHAZELLE_ERROR_POLYGON_TOO_SMALL);
    }

    // Test 2: Concave star triangulation
    std::vector<chazelle::Point> star;
    for (int i = 0; i < 10; ++i) {
        double r = (i % 2 == 0) ? 2.0 : 0.8;
        double angle = i * M_PI / 5.0;
        star.emplace_back(r * std::cos(angle), r * std::sin(angle));
    }

    auto triangles = chazelle::triangulate(star);
    assert(triangles.size() == star.size() - 2);

    double total_tri_area = 0.0;
    for (const auto& t : triangles) {
        assert(t.a < star.size() && t.b < star.size() && t.c < star.size());
        total_tri_area += triangle_area(star[t.a], star[t.b], star[t.c]);
    }
    assert(total_tri_area > 0.0);

    std::cout << "All C++ binding tests passed successfully!" << std::endl;
    return 0;
}
