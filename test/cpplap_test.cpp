#include "cpplap.h"
#include "test.h"
#include <iostream>

using namespace cpplap;

int main()
{
    Vect<float> a(1, 2, 3);
    Vect<float> b(3, 2, 1);
    Vect<float> c = a;

    bool out = true;

    out &= test_exact(a == c, true, "Equality true");
    out &= test_exact(a == b, false, "Equality false");
    out &= test_exact(a * b, 10.0f, "Scalar product");
    out &= test_exact(a * 2, Vect<float>(2, 4, 6), "Vector times scalar");
    out &= test_exact(2 * a, Vect<float>(2, 4, 6), "Scalar times vector");
    out &= test_exact(a + b, Vect<float>(4, 4, 4), "Vector plus vector");
    out &= test_exact(a.norm(), sqrtf(1 * 1 + 2 * 2 + 3 * 3), "First norm");
    out &= test_rel(Vect<float>(1, 1, 1).normalize(), Vect<float>(1 / sqrtf(3), 1 / sqrtf(3), 1 / sqrtf(3)),
                    "Second norm");
    out &=
        test_exact(cross_prod(Vect<float>(0, 1, 2), Vect<float>(0, 0, 1)), Vect<float>(1, 0, 0), "First cross product");
    out &= test_exact(cross_prod(a, b), Vect<float>(-4, 8, -4), "Second cross product");

    Line<float> line(Vect<float>(0, 1, 0), Vect<float>(0, -1, 2));

    out &= test_exact(line.isInLine(Vect<float>(0, 1, 0)), true, "First is in line");
    out &= test_exact(line.isInLine(Vect<float>(0, 0, 2)), true, "Second is in line");
    out &= test_exact(line.isInLine(Vect<float>(0, 1, 2)), false, "Third is in line");

    out &= test_rel(line.distFromRTo(Vect<float>(0, 0, 2)), sqrtf(5.f), "First origin distance");
    out &= test_rel(line.distFromRTo(Vect<float>(0, 2, -2)), -sqrtf(5.f), "Second origin distance");
    out &= test_rel(line.nearestTo(Vect<float>(4, 0, 2)), {0, 0, 2}, "Nearest to");
    out &= test_rel(line.distanceTo(Vect<float>(4, 0, 2)), {4, 0, 0}, "First line point distance");
    out &= test_rel(line.scalarDistanceTo(Vect<float>(-4, -2, 6)), 4.f, "Second line point distance");

    HessePlane<float> hess_plane(Vect<float>(6, 2, 1), Vect<float>(1, 3, -2));

    out &= test_exact(hess_plane.intersectionPointWith(line), Vect<float>(0, 2, -2), "First plane line intersection");

    line = Line<float>(Vect<float>(2, -3, 2), Vect<float>(1, -1, 3));

    ParametricPlane<float> param_plane(Vect<float>(-3, 1, 1), Vect<float>(1, -2, -1), Vect<float>(0, -1, 2));

    out &=
        test_exact(param_plane.intersectionPointWith(line), Vect<float>(-1, 0, -7), "Second plane line intersection");

    hess_plane = HessePlane<float>(Vect<float>(0, 0, 0), Vect<float>(1, 1, 1));

    out &= test_exact(hess_plane.distanceTo(Vect<float>(1, 1, 1)), Vect<float>(1, 1, 1), "First point plane distance");

    out &= test_rel(hess_plane.scalarDistanceTo(Vect<float>(1, 1, 1)), sqrtf(3), "Second point plane distance");

    out &= test_rel(hess_plane.scalarDistanceTo(Vect<float>(-1, -1, -1)), -sqrtf(3), "Third point plane distance");

    out &= test_rel(Vect<float>(6, -6, 7).projectionTo(Vect<float>(-3, 2, -6)), -72. / 49 * Vect<float>(-3, 2, -6),
                    "Projection");

    param_plane = ParametricPlane<float>(Vect<float>(0, 0, 0), Vect<float>(0, 0, 1), Vect<float>(0, 1, 0));

    Plane<float> plane(1, 0, 0, 1);

    out &= test_rel(plane.scalarDistanceTo(Vect<float>(0, -1, -1)), -1.0f, "Fourth plane distance");

    plane = Plane<float>(0, 1, 0, 5);
    out &= test_rel(plane.scalarDistanceTo(Vect<float>(100, 0, -1)), -5.0f, "Fifth plane distance");

    plane = Plane<float>(0, 0, 1, -10);
    out &= test_rel(plane.scalarDistanceTo(Vect<float>(0.123, 123, 0)), 10.0f, "Sixth plane distance");

    plane = Plane<float>(1, 1, 1, 0);
    out &= test_rel(plane.scalarDistanceTo(Vect<float>(2, 2, 2)), sqrtf(12), "Seventh plane distance");

    line = Line<float>(Vect<float>(1, -1, 1), Vect<float>(1, -1, 1));

    out &= test_rel(plane.intersectionPointWith(line), Vect<float>(0, 0, 0), "Third plane line intersection");

    out &= test_rel(plane.intersectionPointWith(line), Vect<float>(0, 0, 0), "Third plane line intersection");

    out &= test_exact(plane.isInPlane(Vect<float>(1, -0.5, -0.5)), true, "First is in plane check");
    out &= test_exact(plane.isInPlane(Vect<float>(1, -0.5, 0.5)), false, "Second is in plane check");

    out &= test_exact(param_plane.isInPlane(Vect<float>(10, 1, 1)), false, "Third is in plane check");

    out &= test_exact(param_plane.checkBounds(Vect<float>(10, 1, 1)), false, "First boundary check");
    out &= test_exact(param_plane.checkBounds(Vect<float>(0, 1, 1)), true, "Second boundary check");
    out &= test_exact(param_plane.checkBounds(Vect<float>(0, 1.1, 1)), false, "Third boundary check");

    out &= test_exact(Matrix<float>(1, 2, 3, 4, 5, 4, 7, 8, 9) * Vect<float>(1, 3, 5), Vect<float>(22, 39, 76),
                      "Matrix times Vector");

    out &= test_rel(Vect<float>(1, 1, 1).rotate(Vect<float>(1, 1, 1), M_PI), Vect<float>(1, 1, 1), "First rotation");
    out &=
        test_rel(Vect<float>(1, 1, 1).rotate(Vect<float>(1, 0, 0), M_PI / 2), Vect<float>(1, -1, 1), "Second rotation");
    out &= test_rel(Vect<float>(1, -1, 0).rotate(Vect<float>(1, 1, 0), M_PI), Vect<float>(-1, 1, 0), "Third rotation");

    out &= test_rel(Vect<float>::SphericalCoords(M_PI, M_PI / 2, 2), Vect<float>(-2, 0, 0), "Spherical coordinates");

    out &= test_rel(Vect<float>::SphericalCoordsAlongX(M_PI / 2, M_PI / 2, 1), Vect<float>(0, 0, 1),
                    "Pherical coords along X");

    out &= test_rel(Vect<float>(-2, 0, 0).phi(), static_cast<float>(M_PI), "Calculate phi");

    out &= test_rel(Vect<float>(-2, 0, 0).theta(), static_cast<float>(M_PI / 2), "Calculate theta");

    if (!out) {
        return -1;
    }

    return 0;
}
