#include "cpplap.h"
#include "units.h"

using namespace cpplap;

int main()
{

    std::cout.precision(4);
    std::cout << std::fixed;
    Vect<double> vect(1, 0, 0);

    std::cout << "Simple vect: " << vect << std::endl;

    vect.rotate({0, 0, 1}, 90.0_deg);
    std::cout << "rotated: " << vect << std::endl;
    std::cout << "phi: " << vect.phi() << std::endl;
    std::cout << "theta: " << vect.theta() << std::endl;

    std::cout << "projection: " << vect.projectionTo({1, 1, 1}) << std::endl;
    std::cout << "cross product: " << vect.cross({0, 0, 1}) << std::endl;
    std::cout << "inner product: " << vect * Vect<double>{2, 2, 2} << std::endl;
    vect *= 5;
    std::cout << "vect * 5: " << vect << std::endl;
    vect.normalize();
    std::cout << "normalized: " << vect << std::endl;
    vect = Vect<double>::SphericalCoords(90_deg, 90_deg, 5);
    std::cout << "From spherical: " << vect << std::endl;

    Vect<double> r{1, 1, 1};
    Vect<double> v(2, 2, 2);

    std::cout << "\nSimple line x_vec = r_vec + t * v_vec" << std::endl;
    Line<double> line(r, v);
    std::cout << "line check valid: " << line.isInLine({0, 0, 0}) << std::endl;
    std::cout << "line check invalid: " << line.isInLine({1, 0, 0}) << std::endl;
    std::cout << "line at -0.5: " << line.at(-0.5) << std::endl;
    std::cout << "line at 2: " << line.at(1) << std::endl;
    std::cout << "line distance: " << line.distanceTo({1, 2, 3}) << std::endl;
    std::cout << "line scalar distance: " << line.scalarDistanceTo({1, 2, 3}) << std::endl;

    std::cout << "\nSimple plane: ax + by + cz = d" << std::endl;
    Plane<double> plane(1, 1, 1, 3);
    std::cout << "intersection: " << plane.intersectionPointWith(line) << std::endl;
    std::cout << "distance vect: " << plane.distanceTo({2, 2, 2}) << std::endl;
    std::cout << "distance scalar: " << plane.scalarDistanceTo({0, 0, 0}) << std::endl;
    std::cout << "check valid: " << plane.isInPlane({2, 1, 0}) << std::endl;
    std::cout << "check invalid: " << plane.isInPlane({2, 1, 1}) << std::endl;

    std::cout << "\nHesse plane: (r_vec - x_vec) * n_vec = 0" << std::endl;
    Vect<double>       n{1, 1, 1};
    HessePlane<double> hesse_plane(r, n);
    std::cout << "intersection: " << hesse_plane.intersectionPointWith(line) << std::endl;
    std::cout << "distance vect:" << hesse_plane.distanceTo({2, 2, 2}) << std::endl;
    std::cout << "distance scalar: " << hesse_plane.scalarDistanceTo({-2, -2, -2}) << std::endl;
    std::cout << "check valid: " << hesse_plane.isInPlane({2, 1, 0}) << std::endl;
    std::cout << "check invalid: " << hesse_plane.isInPlane({2, 1, 1}) << std::endl;

    std::cout << "\nParametric plane: x_vec = r_vec + s*u_vec + t*v_vec " << std::endl;
    Vect<double> u{1, -1, 0};
    v = {0, 1, -1};
    ParametricPlane<double> parametric_plane(r, u, v);
    std::cout << "intersection: " << parametric_plane.intersectionPointWith(line) << std::endl;
    std::cout << "distance vect: " << parametric_plane.distanceTo({2, 2, 2}) << std::endl;
    std::cout << "distance scalar: " << parametric_plane.scalarDistanceTo({-2, -2, -2}) << std::endl;
    std::cout << "check valid: " << parametric_plane.isInPlane({2, 1, 0}) << std::endl;
    std::cout << "check invalid: " << parametric_plane.isInPlane({2, 1, 1}) << std::endl;
    std::cout << "bound check valid: " << parametric_plane.checkBounds({2, 0, 1}) << std::endl;
    std::cout << "bound check invalid: " << parametric_plane.checkBounds({3, -1, 0}) << std::endl;

    return 0;
}