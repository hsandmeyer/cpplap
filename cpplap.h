#pragma once
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace cpplap {

template <class T> struct DefaultPrec {
    static constexpr T value = 1e-5;
};

template <> struct DefaultPrec<double> {
    static constexpr double value = 1e-12;
};

/*
    Equality check of floating point numbers within precision
*/
template <class T> bool cmp_rel(const T a, const T b, const T prec = DefaultPrec<T>::value)
{
    if ((fabs(a) < prec) && (fabs(b) < prec)) {
        return true;
    }
    return (fabs(a - b) / (fabs(a) + fabs(b))) < prec;
}

template <class T> class Matrix;

template <class T> class Vect;

/*
    Simple 3D vector class
*/
template <class T> class Vect {

    T _x0 = 0;
    T _x1 = 0;
    T _x2 = 0;

public:
    Vect() = default;

    Vect(const Vect &a) = default;

    Vect(const T x0, const T x1, const T x2) : _x0(x0), _x1(x1), _x2(x2) {}

    friend std::ostream &operator<<(std::ostream &stream, const Vect<T> &a)
    {
        stream << "[" << a._x0 << "," << a._x1 << "," << a._x2 << "]";
        return stream;
    }

    /* Exact equality */
    friend bool operator==(const Vect<T> &a, const Vect<T> &b)
    {
        return a._x0 == b._x0 && a._x1 == b._x1 && a._x2 == b._x2;
    }

    /**
    Relative equality within prec
    */
    friend bool cmp_rel(const Vect<T> a, const Vect<T> b, const T prec = DefaultPrec<T>::value)
    {
        return cmp_rel(a[0], b[0], prec) && cmp_rel(a[1], b[1], prec) && cmp_rel(a[2], b[2], prec);
    }

    /**
        Inner product
    */
    friend T operator*(const Vect<T> &a, const Vect<T> &b) { return a._x0 * b._x0 + a._x1 * b._x1 + a._x2 * b._x2; }

    /**
        Scalar times vector
    */
    friend Vect<T> operator*(const T a, const Vect<T> &b) { return Vect<T>(a * b._x0, a * b._x1, a * b._x2); }

    /**
        Vector times scalar
    */
    const Vect operator*(const T a) { return Vect(_x0 * a, _x1 * a, _x2 * a); }

    /**
        Vector divided by scalar
    */
    const Vect operator/(const T a) { return Vect(_x0 / a, _x1 / a, _x2 / a); }

    Vect &operator*=(const T a)
    {
        _x0 *= a;
        _x1 *= a;
        _x2 *= a;
        return *this;
    }

    Vect &operator/=(const T a)
    {
        _x0 /= a;
        _x1 /= a;
        _x2 /= a;
        return *this;
    }

    Vect &operator+=(const Vect<T> a)
    {
        _x0 += a._x0;
        _x1 += a._x1;
        _x2 += a._x2;
        return *this;
    }

    Vect &operator-=(const Vect<T> a)
    {
        _x0 -= a._x0;
        _x1 -= a._x1;
        _x2 -= a._x2;
        return *this;
    }

    /**
        Add two vectors
    */
    friend Vect<T> operator+(const Vect<T> &a, const Vect<T> &b)
    {
        return Vect<T>(a._x0 + b._x0, a._x1 + b._x1, a._x2 + b._x2);
    }

    /**
        Subtract two vectors
    */
    friend Vect<T> operator-(const Vect<T> &a, const Vect<T> &b)
    {
        return Vect<T>(a._x0 - b._x0, a._x1 - b._x1, a._x2 - b._x2);
    }

    /**
    Cross product of two vectors
    */
    friend Vect<T> cross_prod(const Vect<T> &a, const Vect<T> &b)
    {
        return Vect<T>(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);
    }

    /**
    Cross product with another vector
    */
    Vect<T> cross(const Vect<T> &a) { return cross_prod(*this, a); }

    /**
    Factory pattern to create a vector from spherical coordinates.
    @param phi angle in xy-plane
    @param theta angle from z-axis
    @param r length of the vector
    */
    static Vect SphericalCoords(const T phi, const T theta, const T r)
    {
        return Vect(r * sin(theta) * cos(phi), r * sin(theta) * sin(phi), r * cos(theta));
    }

    /**
    Factory pattern to create a vector from spherical coordinates along the x-axis
    @param phi angle in yz-plane
    @param theta angle from x-axis
    @param r length of the vector
    */
    static Vect SphericalCoordsAlongX(const T phi, const T theta, const T r)
    {
        return Vect(r * cos(theta), r * sin(theta) * cos(phi), r * sin(theta) * sin(phi));
    }

    /**
    Factory pattern to create a vector from rather unusual coordinates
    @param phi angle in the xz-plane
    @param theta angle in the yz--plane
    @param r length of the vector
    */
    static Vect xyAngularCoords(const T phi, const T theta, const T r = 1)
    {
        Vect ret;
        T    tp = tan(phi);
        T    tt = tan(theta);

        ret._x2 = sqrt(1.0 / (tp * tp + tt * tt + 1));
        ret._x0 = ret._x2 * tp;
        ret._x1 = ret._x2 * tt;
        return r * ret;
    }

    /**
    Get the i-th element.
    */
    T &operator[](const int i)
    {
        switch (i) {
        case 0:
            return _x0;
            break;
        case 1:
            return _x1;
            break;
        case 2:
            return _x2;
            break;
        default:
            throw std::out_of_range("Vect: Index out of range");
        }
    }

    /** Const version of above without reimplementation */
    const T &operator[](const int i) const { return const_cast<T &>((*const_cast<Vect *>(this))[i]); }

    /**
    Project the vector to another vector
    */
    Vect projectionTo(Vect vec) const { return ((*this) * vec) * vec / vec.norm2(); }

    /**
    Rotate the vector around a given axis using given angle
    */
    Vect rotate(Vect axis, T angle)
    {
        axis.normalize();
        *this = axis * (axis * (*this)) + cos(angle) * cross_prod(cross_prod(axis, *this), axis) +
                sin(angle) * cross_prod(axis, *this);
        return *this;
    }

    /**Normalize the vector in place */
    Vect &normalize()
    {
        (*this) = (*this) / norm();
        return *this;
    }

    /** Get the squared Euclidean norm of the vector */
    T norm2() const { return (*this) * (*this); }

    /** Get the Euclidean norm of the vector */
    T norm() const { return sqrt((*this) * (*this)); }

    /** Calculate the angle phi for spherical coordinates*/
    T phi() const { return std::atan2(_x1, _x0); }

    /** Calculate the angle theta for spherical coordinates*/
    T theta() const { return std::acos(_x2 / norm()); }
};

template <class T> class Matrix {

    T _e00;
    T _e01;
    T _e02;
    T _e10;
    T _e11;
    T _e12;
    T _e20;
    T _e21;
    T _e22;

public:
    Matrix(T e00, T e01, T e02, T e10, T e11, T e12, T e20, T e21, T e22)
        : _e00(e00), _e01(e01), _e02(e02), _e10(e10), _e11(e11), _e12(e12), _e20(e20), _e21(e21), _e22(e22){};

    friend Vect<T> operator*(const Matrix<T> &m, const Vect<T> &a)
    {

        T tmp0 = m._e00 * a[0] + m._e01 * a[1] + m._e02 * a[2];
        T tmp1 = m._e10 * a[0] + m._e11 * a[1] + m._e12 * a[2];
        T tmp2 = m._e20 * a[0] + m._e21 * a[1] + m._e22 * a[2];

        return Vect<T>(tmp0, tmp1, tmp2);
    }
};

template <class T = float> class Line {
    /**
    Simple class to describe lines in 3D
    Notation _r + t*_v
    */

    // support vector
    Vect<T> _r;

    // direction vector
    Vect<T> _v;

public:
    Line(const Vect<T> r, const Vect<T> v) : _r(r), _v(v) {}

    Line(const Line &line) = default;

    Vect<T> &getR() { return _r; }

    Vect<T> &getV() { return _v; }

    const Vect<T> &getR() const { return _r; }

    const Vect<T> &getV() const { return _v; }

    /**
    Calculate the resulting vector for a given t
    */
    Vect<T> at(const T t) const { return _r + t * _v; }

    /**
    Check wether a point lies within a the line
    @param vect The vector to be checked
    @param prec The relative precision to account for floating point rounding errors
    */
    bool isInLine(const Vect<T> &vect, T prec = DefaultPrec<T>::value)
    {
        return cross_prod((vect - _r), _v).norm2() < prec;
    }

    /**
    Distance between support vector and another vector. The direction is used to determine the sign of the returned
    distance.
     * If the distance vector points into the opposite direction of the direction vector, we
     * return a negative value.
     */
    T dist(const Vect<T> &vect) const
    {
        Vect<T> dist_vec = (vect - _r);
        T       dist     = dist_vec.norm();
        // Scalar product of vector has negative sign if vectors point in opposite direction.
        // Use scalar product to get the sign of dist
        return copysign(dist, dist_vec * _v);
    }
};

template <class T> class HessePlane {
    /**
    Simple class to describe planes using the hesse normal form
    Notation (x - _r) * _n = 0
    */

    Vect<T> _r;
    Vect<T> _n;

public:
    HessePlane(const Vect<T> r, const Vect<T> n) : _r(r), _n(n) {}

    Vect<T> &getR() { return _r; }

    Vect<T> &getN() { return _n; }

    const Vect<T> &getR() const { return _r; }

    const Vect<T> &getN() const { return _n; }

    friend std::ostream &operator<<(std::ostream &stream, const HessePlane<T> &a)
    {
        stream << "r =\n" << a._r << "n =\n " << a._n;
        return stream;
    }

    /**
    Check whether a vector vec lies within this plane
    */
    bool isInPlane(const Vect<T> vec, const T prec = DefaultPrec<T>::value) { return fabs((vec - _r) * _n) < prec; }

    /**
    Calculate the intersection point of a line and this plane. When the line is parallel the resulting vector will be
    filled with nan or inf
    */
    Vect<T> intersectionPointWith(const Line<T> &line) const
    {
        T t = (getR() * getN() - line.getR() * getN()) / (line.getV() * getN());
        return line.at(t);
    }

    /**
    Calculate the vector of the shortest distance from a given point to this plane
    */
    Vect<T> distanceTo(const Vect<T> &point)
    {
        Line<T> helper_line(point, getN());
        return point - intersectionPointWith(helper_line);
    }

    /**
    Calculate the shortest distance from a given point to this plane. The sign is negative, if the normal vector of the
    plane points to the opposite direction with respect to the given point
    */
    T scalarDistanceTo(const Vect<T> &point)
    {
        Vect<T> dist_vec = distanceTo(point);
        T       dist     = dist_vec.norm();

        // Scalar product of vector has negativ sign if vectors point in opposite direction.
        // Use scalar product to get the sign of dist
        return copysign(dist, dist_vec * getN());
    }
};

template <class T> class ParametricPlane {
    /**
    Simple class to describe a plane using the parametric form
    Notation r + s*u + t*v
    */
    Vect<T> _r;
    Vect<T> _u;
    Vect<T> _v;

public:
    friend std::ostream &operator<<(std::ostream &stream, const ParametricPlane<T> &a)
    {
        stream << "r =\n" << a._r << "u =\n " << a._u << "v =\n " << a._v;
        return stream;
    }

    ParametricPlane(const Vect<T> r, const Vect<T> u, const Vect<T> v) : _r(r), _u(u), _v(v) {}

    Vect<T> &getR() const { return _r; }

    Vect<T> &getU() const { return _u; }

    Vect<T> &getV() const { return _u; }

    /**
    Convert to plane in Hesse normal form
    */
    operator HessePlane<T>() const { return HessePlane<T>(_r, cross_prod(_v, _u)); }

    /**
    Check whether a vector lies within this plane
    */
    bool isInPlane(const Vect<T> vec, const T prec = DefaultPrec<T>::value) const
    {
        HessePlane<T> helper_plane = *this;
        return helper_plane.isInPlane(vec, prec);
    }

    /**
    Check whether a vector lies within this plane and also lies within the parallelogram spanned by _u and _v
    */
    bool checkBounds(const Vect<T> p, const T prec = DefaultPrec<T>::value) const
    {
        Vect<T> diff = _r - p;
        if (!isInPlane(p, prec))
            return false;
        return (diff.projectionTo(_u).norm() <= _u.norm()) & (diff.projectionTo(_v).norm() <= _v.norm());
    }

    /**
    Calculate the intersection point of a line and this plane. When the line is parallel the resulting vector will be
    filled with nan or inf
    */
    Vect<T> intersectionPointWith(const Line<T> &line) const
    {
        HessePlane<T> helper = *this;
        return helper.intersectionPointWith(line);
    }

    /**
    Calculate the vector of the shortest distance from a given point to this plane
    */
    Vect<T> distanceTo(const Vect<T> &point)
    {
        HessePlane<T> helper_plane = *this;
        return helper_plane.distanceTo(point);
    }

    /**
    Calculate the shortest distance from a given point to this plane. The sign is negative, if the normal vector of the
    plane points to the opposite direction with respect to the given point
    */
    T scalarDistanceTo(const Vect<T> &point)
    {
        HessePlane<T> helper_plane = *this;
        return helper_plane.scalarDistanceTo(point);
    }
};

template <class T> class Plane {
    /* Simple class to describe planes using the coordinate form
     Notation a * x + b * y + c * z = d;
     */
    T _a;
    T _b;
    T _c;
    T _d;

public:
    Plane(T a, T b, T c, T d) : _a(a), _b(b), _c(c), _d(d) {}

    /**
    Convert to plane in Hesse normal form
    */
    operator HessePlane<T>() const
    {

        Vect<T> n(_a, _b, _c);
        Vect<T> r;

        if (_a != 0) {
            r[0] = _d / _a;
            r[1] = 0;
            r[2] = 0;
        }
        else if (_b != 0) {
            r[0] = 0;
            r[1] = _d / _b;
            r[2] = 0;
        }
        else {
            r[0] = 0;
            r[1] = 0;
            r[2] = _d / _c;
        }

        return HessePlane<T>(r, n);
    }

    /**
    Check whether a vector lies within this plane
    */
    bool isInPlane(const Vect<T> vec, const T prec = DefaultPrec<T>::value) const
    {
        HessePlane<T> helper_plane = *this;
        return helper_plane.isInPlane(vec, prec);
    }

    /**
    Calculate the intersection point of a line and this plane. When the line is parallel the resulting vector will be
    filled with nan or inf
    */
    Vect<T> intersectionPointWith(const Line<T> &line) const
    {
        HessePlane<T> helper = *this;
        return helper.intersectionPointWith(line);
    }

    /**
    Calculate the vector of the shortest distance from a given point to this plane
    */
    Vect<T> distanceTo(const Vect<T> &point)
    {
        HessePlane<T> helper_plane = *this;
        return helper_plane.distanceTo(point);
    }

    /**
    Calculate the shortest distance from a given point to this plane. The sign is negative, if the normal vector of the
    plane points to the opposite direction with respect to the given point
    */
    T scalarDistanceTo(const Vect<T> &point)
    {
        HessePlane<T> helper_plane = *this;
        return helper_plane.scalarDistanceTo(point);
    }
};

}; // namespace cpplap
