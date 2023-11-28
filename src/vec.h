#pragma once
#include <math.h>

struct Vec3 {
    double e[3];

    double& x()       { return e[0]; }
    double  x() const { return e[0]; }
    double& y()       { return e[1]; }
    double  y() const { return e[1]; }
    double& z()       { return e[2]; }
    double  z() const { return e[2]; }

    double& operator[](int i)       { return e[i]; }
    double  operator[](int i) const { return e[i]; }

    Vec3 operator-() const {
        return {-e[0], -e[1], -e[2]};
    }

    Vec3 operator+(const Vec3 &v) const {
        return {e[0] + v[0], e[1] + v[1], e[2] + v[2]};
    }

    Vec3 operator-(const Vec3 &v) const {
        return {e[0] - v[0], e[1] - v[1], e[2] - v[2]};
    }

    Vec3 operator*(const double d) const {
        return {e[0] * d, e[1] * d, e[2] * d};
    }

    Vec3 operator/(const double d) const {
        return {e[0] / d, e[1] / d, e[2] / d};
    }

    Vec3& operator+=(const Vec3 &v) {
        e[0] += v[0], e[1] += v[1], e[2] += v[2];
        return *this;
    }

    Vec3& operator-=(const Vec3 &v) {
        e[0] -= v[0], e[1] -= v[1], e[2] -= v[2];
        return *this;
    }

    Vec3& operator*=(const double &d) {
        e[0] *= d, e[1] *= d, e[2] *= d;
        return *this;
    }

    Vec3& operator/=(const double &d) {
        e[0] /= d, e[1] /= d, e[2] /= d;
        return *this;
    }

    double operator*(const Vec3 &v) const {
        return e[0] * v[0] + e[1] * v[1] + e[2] * v[2];
    }

    Vec3 operator^(const Vec3 &v) const {
        return { e[1] * v[2] - e[2] * v[1],
                 e[2] * v[0] - e[0] * v[2],
                 e[0] * v[1] - e[1] * v[0]};
    }

    Vec3 operator&(const Vec3 &v) const {
        return {e[0] * v[0], e[1] * v[1], e[2] * v[2]};
    }

    Vec3 operator&=(const Vec3 &v) {
        e[0] *= v[0], e[1] *= v[1], e[2] *= v[2];
        return *this;
    }

    double length_sq() const {
        return e[0] * e[0] + e[1] * e[1] + e[2] * e[2];
    }

    double length() const {
        return sqrt(length_sq());
    }

    Vec3 normalize() const {
        return *this / length();
    }

    Vec3 perpendicular() const {
        if (e[0] != 0)
            return {-e[1], e[0], 0};
        else if (e[1] != 0)
            return {0, -e[2], e[1]};
        else if (e[2] != 0)
            return {e[2], 0, -e[0]};
        else
            return {0, 0, 0};
    }
};

struct Point3 {
    double e[3];

    explicit operator Vec3() const {
        return {e[0], e[1], e[2]};
    }

    double& x()       { return e[0]; }
    double  x() const { return e[0]; }
    double& y()       { return e[1]; }
    double  y() const { return e[1]; }
    double& z()       { return e[2]; }
    double  z() const { return e[2]; }

    double& operator[](int i)       { return e[i]; }
    double  operator[](int i) const { return e[i]; }

    Point3 operator+(const Vec3 &v) const {
        return {e[0] + v[0], e[1] + v[1], e[2] + v[2]}; 
    }

    Point3 operator+=(const Vec3 &v) {
        e[0] += v[0], e[1] += v[1], e[2] += v[2];
        return *this;
    }

    Vec3 operator-(const Point3 &p) const {
        return {e[0] - p[0], e[1] - p[1], e[2] - p[2]};
    }

    double distance_sq(const Point3 &p) const {
        double dx = e[0] - p[0], dy = e[1] - p[1], dz = e[2] - p[2];
        return dx * dx + dy * dy + dz * dz;
    }

    double distance(const Point3 &v) const {
        return sqrt(distance_sq(v));
    }

    template<int elem>
    bool is_less(const Point3 &p) const {
        return e[elem] < p[elem];
    }
};

struct Line {
    Point3 origin;
    Vec3 dir;

    Point3 t(double t) const {
        return origin + (dir * t);
    }
};

inline Point3 operator+(const Vec3& v, const Point3 &p) {
    return p + v;
}

inline Vec3 operator*(const double &d, const Vec3 &v) {
    return v * d;
}

