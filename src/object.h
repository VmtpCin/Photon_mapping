#pragma once
#include <math.h>
#include "vec.h"

constexpr double inf = 1e100;

struct Intersection {
    double t;
    Vec3 normal;

    operator bool() const {
        return t != inf;
    }

    bool operator<(const Intersection &i) const {
        return t < i.t;
    }
};

struct Object {
    double rr[3], ir;

    Object() : rr{0, 0, 1}, ir(1) {};
    Object(const double t_rr[3], double t_ir) : rr{t_rr[0], t_rr[1], t_rr[2]}, ir(t_ir) {}

    virtual Intersection intersect(const Line &l) const {
        return {inf, l.dir};
    };
};

struct Plane : Object {
    Point3 origin;
    Vec3 normal;

    Plane(const Point3 &p, const Vec3 &v) : origin(p), normal(v) {}

    Plane(const Point3 &p, const Vec3 &v, const double t_rr[3], double t_ir)
            : Object(t_rr, t_ir), origin(p), normal(v) { }

    Intersection intersect(const Line &l) const override {
        if (abs(normal * l.dir) < 1e-5)
            return {inf, normal};

        int t = (normal * (origin - l.origin)) / (normal * l.dir);

        return {t > 0 ? t : inf, normal};
    }
};

struct Sphere : Object {
    Point3 center;
    double radius;

    Sphere(const Point3 &p, double r) : center(p), radius(r) {}

    Sphere(const Point3 &p, double r, const double t_rr[3], double t_ir)
            : Object(t_rr, t_ir), center(p), radius(r) { }

    Intersection intersect(const Line &l) const override {
        Vec3 oc = l.origin - center;
        double a = l.dir * l.dir;
        double b = 2 * (oc * l.dir);
        double c = oc * oc - radius * radius;
        
        double delta = b * b - 4 * a * c;

        if (delta < 0)
            return {inf, l.dir};
        
        delta = sqrt(delta);

        double t = -b - delta > 0 ? (-b - delta) / (2 * a)
                                  : (-b + delta) / (2 * a);

        if (t < 1e-5)
            return {inf, l.dir};

        return {t, l.t(t) - center};
    }
};

