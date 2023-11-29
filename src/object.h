#pragma once
#include <math.h>
#include <iostream>
#include "vec.h"

constexpr double t_min = 1e-5;
constexpr double inf = INFINITY;

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
    const Point3 origin;
    const Vec3 normal;

    Plane(const Point3 &p, const Vec3 &v) : origin(p), normal(v) {}

    Plane(const Point3 &p, const Vec3 &v, const double t_rr[3], double t_ir)
            : Object(t_rr, t_ir), origin(p), normal(v) { }

    Intersection intersect(const Line &l) const override {
        if (abs(normal * l.dir) < 1e-5)
            return {inf, normal};

        double t = (normal * (origin - l.origin)) / (normal * l.dir);

        return {t > t_min ? t : inf, normal};
    }
};

struct Sphere : Object {
    const Point3 center;
    const double radius;

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

        double t = -b - delta > (2 * a) * t_min ? (-b - delta) / (2 * a)
                                                : (-b + delta) / (2 * a);

        if (t < 1e-5)
            return {inf, l.dir};

        return {t, l.t(t) - center};
    }
};

struct Triangle : Object {
    const Point3 vertices[3];

    Triangle(const Point3 &p1, const Point3 &p2, const Point3 &p3) : vertices{p1, p2, p3} { }
    Triangle(const Point3 &p1, const Point3 &p2, const Point3 &p3,
             const double t_rr[3], double t_ir) : Object(t_rr, t_ir), vertices{p1, p2, p3} { }

    Intersection intersect(const Line &l) const override {
        const Vec3    oc  = vertices[0] -   l.origin;
        const Vec3 edge1  = vertices[0] - vertices[1];
        const Vec3 edge2  = vertices[0] - vertices[2];
        const Vec3 normal = edge1 ^ edge2;

        constexpr auto signal = [](double d) { return (d > 0) - (d < 0); };
        
        const double holder = l.dir * normal;

        if (abs(holder) < 1e-5)
            return {inf, l.dir};

        const Vec3 vp = oc ^ l.dir;

        const double beta = edge1 * vp;
        if (signal(beta) * signal(holder) < 0)
            return {inf, l.dir};

        const double gamma = edge2 * -vp;
        if (signal(gamma) * signal(holder) < 0)
            return {inf, l.dir};

        if (holder > 0 ? beta + gamma > holder : beta + gamma < holder)
            return {inf, l.dir};

        const double t = (oc * normal) / holder;
        return {t > t_min ? t : inf, normal};
    }
};
