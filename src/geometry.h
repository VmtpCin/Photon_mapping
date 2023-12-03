#pragma once
#include <math.h>
#include <iostream>
#include <vector>
#include "vec.h"

constexpr double t_min = 1e-5;
constexpr double inf = INFINITY;

struct Intersection {
    double t = inf;
    Vec3 normal;

    operator bool() const {
        return t != inf;
    }

    bool operator<(const Intersection &i) const {
        return t < i.t;
    }
};

struct Geometry {
    virtual Intersection intersect(const Line &l) const { return {inf, l.dir}; };
};

struct Object : Geometry {
    Geometry *geometry;
    double rr[3] = {0, 0, 1}, ir = 1;
    Color color = 1;

    ~Object() { delete geometry; }

    Object(Geometry *geo) : geometry(geo) {}
    Object(Geometry *geo, const double t_rr[3], double t_ir)
            : geometry(geo), rr{t_rr[0], t_rr[1], t_rr[2]}, ir(t_ir) {}
    Object(Geometry *geo, const double t_rr[3], double t_ir, Color c)
            : geometry(geo), rr{t_rr[0], t_rr[1], t_rr[2]}, ir(t_ir), color(c) {}

    Intersection intersect(const Line &l) const override { return geometry->intersect(l); };
};

struct Plane : Geometry {
    const Point3 origin;
    const Vec3 normal;

    Plane(const Point3 &p, const Vec3 &v) : origin(p), normal(v) {}

    Intersection intersect(const Line &l) const override;
};

struct Sphere : Geometry {
    const Point3 center;
    const double radius;

    Sphere(const Point3 &p, double r) : center(p), radius(r) {}

    Intersection intersect(const Line &l) const override;
};

struct Triangle : Geometry {
    const Point3 origin;
    const Vec3 edge1, edge2;

    Triangle(const Point3 &p1, const Point3 &p2, const Point3 &p3)
            : origin(p1), edge1(p2 - p1), edge2(p3 - p1) { }

    Intersection intersect(const Line &l) const override;
};

struct Parallelogram : Geometry {
    const Point3 origin;
    const Vec3 edge1, edge2;

    Parallelogram(const Point3 &p1, const Point3 &p2, const Point3 &p3)
                 : origin(p1), edge1(p2 - p1), edge2(p3 - p1) { }

    Intersection intersect(const Line &l) const override;
};
