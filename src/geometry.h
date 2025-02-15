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

    [[nodiscard]] operator bool() const {
        return t != inf;
    }

    bool operator<(const Intersection &i) const {
        return t < i.t;
    }
};

struct Geometry {
    virtual ~Geometry() = default;

    virtual Intersection intersect(const Line &l) const { return {inf, l.dir}; };
};

struct Object : Geometry {
    Geometry *geometry;
    double rr[3] = {0, 0, 0};
    double (*ir)(double);
    Color color = 1.0;

    // ~Object() { if (geometry) delete geometry; }

    Object(Geometry *geo) : geometry(geo) {}
    Object(Geometry *geo, const double t_rr[3])
            : geometry(geo), rr{t_rr[0], t_rr[1], t_rr[2]} {}
    Object(Geometry *geo, const double t_rr[3], Color c)
            : geometry(geo), rr{t_rr[0], t_rr[1], t_rr[2]}, color(c) {}
    Object(Geometry *geo, const double t_rr[3], double (*t_ir)(double))
            : geometry(geo), rr{t_rr[0], t_rr[1], t_rr[2]}, ir(t_ir) {}
    Object(Geometry *geo, const double t_rr[3], double (*t_ir)(double), Color c)
            : geometry(geo), rr{t_rr[0], t_rr[1], t_rr[2]}, ir(t_ir), color(c) {}


    [[nodiscard]] Intersection intersect(const Line &l) const override { return geometry->intersect(l); };
};

struct Plane : Geometry {
    const Point3 origin;
    const Vec3 normal;

    Plane(const Point3 &p, const Vec3 &v) : origin(p), normal(v) {}

    [[nodiscard]] Intersection intersect(const Line &l) const override;
};

struct Sphere : Geometry {
    const Point3 center;
    const double radius;

    Sphere(const Point3 &p, double r) : center(p), radius(r) {}

    [[nodiscard]] Intersection intersect(const Line &l) const override;
};

struct Triangle : Geometry {
    const Point3 origin;
    const Vec3 edge1, edge2;

    Triangle(const Point3 &p1, const Point3 &p2, const Point3 &p3)
            : origin(p1), edge1(p2 - p1), edge2(p3 - p1) { }

    [[nodiscard]] Intersection intersect(const Line &l) const override;
};

struct Parallelogram : Geometry {
    const Point3 origin;
    const Vec3 edge1, edge2;

    Parallelogram(const Point3 &p1, const Point3 &p2, const Point3 &p3)
                 : origin(p1), edge1(p2 - p1), edge2(p3 - p1) { }

    [[nodiscard]] Intersection intersect(const Line &l) const override;
};

struct Bounding_Box : Geometry {
    const Point3 origin;
    const Vec3 diagonal;

    static Point3 generate_origin(const Point3 &p1, const Point3 &p2) {
        return {std::min(p1[0], p2[0]), std::min(p1[1], p2[1]), std::min(p1[2], p2[2])};
    }

    static Vec3 generate_diagonal(const Point3 &p1, const Point3 &p2) {
        const Point3 p = {std::max(p1[0], p2[0]), std::max(p1[1], p2[1]), std::max(p1[2], p2[2])};
        return p - generate_origin(p1, p2);
    }

    Bounding_Box(const Point3 &p1, const Point3 &p2)
                : origin(generate_origin(p1, p2)), diagonal(generate_diagonal(p1, p2)) { }

    [[nodiscard]] Intersection intersect(const Line &l) const override;
};
