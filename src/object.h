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

struct Object {
    double rr[3] = {0, 0, 0}, ir = 1;
    Color color = 1;

    Object() : rr{0, 0, 1}, ir(1) {};
    Object(const double t_rr[3], double t_ir) : rr{t_rr[0], t_rr[1], t_rr[2]}, ir(t_ir) {}
    Object(const double t_rr[3], double t_ir, Color &c)
          : rr{t_rr[0], t_rr[1], t_rr[2]}, ir(t_ir), color(c) {}

    virtual Intersection intersect(const Line &l) const { return {inf, l.dir}; };
};

struct Plane : Object {
    const Point3 origin;
    const Vec3 normal;

    Plane(const Point3 &p, const Vec3 &v) : origin(p), normal(v) {}

    Plane(const Point3 &p, const Vec3 &v, const double t_rr[3], double t_ir)
            : Object(t_rr, t_ir), origin(p), normal(v) { }

    Plane(const Point3 &p, const Vec3 &v, const double t_rr[3], double t_ir, Color c)
            : Object(t_rr, t_ir, c), origin(p), normal(v) { }

    Intersection intersect(const Line &l) const override;
};

struct Sphere : Object {
    const Point3 center;
    const double radius;

    Sphere(const Point3 &p, double r) : center(p), radius(r) {}

    Sphere(const Point3 &p, double r, const double t_rr[3], double t_ir)
            : Object(t_rr, t_ir), center(p), radius(r) { }

    Sphere(const Point3 &p, double r, const double t_rr[3], double t_ir, Color c)
            : Object(t_rr, t_ir, c), center(p), radius(r) { }

    Intersection intersect(const Line &l) const override;
};

struct Triangle : Object {
    const Point3 origin;
    const Vec3 edge1, edge2;

    Triangle(const Point3 &p1, const Point3 &p2, const Point3 &p3)
            : origin(p1), edge1(p1 - p2), edge2(p1 - p3) { }
    Triangle(const Point3 &p1, const Point3 &p2, const Point3 &p3,
             const double t_rr[3], double t_ir) : Object(t_rr, t_ir),
             origin(p1), edge1(p1 - p2), edge2(p1 - p3) { }
    Triangle(const Point3 &p1, const Point3 &p2, const Point3 &p3,
             const double t_rr[3], double t_ir, Color c) : Object(t_rr, t_ir, c),
             origin(p1), edge1(p1 - p2), edge2(p1 - p3) { }

    Intersection intersect(const Line &l) const override;
};

struct Parallelogram : Object {
    const Point3 origin;
    const Vec3 edge1, edge2;

    Parallelogram(const Point3 &p1, const Point3 &p2, const Point3 &p3)
                 : origin(p1), edge1(p1 - p2), edge2(p1 - p3) { }
    Parallelogram(const Point3 &p1, const Point3 &p2, const Point3 &p3,
                  const double t_rr[3], double t_ir) : Object(t_rr, t_ir),
                  origin(p1), edge1(p1 - p2), edge2(p1 - p3) { }
    Parallelogram(const Point3 &p1, const Point3 &p2, const Point3 &p3,
                  const double t_rr[3], double t_ir, Color c) : Object(t_rr, t_ir, c),
                  origin(p1), edge1(p1 - p2), edge2(p1 - p3) { }

    Intersection intersect(const Line &l) const override;
};


extern Point3 bezier_curve(const std::vector<Point3> &points, double t);
extern Point3 bezier_triangle(double s, double t);

extern std::vector<Point3> create_surface_triangle(Point3 (*f)(double, double),
                            int s_step, int t_step, std::vector<Object*> &objs);

extern std::vector<Point3> create_bezier_superfice(const std::vector<std::vector<Point3>> &control,
                                                   int s_step, int t_step, std::vector<Object*> &objs);