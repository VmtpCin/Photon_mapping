#include "object.h"
#include <vector>

Point3 interpolate(const Point3 &p1, const Point3 &p2, double s) {
    return p1 + s * (p2 - p1);
}

Point3 interpolate(const Point3 &p1, const Point3 &p2, const Point3 &p3, double s, double t) {
    return p1 + s * (p2 - p1) + t * (p3 - p1);
}

Intersection Plane::intersect(const Line &l) const {
    const double a = normal * (origin - l.origin);
    const double b = normal * l.dir;

    if (abs(b) < 1e-5 || ((a > 0) ^ (b > 0)))
        return {inf, normal};

    const double t = a / b; 
    return {t > t_min ? t : inf, normal};
}

Intersection Sphere::intersect(const Line &l) const {
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

    if (t < t_min)
        return {inf, l.dir};

    return {t, l.t(t) - center};
}

Intersection Triangle::intersect(const Line &l) const {
    const Vec3      oc  = origin - l.origin;
    const Vec3   normal =  edge1 ^ edge2;
    const double holder =  l.dir * normal;

    if (abs(holder) < 1e-5)
        return {inf, l.dir};

    const Vec3 vp = oc ^ l.dir;
    const double b = edge1 * vp;
    const double c = edge2 * -vp;
    const double t_n = oc * normal;

    if (holder > 0 ? b < 0 || c < 0 || b + c > holder || t_n < 0
                   : b > 0 || c > 0 || b + c < holder || t_n > 0)
        return {inf, l.dir};

    const double t = t_n / holder;
    return {t > t_min ? t : inf, normal};
}

Intersection Parallelogram::intersect(const Line &l) const {
    const Vec3      oc  = origin - l.origin;
    const Vec3   normal =  edge1 ^ edge2;
    const double holder =  l.dir * normal;

    if (abs(holder) < 1e-5)
        return {inf, l.dir};

    const Vec3 vp = oc ^ l.dir;
    const double b = edge1 * vp;
    const double c = edge2 * -vp;
    const double t_n = oc * normal;

    if (holder > 0 ? 0 > b || b > holder || 0 > c || c > holder || t_n < 0
                   : 0 < b || b < holder || 0 < c || c < holder || t_n > 0)
        return {inf, l.dir};

    const double t = t_n / holder;
    return {t > t_min ? t : inf, normal};
}

Point3 bezier_curve(const std::vector<Point3> &points, double t) {
    std::vector<Point3> control;

    control.assign(points.begin(), points.end());

    while (control.size() > 1) {
        for (int i = 0; i < control.size() - 1; i++)
            control[i] = interpolate(control[i], control[i + 1], t);

        control.pop_back();
    }

    return control[0];
}

Point3 bezier_triangle(const std::vector<std::vector<Point3>> &points, double s, double t) {
    std::vector<std::vector<Point3>> control;

    for (const auto &p : points)
        control.push_back(std::vector<Point3>(p.begin(), p.end()));

    while (control.size() > 1) {
        for (int i = 0; i < control.size() - 1; i++)
            for (int j = 0; j < control[i].size(); j++)
                control[i][j] = interpolate(control[i][j], control[i + 1][j], control[i + 1][j + 1], s, t);

        control.pop_back();
    }

    return control[0][0];
}

Point3 bezier_triangle(double s, double t) {
    std::vector<std::vector<Point3>> v;
    v.push_back({{1, 0, 0}});
    v.push_back({{1, 1, 0}, {1, 1, 1}});
    v.push_back({{1, 0, 2}, {1, 2, 2}, {1, 2, 3}});

    return bezier_triangle(v, s, t);
}

std::vector<Point3> create_surface_triangle(Point3 (*f)(double, double),
                    int s_step, int t_step, std::vector<Object*> &objs) {
    std::vector<Point3> curve;

    for (int it_s = 0; it_s <= s_step; ++it_s) {
        for (int it_t = 0; it_t <= t_step; ++it_t) {
            double s = double(it_s) / s_step;
            double t = double(it_t) / t_step;

            if (s + t > 1)
                break;

            curve.push_back(f(s, t));
        }
    }

    return curve;
}

std::vector<Point3> create_bezier_superfice(const std::vector<std::vector<Point3>> &control,
                                            int s_step, int t_step, std::vector<Object*> &objs) {
    std::vector<Point3> curve;

    for (int it_s; it_s <= s_step; ++it_s) {
        std::vector<Point3> next_it;
        double s = double(it_s) / s_step;

        for (auto &c : control)
            next_it.push_back(bezier_curve(c, s));

        for (int it_t = 0; it_t <= t_step; ++it_t)
            curve.push_back(bezier_curve(next_it, double(it_t) / t_step));
    }

    return curve;
}
