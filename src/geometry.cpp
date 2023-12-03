#include "geometry.h"
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
    double ax2 = 2 * (l.dir * l.dir);
    double b   = 2 * (   oc * l.dir);
    double c   = oc * oc - radius * radius;
    
    double delta = b * b - 2 * ax2 * c;

    if (delta < 0)
        return {inf, l.dir};
    
    delta = sqrt(delta);

    double t = -b - delta > ax2 * t_min ? (-b - delta) / ax2
                                        : (-b + delta) / ax2;

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

    const Vec3 vp = l.dir ^ oc;
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

    if (holder > 0 ? b < 0 || b > holder || c < 0 || c > holder || t_n < 0
                   : b > 0 || b < holder || c > 0 || c < holder || t_n > 0)
        return {inf, l.dir};

    const double t = t_n / holder;
    return {t > t_min ? t : inf, normal};
}
