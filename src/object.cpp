#include "object.h"

Intersection Plane::intersect(const Line &l) const {
    if (abs(normal * l.dir) < 1e-5)
        return {inf, normal};

    double t = (normal * (origin - l.origin)) / (normal * l.dir);

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
    const double t = oc * normal;

    if (holder > 0 ? b < 0 || c < 0 || b + c > holder || t < 0 || t < t_min * holder
                   : b > 0 || c > 0 || b + c < holder || t > 0 || t > t_min * holder)
        return {inf, l.dir};

    return {t / holder, normal};
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
    const double t = oc * normal;

    if (holder > 0 ? 0 > b || b > holder || 0 > c || c > holder || t < 0 || t < t_min * holder
                   : 0 < b || b < holder || 0 < c || c < holder || t > 0 || t > t_min * holder)
        return {inf, l.dir};

    return {abs(t) > t_min * abs(holder) ? t / holder : inf, normal};
}
