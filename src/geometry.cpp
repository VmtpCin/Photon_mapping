#include "geometry.h"
#include <vector>

Point3 interpolate(const Point3 &p1, const Point3 &p2, double s) {
    return p1 + s * (p2 - p1);
}

Point3 interpolate(const Point3 &p1, const Point3 &p2, const Point3 &p3, double s, double t) {
    return p1 + s * (p2 - p1) + t * (p3 - p1);
}

Intersection Plane::intersect(const Line &l) const {
    constexpr double epsilon = 1e-5;
    const double a = normal * (origin - l.origin);
    const double b = normal * l.dir;

    if (abs(b) < epsilon || ((a > 0) != (b > 0)))
        return {inf, normal};

    const double t = a / b; 
    return {t > t_min ? t : inf, b < 0 ? normal : -normal};
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
    constexpr double epsilon = 1e-5;

    const Vec3      oc  = origin - l.origin;
    const Vec3   normal =  edge1 ^ edge2;
    const double holder =  l.dir * normal;

    if (abs(holder) < epsilon)
        return {inf, l.dir};

    const Vec3 vp = l.dir ^ oc;
    const double b = edge1 * vp;
    const double c = edge2 * -vp;
    const double t_n = oc * normal;

    if (holder > 0 ? b < 0 || c < 0 || b + c > holder || t_n < 0
                   : b > 0 || c > 0 || b + c < holder || t_n > 0)
        return {inf, l.dir};

    const double t = t_n / holder;
    return {t > t_min ? t : inf, holder < 0 ? normal : -normal};
}

Intersection Parallelogram::intersect(const Line &l) const {
    const Vec3      oc  = origin - l.origin;
    const Vec3   normal =  edge1 ^ edge2;
    const double holder =  l.dir * normal;

    if (abs(holder) < 1e-5)
        return {inf, l.dir};

    const Vec3 vp = l.dir ^ oc;
    const double b = edge1 * vp;
    const double c = edge2 * -vp;
    const double t_n = oc * normal;

    if (holder > 0 ? b < 0 || b > holder || c < 0 || c > holder || t_n < 0
                   : b > 0 || b < holder || c > 0 || c < holder || t_n > 0)
        return {inf, l.dir};

    const double t = t_n / holder;
    return {t > t_min ? t : inf, holder < 0 ? normal : -normal};
}

Intersection Bounding_Box::intersect(const Line &l) const {
    const Vec3 oc1 = origin - l.origin;
    const Vec3 oc2 = oc1 + diagonal;

    const double h1 = l.dir.x(), h2 = l.dir.y(), h3 = l.dir.z();
    const Vec3 vp1 = oc1 ^ l.dir, vp2 = oc2 ^ l.dir;

    const double a1 = diagonal.x() * vp1.x(), a2 = diagonal.x() * vp2.x();
    const double b1 = diagonal.y() * vp1.y(), b2 = diagonal.y() * vp2.y();
    const double c1 = diagonal.z() * vp1.z(), c2 = diagonal.z() * vp2.z();

    constexpr auto check1 = [](double holder, double b, double c) {
        return holder > 0 ? 0 <= b && b <= holder && 0 <= c && c <= holder
                          : 0 >= b && b >= holder && 0 >= c && c >= holder;
    };

    constexpr auto check2 = [](double holder, double t) {
        return holder > 0 ? t > 0 : t < 0;
    };

    Intersection inter;

    if (std::abs(h1) > 1e-5 && check1(h1, -b1, c1) && check2(h1, oc1.x()))
        inter = std::min(inter, {oc1.x() / h1, { 1, 0, 0}});

    if (std::abs(h1) > 1e-5 && check1(h1, b2, -c2) && check2(h1, oc2.x()))
        inter = std::min(inter, {oc2.x() / h1, {-1, 0, 0}});

    if (std::abs(h2) > 1e-5 && check1(h2, -c1, a1) && check2(h2, oc1.y()))
        inter = std::min(inter, {oc1.y() / h2, {0,  1, 0}});

    if (std::abs(h2) > 1e-5 && check1(h2, c2, -a2) && check2(h2, oc2.y()))
        inter = std::min(inter, {oc2.y() / h2, {0, -1, 0}});

    if (std::abs(h3) > 1e-5 && check1(h3, -a1, b1) && check2(h3, oc1.z()))
        inter = std::min(inter, {oc1.z() / h3, {0, 0,  1}});

    if (std::abs(h3) > 1e-5 && check1(h3, a2, -b2) && check2(h3, oc2.z()))
        inter = std::min(inter, {oc2.z() / h3, {0, 0, -1}});

    return inter;
}
