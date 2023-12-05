#pragma once
#include "vec.h"
#include <algorithm>
#include <vector>

constexpr double k_kdtree = 1;

struct KDTree : std::vector<Photon> {
private:
    template <int depth>
    void sort(size_t b, size_t e) {
        if (e - b <= 1)
            return;

        constexpr auto comp = [](const Photon &p1, const Photon &p2) {
            return p1.point[depth] < p2.point[depth];
        };

        const auto &v_b = begin();
        const auto middle = (b + e) / 2;
        std::nth_element(v_b + b, v_b + middle, v_b + e, comp);

        sort<(depth + 1) % 3>(b, middle);
        sort<(depth + 1) % 3>(middle + 1, e);
    }

    template <int depth>
    void find_all_near(const Point3 &p, double r_sq, size_t b, size_t e, KDTree &kdt) const {
        if (b == e)
            return;

        const size_t m = (b + e) / 2;
        const Point3 &cur_p = (*this)[m].point;

        if (p.distance_sq(cur_p) <= r_sq)
            kdt.push_back((*this)[m]);

        const double dist = p[depth] - cur_p[depth];

        if (dist < 0) {
            find_all_near<(depth + 1) % 3>(p, r_sq, b, m, kdt);

            if (dist * dist <= r_sq)
                find_all_near<(depth + 1) % 3>(p, r_sq, m + 1, e, kdt);
        } else {
            find_all_near<(depth + 1) % 3>(p, r_sq, m + 1, e, kdt);

            if (dist * dist <= r_sq)
                find_all_near<(depth + 1) % 3>(p, r_sq, b, m, kdt);
        }
    }

    template <int depth>
    Color get_intensity(const Point3 &p, const Vec3 &n, double r_sq, size_t b, size_t e) const {
        if (b == e)
            return 0;

        Color intensity = 0;
        const size_t m = (b + e) / 2;
        const Point3 &cur_p = (*this)[m].point;

        if (p.distance_sq(cur_p) <= r_sq)
            intensity += (*this)[m].I * std::abs(n * (*this)[m].dir); // * (1 - sqrt(p.distance_sq(cur_p)/(k_kdtree * r_sq)));

        const double dist = p[depth] - cur_p[depth];

        if (dist < 0) {
            intensity += get_intensity<(depth + 1) % 3>(p, n, r_sq, b, m);

            if (dist * dist <= r_sq)
                intensity += get_intensity<(depth + 1) % 3>(p, n, r_sq, m + 1, e);
        } else {
            intensity += get_intensity<(depth + 1) % 3>(p, n, r_sq, m + 1, e);

            if (dist * dist <= r_sq)
                intensity += get_intensity<(depth + 1) % 3>(p, n, r_sq, b, m);
        }

        return intensity;
    }

public:
    void sort() { shrink_to_fit(); return sort<0>(0, size()); }

    KDTree find_all_near(const Point3 &p, double radius) const {
        KDTree photons;
        find_all_near<0>(p, radius * radius, 0, size(), photons);
        return photons;
    }

    Color get_intensity(const Point3 &p, const Vec3 &n, double radius) const {
        return get_intensity<0>(p, n, radius * radius, 0, size()) / (M_PI * radius * radius);
    }
};
