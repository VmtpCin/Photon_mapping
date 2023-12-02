#pragma once
#include "vec.h"
#include <algorithm>
#include <vector>

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

        size_t m = (b + e) / 2;
        const Point3 &cur_p = (*this)[m].point;

        if (p.distance_sq(cur_p) <= r_sq)
            kdt.push_back((*this)[m]);

        const double dist = p[depth] - cur_p[depth];

        if (p[depth] < cur_p[depth]) {
            find_all_near<(depth + 1) % 3>(p, r_sq, b, m, kdt);

            if (dist * dist <= r_sq)
                find_all_near<(depth + 1) % 3>(p, r_sq, m + 1, e, kdt);
        } else {
            find_all_near<(depth + 1) % 3>(p, r_sq, m + 1, e, kdt);

            if (dist <= r_sq)
                find_all_near<(depth + 1) % 3>(p, r_sq, b, m, kdt);
        }
    }

public:
    void sort() { return sort<0>(0, size()); }

    KDTree find_all_near(const Point3 &p, double radius) const {
        KDTree photons;
        find_all_near<0>(p, radius * radius, 0, size(), photons);
        return photons;
    }
};
