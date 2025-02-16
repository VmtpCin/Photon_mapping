#pragma once
#include "vec.h"
#include "geometry.h"
#include <algorithm>
#include <vector>

struct Photon {
    Point3 point;
    const Object *obj;
    Vec3 dir;
    Color I;
};

struct Light {
    Point3 pos;
    Color Intensity;
    int quantity;
};

struct KDTree : std::vector<Photon> {
private:
    static constexpr float k_kdtree = 2;

    template <int depth>
    void sort(size_t b, size_t e) {
        if (e - b <= 1)
            return;

        constexpr auto comp = [](const Photon &p1, const Photon &p2) {
            return p1.point[depth] < p2.point[depth];
        };

        const auto &v_b = begin();
        const size_t middle = (b + e) / 2;
        std::nth_element(v_b + b, v_b + middle, v_b + e, comp);

        sort<(depth + 1) % 3>(b, middle);
        sort<(depth + 1) % 3>(middle + 1, e);
    }

    template <int depth>
    void find_all_near(const Point3 &p, const Object *obj, double r_sq, size_t b, size_t e, KDTree &kdt) const {
        if (b == e)
            return;

        const size_t m = (b + e) / 2;
        const Point3 &cur_p = (*this)[m].point;

        if (p.distance_sq(cur_p) <= r_sq)
            kdt.push_back((*this)[m]);

        const double dist = p[depth] - cur_p[depth];

        if (dist <= 0 || dist * dist <= r_sq)
            find_all_near<(depth + 1) % 3>(p, obj, r_sq, b, m, kdt);
        
        if (dist >= 0 || dist * dist <= r_sq)
            find_all_near<(depth + 1) % 3>(p, obj, r_sq, m + 1, e, kdt);
    }

    template <int depth>
    Color get_intensity(const Point3 &p, const Object *obj, const Vec3 &n, double r_sq, size_t b, size_t e) const {
        if (b == e)
            return 0;

        Color intensity = 0;
        const size_t m = (b + e) / 2;
        const Point3 &cur_p = (*this)[m].point;

        if (p.distance_sq(cur_p) <= r_sq) {
            const float d = -(n * ((*this)[m].dir));

            if ((*this)[m].obj == obj && d > 0)
                intensity += (*this)[m].I * d * (1 - sqrt(p.distance_sq(cur_p)/(k_kdtree * r_sq)));
        }

        const float dist = p[depth] - cur_p[depth];

        if (dist <= 0 || dist * dist <= r_sq)
            intensity += get_intensity<(depth + 1) % 3>(p, obj, n, r_sq, b, m);
        
        if (dist >= 0 || dist * dist <= r_sq)
            intensity += get_intensity<(depth + 1) % 3>(p, obj, n, r_sq, m + 1, e);

        return intensity;
    }

public:
    void combine_close_photons(double threshold) {
        if (empty()) return;

        std::vector<Photon> combined;
        combined.reserve(size());

        const double threshold_sq = threshold * threshold;
        Photon current = (*this)[0];

        for (size_t i = 1; i < size(); ++i) {
            Photon &next = (*this)[i];

            if (current.point.distance_sq(next.point) <= threshold_sq) {
                current.I += next.I;
            } else {
                combined.push_back(current);  // Keep the merged photon
                current = next;               // Start a new group
            }
        }

        std::cout << combined.size() << std::endl << size() << std::endl;
    }

    void sort() { shrink_to_fit(); sort<0>(0, size()); }

    KDTree find_all_near(const Point3 &p, const Object *obj, double radius) const {
        KDTree photons;
        find_all_near<0>(p, obj, radius * radius, 0, size(), photons);
        return photons;
    }

    Color get_intensity(const Point3 &p, const Object *obj, const Vec3 &n, double radius) const {
        return 3 * get_intensity<0>(p, obj, n.normalize(), radius * radius, 0, size()) / (M_PI * radius * radius);
    }
};
