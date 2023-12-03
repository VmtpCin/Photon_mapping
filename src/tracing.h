#pragma once
#include "kdtree.h"
#include "vec.h"
#include "object.h"
#include <vector>
#include <random>
#include <string>

extern const std::string path;

struct Rand {
private:
    std::random_device rd;
    std::mt19937 gene;
    std::uniform_real_distribution<double> dis;

public:
    Rand(double first, double last) : gene(rd()), dis(first, last) { }

    double gen() {
        return dis(gene);
    }
};

struct Camera {
    Point3 origin;
    Vec3 oa, up;
    int hres, vres;
    double tamx, tamy;
    Vec3 vec_initial, desl_h, desl_v;
    mutable std::vector<std::vector<Color>> grid;

    Camera(const Point3 &o, const Point3 &t, const Vec3 &u, int hr,
           int vr, double tx, double ty, double dist = 1) : origin(o),
           hres(hr), vres(vr), tamx(tx), tamy(ty) {

            oa = (t - o ).normalize();
        Vec3 b = (u ^ oa).normalize();
            up = (b ^ oa).normalize();


        desl_h = ((2 * tamx) / (hres - 1)) * b;
        desl_v = ((2 * tamy) / (vres - 1)) * up;

        vec_initial = oa * dist - tamx * b - tamy * up;

        grid.assign(hres, std::vector<Color>(vres, Color()));
    }

    Vec3 pixel_ray(int i, int j) const {
        return vec_initial + i * desl_h + j * desl_v;
    }
};

extern void raycast(const Camera &cam, const std::vector<Object*> &objs);
extern KDTree emit_photons(const Point3 &p, int num, double power,
                           const std::vector<Object*> &objs);
extern void starfield_projection(const Camera &cam, const KDTree &kdt);
extern void visualize_photomap(const Camera &cam, const std::vector<Object*> &objs, const KDTree &kdt);
