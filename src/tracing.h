#pragma once
#include "kdtree.h"
#include "vec.h"
#include "geometry.h"
#include <fstream>
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
        Vec3 h = (2 * tamx * i) * (oa ^ up) / (hres - 1);
        Vec3 v = (2 * tamy * j) *    up     / (vres - 1);
        return vec_initial + h + v;
    }

    void print() const {
        constexpr double max_v = 1;
        std::ofstream outFile(path);
        outFile << "P3\n" << hres << " " << vres << "\n255\n";

        // for (const auto &vs : grid)
        //     for (const auto &v : vs)
        //         max_v = std::max({max_v, v.R, v.G, v.B});

        for (int j = 0; j < hres; ++j)
            for (int i = 0; i < vres; ++i)
                outFile << std::clamp(int(255 * grid[i][j].R / max_v), 0, 255) << " "
                        << std::clamp(int(255 * grid[i][j].G / max_v), 0, 255) << " "
                        << std::clamp(int(255 * grid[i][j].B / max_v), 0, 255) << std::endl;

        outFile.close();
    }
};

extern void simplecast(const Camera &cam, const std::vector<Object*> &objs);
extern KDTree emit_photons(const Point3 &p, int num, double power,
                           const std::vector<Object*> &objs);
extern KDTree emit_photons_th(const Point3 &p, int num, double power,
                              const std::vector<Object*> &objs, int threads);

extern void starfield_projection(const Camera &cam, const KDTree &kdt);
extern void visualize_radiance(const Camera &cam, const std::vector<Object*> &objs, const KDTree &kdt);

extern void raycast(const Camera &cam, const std::vector<Object*> &objs, const KDTree &kdt);
extern void raycast_th(const Camera &cam, const std::vector<Object*> &objs, const KDTree &kdt, int threads);