#pragma once
#include "kdtree.h"
#include "vec.h"
#include "geometry.h"
#include <fstream>
#include <vector>
#include <random>
#include <string>
#include <iostream>

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
           int vr, double dist = 1, double tx = 1, double ty = 1) : origin(o),
           hres(hr), vres(vr), tamx(tx), tamy(ty) {

            oa = (t - o ).normalize();
        Vec3 b = (u ^ oa).normalize();
            up = (b ^ oa).normalize();


        desl_h = (tamx / (hres - 1)) * b;
        desl_v = (tamy / (vres - 1)) * up;

        vec_initial = oa * dist - tamx * b / 2 - tamy * up / 2;

        grid.assign(hres, std::vector<Color>(vres, Color()));
    }

    Vec3 pixel_ray(int i, int j) const {
        Vec3 h = (tamx * i) * (oa ^ up) / (hres - 1);
        Vec3 v = (tamy * j) *    up     / (vres - 1);
        return vec_initial + h + v;
    }

    void print() const {
        constexpr double max_v = 1;
        std::ofstream outFile(path);
        outFile << "P3\n" << hres << " " << vres << "\n255\n";

        for (int j = 0; j < vres; ++j)
            for (int i = 0; i < hres; ++i) {
                Color3 c(grid[i][j]);

                outFile << std::clamp(int(255 * c.R / max_v), 0, 255) << " "
                        << std::clamp(int(255 * c.G / max_v), 0, 255) << " "
                        << std::clamp(int(255 * c.B / max_v), 0, 255) << std::endl;
            }

        outFile.close();
    }
};

extern void simplecast(const Camera &cam, const std::vector<const Object*> &objs);
extern std::pair<KDTree, KDTree> emit_photons(const std::vector<Light> &lights,
                                              const std::vector<const Object*> &objs);
extern std::pair<KDTree, KDTree> emit_photons_th(const std::vector<Light> &lights,
                                                 const std::vector<const Object*> &objs, int threads);

extern void starfield_projection(const Camera &cam, const KDTree &kdt);
extern void visualize_radiance(const Camera &cam, const std::vector<const Object*> &objs, const KDTree &kdt);

extern void raycast(const Camera &cam, const std::vector<const Object*> &objs, const std::vector<Light> &lights, const KDTree &kdt, const KDTree &kdt_refraction);
extern void raycast_th(const Camera &cam, const std::vector<const Object*> &objs, const std::vector<Light> &lights, const KDTree &kdt, const KDTree &kdt_refraction, int threads);