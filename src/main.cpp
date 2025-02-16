#include <bits/stdc++.h>
#include "kdtree.h"
#include "tracing.h"
#include "vec.h"
#include "geometry.h"
#include "surface.h"
#include "reader.h"
using namespace std;

#ifdef MY_PATH
const std::string path = MY_PATH;
#else
const std::string path = "out.ppm";
#endif

int main(int argc, char *argv[]) {
    clock_t start_time = clock();

    setupColor();

    vector<const Object*> objects;
    double trans[] = {0, 0, 1};
    double reflx[] = {0, 1, 0};
    double opaco[] = {1, 0, 0};

    constexpr auto n_glass = [](double wavelength) {
        double n = 0.9;
        double wl_sq = wavelength * wavelength;

        n += (1.03961212  * wl_sq) / (wl_sq - 600069867e-3);
        n += (0.231792344 * wl_sq) / (wl_sq - 2.00179144e-2);
        n += (1.01046945  * wl_sq) / (wl_sq - 1.03560653e2);

        return n;
    };

    constexpr auto n_water = [](double wavelength) {
        double n = 1.0;
        double wl_sq = wavelength * wavelength;
    
        n += (0.5675888 * wl_sq) / (wl_sq - 0.050263605);
        n += (0.1724913 * wl_sq) / (wl_sq - 0.0664945);
        n += (0.0205836 * wl_sq) / (wl_sq - 0.1531933);
    
        return std::sqrt(n);
    };

    constexpr auto n_air = [](double wavelength) {
        return 1.0;
    };

    //objects.push_back(new Object(new Parallelogram({v1x, v1y, v1z}, {v2x, v2y, v2z}, {v3x,  v3y, v3z}), tipo do objeto, cor rgb)); // right
    objects.push_back(new Object(new Parallelogram({0, -1, -1}, {4, -1, -1}, {0,  1, -1}), opaco, Color(Color3(1, 0, 0)))); // right
    objects.push_back(new Object(new Parallelogram({0, -1,  1}, {0,  1,  1}, {4, -1,  1}), opaco, Color(Color3(0, 1, 0)))); // left
    objects.push_back(new Object(new Parallelogram({4, -1, -1}, {4, -1,  1}, {4,  1, -1}), opaco, Color(Color3(0, 0, 1)))); // back
    objects.push_back(new Object(new Parallelogram({0, -1, -1}, {4, -1, -1}, {0, -1,  1}), opaco, Color(Color3(1, 1, 1)))); // down
    objects.push_back(new Object(new Parallelogram({0,  1, -1}, {4,  1, -1}, {0,  1,  1}), opaco, Color(Color3(1, 1, 1)))); // up

    objects.push_back(new Object(new Sphere({2, -0.55, 0}, 0.4), trans, n_glass));

    // std::vector<std::vector<Point3>> bezier_control;
    // bezier_control.push_back(std::vector<Point3>({{1.5, -1, -1}, {2, 0.5, -1}, {2.5, -1, -1}, {3, 1, -1}}));
    // bezier_control.push_back(std::vector<Point3>({{1.5, -1,  1}, {2, 0.5,  1}, {2.5, -1,  1}, {3, 1, 1}}));
    // create_bezier_superfice(bezier_control, 100, 1, objects);

    std::vector<Light> lights;
    lights.push_back({{2, 0.3, 0}, 250 * Color(Color3(1, 1, 1)), int(1e5)});

    const Camera cam({-1, 0, 0}, {1, 0, 0}, {0, 1, 0}, 500, 500);

    auto [kdt, kdt_refraction] = emit_photons_th(lights, objects, 12);
    kdt.sort();
    kdt_refraction.sort();

    // simplecast(cam, objects);

    // starfield_projection(cam, kdt_refraction);
    // visualize_radiance(cam, objects, kdt);
    // raycast(cam, objects, lights, kdt, kdt_refraction);
    raycast_th(cam, objects, lights, kdt, kdt_refraction, 12);

    cam.print();

    printf("It took %fs\n", (clock() - start_time)/1000.0);

    return 0;
}
