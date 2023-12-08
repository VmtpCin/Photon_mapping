#include <bits/stdc++.h>
#include "kdtree.h"
#include "tracing.h"
#include "vec.h"
#include "geometry.h"
#include "surface.h"
using namespace std;

#ifdef MY_PATH
const std::string path = MY_PATH;
#else
const std::string path = "out.ppm";
#endif

int main() {
    clock_t start_time = clock();

    vector<Object*> objects;
    double trans[] = {0, 0.15, 0.8};
    double reflx[] = {0, 1, 0};
    double opaco[] = {1, 0, 0};

    objects.push_back(new Object(new Parallelogram({0, -1, -1}, {4, -1, -1}, {0,  1, -1}), opaco, {1, 0, 0})); // right
    objects.push_back(new Object(new Parallelogram({0, -1,  1}, {0,  1,  1}, {4, -1,  1}), opaco, {0, 1, 0})); // left
    objects.push_back(new Object(new Parallelogram({4, -1, -1}, {4, -1,  1}, {4,  1, -1}), opaco, {0, 0, 1})); // back
    objects.push_back(new Object(new Parallelogram({0, -1, -1}, {4, -1, -1}, {0, -1,  1}), opaco)); // down
    objects.push_back(new Object(new Parallelogram({0,  1, -1}, {4,  1, -1}, {0,  1,  1}), opaco)); // up

    // objects.push_back(new Object(new Sphere({3, -0.6,  0.5}, 0.4), reflx));
    // objects.push_back(new Object(new Sphere({2.5, -0.6, -0.5}, 0.4), trans, 3));

    // std::vector<std::vector<Point3>> bezier_control;
    // bezier_control.push_back(std::vector<Point3>({{1.5, -1, -1}, {2, 0.5, -1}, {2.5, -1, -1}, {3, 1, -1}}));
    // bezier_control.push_back(std::vector<Point3>({{1.5, -1,  1}, {2, 0.5,  1}, {2.5, -1,  1}, {3, 1, 1}}));
    // objects.push_back(new Object(new Sphere({2, 0.3, 0}, 0.2), opaco, {1, 1, 0}));
    // objects.push_back(new Object(new Sphere({3.4, -0.6, 0}, 0.4), opaco, {1, 0, 1}));

    // create_bezier_superfice(bezier_control, 100, 1, objects);
    // printf("%zu\n", objects.size());

    const Camera cam({-1, 0, 0}, {1, 0, 0}, {0, 1, 0}, 500, 500, 0.5, 0.5);

    KDTree kdt = emit_photons_th({2, 0, 0}, 1e6, 30, objects, 12);
    kdt.sort();

    // simplecast(cam, objects);

    // starfield_projection(cam, kdt);
    // visualize_radiance(cam, objects, kdt);
    // raycast(cam, objects, kdt);
    raycast_th(cam, objects, kdt, 12);

    cam.print();

    printf("It took %fs\n", (clock() - start_time)/1000.0);
}
