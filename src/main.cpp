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
    double trans[] = {0, 0, 1};
    double reflx[] = {0, 1, 0};
    double opaco[] = {1, 0, 0};

    objects.push_back(new Object(new Parallelogram({0, -1, -1}, {4, -1, -1}, {0,  1, -1}), opaco, {1, 0, 0})); // right
    objects.push_back(new Object(new Parallelogram({0, -1,  1}, {4, -1,  1}, {0,  1,  1}), opaco, {0, 1, 0})); // left
    objects.push_back(new Object(new Parallelogram({4, -1, -1}, {4, -1,  1}, {4,  1, -1}))); // back
    objects.push_back(new Object(new Parallelogram({0, -1, -1}, {4, -1, -1}, {0, -1,  1}))); // down
    objects.push_back(new Object(new Parallelogram({0,  1, -1}, {4,  1, -1}, {0,  1,  1}))); // up

    objects.push_back(new Object(new Sphere({3, -0.6,  0.5}, 0.4), reflx));
    objects.push_back(new Object(new Sphere({2.5, -0.6, -0.5}, 0.4), trans, 3));

    const Camera cam({-1, 0, 0}, {1, 0, 0}, {0, 1, 0}, 500, 500, 0.5, 0.5);

    simplecast(cam, objects);

    KDTree kdt = emit_photons({2.5, 0.75, 0}, 1e6, 1e5, objects);
    kdt.sort();

    // starfield_projection(cam, kdt);
    // visualize_radiance(cam, objects, kdt);
    raycast(cam, objects, kdt);

    cam.print();

    printf("It took %fs\n", (clock() - start_time)/1000.0);
}
