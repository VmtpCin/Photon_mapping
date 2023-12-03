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
    double trans[] = {0, 0, 0};
    double opaco[] = {1, 0, 0};

    objects.push_back(new Object(new Sphere({3, 0, 0}, 0.5), trans, 0.0, {0, 1, 0}));
    objects.push_back(new Object(new Plane({0, -1, 0}, {0, 1, 0}), opaco, 0.0, {0, 1, 0}));

    // objects.push_back(new Object(new Triangle({4, 0, 0}, {5, 0, 1}, {5, 1, 0})));

    const Camera cam({-1, 0, 0}, {1, 0, 0}, {0, 1, 0}, 500, 500, 0.5, 0.5);

    // raycast(cam, objects);

    KDTree kdt = emit_photons({3, 3, 3}, 1e6, 1e6, objects);
    kdt.sort();

    // starfield_projection(cam, kdt);
    visualize_photomap(cam, objects, kdt);

    cam.print();

    printf("It took %fs\n", (clock() - start_time)/1000.0);
}
