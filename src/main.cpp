#include <bits/stdc++.h>
#include "kdtree.h"
#include "tracing.h"
#include "vec.h"
#include "object.h"
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

    Sphere sp({3, 1, 0}, 0.5, trans, 1.5);
    objects.push_back(&sp);

    Plane pl({0, -1, 0}, {0, 1, 0}, opaco, 0);
    objects.push_back(&pl);

    // Parallelogram tr({4, 0, 0}, {5, 0, 1}, {5, 1, 0}, trans, 1);
    // objects.push_back(&tr);

    const Camera cam({-1, 0, 0}, {1, 0, 0}, {0, 1, 0}, 500, 500, 0.5, 0.5);

    // raycast(cam, objects);

    auto photons = emit_photons({3, 3, 0}, 1e7, objects);
    starfield_projection(cam, photons);

    KDTree kd(photons);

    printf("It took %fs\n", (clock() - start_time)/1000.0);
}
