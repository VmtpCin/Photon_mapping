#include <bits/stdc++.h>
#include "render.h"
#include "vec.h"
#include "object.h"
using namespace std;


int main() {
    vector<Object*> objects;
    double trans[] = {0, 0, 1};
    double opaco[] = {0, 0, 0};

    Sphere sp({3, 1, 0}, 0.25, trans, 1.5);
    objects.push_back(&sp);

    Plane pl({0, -1, 0}, {0, 1, 0}, opaco, 0);
    objects.push_back(&pl);

    const Camera cam({-1, 0, 0}, {1, 0, 0}, {0, 1, 0}, 500, 500, 0.5, 0.5);

    // raycast(cam, objects);
    auto photons = emit_photons({3, 3, 0}, 1e6, cam.up, objects);

    starfield_projection(cam, photons);

    // printf("%zu end!", photons.size());
}