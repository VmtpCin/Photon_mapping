#pragma once
#include "vec.h"
#include "render.h"
#include <algorithm>
#include <vector>

struct KDTree {
    struct Node {
        Point3 p;
        Node *right, *left;        
    };

    Node *root;

    template<int depth>
    Node* insert(std::vector<Photon>::iterator photon_b,
                 std::vector<Photon>::iterator photon_e) {
        if (photon_b == photon_e)
            return nullptr;
        if (photon_e - photon_b == 1)
            return &Node({points[0], nullptr, nullptr});

        const auto points_m = (photon_b + photon_e) / 2;
        std::nth_element(photon_b, points_m, photon_e, Photon::is_less<depth>());

        Node n;
        n.p = *points_m;
        n.left  = insert<(depth + 1) % 3>(photon_b, points_m);
        n.right = insert<(depth + 1) % 3>(points_m + 1, photon_e);

        return &n;
    }

    void insert(std::vector<Photon> &v) {
        root = insert<0>(v.begin(), v.end());
    }
};
