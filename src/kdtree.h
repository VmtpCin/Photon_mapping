#pragma once
#include "vec.h"
#include "tracing.h"
#include <algorithm>
#include <vector>
#include <iostream>

struct KDTree {
    struct Node {
        Photon p;
        Node *left, *right;        
    };

    Node *root;

    KDTree() {}
    KDTree(std::vector<Photon> &v) {
        insert(v);
    }

    template<int depth>
    Node* insert(std::vector<Photon> &v, size_t begin, size_t end) {
        if (begin == end)
            return nullptr;
        if (end - begin == 1)
            return new Node({v[begin], nullptr, nullptr});

        const auto middle = (begin + end) / 2;

        constexpr auto comparator = [](const Photon &p1, const Photon &p2) {
            return p1.point[depth] < p2.point[depth];
        };

        const auto &v_b = v.begin();
        std::nth_element(v_b + begin, v_b + middle, v_b + end, comparator);
        Node *n = new Node({v[middle],
                            insert<(depth + 1) % 3>(v, begin, middle),
                            insert<(depth + 1) % 3>(v, middle + 1, end)});

        return n;
    }

    void insert(std::vector<Photon> &v) {
        root = insert<0>(v, 0, v.size());
    }
};
