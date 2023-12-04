#include "tracing.h"
#include "geometry.h"
#include <iostream>

void simplecast(const Camera &cam, const std::vector<Object*> &objs) {
    for (int i = 0; i < cam.vres; ++i)
        for (int j = 0; j < cam.hres; ++j) {
            Line l{cam.origin, cam.pixel_ray(i, j)};
            Intersection inter;
            Object *obj;

            for (const auto &o : objs) {
                auto temp = o->intersect(l);
                if (temp < inter)
                    inter = temp, obj = o;
            }

            if (inter) {
                Vec3 n = inter.normal.normalize();
                Vec3 d = l.dir.normalize();

                cam.grid[i][j] = obj->color * std::abs(n * d);
            } else {
                cam.grid[i][j] = 0;
            }
        }
}

Color castray(const Line &l, const std::vector<Object*> &objs, const KDTree &kdt, int depth) {
    if (depth > 5)
        return 0;

    Intersection inter;
    Object *obj = nullptr;
    Color result = 0;

    for (const auto &o : objs) {
        auto temp = o->intersect(l);
        if (temp < inter)
            inter = temp, obj = o;
    }

    if (inter) {
        Point3 p = l.t(inter.t);
        Vec3 &normal = inter.normal;

        if (obj->rr[1] > 0) {
            Vec3 n_dir = l.dir - 2 * normal * (normal * l.dir) / normal.length_sq();
            result += obj->rr[1] * castray({p, n_dir}, objs, kdt, depth + 1);
        }

        if (obj->rr[2] > 0) {
            Vec3 n = normal.normalize();
            Vec3 dir = l.dir.normalize();
            double cosI = -dir * n;
            double eta = 1.0/obj->ir;

            if (cosI < 0) {
                cosI *= -1;
                n *= -1;
                eta = 1/eta;
            }

            double temp = 1 - (eta * eta) * (1 - cosI * cosI);

            if (temp > 0) {
                Vec3 n_dir = eta * dir - (sqrt(temp) - cosI * eta) * n;
                result += obj->rr[2] * castray({p, n_dir}, objs, kdt, depth + 1);
            } else {
                Vec3 n_dir = dir - 2 * n * (n * dir);
                result += obj->rr[2] * castray({p, n_dir}, objs, kdt, depth + 1);
            }
        }

        if (obj->rr[1] + obj->rr[2] < 1)
            result += (1 - (obj->rr[1] + obj->rr[2]))
                    * kdt.get_intensity(l.t(inter.t), inter.normal.normalize(), 0.15);
    }

    return result;
}

void raycast(const Camera &cam, const std::vector<Object*> &objs, const KDTree &kdt) {
    for (int i = 0; i < cam.vres; ++i)
        for (int j = 0; j < cam.hres; ++j)
            cam.grid[i][j] = castray({cam.origin, cam.pixel_ray(i, j)}, objs, kdt, 0);
}
