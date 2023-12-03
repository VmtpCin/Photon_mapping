#include "tracing.h"
#include "object.h"

void raycast(const Camera &cam, const std::vector<Object*> &objs) {
    for (int i = 0; i < cam.vres; ++i)
        for (int j = 0; j < cam.hres; ++j) {
            Line l{cam.origin, cam.pixel_ray(i, j)};
            Intersection inter;

            for (const auto &obj : objs) {
                auto temp = obj->intersect(l);
                if (temp < inter)
                    inter = temp;
            }

            bool b_inter = inter.t < 1e100;
            short color = b_inter ? 255 * abs(   inter.normal.normalize()
                                               * l.dir.normalize())
                                          : 0;

            cam.grid[i][j] = color;
        }
}
