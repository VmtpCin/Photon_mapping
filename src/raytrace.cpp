#include "tracing.h"
#include "object.h"
#include <fstream>

void raycast(const Camera &cam, const std::vector<Object*> &objs) {
    std::ofstream outFile(path);

    outFile << "P3\n" << cam.hres << " " << cam.vres << "\n255\n";
    for (int i = 0; i < cam.vres; ++i)
        for (int j = 0; j < cam.hres; ++j) {
            Line l{cam.origin, cam.pixel_ray(j, i)};
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

            outFile << color << " "
                    << color << " "
                    << color << std::endl;
        }

    outFile.close();
}
