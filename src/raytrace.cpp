#include "render.h"
#include "object.h"
#include <fstream>

void raycast(const Camera &cam, const std::vector<Object*> &objs) {
    std::ofstream outFile(path);

    outFile << "P3\n" << cam.hres << " " << cam.vres << "\n255\n";
    for (int i = 0; i < cam.vres; ++i)
        for (int j = 0; j < cam.hres; ++j) {
            Line l{cam.origin, cam.pixel_ray(j, i)};
            Intersection inter = {1e100, {0, 0, 0}};

            for (const auto &obj : objs) {
                auto temp = obj->intersect(l);
                if (temp < inter)
                    inter = temp;
            }

            bool color = inter.t < 1e100;
            outFile << 255 * color << " "
                    << 255 * color << " "
                    << 255 * color << std::endl;
        }

    outFile.close();
}
