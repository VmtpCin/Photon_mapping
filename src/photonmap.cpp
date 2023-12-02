#include "tracing.h"
#include "kdtree.h"
#include "object.h"
#include <fstream>

void russian_rolette(const Line &l, double ir, const std::vector<Object*> &objs,
                     KDTree &kdt) {
    Intersection inter({1e100, {0, 0, 0}});
    Object *obj = nullptr;

    for (const auto &o : objs) {
        auto temp = o->intersect(l);
        if (temp < inter)
            inter = temp, obj = o;
    }

    if (obj) {
        Rand r(0, 1);
        Point3 p = l.t(inter.t);
        Vec3 &normal = inter.normal;
        double dice = r.gen();

        if (dice < obj->rr[0]) { // difusao
            double r1 = r.gen(), r2 = r.gen();
            r1 = sqrt(r1);

            const Vec3 w = normal.normalize();
            const Vec3 v = w.perpendicular().normalize();
            const Vec3 u = v ^ w;

            double teta = acos(r1);
            double phi  = 2 * M_PI * r2;
            Vec3  n_dir = sin(teta) * cos(phi) * u + sin(teta) * sin(phi) * v + r1 * w;
 
            Line n_l({p, n_dir});

            kdt.push_back({p, l.dir, 0});

            russian_rolette(n_l, ir, objs, kdt);
        } else if (dice < obj->rr[0] + obj->rr[1]) { // reflexao
            Vec3 n_dir = 2 * normal * (normal * -l.dir) - (-l.dir);
            russian_rolette({p, n_dir}, ir, objs, kdt);
        } else if (dice < obj->rr[0] + obj->rr[1] + obj->rr[2]) { // refracao
            Vec3 wo = -l.dir.normalize();
            Vec3  n = normal;
            double cos_teta_i = wo * n;
            double eta = obj->ir;

            if (cos_teta_i < 0) {
                cos_teta_i *= -1;
                n *= -1;
                eta = 1/eta;
            }

            double temp = 1 - (1 - cos_teta_i * cos_teta_i) / (eta * eta);

            if (temp > 0) {
                double cos_teta2 = sqrt(temp);
                Vec3 n_dir = l.dir / eta - (cos_teta2 - cos_teta_i / eta) * n;
                russian_rolette({p, n_dir}, obj->ir, objs, kdt);
            }
        } else { // absorcao
            kdt.push_back({p, l.dir, 0});
        }
    }
}

KDTree emit_photons(const Point3 &p, int num, const std::vector<Object*> &objs) {
    KDTree kdt;

    Line l;
    l.origin = p;
    Rand r(-1, 1);

    for (int i = 0; i < num; ++i) {
        do
            l.dir = {r.gen(), r.gen(), r.gen()};
        while (l.dir.length_sq() > 1);

        russian_rolette(l, 1, objs, kdt);
    }

    return kdt;
}

void starfield_projection(const Camera &cam, const KDTree &photons) {
    std::ofstream outFile(path);

    bool grid[500][500];

    for (int i = 0; i < 500; ++i)
        for (int j = 0; j < 500; ++j)
            grid[i][j] = false;

    Plane pl_proj(cam.vec_initial + cam.origin, cam.oa);
    Line t_l({cam.origin, cam.vec_initial});
    Point3 p_i = t_l.t(pl_proj.intersect(t_l).t);

    for (const auto &photon : photons) {
        Line l = {photon.point, cam.origin - photon.point};
        auto inter = pl_proj.intersect(l);

        if (inter && l.dir * cam.oa < 0) {
            Vec3 v = l.t(inter.t) - p_i;
            int x = round(v * cam.desl_v / cam.desl_v.length_sq());
            int y = round(v * cam.desl_h / cam.desl_h.length_sq());

            if (0 <= x && x < cam.hres && 0 <= y && y < cam.vres)
                grid[x][y] = true;
        }
    }

    outFile << "P3\n" << cam.hres << " " << cam.vres << "\n255\n";
    for (int i = 0; i < cam.vres; ++i) {
        for (int j = 0; j < cam.hres; ++j)
            outFile << 255 * grid[i][j] << " "
                    << 255 * grid[i][j] << " "
                    << 255 * grid[i][j] << std::endl;
    }

    outFile.close();
}
