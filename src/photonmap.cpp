#include "render.h"
#include "object.h"
#include <fstream>

void russian_rolette(const Line &l, const Vec3 &up, double ir,
                     const std::vector<Object*> &objs, std::vector<Photon> &list) {
    Intersection inter({1e100, {0, 0, 0}});
    Object *obj = nullptr;

    for (const auto &o : objs) {
        auto temp = o->intersect(l);
        if (temp < inter)
            inter = temp, obj = o;
    }

    if (obj != nullptr) {
        Rand r(0, 1);
        Point3 p = l.t(inter.t);
        Vec3 &normal = inter.normal;
        double dice = r.gen();

        if (dice < obj->rr[0]) { // difusao
            double r1 = r.gen(), r2 = r.gen();

            const Vec3 w = (normal ^ up).normalize();
            const Vec3 v = (w ^ up).normalize();
            const Vec3 u = (v ^  w).normalize();

            double teta = 1 / cos(sqrt(r1));
            double phi  = 2 * M_PI * r2;
            Vec3     sp = sin(teta) * cos(phi) * u + sin(teta) * sin(phi) * v + cos(teta) * w;
            Vec3  n_dir = sp[0] * u + sp[1] * v + sp[2] * w;
 
            Line n_l({p, n_dir});

            russian_rolette(n_l, up, ir, objs, list);
        } else if (dice < obj->rr[0] + obj->rr[1]) { // reflexao
            Vec3 n_dir = 2 * normal * (normal * -l.dir) - (-l.dir);
            russian_rolette({p, n_dir}, up, ir, objs, list);
        } else if (dice < obj->rr[0] + obj->rr[1]) { // refracao
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
                russian_rolette({p, n_dir}, up, obj->ir, objs, list);
            }
        } else { // absorcao
            list.push_back({p, l.dir, 0});
        }
    }
}

std::vector<Photon> emit_photons(const Point3 &p, int num,
                    const Vec3 &up, const std::vector<Object*> &objs) {
    std::vector<Photon> list;

    Line l;
    l.origin = p;
    Rand r(-1, 1);

    for (int i = 0; i < num; ++i) {
        do
            l.dir = {r.gen(), r.gen(), r.gen()};
        while (l.dir.length_sq() > 1);

        russian_rolette(l, up, 1, objs, list);
    }

    return list;
}

void starfield_projection(const Camera &cam, const std::vector<Photon> &photons) {
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
            double v_del_h = v * cam.desl_h / cam.desl_h.length_sq();
            double v_del_v = v * cam.desl_v / cam.desl_v.length_sq();

            int x = round(v_del_v);
            int y = round(v_del_h);

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