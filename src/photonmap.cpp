#include "tracing.h"
#include "kdtree.h"
#include "geometry.h"

void russian_rolette(const Line &l, double ir, Color intensity,
                     const std::vector<Object*> &objs, KDTree &kdt) {
    Intersection inter;
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
        Color new_intensity = intensity & obj->color;

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

            kdt.push_back({p, l.dir, intensity / 2});

            russian_rolette(n_l, ir, new_intensity / 2, objs, kdt);
        } else if (dice < obj->rr[0] + obj->rr[1]) { // reflexao
            Vec3 n_dir = 2 * normal * (normal * -l.dir) - (-l.dir);
            russian_rolette({p, n_dir}, ir, new_intensity, objs, kdt);
        } else if (dice < obj->rr[0] + obj->rr[1] + obj->rr[2]) { // refracao
            Vec3 n = normal.normalize();
            double cos_teta_i = -l.dir.normalize() * n;
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
                russian_rolette({p, n_dir}, obj->ir, new_intensity, objs, kdt);
            } else {
                Vec3 n_dir = 2 * normal * (normal * -l.dir) - (-l.dir);
                russian_rolette({p, n_dir}, ir, new_intensity, objs, kdt);
            }
        } else { // absorcao
            kdt.push_back({p, l.dir, intensity});
        }
    }
}

KDTree emit_photons(const Point3 &p, int num, double power, const std::vector<Object*> &objs) {
    KDTree kdt;

    Line l;
    l.origin = p;
    Rand r(-1, 1);
    Color intensity = power / num;

    for (int i = 0; i < num; ++i) {
        do
            l.dir = {r.gen(), r.gen(), r.gen()};
        while (l.dir.length_sq() > 1);

        russian_rolette(l, 1, intensity, objs, kdt);
    }

    return kdt;
}

void visualize_photomap(const Camera &cam, const std::vector<Object*> &objs, const KDTree &kdt) {
    for (int i = 0; i < cam.vres; ++i)
        for (int j = 0; j < cam.hres; ++j) {
            Line l{cam.origin, cam.pixel_ray(i, j)};
            Intersection inter;

            for (const auto &obj : objs) {
                auto temp = obj->intersect(l);
                if (temp < inter)
                    inter = temp;
            }

            cam.grid[i][j] = inter ? kdt.get_intensity(l.t(inter.t), 0.2) : Color(0);
        }
}

void starfield_projection(const Camera &cam, const KDTree &photons) {
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
                cam.grid[x][y] = 255;
        }
    }
}
