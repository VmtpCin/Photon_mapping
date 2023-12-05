#include "tracing.h"
#include "kdtree.h"
#include "geometry.h"
#include <mutex>
#include <atomic>
#include <thread>

void russian_rolette(const Line &l, double ir, Color intensity,
                     const std::vector<Object*> &objs, KDTree &kdt, int depth) {
    if (depth > 5)
        return;

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

            kdt.push_back({p, l.dir.normalize(), new_intensity / 2});

            russian_rolette(n_l, ir, new_intensity / 2, objs, kdt, depth + 1);
        } else if (dice < obj->rr[0] + obj->rr[1]) { // reflexao
            Vec3 n_dir = l.dir - 2 * normal * (normal * l.dir) / normal.length_sq();
            russian_rolette({p, n_dir}, ir, new_intensity, objs, kdt, depth + 1);
        } else if (dice < obj->rr[0] + obj->rr[1] + obj->rr[2]) { // refracao
            Vec3 n = normal.normalize();
            Vec3 dir = l.dir.normalize();
            double cosI = -dir * n;
            double eta = 1/obj->ir;

            if (cosI < 0) {
                cosI *= -1;
                n *= -1;
                eta = 1/eta;
            }

            double temp = 1 - (eta * eta) * (1 - cosI * cosI);

            if (temp > 0) {
                Vec3 n_dir = eta * dir - (sqrt(temp) - cosI * eta) * n;
                russian_rolette({p, n_dir}, obj->ir, new_intensity, objs, kdt, depth + 1);
            } else {
                Vec3 n_dir = dir - 2 * n * (n * dir);
                russian_rolette({p, n_dir}, ir, new_intensity, objs, kdt, depth + 1);
            }
        } else { // absorcao
            kdt.push_back({p, l.dir.normalize(), intensity});
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

        russian_rolette(l, 1.0, intensity, objs, kdt, 0);
    }

    return kdt;
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
            int x = round(v * cam.desl_h / cam.desl_h.length_sq());
            int y = round(v * cam.desl_v / cam.desl_v.length_sq());

            if (0 <= x && x < cam.hres && 0 <= y && y < cam.vres)
                cam.grid[x][y] = 1e6;
        }
    }
}

void visualize_radiance(const Camera &cam, const std::vector<Object*> &objs, const KDTree &kdt) {
    for (int i = 0; i < cam.vres; ++i)
        for (int j = 0; j < cam.hres; ++j) {
            Line l{cam.origin, cam.pixel_ray(i, j)};
            Intersection inter;

            for (const auto &obj : objs) {
                auto temp = obj->intersect(l);
                if (temp < inter)
                    inter = temp;
            }

            cam.grid[i][j] = inter ? kdt.get_intensity(l.t(inter.t), inter.normal.normalize(), 0.2) : Color(0);
        }
}

std::atomic_int32_t th_cnt;
std::vector<KDTree> th_kdts;

void emit_photons_th_aux(const Point3 &p, int num, Color intensity, const std::vector<Object*> &objs, int idx) {
    Line l;
    l.origin = p;
    Rand r(-1, 1);

    while (th_cnt.fetch_add(1, std::memory_order_acq_rel) < num) {
        do
            l.dir = {r.gen(), r.gen(), r.gen()};
        while (l.dir.length_sq() > 1);

        russian_rolette(l, 1.0, intensity, objs, th_kdts[idx], 0);
    }
}

KDTree emit_photons_th(const Point3 &p, int num, double power, const std::vector<Object*> &objs, int threads) {
    KDTree kdt;
    std::thread th[threads - 1];
    Color intensity = power / num;

    th_kdts.assign(threads, KDTree());
    th_cnt.store(0, std::memory_order_release);

    for (int i = 0; i < threads - 1; ++i)
        th[i] = std::thread(emit_photons_th_aux, p, num, intensity, objs, i);

    emit_photons_th_aux(p, num, intensity, objs, threads - 1);

    for (auto &t : th)
        t.join();

    for (int i = 0; i < threads; ++i) {
        kdt.insert(kdt.end(), th_kdts[i].begin(), th_kdts[i].end());
        th_kdts[i].clear();
    }

    th_kdts.clear();

    return kdt;
}
