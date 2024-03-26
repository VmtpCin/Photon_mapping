#include "tracing.h"
#include "kdtree.h"
#include "geometry.h"
#include <mutex>
#include <atomic>
#include <thread>

void russian_rolette(bool refraction, const Line &l, Color intensity, std::vector<const Object*> &ir,
                     const std::vector<const Object*> &objs, KDTree &kdt, KDTree &kdt_refraction, int depth) {
    if (depth > 5 || !intensity.not_empty())
        return;

    Intersection inter;
    const Object *obj = nullptr;

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

            kdt.push_back({p, obj, l.dir, new_intensity / 2});

            russian_rolette(refraction, n_l, new_intensity / 2, ir, objs, kdt, kdt_refraction, depth + 1);
        } else if (dice < obj->rr[0] + obj->rr[1]) { // reflexao
            Vec3 n_dir = l.dir - 2 * normal * (normal * l.dir) / normal.length_sq();
            russian_rolette(refraction, {p, n_dir}, new_intensity, ir, objs, kdt, kdt_refraction, depth + 1);
        } else if (dice < obj->rr[0] + obj->rr[1] + obj->rr[2]) { // refracao
            Color white = std::min({intensity.R, intensity.G, intensity.B});
            Color dominant = intensity - white;

            Vec3 n = normal.normalize();
            double cosI = -l.dir * n;
            bool getting_in;

            if (cosI < 0) {
                getting_in = false;
                cosI *= -1;
                n *= -1;
            } else {
                getting_in = true;
            }

            if (dominant.not_empty()) {
                double eta;
                double wavelength = dominant.to_wavelength();

                if (getting_in)
                    eta = ir.back()->ir(wavelength)/obj->ir(wavelength);
                else
                    eta = obj->ir(wavelength)/ir[ir.size() - 2]->ir(wavelength);

                double temp = 1 - (eta * eta) * (1 - cosI * cosI);

                if (temp > 0) {
                    if (getting_in) ir.push_back(obj);
                    else ir.pop_back();

                    Vec3 n_dir = eta * l.dir - (sqrt(temp) - cosI * eta) * n;
                    russian_rolette(true, {p, n_dir}, dominant, ir, objs, kdt, kdt_refraction, depth + 1);

                    if (getting_in) ir.pop_back();
                    else ir.push_back(obj);
                } else {
                    Vec3 n_dir = l.dir - 2 * n * (n * l.dir);
                    russian_rolette(refraction, {p, n_dir}, dominant, ir, objs, kdt, kdt_refraction, depth + 1);
                }
            }

            if (white.not_empty()) {
                Rand r1(0, 2.99);
                Rand r2(0, 1);

                constexpr int iterations = 10;
                for (int i = 0; i < iterations; ++i) {
                    int c1, c2;

                    do
                        c1 = r1.gen(), c2 = r1.gen();
                    while (c1 == c2);

                    Color c = {c1 == 0 ? white.R : 0, c1 == 1 ? white.G : 0, c1 == 2 ? white.B : 0};
                    c += r2.gen() * Color({c2 == 0 ? white.R : 0, c2 == 1 ? white.G : 0, c2 == 2 ? white.B : 0});

                    double eta;
                    double wavelength = c.to_wavelength();

                    if (getting_in)
                        eta = ir.back()->ir(wavelength)/obj->ir(wavelength);
                    else
                        eta = obj->ir(wavelength)/ir[ir.size() - 2]->ir(wavelength);

                    double temp = 1 - (eta * eta) * (1 - cosI * cosI);

                    if (temp > 0) {
                        if (getting_in) ir.push_back(obj);
                        else ir.pop_back();

                        Vec3 n_dir = eta * l.dir - (sqrt(temp) - cosI * eta) * n;
                        russian_rolette(true, {p, n_dir}, 2 * c / iterations, ir, objs, kdt, kdt_refraction, depth + 1);

                        if (getting_in) ir.pop_back();
                        else ir.push_back(obj);
                    } else {
                        Vec3 n_dir = l.dir - 2 * n * (n * l.dir);
                        russian_rolette(refraction, {p, n_dir}, 2 * c / iterations, ir, objs, kdt, kdt_refraction, depth + 1);
                    }
                }
            }
        } else { // absorcao
            if (!refraction)  kdt.push_back({p, obj, l.dir, new_intensity});
            else   kdt_refraction.push_back({p, obj, l.dir, new_intensity});
        }
    }
}

std::pair<KDTree, KDTree> emit_photons(const std::vector<Light> &lights, const std::vector<const Object*> &objs) {
    KDTree kdt, kdt_refraction;

    for (const auto& light : lights) {
        Line l;
        l.origin = light.pos;
        Rand r(-1, 1);
        Light new_light = {light.pos, light.Intensity / light.quantity, light.quantity};
        std::vector<const Object*> ir;

        double rr[3] = {0, 0, 0};
        constexpr auto n_air = [](double wavelength) { return 1.0; };

        ir.push_back(new Object(nullptr, rr, n_air));

        for (int i = 0; i < light.quantity; ++i) {
            do
                l.dir = {r.gen(), r.gen(), r.gen()};
            while (l.dir.length_sq() > 1);

            l.dir = l.dir.normalize();

            russian_rolette(false, l, new_light.Intensity, ir, objs, kdt, kdt_refraction, 0);
        }

        delete ir[0];
    }
    return {kdt, kdt_refraction};
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
                cam.grid[x][y] += 3000 * photon.I;
        }
    }
}

void visualize_radiance(const Camera &cam, const std::vector<const Object*> &objs, const KDTree &kdt) {
    for (int i = 0; i < cam.vres; ++i)
        for (int j = 0; j < cam.hres; ++j) {
            Line l{cam.origin, cam.pixel_ray(i, j)};
            Intersection inter;
            const Object *obj = nullptr;

            for (const auto &o : objs) {
                auto temp = o->intersect(l);
                if (temp < inter)
                    inter = temp, obj = o;
            }

            cam.grid[i][j] = inter ? kdt.get_intensity(l.t(inter.t), obj, inter.normal.normalize(), 0.2) : Color(0);
        }
}

std::atomic_int32_t th_cnt;
std::vector<KDTree> th_kdts;
std::vector<KDTree> th_kdts_refraction;

void emit_photons_th_aux(const Light &light, const std::vector<const Object*> &objs, int idx) {
    Line l;
    l.origin = light.pos;
    Rand r(-1, 1);
    std::vector<const Object*> ir;

    double rr[3] = {0, 0, 0};
    constexpr auto n_air = [](double wavelength) { return 1.0; };

    ir.push_back(new Object(nullptr, rr, n_air));


    while (th_cnt.fetch_add(1, std::memory_order_acq_rel) < light.quantity) {
        do
            l.dir = {r.gen(), r.gen(), r.gen()};
        while (l.dir.length_sq() > 1);

        l.dir = l.dir.normalize();

        russian_rolette(false, l, light.Intensity, ir, objs, th_kdts[idx], th_kdts_refraction[idx], 0);
    }

    delete ir[0];
}

std::pair<KDTree, KDTree> emit_photons_th(const std::vector<Light> &lights, const std::vector<const Object*> &objs, int threads) {
    KDTree kdt, kdt_refraction;
    std::thread th[threads - 1];

    th_kdts.assign(threads, KDTree());
    th_kdts_refraction.assign(threads, KDTree());

    for (const auto &light : lights) {
        Light new_light = {light.pos, light.Intensity / light.quantity, light.quantity};
        th_cnt.store(0, std::memory_order_release);

        for (int i = 0; i < threads - 1; ++i)
            th[i] = std::thread(emit_photons_th_aux, new_light, objs, i);

        emit_photons_th_aux(new_light, objs, threads - 1);

        for (auto &t : th)
            t.join();
    }

    for (int i = 0; i < threads; ++i) {
        kdt.insert(kdt.end(), th_kdts[i].begin(), th_kdts[i].end());
        th_kdts[i].clear();

        kdt_refraction.insert(kdt_refraction.end(), th_kdts_refraction[i].begin(), th_kdts_refraction[i].end());
        th_kdts_refraction[i].clear();
    }

    th_kdts.clear();
    th_kdts_refraction.clear();

    return {kdt, kdt_refraction};
}
