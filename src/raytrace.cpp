#include "tracing.h"
#include "geometry.h"
#include <atomic>
#include <iostream>
#include <thread>

void simplecast(const Camera &cam, const std::vector<const Object*> &objs) {
    for (int i = 0; i < cam.vres; ++i)
        for (int j = 0; j < cam.hres; ++j) {
            Line l{cam.origin, cam.pixel_ray(i, j).normalize()};
            Intersection inter;
            const Object *obj;

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

Color castray(const Line &l, const std::vector<const Object*> &objs, const std::vector<Light> &lights, const KDTree &kdt,
              const KDTree &kdt_refraction, std::vector<const Object*> &ir, double wavelength = 0, int depth = 0) {
    if (depth > 5)
        return 0;

    Intersection inter;
    const Object *obj;
    Color result = 0;

    for (const auto &o : objs) {
        auto temp = o->intersect(l);
        if (temp < inter)
            inter = temp, obj = o;
    }

    if (inter) {
        Point3 p = l.t(inter.t);
        Vec3 &normal = inter.normal;

        for (const auto &light : lights) {
            const Line l_l = {p, light.pos - p};

            for (const auto &o : objs) {
                auto temp = o->intersect(l);
                if (temp < inter)
                    inter = temp, obj = o;
            }

            if ((!inter || inter > 1) && (obj->rr[0] || obj->rr[1])) {
                Color c = light.Intensity / (M_PI * l_l.dir.length_sq());

                if (obj->rr[0])
                    result += c * obj->rr[0] * (l_l.dir.normalize() * normal);

                // if (obj->rr[1]) {
                //     Vec3 n_dir = l.dir - 2 * normal * (normal * l.dir);
                //     result += c * obj->rr[1] * (l_l.dir.normalize() * n_dir);
                // }
            }
        }

        // reflexão
        if (obj->rr[1] > 0) {
            Vec3 n_dir = l.dir - 2 * normal * (normal * l.dir) / normal.length_sq();
            result += obj->rr[1] * castray({p, n_dir}, objs, lights, kdt, kdt_refraction, ir, wavelength, depth + 1);
        }

        // refração
        if (obj->rr[2] > 0) {            
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

            if (wavelength) {
                double eta;
                if (getting_in)
                    eta = ir.back()->ir(wavelength)/obj->ir(wavelength);
                else
                    eta = obj->ir(wavelength)/ir[ir.size() - 2]->ir(wavelength);


                double temp = 1 - (eta * eta) * (1 - cosI * cosI);

                if (temp > 0) {
                    if (getting_in) ir.push_back(obj);
                    else ir.pop_back();

                    Vec3 n_dir = eta * l.dir - (sqrt(temp) - cosI * eta) * n;
                    result += obj->rr[2] * castray({p, n_dir}, objs, lights, kdt, kdt_refraction, ir, wavelength, depth + 1);

                    if (getting_in) ir.pop_back();
                    else ir.push_back(obj);
                } else {
                    Vec3 n_dir = l.dir - 2 * n * (n * l.dir);
                    result += obj->rr[2] * castray({p, n_dir}, objs, lights, kdt, kdt_refraction, ir, wavelength, depth + 1);
                }
            } else {
                constexpr int iterations = 2;
                for (int i = 0; i < iterations; ++i) {
                    double eta;
                    double wl = iterations == 1 ? 565 : 380 + 370 * i / (iterations - 1);

                    if (getting_in)
                        eta = ir.back()->ir(wl)/obj->ir(wl);
                    else
                        eta = obj->ir(wl)/ir[ir.size() - 2]->ir(wl);

                    double temp = 1 - (eta * eta) * (1 - cosI * cosI);

                    if (temp > 0) {
                        if (getting_in) ir.push_back(obj);
                        else ir.pop_back();

                        Vec3 n_dir = eta * l.dir - (sqrt(temp) - cosI * eta) * n;
                        result += obj->rr[2] * castray({p, n_dir}, objs, lights, kdt, kdt_refraction, ir, wl, depth + 1) / iterations;

                        if (getting_in) ir.pop_back();
                        else ir.push_back(obj);
                    } else {
                        Vec3 n_dir = l.dir - 2 * n * (n * l.dir);
                        result += obj->rr[2] * castray({p, n_dir}, objs, lights, kdt, kdt_refraction, ir, wavelength, depth + 1) / iterations;
                    }
                }
            }
        }

        // difusão
        if (obj->rr[1] + obj->rr[2] < 1)
            result += (1 - (obj->rr[1] + obj->rr[2]))
                    * (              kdt.get_intensity(l.t(inter.t), obj, inter.normal.normalize(), 0.15)
                        + kdt_refraction.get_intensity(l.t(inter.t), obj, inter.normal.normalize(), 0.04));
    }

    return result;
}

void raycast(const Camera &cam, const std::vector<const Object*> &objs, const std::vector<Light> &lights, const KDTree &kdt, const KDTree &kdt_refraction) {
    double rr[3] = {0, 0, 0};
    constexpr auto n_air = [](double wavelength) { return 1.0; };
    std::vector<const Object*> ir;
    ir.push_back(new Object(nullptr, rr, n_air));

    for (int i = 0; i < cam.vres; ++i)
        for (int j = 0; j < cam.hres; ++j)
            cam.grid[i][j] = castray({cam.origin, cam.pixel_ray(i, j).normalize()}, objs, lights, kdt, kdt_refraction, ir);
}

std::atomic_int32_t it_i;
std::vector<std::vector<Color>> th_rc_aux; 

void raycast_th_aux(const Camera &cam, const std::vector<const Object*> &objs, const std::vector<Light> &lights, const KDTree &kdt, const KDTree &kdt_refraction) {
    int i;
    double rr[3] = {0, 0, 0};
    constexpr auto n_air = [](double wavelength) { return 1.0; };
    std::vector<const Object*> ir;
    ir.push_back(new Object(nullptr, rr, n_air));

    while((i = it_i.fetch_add(1, std::memory_order_acq_rel)) < cam.vres)
        for (int j = 0; j < cam.hres; ++j)
            th_rc_aux[i][j] = castray({cam.origin, cam.pixel_ray(i, j).normalize()}, objs, lights, kdt, kdt_refraction, ir);
}

void raycast_th(const Camera &cam, const std::vector<const Object*> &objs, const std::vector<Light> &lights, const KDTree &kdt, const KDTree &kdt_refraction, int threads) {
    th_rc_aux.swap(cam.grid);

    std::thread th[threads - 1];

    it_i.store(0, std::memory_order_release);

    for (int i = 0; i < threads - 1; ++i)
        th[i] = std::thread(raycast_th_aux, cam, objs, lights, kdt, kdt_refraction);

    raycast_th_aux(cam, objs, lights, kdt, kdt_refraction);

    for (auto &t : th)
        t.join();

    th_rc_aux.swap(cam.grid);
}
