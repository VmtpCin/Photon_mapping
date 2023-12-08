#include "surface.h"

Point3 bezier_curve(const std::vector<Point3> &points, double t) {
    std::vector<Point3> control;

    control.assign(points.begin(), points.end());

    while (control.size() > 1) {
        for (int i = 0; i < control.size() - 1; i++)
            control[i] = interpolate(control[i], control[i + 1], t);

        control.pop_back();
    }

    return control[0];
}

Point3 bezier_triangle(const std::vector<std::vector<Point3>> &points, double s, double t) {
    std::vector<std::vector<Point3>> control;

    for (const auto &p : points)
        control.push_back(std::vector<Point3>(p.begin(), p.end()));

    while (control.size() > 1) {
        for (int i = 0; i < control.size() - 1; i++)
            for (int j = 0; j < control[i].size(); j++)
                control[i][j] = interpolate(control[i][j], control[i + 1][j], control[i + 1][j + 1], s, t);

        control.pop_back();
    }

    return control[0][0];
}

Point3 bezier_triangle(double s, double t) {
    std::vector<std::vector<Point3>> v;
    v.push_back({{1, 0, 0}});
    v.push_back({{1, 1, 0}, {1, 1, 1}});
    v.push_back({{1, 0, 2}, {1, 2, 2}, {1, 2, 3}});

    return bezier_triangle(v, s, t);
}

std::vector<Point3> create_surface_triangle(Point3 (*f)(double, double),
                    int s_step, int t_step, std::vector<Object*> &objs) {
    std::vector<Point3> prev_curve, curve;

    for (int it_s = 0; it_s <= s_step; ++it_s) {
        for (int it_t = 0; it_t <= t_step; ++it_t) {
            double s = double(it_s) / s_step;
            double t = double(it_t) / t_step;

            if (s + t > 1)
                break;

            curve.push_back(f(s, t));
        }

        const double trans[3] = {0, 0.2, 0.8};
        for (int i = 0; i + 1 < prev_curve.size(); ++i) {
            objs.push_back(new Object(new Triangle(prev_curve[i],          curve[i], curve[i + 1]), trans, 1.3));
            objs.push_back(new Object(new Triangle(prev_curve[i], prev_curve[i + 1], curve[i + 1]), trans, 1.3));
        }

        prev_curve.swap(curve);
        curve.clear();
    }

    return curve;
}

std::vector<Point3> create_bezier_superfice(const std::vector<std::vector<Point3>> &control,
                                            int s_step, int t_step, std::vector<Object*> &objs) {
    std::vector<Point3> prev_curve, curve;

    for (int it_s; it_s <= s_step; ++it_s) {
        std::vector<Point3> next_it;
        double s = double(it_s) / s_step;

        for (auto &c : control)
            next_it.push_back(bezier_curve(c, s));

        for (int it_t = 0; it_t <= t_step; ++it_t)
            curve.push_back(bezier_curve(next_it, double(it_t) / t_step));

        const double trans[3] = {0, 0.2, 0.8};
        for (int i = 0; i + 1 < prev_curve.size(); ++i) {
            objs.push_back(new Object(new Triangle(prev_curve[i],          curve[i], curve[i + 1]), trans, 1.3));
            objs.push_back(new Object(new Triangle(prev_curve[i], prev_curve[i + 1], curve[i + 1]), trans, 1.3));
        }

        prev_curve.swap(curve);
        curve.clear();
    }

    return curve;
}
