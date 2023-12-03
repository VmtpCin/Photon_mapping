#include "geometry.h"

extern Point3 bezier_curve(const std::vector<Point3> &points, double t);
extern Point3 bezier_triangle(double s, double t);

extern std::vector<Point3> create_surface_triangle(Point3 (*f)(double, double),
                            int s_step, int t_step, std::vector<Object*> &objs);

extern std::vector<Point3> create_bezier_superfice(const std::vector<std::vector<Point3>> &control,
                                                   int s_step, int t_step, std::vector<Object*> &objs);
