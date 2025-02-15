#include "matrix.h"
#include "vec.h"

Matrix<float> Color::C;

void setupColor() {
    Matrix<float> M(3, CN);

    for (size_t i = 0; i < CN; ++i) {
        Color3 c(Color::get_wavelength(i));

        M[0][i] = c.R;
        M[1][i] = c.G;
        M[2][i] = c.B;
    }

    Matrix<float> M_t = M.transpose();

    Color::C = M_t * (M * M_t).inverse();
}

Color::Color(const Color3 &rgb) {
    constexpr auto calculate = [&](Matrix<float> m) {
        Matrix<float> Color(3, 1);
        for (size_t i = 0; i < CN; ++i) {
            if (m[i][0] < 0) continue;
            Color3 c(get_wavelength(i));
    
            Color[0][0] += m[i][0] * c.R;
            Color[1][0] += m[i][0] * c.G;
            Color[2][0] += m[i][0] * c.B;
        }

        return Color;
    };

    const Matrix<float> goal({{rgb.R}, {rgb.G}, {rgb.B}});    
    Matrix<float> M = C * goal, newgoal;

    size_t cnt = 0;
    do {
        newgoal = goal - calculate(M);
        M += C * newgoal;

        for (size_t i = 0; i < CN; ++ i)
            if (M[i][0] < 0)
                M[i][0] = 0;

    } while (++cnt < 1e4 && newgoal.sumValuesAbsolute() >= 5e-3);

    for (size_t i = 0; i < CN; ++i)
        channels[i] = M[i][0];
}
