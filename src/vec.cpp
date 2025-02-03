#include "vec.h"

Color::Color(const Color3 &rgb) {
    for (size_t i = 0; i < CN; ++i) {
        float wl = Color3::get_wavelength(i);
        const Color3 mappedRGB(wl);

        (*this)[i] = rgb.R * mappedRGB.R + rgb.G * mappedRGB.G + rgb.B * mappedRGB.B;
    }

    float sum = 0.0;
    for (float value : (*this))
        sum += value;

    if (sum > 1e-5)
        for (float& value : (*this))
            value /= sum;
}
