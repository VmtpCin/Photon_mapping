#pragma once
#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <string>
#include <math.h>

constexpr size_t CN = 32;

struct Vec3 {
    double e[3];

    double& x()       { return e[0]; }
    double  x() const { return e[0]; }
    double& y()       { return e[1]; }
    double  y() const { return e[1]; }
    double& z()       { return e[2]; }
    double  z() const { return e[2]; }

    double& operator[](int i)       { return e[i]; }
    double  operator[](int i) const { return e[i]; }

    operator std::string() const {
        return     "(" + std::to_string(x())
                + ", " + std::to_string(y())
                + ", " + std::to_string(z()) + ")";
    }

    friend std::ostream& operator<<(std::ostream &os, const Vec3 &v) {
        return os << std::string(v);
    }

    Vec3 operator-() const {
        return {-e[0], -e[1], -e[2]};
    }

    Vec3 operator+(const Vec3 &v) const {
        return {e[0] + v[0], e[1] + v[1], e[2] + v[2]};
    }

    Vec3 operator-(const Vec3 &v) const {
        return {e[0] - v[0], e[1] - v[1], e[2] - v[2]};
    }

    Vec3 operator*(const double d) const {
        return {e[0] * d, e[1] * d, e[2] * d};
    }

    Vec3 operator/(const double d) const {
        return {e[0] / d, e[1] / d, e[2] / d};
    }

    Vec3& operator+=(const Vec3 &v) {
        e[0] += v[0], e[1] += v[1], e[2] += v[2];
        return *this;
    }

    Vec3& operator-=(const Vec3 &v) {
        e[0] -= v[0], e[1] -= v[1], e[2] -= v[2];
        return *this;
    }

    Vec3& operator*=(const double &d) {
        e[0] *= d, e[1] *= d, e[2] *= d;
        return *this;
    }

    Vec3& operator/=(const double &d) {
        e[0] /= d, e[1] /= d, e[2] /= d;
        return *this;
    }

    double operator*(const Vec3 &v) const {
        return e[0] * v[0] + e[1] * v[1] + e[2] * v[2];
    }

    Vec3 operator^(const Vec3 &v) const {
        return { e[1] * v[2] - e[2] * v[1],
                 e[2] * v[0] - e[0] * v[2],
                 e[0] * v[1] - e[1] * v[0]};
    }

    Vec3 operator&(const Vec3 &v) const {
        return {e[0] * v[0], e[1] * v[1], e[2] * v[2]};
    }

    Vec3& operator&=(const Vec3 &v) {
        e[0] *= v[0], e[1] *= v[1], e[2] *= v[2];
        return *this;
    }

    double length_sq() const {
        return e[0] * e[0] + e[1] * e[1] + e[2] * e[2];
    }

    double length() const {
        return sqrt(length_sq());
    }

    Vec3 normalize() const {
        return *this / length();
    }

    Vec3 perpendicular() const {
        if (e[0] != 0)
            return {-e[1], e[0], 0};
        else if (e[1] != 0)
            return {0, -e[2], e[1]};
        else if (e[2] != 0)
            return {e[2], 0, -e[0]};
        else
            return {0, 0, 0};
    }
};

struct Point3 {
    double e[3];

    explicit operator Vec3() const {
        return {e[0], e[1], e[2]};
    }

    double& x()       { return e[0]; }
    double  x() const { return e[0]; }
    double& y()       { return e[1]; }
    double  y() const { return e[1]; }
    double& z()       { return e[2]; }
    double  z() const { return e[2]; }

    double& operator[](int i)       { return e[i]; }
    double  operator[](int i) const { return e[i]; }

    operator std::string() const {
        return     "(" + std::to_string(x())
                + ", " + std::to_string(y())
                + ", " + std::to_string(z()) + ")";
    }

    friend std::ostream& operator<<(std::ostream &os, const Point3 &p) {
        return os << std::string(p);
    }

    Point3 operator+(const Vec3 &v) const {
        return {e[0] + v[0], e[1] + v[1], e[2] + v[2]}; 
    }

    Point3& operator+=(const Vec3 &v) {
        e[0] += v[0], e[1] += v[1], e[2] += v[2];
        return *this;
    }

    Vec3 operator-(const Point3 &p) const {
        return {e[0] - p[0], e[1] - p[1], e[2] - p[2]};
    }

    double distance_sq(const Point3 &p) const {
        double dx = e[0] - p[0], dy = e[1] - p[1], dz = e[2] - p[2];
        return dx * dx + dy * dy + dz * dz;
    }

    double distance(const Point3 &p) const {
        return sqrt(distance_sq(p));
    }
};

struct Line {
    Point3 origin;
    Vec3 dir;

    Point3 t(double t) const {
        return origin + (dir * t);
    }
};

struct Color {
    float channels[CN];

    float* begin() { return channels; }
    const float* begin() const { return channels; }

    float* end() { return channels + CN; }
    const float* end() const { return channels + CN; }

    Color(int i = 0) {
        for (float &c : channels)
            c = static_cast<float>(i);
    }

    Color(std::initializer_list<float> values) {
        if (values.size() == 1) {
            for (float &c : channels)
                c = *values.begin();

            return;
        }

        if (values.size() != CN)
            throw std::invalid_argument("Number of arguments (" + std::to_string(values.size()) + 
                                        ") must match the size of the Color (" + std::to_string(CN) + ").");

        size_t i = 0;
        for (const float &v : values)
            channels[i++] = v;
    }

    explicit Color(const struct Color3& c);

    float& operator[](size_t i)       { return channels[i]; }
    float  operator[](size_t i) const { return channels[i]; }

    Color operator&(const Color &c) const {
        Color result;

        std::transform(channels, channels + CN, c.channels, result.channels,
                   [](float a, float b) { return a * b; });

        return result;
    }

    Color operator+(const Color &c) const {
        Color result;

        std::transform(channels, channels + CN, c.channels, result.channels,
                   [](float a, float b) { return a + b; });

        return result;
    }

    Color operator-(const Color &c) const {
        Color result;

        std::transform(channels, channels + CN, c.channels, result.channels,
                   [](float a, float b) { return a - b; });

        return result;
    }

    Color operator*(const float d) const {
        Color result;

        for (size_t i = 0; i < CN; ++i)
            result[i] = (*this)[i] * d;

        return result;
    }

    Color operator/(const double d) const {
        Color result;

        for (size_t i = 0; i < CN; ++i)
            result[i] = (*this)[i] / d;

        return result;
    }

    Color& operator&=(const Color& c) {
        for (size_t i = 0; i < CN; ++i)
            (*this)[i] *= c[i];
        return *this;
    }

    Color& operator+=(const Color &c) {
        for (size_t i = 0; i < CN; ++i)
            (*this)[i] += c[i];
        return *this;
    }

    Color& operator*=(const double d) {
        for (float &c : channels)
            c *= d;
        return *this;
    }

    Color& operator/=(const double d) {
        for (float &c : channels)
            c /= d;
        return *this;
    }

    bool not_empty() const {
        for (const float &c : channels)
            if (c > 1e-10)
                return true;

        return false;
    }
};

struct Color3 {
    float R = 0, G = 0, B = 0;

    Color3(float r, float g, float b) : R(r), G(g), B(b) {}

    Color3 (float wl) {
        if (wl >= 380 && wl < 440) {
            R = -(wl - 440) / (440 - 380);
            G = 0.0;
            B = 1.0;
        } else if (wl >= 440 && wl < 490) {
            R = 0.0;
            G = (wl - 440) / (490 - 440);
            B = 1.0;
        } else if (wl >= 490 && wl < 510) {
            R = 0.0;
            G = 1.0;
            B = -(wl - 510) / (510 - 490);
        } else if (wl >= 510 && wl < 580) {
            R = (wl - 510) / (580 - 510);
            G = 1.0;
            B = 0.0;
        } else if (wl >= 580 && wl < 645) {
            R = 1.0;
            G = -(wl - 645) / (645 - 580);
            B = 0.0;
        } else if (wl >= 645 && wl <= 750) {
            R = 1.0;
            G = 0.0;
            B = 0.0;
        }

        float factor = 1.0;
        if (wl >= 380 && wl < 420)
            factor = 0.3 + 0.7 * (wl - 380) / (420 - 380);
        else if (wl >= 645 && wl <= 750)
            factor = 0.3 + 0.7 * (750 - wl) / (750 - 645);

        *this *= factor;
    }

    explicit Color3(const Color &c) {
        for (size_t i = 0; i < CN; ++i) {
            if (c[i] <= 1e-5) continue;

            const float wl = get_wavelength(i);
            *this += Color3(wl) * c[i];
        }

        R /= 0.510f;
        G /= 0.523f;
        B /= 0.326f;
    }

    static constexpr float get_wavelength(size_t idx) {
        constexpr float minWavelength = 380.0;
        constexpr float maxWavelength = 750.0;
        if (idx >= CN)
            throw std::out_of_range("Index out of range");

        constexpr float segmentSize = (maxWavelength - minWavelength) / CN;

        return minWavelength + idx * segmentSize;
    }

    Color3 operator+(const Color3 &c) const {
        return {R + c.R, G + c.G, B + c.B};
    }

    Color3 operator*(float f) const {
        return {R * f, G * f, B * f};
    }

    Color3& operator+=(const Color3 &c) {
        R += c.R, G += c.G, B += c.B;
        return *this;
    }

    Color3& operator*=(float f) {
        R *= f, G *= f, B *= f;
        return *this;
    }
};

inline Point3 operator+(const Vec3& v, const Point3 &p) {
    return p + v;
}

inline Vec3 operator*(const double &d, const Vec3 &v) {
    return v * d;
}

inline Color operator*(const double &d, const Color &c) {
    return c * d;
}

extern Point3 interpolate(const Point3 &p1, const Point3 &p2, const Point3 &p3, double s, double t);
extern Point3 interpolate(const Point3 &p1, const Point3 &p2, double s);

