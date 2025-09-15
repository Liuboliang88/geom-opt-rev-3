#pragma once

#include <string>
#include <cstddef>

// Vec3d is a class for representing
// a vector or a point in 3D space.
class Vec3f64 {
public:
    double x, y, z;

    explicit Vec3f64();
    explicit Vec3f64(const std::string &input);
    explicit Vec3f64(double x, double y, double z);

    double len() const;
    double len2() const;
    std::string to_string() const;

    /* access operator */
    double  operator[](std::size_t i) const;
    double& operator[](std::size_t i);

    /* runtime range check */
    double  at(std::size_t i) const;
    double& at(std::size_t i);

    Vec3f64&  operator+=(const Vec3f64 &vec);
    Vec3f64&  operator-=(const Vec3f64 &vec);
    Vec3f64&  operator*=(const Vec3f64 &vec);
    Vec3f64&  operator/=(const Vec3f64 &vec);
    Vec3f64&  operator+=(double val);
    Vec3f64&  operator-=(double val);
    Vec3f64&  operator*=(double val);
    Vec3f64&  operator/=(double val);
};

Vec3f64  operator+(const Vec3f64 &vec); // positive
Vec3f64  operator-(const Vec3f64 &vec); // negative
Vec3f64  operator+(const Vec3f64 &vec1, const Vec3f64 &vec2);
Vec3f64  operator-(const Vec3f64 &vec1, const Vec3f64 &vec2);
Vec3f64  operator*(const Vec3f64 &vec1, const Vec3f64 &vec2);
Vec3f64  operator/(const Vec3f64 &vec1, const Vec3f64 &vec2);
Vec3f64  operator+(const Vec3f64 &vec, double val);
Vec3f64  operator-(const Vec3f64 &vec, double val);
Vec3f64  operator*(const Vec3f64 &vec, double val);
Vec3f64  operator/(const Vec3f64 &vec, double val);
Vec3f64  operator+(double val, const Vec3f64 &vec);
Vec3f64  operator-(double val, const Vec3f64 &vec);
Vec3f64  operator*(double val, const Vec3f64 &vec);
Vec3f64  operator/(double val, const Vec3f64 &vec);

/* dot product */
double  dot(const Vec3f64 &vec1, const Vec3f64 &vec2);

/* cross product */
Vec3f64 crs(const Vec3f64 &vec1, const Vec3f64 &vec2);
