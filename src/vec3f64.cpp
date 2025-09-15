#include "vec3f64.hpp"
#include <string>
#include <sstream>
#include <cstdio>
#include <cmath>
#include <cassert>
#include <cstdlib>
#include <iostream>

Vec3f64::Vec3f64() : x(0.0), y(0.0), z(0.0) {}

Vec3f64::Vec3f64(double x, double y, double z)
: x(x), y(y), z(z) {}

Vec3f64::Vec3f64(const std::string &input) {
    std::istringstream iss(input);
    iss >> x >> y >> z;
}

double Vec3f64::len() const
{ return std::sqrt(x*x + y*y + z*z); }

double Vec3f64::len2() const
{ return x*x + y*y + z*z; }

std::string Vec3f64::to_string() const {
    char buff[128];
    std::sprintf(buff, " %16.8lf %16.8lf %16.8lf", x, y, z);
    return std::string(buff);
}

/* access operator */
double  Vec3f64::operator[](std::size_t i) const {
    assert(i < 3);
    switch (i) {
        case 0: return x;
        case 1: return y;
        case 2: return z;
        default: std::abort(); // shouldn't be reached
    }

    return 0.0;
}

double& Vec3f64::operator[](std::size_t i) {
    assert(i < 3);
    switch (i) {
        case 0: return x;
        case 1: return y;
        case 2: return z;
        default: std::abort(); // shouldn't be reached
    }

    return x;
}

/* same as operator[] */
double Vec3f64::at(std::size_t i) const {
    if (i >= 3) {
        std::cerr << "Vec3f64::at(std::size_t i) error:" << std::endl;
        std::cerr << "index out of range: i = " << i << std::endl;
        std::abort();
    }

    return (*this)[i];
}

double& Vec3f64::at(std::size_t i) {
    if (i >= 3) {
        std::cerr << "Vec3f64::at(std::size_t i) error:" << std::endl;
        std::cerr << "index out of range: i = " << i << std::endl;
        std::abort();
    }

    return (*this)[i];
}

Vec3f64& Vec3f64::operator+=(const Vec3f64 &vec)
{ x += vec.x; y += vec.y; z += vec.z; return *this; }

Vec3f64& Vec3f64::operator-=(const Vec3f64 &vec)
{ x -= vec.x; y -= vec.y; z -= vec.z; return *this; }

Vec3f64& Vec3f64::operator*=(const Vec3f64 &vec)
{ x *= vec.x; y *= vec.y; z *= vec.z; return *this; }

Vec3f64& Vec3f64::operator/=(const Vec3f64 &vec)
{ x /= vec.x; y /= vec.y; z /= vec.z; return *this; }

Vec3f64&  Vec3f64::operator+=(double val)
{ x += val; y += val; z += val; return *this; }

Vec3f64&  Vec3f64::operator-=(double val)
{ x -= val; y -= val; z -= val; return *this; }

Vec3f64&  Vec3f64::operator*=(double val)
{ x *= val; y *= val; z *= val; return *this; }

Vec3f64&  Vec3f64::operator/=(double val)
{ x /= val; y /= val; z /= val; return *this; }


/*  operators */
Vec3f64 operator+(const Vec3f64 &vec)
{ return Vec3f64(+vec.x, +vec.y, +vec.z); }

Vec3f64 operator-(const Vec3f64 &vec)
{ return Vec3f64(-vec.x, -vec.y, -vec.z); }


Vec3f64  operator+(const Vec3f64 &vec1, const Vec3f64 &vec2)
{ return Vec3f64(vec1.x + vec2.x, vec1.y + vec2.y, vec1.z + vec2.z); }

Vec3f64  operator-(const Vec3f64 &vec1, const Vec3f64 &vec2)
{ return Vec3f64(vec1.x - vec2.x, vec1.y - vec2.y, vec1.z - vec2.z); }

Vec3f64  operator*(const Vec3f64 &vec1, const Vec3f64 &vec2)
{ return Vec3f64(vec1.x * vec2.x, vec1.y * vec2.y, vec1.z * vec2.z); }

Vec3f64  operator/(const Vec3f64 &vec1, const Vec3f64 &vec2)
{ return Vec3f64(vec1.x / vec2.x, vec1.y / vec2.y, vec1.z / vec2.z); }

Vec3f64  operator+(const Vec3f64 &vec, double val)
{ return Vec3f64(vec.x + val, vec.y + val, vec.z + val); }

Vec3f64  operator-(const Vec3f64 &vec, double val)
{ return Vec3f64(vec.x - val, vec.y - val, vec.z - val); }

Vec3f64  operator*(const Vec3f64 &vec, double val)
{ return Vec3f64(vec.x * val, vec.y * val, vec.z * val); }

Vec3f64  operator/(const Vec3f64 &vec, double val)
{ return Vec3f64(vec.x / val, vec.y / val, vec.z / val); }

Vec3f64  operator+(double val, const Vec3f64 &vec)
{ return Vec3f64(val + vec.x, val + vec.y, val + vec.z); }

Vec3f64  operator-(double val, const Vec3f64 &vec)
{ return Vec3f64(val - vec.x, val - vec.y, val - vec.z); }

Vec3f64  operator*(double val, const Vec3f64 &vec)
{ return Vec3f64(val * vec.x, val * vec.y, val * vec.z); }

Vec3f64  operator/(double val, const Vec3f64 &vec)
{ return Vec3f64(val / vec.x, val / vec.y, val / vec.z); }

/* dot product */
double dot(const Vec3f64 &vec1, const Vec3f64 &vec2)
{ return vec1.x * vec2.x + vec1.y * vec2.y + vec1.z * vec2.z; }

/* cross product */
Vec3f64 crs(const Vec3f64 &vec1, const Vec3f64 &vec2) {
    return Vec3f64(
        vec1.y * vec2.z - vec1.z * vec2.y,
        vec1.z * vec2.x - vec1.x * vec2.z,
        vec1.x * vec2.y - vec1.y * vec2.x
    );
}

