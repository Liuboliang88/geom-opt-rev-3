#pragma once

#include "matrix.hpp"

class Hessian {
public:
    Hessian() = delete;
    Hessian(std::size_t dim): dim(dim), stp(0) {}
    ~Hessian() {}

    Matrix hessian() const;

    bool update_hessian_bfgs(const Matrix& x, const Matrix& g);

    // GDIIS: https://doi.org/10.1039/B108658H
    bool update_hessian_gdiis(const Matrix& x, const Matrix& g);

private:
    std::size_t dim = 0;
    std::size_t stp = 0;
    
    Matrix x0, g0, h0;
    Matrix x1, g1, h1;

    void _update_hessian_init(const Matrix& x, const Matrix& g);

    Matrix _delta_hessian_bfgs(const Matrix& sk, const Matrix& yk, const Matrix Hk) const;
    Matrix _delta_hessian_sr1(const Matrix& sk, const Matrix& yk, const Matrix Hk) const;
    Matrix _delta_hessian_psb(const Matrix& sk, const Matrix& yk, const Matrix Hk) const;
    double _phi_bofill(const Matrix& sk, const Matrix& yk, const Matrix Hk) const;
};
