#include "hessian.hpp"
#include <iostream>
#include <cstdlib>
#include <cassert>
#include <cmath>

Matrix Hessian::hessian() const {
    if (stp == 0 || h1.size() == 0) {
        std::cerr << "Matrix Hessian::hessian() const" << std::endl;
        std::cerr << "error: hessian matrix is not initialized." << std::endl;
        std::abort();
    }

    return h1;
}

bool Hessian::update_hessian_bfgs(const Matrix& x, const Matrix& g) {
    // check dimension
    if(x.rows() != dim || x.cols() != 1) {
        std::cerr << "bool Hessian::update_hessian_bfgs(const Matrix& x, const Matrix& g)" << std::endl;
        std::cerr << "error: x.rows() = " << x.rows() << ", x.cols() = " << x.cols() << ". dim = " << dim << std::endl;
        std::abort();
    }

    if(g.rows() != dim || g.cols() != 1) {
        std::cerr << "bool Hessian::update_hessian_bfgs(const Matrix& x, const Matrix& g)" << std::endl;
        std::cerr << "error: g.rows() = " << g.rows() << ", g.cols() = " << g.cols() << ". dim = " << dim << std::endl;
        std::abort();
    }

    if (stp == 0) {
        _update_hessian_init(x, g);
        ++stp;
        return true;
    }

    x1 = x;
    g1 = g;

    Matrix sk = x1 - x0;
    Matrix yk = g1 - g0;

    // for geom opt, 1e-4 is ok
    if (sk.norm() < 1e-4) {
        assert(yk.norm() < 1e-2);
        std::cout << "bool Hessian::update_hessian_bfgs(const Matrix& x, const Matrix& g)" << std::endl;
        std::cout << "warning: points are too close. Hessian is not updated." << std::endl;

        ++stp;
        return false;
    }


    h1 = h0 + _delta_hessian_bfgs(sk, yk, h0);

    x0 = x1;
    g0 = g1;
    h0 = h1;

    ++stp;
    return true;
}

bool Hessian::update_hessian_gdiis(const Matrix& x, const Matrix& g) {
    // check dimension
    if(x.rows() != dim || x.cols() != 1) {
        std::cerr << "bool Hessian::update_hessian_bofill(const Matrix& x, const Matrix& g)" << std::endl;
        std::cerr << "error: x.rows() = " << x.rows() << ", x.cols() = " << x.cols() << ". dim = " << dim << std::endl;
        std::abort();
    }

    if(g.rows() != dim || g.cols() != 1) {
        std::cerr << "bool Hessian::update_hessian_bofill(const Matrix& x, const Matrix& g)" << std::endl;
        std::cerr << "error: g.rows() = " << g.rows() << ", g.cols() = " << g.cols() << ". dim = " << dim << std::endl;
        std::abort();
    }

    if (stp == 0) {
        _update_hessian_init(x, g);
        ++stp;
        return true;
    }

    x1 = x;
    g1 = g;

    Matrix sk = x1 - x0;
    Matrix yk = g1 - g0;

    // for geom opt, 1e-4 is ok
    if (sk.norm() < 1e-4) {
        assert(yk.norm() < 1e-2);
        std::cout << "bool Hessian::update_hessian_gdiis(const Matrix& x, const Matrix& g)" << std::endl;
        std::cout << "warning: points are too close. Hessian is not updated." << std::endl;

        ++stp;
        return false;
    }

    Matrix dH_BFGS = _delta_hessian_bfgs(sk, yk, h0);
    Matrix dH_SR1 = _delta_hessian_sr1(sk, yk, h0);
    double phi_Bofill = _phi_bofill(sk, yk, h0);
    double phi = std::sqrt(phi_Bofill);

    h1 = h0 + phi * dH_SR1 + (1.0 - phi) * dH_BFGS;

    x0 = x1;
    g0 = g1;
    h0 = h1;

    ++stp;
    return true;
}

void Hessian::_update_hessian_init(const Matrix& x, const Matrix& g) {
    x0 = x1 = x;
    g0 = g1 = g;
    h0 = h1 = Matrix(dim, dim, 0.0);

    for (std::size_t i = 0; i < dim; ++i) {
        h0(i,i) = h1(i,i) = 1.0;
    }
}

Matrix Hessian::_delta_hessian_bfgs(const Matrix& sk, const Matrix& yk, const Matrix Hk) const {
    double yts = (yk.trans() % sk)(0);
    double stHs = (sk.trans() % Hk % sk)(0);
    Matrix Hksk = Hk % sk;
    return (yk % yk.trans()) / yts - (Hksk % Hksk.trans()) / stHs;
}

Matrix Hessian::_delta_hessian_sr1(const Matrix& sk, const Matrix& yk, const Matrix Hk) const {
    Matrix v = Hk % sk - yk;
    double vts = (v.trans() % sk)(0);
    return - v % v.trans() / vts;
}

Matrix Hessian::_delta_hessian_psb(const Matrix& sk, const Matrix& yk, const Matrix Hk) const {
    Matrix v = Hk % sk - yk;
    double vts = (v.trans() % sk)(0);
    double sts = (sk.trans() % sk)(0);

    Matrix sst = sk % sk.trans();
    Matrix vst = v % sk.trans();
    return vts * sst / (sts * sts) - (vst + vst.trans()) / sts;
}

double Hessian::_phi_bofill(const Matrix& sk, const Matrix& yk, const Matrix Hk) const{
    Matrix v = Hk % sk - yk;
    double vts = (v.trans() % sk)(0);
    double vtv = (v.trans() % v)(0);
    double sts = (sk.trans() % sk)(0);
    return vts * vts / (vtv * sts);
}