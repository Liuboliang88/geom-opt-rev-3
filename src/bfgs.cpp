#include "bfgs.hpp"
#include "matrix.hpp"
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cassert>

Matrix bfgs_step(const Matrix& x, const Matrix& g, const Matrix& H) {
    return x - H.inver() % g;
}

Matrix bfgs_step_modified(const Matrix& x, const Matrix& g, const Matrix& H) {
    SelfAdjointEigenSolver es(H);
    Matrix ee = es.eigen_val();
    Matrix ev = es.eigen_vec();

    Matrix invE(H.rows(), H.cols());
    for (std::size_t i = 0; i < H.rows(); ++i) {
        if (std::abs(ee(i) < 1e-6)) {
            std::cout << "ee too small." << std::endl;
            std::abort();
        }

        invE(i, i) = 1.0 / ee(i);
        if (invE(i,i) < 0.0) invE(i, i) *= -1.0;
    }

    Matrix invH = ev % invE % ev.trans();
    return  x - invH % g;
}