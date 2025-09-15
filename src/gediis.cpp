#include "gediis.hpp"
#include <cassert>
#include <cmath>
#include <iostream>
#include <cstdlib>

GEDIIS::GEDIIS(std::size_t dim, std::size_t minSpace, std::size_t maxSpace)
: dim(dim), minSpace(minSpace), maxSpace(maxSpace) {}

bool GEDIIS::update_diis_space(const Matrix& x, const Matrix& g, double E) {
    //TODO: check

    xHist.push_back(x);
    gHist.push_back(g);
    EHist.push_back(E);

    assert(xHist.size() == gHist.size());
    assert(xHist.size() == EHist.size());

    if (xHist.size() >= 2) {
        Matrix B = make_B_matrix();
        if (std::abs(B.det()) < 1e-5) {
            xHist.pop_back();
            gHist.pop_back();
            EHist.pop_back();

            return false;
        }
    }


    if (xHist.size() > maxSpace) {
        xHist.erase(xHist.begin());
        gHist.erase(gHist.begin());
        EHist.erase(EHist.begin());
    }

    return true;
}


bool GEDIIS::diis_space_already() const {
    assert(xHist.size() == gHist.size());
    assert(xHist.size() == EHist.size());
    return xHist.size() >= minSpace;
}


Matrix GEDIIS::extrapolate_new_vector() const {
    assert(xHist.size() == gHist.size());
    assert(xHist.size() == EHist.size());
    assert(xHist.size() >= minSpace);

    std::size_t N = xHist.size();

    Matrix B = make_B_matrix();
    Matrix L = make_L_matrix();

    double absDet = std::abs(B.det());
    if (absDet < 1.0e-6) {
        std::cout << "det(B) is too small." << std::endl;
        std::cout << "det(B) = " << absDet <<  std::endl;
        std::abort();
    }

    Matrix c = B.solve(L);

    assert(c.rows() == N+1);
    assert(c.cols() == 1);

    double cSum = 0.0;
    for (std::size_t i = 0; i < N; ++i) {
        cSum += c(i);
    }
    std::cout << "c sum = " << cSum << std::endl;

    Matrix xnew(dim, 1, 0.0);
    for(std::size_t i = 0; i < N; ++i) {
        xnew += xHist[i] * c(i) / cSum;
    }
    
    return xnew;
}


Matrix GEDIIS::make_B_matrix() const {
    assert(xHist.size() == gHist.size());
    assert(xHist.size() == EHist.size());
    assert(xHist.size() >= 2);

    std::size_t N = xHist.size();
    Matrix B(N+1, N+1);

    for(std::size_t i = 0; i < N; ++i) {
    for(std::size_t j = 0; j < N; ++j) {
        B(i,j) = ((gHist[i] - gHist[j]).trans() 
                        % (xHist[i] - xHist[j]))(0);
    }}

    for (std::size_t i = 0; i < N; ++i) {
        B(i, N) = B(N, i) = 1.0;
    }

    B(N, N) = 0.0;

    return B;
}

Matrix GEDIIS::make_L_matrix() const {
    assert(xHist.size() == gHist.size());
    assert(xHist.size() == EHist.size());
    assert(xHist.size() >= minSpace);

    std::size_t N = xHist.size();

    Matrix L(N+1, 1);
    for (std::size_t i = 0; i < N; ++i) {
        L(i) = EHist[i];
    }

    L(N) = 1.0;

    return L;
}
