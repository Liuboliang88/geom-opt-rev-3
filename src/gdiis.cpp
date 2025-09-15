#include "gdiis.hpp"
#include <cassert>

GDIIS::GDIIS(std::size_t dim, std::size_t minSpace, std::size_t maxSpace)
: dim(dim), minSpace(minSpace), maxSpace(maxSpace) {}

void GDIIS::update_diis_space(const Matrix& x, const Matrix& e) {
    //TODO: check

    xHist.push_back(x);
    eHist.push_back(e);

    assert(xHist.size() == eHist.size());
    if (xHist.size() > maxSpace) {
        xHist.erase(xHist.begin());
        eHist.erase(eHist.begin());
    }
}


bool GDIIS::diis_space_already() const {
    assert(xHist.size() == eHist.size());
    return xHist.size() >= minSpace;
}


Matrix GDIIS::extrapolate_new_vector() const {
    assert(xHist.size() == eHist.size());
    assert(xHist.size() >= minSpace);

    std::size_t N = xHist.size();
    Matrix B(N+1, N+1);

    for(std::size_t i = 0; i < N; ++i) {
    for(std::size_t j = 0; j < N; ++j) {
        B(i,j) = (eHist[i].trans() % eHist[j])(0);
    }}

    for (std::size_t i = 0; i < N; ++i) {
        B(i, N) = B(N, i) = 1.0;
    }

    B(N, N) = 0.0;

    Matrix L(N+1, 1);
    L(N) = 1.0;


    Matrix c = B.solve(L);

    assert(c.rows() == N+1);
    assert(c.cols() == 1);

    Matrix xnew = Matrix(dim, 1, 0.0);
    for(std::size_t i = 0; i < N; ++i) {
        xnew += xHist[i] * c(i);
    }
    
    return xnew;
}