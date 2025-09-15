#include "rfo.hpp"
#include <cassert>
#include <iostream>
#include <cmath>
#include <limits>
#include <vector>
#include "matrix.hpp"


Matrix rfo_step(const Matrix& x, const Matrix& g, const Matrix& H) {
    size_t N = x.size();
    assert(x.rows() == N && x.cols() == 1);
    assert(g.rows() == N && g.cols() == 1);
    assert(H.rows() == N && H.cols() == N);

    // 1) 构造 (N+1)x(N+1) 扩展矩阵 [ H  g ; g^T  0 ]  RF矩阵
    Matrix RFOmat(N+1, N+1, 0.0);
    for (std::size_t i = 0; i < N; ++i) {
    for (std::size_t j = 0; j < N; ++j) {
        RFOmat(i, j) = H(i, j);
    }}

    for (std::size_t i = 0; i < N; ++i) {
        RFOmat(i, N) = g(i, 0);
        RFOmat(N, i) = g(i, 0);
    }

    RFOmat(N, N) = 0.0;

    SelfAdjointEigenSolver es(RFOmat);
    Matrix ev = es.eigen_vec();
    Matrix bestVec(N+1, 1);

    for (std::size_t i = 0; i < N+1; ++i) {
        bestVec(i) = ev(i, 0);
    }

    if (std::abs(bestVec(N)) < 1e-8) {
        std::cerr << "❌ RFO failed: bestVec(N) is too small.\n";
        std::abort();
    }

    // 4) 构造步长方向 p = v_head / v_tail
    Matrix p(N, 1);
    double norm2 = 0.0;
    for (size_t i = 0; i < N; ++i) {
        p(i) = bestVec(i) / bestVec(N);
        // norm2 += p(i) * p(i);
    }
    // double step_norm = std::sqrt(norm2);

    return x + p;
}