#pragma once

#include "matrix.hpp"
#include <cstddef>

class GDIIS {
public:
    GDIIS(std::size_t dim, std::size_t minSpace = 4, 
                           std::size_t maxSpace = 10);

    void update_diis_space(const Matrix& x, const Matrix& e);
    bool diis_space_already() const;
    Matrix extrapolate_new_vector() const;

private:
    std::size_t dim = 0;
    std::size_t minSpace = 0;
    std::size_t maxSpace = 0;
    std::vector<Matrix> xHist;
    std::vector<Matrix> eHist;
};