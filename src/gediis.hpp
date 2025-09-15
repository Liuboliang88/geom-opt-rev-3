#pragma once

#include "matrix.hpp"
#include <cstddef>

class GEDIIS {
public:
    GEDIIS(std::size_t dim, std::size_t minSpace = 4, 
                            std::size_t maxSpace = 10);

    bool update_diis_space(const Matrix& x, const Matrix& g, double E);
    bool diis_space_already() const;
    Matrix extrapolate_new_vector() const;

private:
    std::size_t dim = 0;
    std::size_t minSpace = 0;
    std::size_t maxSpace = 0;
    std::vector<Matrix> xHist;
    std::vector<Matrix> gHist;
    std::vector<double> EHist;

    Matrix make_B_matrix() const;
    Matrix make_L_matrix() const;
};