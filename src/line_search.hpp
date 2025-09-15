#pragma once

#include "matrix.hpp"
#include <vector>

double  fit_quartic_function(const Matrix& x0, const Matrix& g0, double E0, 
                             const Matrix& x1, const Matrix& g1, double E1, std::vector<double>& coeff);