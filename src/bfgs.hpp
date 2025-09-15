#pragma once

#include "matrix.hpp"

Matrix bfgs_step(const Matrix& x, const Matrix& g, const Matrix& H);
Matrix bfgs_step_modified(const Matrix& x, const Matrix& g, const Matrix& H);