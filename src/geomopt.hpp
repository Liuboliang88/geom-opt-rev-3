#pragma once

#include "hessian.hpp"
#include "cell.hpp"
#include "rfo.hpp"
#include <string>

class GeomOpt {
public:
    GeomOpt(Cell& cell);

    // Gaussian 4 条收敛判据
    bool run(double gMaxTol = 4.5e-4,  // Max grad
             double gRMSTol = 3.0e-4,  // RMS grad
             double dMaxTol = 1.8e-3,  // Max disp
             double dRMSTol = 1.2e-3,  // RMS disp
             std::size_t maxCycle = 200);

private:
    Cell& cell;
    Hessian hess;

    std::vector<Matrix> xHist;
    std::vector<Matrix> gHist;
    std::vector<double> EHist;

};