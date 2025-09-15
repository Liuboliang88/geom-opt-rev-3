#include "line_search.hpp"
#include <cassert>
#include <iostream>
#include <cstdlib>
#include <cmath>

double fit_quartic_function(const Matrix& x0, const Matrix& g0, double E0, 
                            const Matrix& x1, const Matrix& g1, double E1, std::vector<double>& coeff) {
    std::size_t N = x0.size();

    assert(x0.rows() == N && x0.cols() == 1);
    assert(x1.rows() == N && x1.cols() == 1);
    assert(g0.rows() == N && g0.cols() == 1);
    assert(g1.rows() == N && g1.cols() == 1);

    const Matrix d = x1 - x0;

    double g0d = (g0.trans() % d)(0);
    double g1d = (g1.trans() % d)(0);

    double a0 = E0;
    double a1 = g0d;

    double c1 = E1 - E0 - g0d;
    double c2 = g1d - g0d;

    // -60. c1^2 + 27. c1 c2 - 1. c2^2
    double inSqrt = -60.0*c1*c1 + 27.0*c1*c2 - c2*c2;
    std::cerr << "in sqrt = " << inSqrt << std::endl;

    if (inSqrt < -1e-2) {
        std::cerr << "quartic fit failed." << std::endl;
        return -100.0;
    } 
    
    if (inSqrt < 0.0) {
        std::cerr << "in sqrt = " << inSqrt << std::endl;
        inSqrt = 0.0;
    }

    // a2 -> 0.00613497 (-12. c1 + 19. c2 - 11.3137 Sqrt[-60. c1^2 + 27. c1 c2 - 1. c2^2])
    // a3 -> 0.0122699 (80. c1 - 18. c2 -  1.41421 Sqrt[-60. c1^2 + 27. c1 c2 - 1. c2^2])
    // a4 -> 0.00613497 (15. c1 + 17. c2 + 14.1421 Sqrt[-60. c1^2 + 27. c1 c2 - 1. c2^2])

    double sqt = std::sqrt(inSqrt);
    double a2 = 0.00613497 * (-12.0 * c1 + 19.0 * c2 - 11.3137 * sqt);
    double a3 = 0.01226990 * ( 80.0 * c1 - 18.0 * c2 - 1.41421 * sqt);
    double a4 = 0.00613497 * ( 15.0 * c1 + 17.0 * c2 + 14.1421 * sqt);

    if (a4 < 0.0) {
        std::cerr << "a4 = " << a4 << std::endl;
        return -100.0;
    }

    double alpha = -0.25 * a3 / a4;
    coeff = {a0, a1, a2, a3, a4};

    auto first_deriv = [&](double x) {
        return a1 + 2.0*a2 * x + 3.0*a3 * x*x + 4.0*a4 * x*x*x;
    };

    std::cout << "a4 = " << a4 << std::endl;
    std::cout << "frist derviative = " << first_deriv(alpha) << std::endl;

    return alpha;
}