#include "line_search.hpp"
#include <cassert>
#include <iostream>
#include <cstdlib>
#include <cmath>

double fit_quartic_function(const Matrix& x0, const Matrix& g0, double E0, 
                            const Matrix& x1, const Matrix& g1, double E1, std::vector<double>& coeff) {
    std::size_t N = x0.size();

    // 用更安全的检查替代assert
    if (x0.rows() != N || x0.cols() != 1 || x1.rows() != N || x1.cols() != 1 ||
        g0.rows() != N || g0.cols() != 1 || g1.rows() != N || g1.cols() != 1) {
        std::cerr << "Error: Matrix dimension check failed in quartic fit" << std::endl;
        return -100.0;
    }

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
        std::cerr << "quartic fit failed: discriminant too negative." << std::endl;
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
        std::cerr << "quartic fit failed: a4 = " << a4 << " < 0 (no minimum)" << std::endl;
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

// 三次拟合函数（高斯风格）
double fit_cubic_function(const Matrix& x0, const Matrix& g0, double E0, 
                          const Matrix& x1, const Matrix& g1, double E1, std::vector<double>& coeff) {
    std::size_t N = x0.size();
    
    // 用更安全的检查替代assert
    if (x0.rows() != N || x0.cols() != 1 || x1.rows() != N || x1.cols() != 1 ||
        g0.rows() != N || g0.cols() != 1 || g1.rows() != N || g1.cols() != 1) {
        std::cerr << "Error: Matrix dimension check failed in cubic fit" << std::endl;
        return -100.0;
    }
    
    const Matrix d = x1 - x0;
    double g0d = (g0.trans() % d)(0);
    double g1d = (g1.trans() % d)(0);
    
    // 三次多项式：f(α) = a0 + a1*α + a2*α² + a3*α³
    // 约束条件：f(0)=E0, f'(0)=g0d, f(1)=E1, f'(1)=g1d
    double a0 = E0;
    double a1 = g0d;
    
    // 解2x2线性方程组求a2, a3
    // a2 + a3 = E1 - E0 - g0d
    // 2*a2 + 3*a3 = g1d - g0d
    double b1 = E1 - E0 - g0d;
    double b2 = g1d - g0d;
    
    double det = 1.0;  // 1*3 - 1*2 = 1
    double a3 = b2 - 2.0*b1;
    double a2 = 3.0*b1 - b2;
    
    // 检查是否有最小值
    // f'(α) = a1 + 2*a2*α + 3*a3*α²
    if (std::abs(a3) < 1e-10) {
        // 退化为二次函数
        if (a2 <= 0) {
            std::cerr << "Cubic fit failed: no minimum" << std::endl;
            return -100.0;
        }
        double alpha = -a1 / (2.0*a2);
        coeff = {a0, a1, a2, 0.0};
        return alpha;
    }
    
    // 求导数为零的点
    double discriminant = 4.0*a2*a2 - 12.0*a3*a1;
    if (discriminant < 0) {
        std::cerr << "Cubic fit failed: no real critical points" << std::endl;
        return -100.0;
    }
    
    double sqrt_disc = std::sqrt(discriminant);
    double alpha1 = (-2.0*a2 + sqrt_disc) / (6.0*a3);
    double alpha2 = (-2.0*a2 - sqrt_disc) / (6.0*a3);
    
    // 选择二阶导数为正的根（最小值点）
    double alpha_min = -100.0;
    double f2_alpha1 = 2.0*a2 + 6.0*a3*alpha1;
    double f2_alpha2 = 2.0*a2 + 6.0*a3*alpha2;
    
    if (f2_alpha1 > 0 && alpha1 >= 0.0 && alpha1 <= 1.5) {
        alpha_min = alpha1;
    } else if (f2_alpha2 > 0 && alpha2 >= 0.0 && alpha2 <= 1.5) {
        alpha_min = alpha2;
    }
    
    if (alpha_min == -100.0) {
        std::cerr << "Cubic fit failed: no valid minimum in range" << std::endl;
        return -100.0;
    }
    
    coeff = {a0, a1, a2, a3};
    std::cout << "Cubic fit successful: alpha = " << alpha_min << std::endl;
    return alpha_min;
}

// 计算多项式在给定参数处的能量值
double evaluate_polynomial_energy(double alpha, const std::vector<double>& coeff) {
    if (coeff.size() == 4) {
        // 三次多项式
        double a0 = coeff[0], a1 = coeff[1], a2 = coeff[2], a3 = coeff[3];
        return a0 + a1*alpha + a2*alpha*alpha + a3*alpha*alpha*alpha;
    } else if (coeff.size() == 5) {
        // 四次多项式
        double a0 = coeff[0], a1 = coeff[1], a2 = coeff[2], a3 = coeff[3], a4 = coeff[4];
        return a0 + a1*alpha + a2*alpha*alpha + a3*alpha*alpha*alpha + a4*alpha*alpha*alpha*alpha;
    } else {
        std::cerr << "Error: coefficient vector size must be 4 (cubic) or 5 (quartic)" << std::endl;
        return 0.0;
    }
}

// 保持向后兼容的四次多项式能量计算
double evaluate_quartic_energy(double alpha, const std::vector<double>& coeff) {
    return evaluate_polynomial_energy(alpha, coeff);
}

// 计算多项式在给定参数处的梯度
Matrix evaluate_polynomial_gradient(double alpha, const Matrix& x0, const Matrix& x1, 
                                    const Matrix& g0, const Matrix& g1, const std::vector<double>& coeff) {
    // 使用线性插值计算梯度（这是一个近似）
    Matrix g_interp = (1.0 - alpha) * g0 + alpha * g1;
    return g_interp;
}

// 保持向后兼容的四次多项式梯度计算
Matrix evaluate_quartic_gradient(double alpha, const Matrix& x0, const Matrix& x1, 
                                 const Matrix& g0, const Matrix& g1, const std::vector<double>& coeff) {
    return evaluate_polynomial_gradient(alpha, x0, x1, g0, g1, coeff);
}

// 高斯风格的智能line search
double gaussian_style_line_search(const Matrix& x0, const Matrix& g0, double E0, 
                                 const Matrix& x1, const Matrix& g1, double E1, 
                                 std::vector<double>& coeff, double& best_energy) {
    
    std::cout << "=== Gaussian-style line search ===" << std::endl;
    
    // 安全检查：确保输入矩阵维度一致
    if (x0.size() != x1.size() || g0.size() != g1.size() || 
        x0.size() != g0.size() || x1.size() != g1.size()) {
        std::cerr << "Error: Matrix dimension mismatch in line search" << std::endl;
        best_energy = E0;
        return -100.0;
    }
    
    // 检查能量值是否合理
    if (std::isnan(E0) || std::isnan(E1) || std::isinf(E0) || std::isinf(E1)) {
        std::cerr << "Error: Invalid energy values in line search" << std::endl;
        best_energy = E0;
        return -100.0;
    }
    
    // 确定哪个点是最佳点（能量最低）
    bool current_is_best = (E0 < E1);
    std::cout << "Current point is " << (current_is_best ? "BEST" : "NOT best") 
              << " (E0=" << E0 << ", E1=" << E1 << ")" << std::endl;
    
    // 步骤1：尝试约束四次拟合（高斯的主要策略）
    std::cout << "Step 1: Trying constrained quartic fit..." << std::endl;
    double alpha_quartic = fit_quartic_function(x0, g0, E0, x1, g1, E1, coeff);
    
    if (alpha_quartic != -100.0) {
        double E_quartic = evaluate_polynomial_energy(alpha_quartic, coeff);
        
        // 高斯的接受性检查
        bool acceptable = false;
        if (current_is_best) {
            // 如果当前点是最佳点，任何四次步长都可接受
            acceptable = true;
            std::cout << "Quartic step acceptable: current is best point" << std::endl;
        } else {
            // 如果当前点不是最佳点，必须在两点之间
            if (alpha_quartic >= 0.0 && alpha_quartic <= 1.0) {
                acceptable = true;
                std::cout << "Quartic step acceptable: between current and best" << std::endl;
            } else {
                std::cout << "Quartic step rejected: not between points (alpha=" << alpha_quartic << ")" << std::endl;
            }
        }
        
        if (acceptable && E_quartic < E0 - 1e-7) {
            best_energy = E_quartic;
            std::cout << "✅ Quartic fit successful: alpha = " << alpha_quartic 
                      << ", energy improvement = " << (E0 - E_quartic) << std::endl;
            return alpha_quartic;
        }
    }
    
    // 步骤2：尝试三次拟合
    std::cout << "Step 2: Quartic failed, trying cubic fit..." << std::endl;
    double alpha_cubic = fit_cubic_function(x0, g0, E0, x1, g1, E1, coeff);
    
    if (alpha_cubic != -100.0) {
        double E_cubic = evaluate_polynomial_energy(alpha_cubic, coeff);
        
        // 高斯的三次步长接受性检查
        bool acceptable = false;
        if (alpha_cubic >= 0.0 && alpha_cubic <= 1.0) {
            // 在两点之间
            acceptable = true;
            std::cout << "Cubic step acceptable: between points" << std::endl;
        } else if (current_is_best && alpha_cubic <= 1.0) {
            // 当前点是最佳点，且步长不大于前一步
            acceptable = true;
            std::cout << "Cubic step acceptable: current is best and step ≤ previous" << std::endl;
        } else {
            std::cout << "Cubic step rejected: alpha=" << alpha_cubic << std::endl;
        }
        
        if (acceptable && E_cubic < E0 - 1e-7) {
            best_energy = E_cubic;
            std::cout << "✅ Cubic fit successful: alpha = " << alpha_cubic 
                      << ", energy improvement = " << (E0 - E_cubic) << std::endl;
            return alpha_cubic;
        }
    }
    
    // 步骤3：所有拟合都失败的情况
    std::cout << "Step 3: All fits failed, applying fallback strategy..." << std::endl;
    
    if (current_is_best) {
        // 如果当前点是最佳点，不采取任何线搜索步长
        std::cout << "❌ No line search step: current point is best and all fits failed" << std::endl;
        best_energy = E0;
        return -100.0;
    } else {
        // 如果当前点不是最佳点，取中点
        double alpha_midpoint = 0.5;
        double E_midpoint = 0.5 * E0 + 0.5 * E1;  // 线性插值估计
        
        std::cout << "✅ Fallback to midpoint: alpha = " << alpha_midpoint 
                  << ", estimated energy = " << E_midpoint << std::endl;
        
        best_energy = E_midpoint;
        // 创建线性插值的系数向量
        coeff = {E0, E1 - E0, 0.0, 0.0};
        return alpha_midpoint;
    }
}

// 保持向后兼容的智能line search
double smart_line_search(const Matrix& x0, const Matrix& g0, double E0, 
                        const Matrix& x1, const Matrix& g1, double E1, 
                        std::vector<double>& coeff, double& best_energy) {
    return gaussian_style_line_search(x0, g0, E0, x1, g1, E1, coeff, best_energy);
}


hpp文件

#pragma once

#include "matrix.hpp"
#include <vector>

double  fit_quartic_function(const Matrix& x0, const Matrix& g0, double E0, 
                             const Matrix& x1, const Matrix& g1, double E1, std::vector<double>& coeff);

// 三次拟合函数（高斯风格）
double  fit_cubic_function(const Matrix& x0, const Matrix& g0, double E0, 
                           const Matrix& x1, const Matrix& g1, double E1, std::vector<double>& coeff);

// 计算多项式在给定参数处的能量值
double  evaluate_polynomial_energy(double alpha, const std::vector<double>& coeff);

// 计算多项式在给定参数处的梯度
Matrix  evaluate_polynomial_gradient(double alpha, const Matrix& x0, const Matrix& x1, 
                                     const Matrix& g0, const Matrix& g1, const std::vector<double>& coeff);

// 高斯风格的智能line search
double  gaussian_style_line_search(const Matrix& x0, const Matrix& g0, double E0, 
                                  const Matrix& x1, const Matrix& g1, double E1, 
                                  std::vector<double>& coeff, double& best_energy);

// 保持向后兼容
double  evaluate_quartic_energy(double alpha, const std::vector<double>& coeff);
Matrix  evaluate_quartic_gradient(double alpha, const Matrix& x0, const Matrix& x1, 
                                  const Matrix& g0, const Matrix& g1, const std::vector<double>& coeff);
double  smart_line_search(const Matrix& x0, const Matrix& g0, double E0, 
                         const Matrix& x1, const Matrix& g1, double E1, 
                         std::vector<double>& coeff, double& best_energy);
