// File: src/geomopt.cpp
#include "geomopt.hpp"
#include "matrix.hpp"
#include "rfo.hpp"
#include "bfgs.hpp"
#include "gdiis.hpp"
#include "gediis.hpp"
#include "line_search.hpp"

#include <cstdio>
#include <iostream>
#include <cmath>

// --- Gaussian 收敛检查工具函数 ---
static bool check_convergence(double gMax, double gRMS,
                              double dMax, double dRMS,
                              double gMaxTol, double gRMSTol,
                              double dMaxTol, double dRMSTol) {
    return (gMax < gMaxTol) && (gRMS < gRMSTol) &&
           (dMax < dMaxTol) && (dRMS < dRMSTol);
}

// --- 计算梯度最大值和 RMS ---
static double max_grad(const Matrix& g) {
    double gMax = 0.0;
    for (std::size_t i = 0; i < g.rows(); ++i) {
        double g2 = g(i,0)*g(i,0) + g(i,1)*g(i,1) + g(i,2)*g(i,2);
        gMax = std::max(gMax, std::sqrt(g2));
    }
    return gMax;
}
static double rms_grad(const Matrix& g) {
    double gRMS = 0.0;
    for (std::size_t i = 0; i < g.rows(); ++i) {
        gRMS += g(i,0)*g(i,0) + g(i,1)*g(i,1) + g(i,2)*g(i,2);
    }
    return std::sqrt(gRMS / g.rows());
}

// --- 计算位移最大值和 RMS ---
static double max_disp(const Matrix& dx) {
    double dMax = 0.0;
    for (std::size_t i = 0; i < dx.rows(); ++i) {
        double d2 = dx(i,0)*dx(i,0) + dx(i,1)*dx(i,1) + dx(i,2)*dx(i,2);
        dMax = std::max(dMax, std::sqrt(d2));
    }
    return dMax;
}
static double rms_disp(const Matrix& dx) {
    double dRMS = 0.0;
    for (std::size_t i = 0; i < dx.rows(); ++i) {
        dRMS += dx(i,0)*dx(i,0) + dx(i,1)*dx(i,1) + dx(i,2)*dx(i,2);
    }
    return std::sqrt(dRMS / dx.rows());
}


// ---------- GeomOpt ----------
GeomOpt::GeomOpt(Cell& cell) : cell(cell), hess(cell.coords.size()) {}

bool GeomOpt::run(double maxGrad, 
                  double rmsGrad, 
                  double maxDisp,
                  double rmsDisp,
                  std::size_t maxCycle) {
    const std::size_t nAtom = cell.num_atoms();

    GDIIS  gdiis (3 * nAtom);
    GEDIIS gediis(3 * nAtom);

    for (std::size_t cycle = 0; cycle < maxCycle; ++cycle) {
        std::printf("cycle = %zu    ", cycle);

        // === 1) 从 Gaussian 得到本轮梯度与能量 ===
        // 说明：calculate_gradient_from_gauss() 会写 gauss.gjf 并调用 g16 生成 gauss.log
        Matrix x = cell.coords;                          // 当前坐标（nAtom x 3）
        Matrix g = cell.calculate_gradient_from_gauss(); // 本轮梯度（会触发 g16）
        double E = cell.get_energy_from_gauss();        // 从同一个 gauss.log 里取能量


        Matrix xv = x.reshape(3*nAtom, 1);
        Matrix gv = g.reshape(3*nAtom, 1);

        xHist.push_back(xv);
        gHist.push_back(gv);
        EHist.push_back(E);

        // line search
        if (xHist.size() > 2) {
            std::size_t N = xHist.size();
            Matrix x0 = xHist[N-1];
            Matrix x1 = xHist[N-2];
            Matrix g0 = gHist[N-1];
            Matrix g1 = gHist[N-2];

            double E0 = EHist[N-1];
            double E1 = EHist[N-2];

            std::vector<double> coeff;
            double alpha = fit_quartic_function(x0, g0, E0, x1, g1, E1, coeff);
            std::cout << "alpha = " << alpha << std::endl;

            if (alpha > -1.0 && alpha < 2) {
                xv = (1.0-alpha) * x0 + alpha * x1;
                gv = (1.0-alpha) * g0 + alpha * g1;
                std::cout << "line search is used." << std::endl;
            }else  {
                 std::cout << "line search is not used" << std::endl;    
            }

           // TODO: change last history by fitted coord, grad and energy?
        }


        // === 2) 立即做“高斯 4 条收敛判据”检查（每一轮都检查）===
        // 只要 gauss.log 上显示四个 YES，就按 Gaussian 的判据提前收敛退出
        if (cell.is_converged_from_gauss("gauss.log")) {
            std::cout << "✅ Gaussian convergence table: ALL 4 criteria = YES (cycle "
                      << cycle << ").\n";
            // 同时也打印一下本轮梯度概览，方便定位
            const double gMax_now = max_grad(g);
            const double gRMS_now = rms_grad(g);
            std::printf("max grad = %.6lf, rms grad = %.6lf\n", gMax_now, gRMS_now);
            std::cout << "opt converged.\n";
            return true;
        }

        // === 3) 传统 max/rms 阈值（你原来的自定义阈值逻辑，保留）===
        const double gMax = max_grad(g);
        const double gRMS = rms_grad(g);
        std::printf("max grad = %.6lf, rms grad = %.6lf\n", gMax, gRMS);

        // if (gMax < maxGrad && gRMS < rmsGrad) {
        //     std::cout << "opt converged (custom threshold).\n";
        //     return true;
        // }

        hess.update_hessian_gdiis(xv, gv);
        Matrix H = hess.hessian();

        // // === 5) 选择步子（保留你原有的 RFO / GEDIIS 混合逻辑）===
        // // GDIIS：用 H^{-1} g 作为误差向量进入空间
        // Matrix e = H.inver() % gv;
        // gdiis.update_diis_space(xv, e);

        // // GEDIIS：能量 + 梯度进入 DIIS 空间
        // bool addSucc = gediis.update_diis_space(xv, gv * 100.0, E);
        // std::cout << "gediis.update_diis_space = " << (addSucc ? 1 : 0) << std::endl;

        Matrix xnew;
        // if (gediis.diis_space_already() && addSucc && gRMS < 0.01) {
        //     std::cout << "gediis step" << std::endl;
        //     xnew = gediis.extrapolate_new_vector();
        // } else {
        //     std::cout << "rfo step" << std::endl;
        //     xnew = rfo_step(xv, gv, H);
        // }

        std::cout << "rfo step" << std::endl;
        xnew = rfo_step(xv, gv, H);

        // === 6) trust 半径控制（与你原本一致）===
        Matrix coordNew = xnew.reshape(nAtom, 3);
        Matrix dcoord   = coordNew - cell.coords;

        const double BohrPerAngstrom = 0.529177210544;
        const double trust = 0.2 * BohrPerAngstrom; // Bohr
        double dmax = 0.0;
        for (std::size_t i = 0; i < nAtom; ++i) {
            const double dr = std::sqrt(  dcoord(i,0)*dcoord(i,0)
                                        + dcoord(i,1)*dcoord(i,1)
                                        + dcoord(i,2)*dcoord(i,2)  );
            dmax = std::max(dmax, dr);
        }

        if (dmax > trust) {
            const double factor = trust / dmax;
            std::printf("trust factor = %.3lf\n", factor);
            dcoord *= factor;
        }

        // === 7) 更新坐标，进入下一轮 ===
        cell.coords += dcoord;

        std::cout << std::endl;
    }

    return false; // 未达成收敛
}
