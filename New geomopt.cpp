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

        // 将真实的Gaussian计算结果添加到历史记录
        // 重要：历史记录应该只包含真实的量化计算结果，不包含插值/拟合数据
        xHist.push_back(xv);
        gHist.push_back(gv);
        EHist.push_back(E);
        // === 改进的自适应线搜索策略 ===
        const double gMax_current = max_grad(g);
        bool should_use_line_search = false;
        std::string ls_reason = "";
        
        if (xHist.size() > 1) {
            double E_curr = EHist[EHist.size()-1];
            double E_prev = EHist[EHist.size()-2];
            
            // 策略1：早期积极干预 (前25步)
            if (cycle < 25) {
                if (gMax_current > 0.012) {  // 重新调高阈值，减少无效触发
                    should_use_line_search = true;
                    ls_reason = "Early high gradient guidance";
                }
                // 早期能量停滞检测 - 更严格
                else if (xHist.size() >= 3 && gMax_current > 0.008) {
                    double recent_improve = std::abs(E_curr - E_prev);
                    double prev_improve = std::abs(EHist[EHist.size()-2] - EHist[EHist.size()-3]);
                    if (recent_improve < 0.3 * prev_improve && prev_improve > 5e-6) {
                        should_use_line_search = true;
                        ls_reason = "Early convergence stagnation";
                    }
                }
            }
            
            // 策略2：中期救援和振荡处理 (30-100步)
            else if (cycle < 100) {
                if (E_curr > E_prev + 3e-6 && gMax_current > 0.003) {  // 更敏感的救援
                    should_use_line_search = true;
                    ls_reason = "Mid-stage energy increase rescue";
                }
                // 梯度振荡检测
                else if (xHist.size() >= 4 && gMax_current > 0.004) {
                    // 简化的振荡检测：检查最近4步的梯度变化
                    std::vector<double> recent_grads;
                    for (int i = std::max(0, (int)gHist.size()-4); i < gHist.size(); ++i) {
                        recent_grads.push_back(max_grad(gHist[i].reshape(cell.num_atoms(), 3)));
                    }
                    if (recent_grads.size() >= 3) {
                        bool has_oscillation = false;
                        for (int i = 1; i < recent_grads.size()-1; ++i) {
                            if ((recent_grads[i] > recent_grads[i-1] && recent_grads[i] > recent_grads[i+1]) ||
                                (recent_grads[i] < recent_grads[i-1] && recent_grads[i] < recent_grads[i+1])) {
                                has_oscillation = true;
                                break;
                            }
                        }
                        if (has_oscillation) {
                            should_use_line_search = true;
                            ls_reason = "Mid-stage gradient oscillation";
                        }
                    }
                }
            }
            
            // 策略3：后期精细救援 (100+步)
            else {
                if (E_curr > E_prev + 8e-7 && gMax_current > 0.001) {
                    should_use_line_search = true;
                    ls_reason = "Late-stage precision rescue";
                }
            }
        }
        
        if (should_use_line_search) {
            std::cout << "🔍 Line search triggered: " << ls_reason << std::endl;
            
            std::size_t N = xHist.size();
            double E_current = EHist[N-1];
            
            // 智能参考点选择：扩大搜索范围
            std::size_t search_range = std::min(N-1, std::size_t(8));  // 增加到8步
            if (cycle < 20) search_range = std::min(N-1, std::size_t(3));  // 早期保守
            
            std::size_t best_ref_idx = N-2;
            double E_best_ref = EHist[N-2];
            
            for (std::size_t i = N-1-search_range; i < N-1; ++i) {
                if (EHist[i] < E_best_ref) {
                    E_best_ref = EHist[i];
                    best_ref_idx = i;
                }
            }
            
            Matrix x_ref = xHist[best_ref_idx];
            Matrix g_ref = gHist[best_ref_idx];
            double E_ref = E_best_ref;
            
            if (best_ref_idx != N-2) {
                std::cout << "Using optimal reference from step " << best_ref_idx 
                          << " (E=" << E_ref << ")" << std::endl;
            }
            
            // 多候选点策略
            std::vector<double> candidate_alphas;
            if (E_current > E_ref) {
                // 能量上升：倾向参考点
                candidate_alphas = {0.3, 0.5, 0.7, 0.85};
            } else {
                // 能量下降：保守插值
                candidate_alphas = {0.2, 0.4, 0.6};
            }
            
            // 根据梯度和阶段调整候选点
            if (gMax_current > 0.02) {
                candidate_alphas.push_back(0.9);  // 高梯度时更激进
            }
            if (cycle > 50) {
                candidate_alphas.insert(candidate_alphas.begin(), 0.1);  // 后期精细调整
            }
            
            // 评估候选点并选择最佳
            double best_alpha = -1.0;
            double best_improvement = -1e10;
            
            // 自适应改善阈值 - 更严格，减少无效线搜索
            double min_improvement;
            if (cycle < 15) {
                min_improvement = 2e-6;      // 早期更严格，确保有意义的改善
            } else if (cycle < 60) {
                min_improvement = 8e-7;      // 中期中等
            } else {
                min_improvement = 3e-7;      // 后期适度严格
            }
            
            for (double alpha : candidate_alphas) {
                // 简化能量估计：线性插值 + 二次修正
                double E_est = (1.0 - alpha) * E_current + alpha * E_ref;
                
                // 添加基于梯度的二次修正
                Matrix dx = x_ref - xHist[N-1];
                double grad_correction = 0.0;
                for (std::size_t i = 0; i < dx.size(); ++i) {
                    grad_correction += gHist[N-1](i) * dx(i);
                }
                E_est += alpha * (1.0 - alpha) * grad_correction * 0.1;  // 小的二次修正
                
                double improvement = E_current - E_est;
                
                if (improvement > min_improvement && alpha >= 0.0 && alpha <= 1.0) {
                    if (improvement > best_improvement) {
                        best_improvement = improvement;
                        best_alpha = alpha;
                    }
                }
            }
            
            // 应用最佳线搜索步长
            if (best_alpha > 0.0) {
                Matrix xv_fitted = (1.0 - best_alpha) * xHist[N-1] + best_alpha * x_ref;
                Matrix gv_fitted = (1.0 - best_alpha) * gHist[N-1] + best_alpha * g_ref;
                
                // 梯度合理性检查
                double grad_norm_fitted = 0.0, grad_norm_current = 0.0;
                for (std::size_t i = 0; i < gv_fitted.size(); ++i) {
                    grad_norm_fitted += gv_fitted(i) * gv_fitted(i);
                    grad_norm_current += gHist[N-1](i) * gHist[N-1](i);
                }
                grad_norm_fitted = std::sqrt(grad_norm_fitted);
                grad_norm_current = std::sqrt(grad_norm_current);
                
                if (grad_norm_fitted < 1.5 * grad_norm_current) {  // 放宽梯度检查
                    xv = xv_fitted;
                    gv = gv_fitted;
                    
                    std::cout << "✅ Enhanced line search applied: α=" << best_alpha 
                              << ", ΔE=" << best_improvement 
                              << ", ref_step=" << best_ref_idx << std::endl;
                } else {
                    std::cout << "❌ Line search rejected: gradient norm constraint" << std::endl;
                }
            } else {
                std::cout << "❌ Line search failed: no acceptable candidate found" << std::endl;
            }
            
        } else if (xHist.size() > 1) {
            std::cout << "Line search skipped: gMax=" << gMax_current 
                      << " (adaptive thresholds by stage)" << std::endl;
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

HPP文件

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
             std::size_t maxCycle = 400);

private:
    Cell& cell;
    Hessian hess;

    std::vector<Matrix> xHist;
    std::vector<Matrix> gHist;
    std::vector<double> EHist;

};

