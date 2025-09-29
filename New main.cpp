#include <iostream>
#include <string>
#include <algorithm>
#include <fstream>
#include <cstdio>
#include "cell.hpp"
#include "geomopt.hpp"
#include "rfo.hpp"

int main(int argc, char* argv[]) {
    if (argc < 2 || argc > 3) {
        std::cerr << "Usage: opt POSCAR [optimization_method]\n";
        std::cerr << "optimization_method options:\n";
        std::cerr << "  rfo          - RFO optimization (default)\n";
        std::cerr << "  rfo-linesearch - RFO with line search\n";
        std::cerr << "  gaussian     - Pure Gaussian optimization\n";
        return 1;
    }

    Cell cell = Cell::from_poscar(argv[1]);
    
    // 保存初始构型用于调试
    std::cout << "=== 初始构型信息 ===" << std::endl;
    std::cout << "原子数: " << cell.num_atoms() << std::endl;
    std::cout << "原子符号: ";
    for (const auto& symbol : cell.symbols) {
        std::cout << symbol << " ";
    }
    std::cout << std::endl;
    
    // 保存初始构型为POSCAR格式
    std::ofstream initial_poscar("initial_structure.poscar");
    initial_poscar << cell.to_poscar();
    initial_poscar.close();
    std::cout << "初始构型已保存到: initial_structure.poscar" << std::endl;
    
    // 保存初始构型为Gaussian gjf格式
    std::ofstream initial_gjf("initial_structure.gjf");
    initial_gjf << "%chk=initial.chk\n"
               << "# hf/6-31G(d) opt\n\n"
               << "Initial structure from POSCAR\n\n"
               << "0 1\n";
    
    Matrix coords_unit_A = cell.coords / 0.529177210544;  // 转换为埃单位
    for (size_t i = 0; i < coords_unit_A.rows(); ++i) {
        char buff[128];
        std::sprintf(buff, "  %16.8lf  %16.8lf  %16.8lf", 
                     coords_unit_A(i, 0), coords_unit_A(i, 1), coords_unit_A(i, 2));
        initial_gjf << cell.symbols[i] << "  " << std::string(buff) << std::endl;
    }
    initial_gjf << std::endl;
    initial_gjf.close();
    std::cout << "初始构型已保存到: initial_structure.gjf" << std::endl;
    std::cout << "========================" << std::endl;

    // 确定优化方法
    std::string opt_method = "rfo";  // 默认使用RFO
    if (argc == 3) {
        opt_method = argv[2];
    }
    
    std::cout << "🔧 优化方法: " << opt_method << std::endl;
    std::cout << "========================" << std::endl;
    
    bool success = false;
    
    if (opt_method == "gaussian") {
        // 纯Gaussian优化
        std::cout << "🚀 使用纯Gaussian几何优化" << std::endl;
        success = cell.pure_gaussian_optimize("hf/6-31G(d)", 200);
        
        if (success) {
            std::cout << "✅ Gaussian优化成功完成" << std::endl;
            
            // 保存优化后的构型
            std::ofstream final_poscar("final_structure_gaussian.poscar");
            final_poscar << cell.to_poscar();
            final_poscar.close();
            std::cout << "优化后构型已保存到: final_structure_gaussian.poscar" << std::endl;
        } else {
            std::cout << "❌ Gaussian优化失败" << std::endl;
        }
        
    } else if (opt_method == "rfo" || opt_method == "rfo-linesearch") {
        // RFO优化（带或不带线搜索）
        if (opt_method == "rfo-linesearch") {
            std::cout << "🔄 使用RFO + 线搜索优化" << std::endl;
        } else {
            std::cout << "🔄 使用RFO优化" << std::endl;
            // 注意：当前代码中线搜索是默认启用的，如果要纯RFO需要修改geomopt.cpp
        }
        
        GeomOpt opt(cell);
        success = opt.run();

        if (cell.is_converged_from_gauss("gauss.log")) {
            std::cout << "✅ Gaussian convergence table: ALL 4 criteria = YES.\n";
            
            // 保存优化后的构型
            std::ofstream final_poscar("final_structure_rfo.poscar");
            final_poscar << cell.to_poscar();
            final_poscar.close();
            std::cout << "优化后构型已保存到: final_structure_rfo.poscar" << std::endl;
        } else {
            std::cout << "⚠️ Not converged per Gaussian criteria.\n";
        }
        
    } else {
        std::cerr << "❌ 未知的优化方法: " << opt_method << std::endl;
        std::cerr << "支持的方法: rfo, rfo-linesearch, gaussian" << std::endl;
        return 1;
    }
    
    return success ? 0 : 2;

}
