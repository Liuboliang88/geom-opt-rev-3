// File: cell.hpp
#pragma once

#include <string>
#include <vector>
#include "vec3f64.hpp"
#include "matrix.hpp"

class Cell {
public:
    std::vector<std::string> species;   // 元素符号顺序
    std::vector<std::size_t> numAtom;   // 每个元素的原子数
    
    std::vector<std::string> symbols;   // 每个原子对应的元素符号
    
    Matrix lattice;                     // 3x3 矩阵
    Matrix coords;                      // 原子笛卡尔坐标

    Cell();  // 默认构造函数
    static Cell from_poscar(const std::string& filename);

    std::size_t num_atoms() const;

    std::string to_geometry() const;
    std::string to_poscar() const;
    Matrix calculate_gradient_from_vasp() const;
    Matrix calculate_gradient_from_gauss() const;
    
    // 新增：从Gaussian输出中解析能量（Hartree单位）
    double get_energy_from_gauss() const;

    bool is_converged_from_gauss(const std::string& log_file = "gauss.log") const;
};
