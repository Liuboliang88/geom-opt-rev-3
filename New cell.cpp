// File: cell.cpp
#include "cell.hpp"
#include <fstream>
#include <sstream>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <cassert>

static const double BohrPerAngstrom = 0.529177210544;

static std::string to_string(const Matrix& m) {
    char buff[128];
    std::ostringstream oss;

    for (std::size_t i = 0; i < m.rows(); ++i) {
        for (std::size_t j = 0; j < m.cols(); ++j) {
            std::sprintf(buff, "  %16.8lf", m(i, j));
            oss << std::string(buff);
        }
        oss << std::endl;
    }
    
    return oss.str();
}


// ===== 构造函数 =====
Cell::Cell() {
    std::system(
      "rm -rf "
      "CHG CHGCAR WAVECAR OUTCAR "
      "OSZICAR vasprun.xml vasp.log "
      "DOSCAR EIGENVAL IBZKPT XDATCAR "
      "PCDAT CONTCAR REPORT force.* "
    );
}
// ===== 从POSCAR读取晶胞 =====
Cell Cell::from_poscar(const std::string& filename) {
    Cell cell;
    std::ifstream fin(filename);
    if (!fin) {
        std::cerr << "❌ Failed to open POSCAR file: " << filename << "\n";
        std::abort();
    }
    // 1) header + scale
    std::string line;
    std::getline(fin, line);  // comment

    std::getline(fin, line);  // scale factor
    double scale = std::stod(line);

    //  2) lattice vectors
    cell.lattice = Matrix(3, 3);
    for (std::size_t i = 0; i < 3; ++i) {
        std::getline(fin, line);
        std::istringstream iss(line);
        iss >> cell.lattice(i, 0) >> cell.lattice(i, 1) >> cell.lattice(i, 2);
    }

    cell.lattice *= scale * BohrPerAngstrom;


    //  3) species & counts
    std::getline(fin, line);
    std::istringstream iss_species(line);
    std::string sp;
    while (iss_species >> sp) cell.species.push_back(sp);

    
    std::getline(fin, line);
    std::istringstream iss_num(line);
    int count;
    while (iss_num >> count) cell.numAtom.push_back(count);

    // 4) build symbol list
    for (std::size_t i = 0; i < cell.species.size(); ++i) {
        for (std::size_t j = 0; j < cell.numAtom[i]; ++j) {
            cell.symbols.push_back(cell.species[i]);
        }
    }


    // 5) detect "Direct" or "Cartesian" line
    bool is_direct = false;
    while (std::getline(fin, line)) {
        if (line.find("Direct") != std::string::npos) {
            is_direct = true;
            break;
        }
        if (line.find("Cartesian") != std::string::npos) {
            is_direct = false;
            break;
        }
    }

    std::size_t nAtom = cell.symbols.size();

    // 6) read coordinates, 并在 fraction 模式下转换到 Cartesian
    cell.coords = Matrix(nAtom, 3);
    for (std::size_t i = 0; i < cell.coords.rows(); ++i) {
        fin >> cell.coords(i, 0) >> cell.coords(i, 1) >> cell.coords(i, 2);

        if (is_direct) {
            Matrix frac(1, 3);
            frac(0, 0) = cell.coords(i, 0);
            frac(0, 1) = cell.coords(i, 1);
            frac(0, 2) = cell.coords(i, 2);

            Matrix cart = frac % cell.lattice;

            cell.coords(i, 0) = cart(0, 0);
            cell.coords(i, 1) = cart(0, 1);
            cell.coords(i, 2) = cart(0, 2);
        }
    }

    for (std::size_t i = 0; i < cell.coords.size(); ++i) {
        if (std::isnan(cell.coords(i))) {
            std::cerr << "❌ 检测到NaN坐标,终止优化！" << std::endl;
            std::abort();
        }
    }

    cell.coords *= BohrPerAngstrom;

    return cell;
}

std::size_t Cell::num_atoms() const {
    assert(symbols.size() == coords.rows());
    return symbols.size();
}


std::string Cell::to_geometry() const {
    char buff[128];
    std::ostringstream oss;
    for (size_t i = 0; i < coords.size(); ++i) {
        std::sprintf(buff, "  %16.8lf  %16.8lf  %16.8lf", 
                     coords(i, 0), coords(i, 1), coords(i, 2));
        oss << symbols[i] << std::string(buff) << std::endl;
    }
    return oss.str();
}

std::string Cell::to_poscar() const {
    std::ostringstream oss;
    oss << "VASP cell\n1.0\n";
    oss << to_string(lattice / BohrPerAngstrom);

    for (const auto& sp : species) {
        oss << sp << " ";
    }
    oss << std::endl;

    for (const auto& num : numAtom) {
        oss << num << " ";
    }

    oss << std::endl;
    oss << "Cartesian" << std::endl;
    oss << to_string(coords / BohrPerAngstrom);
    return oss.str();
}

// ===== VASP梯度计算 =====
Matrix Cell::calculate_gradient_from_vasp() const {
    const std::string job = "force";
    // 写 POSCAR
    std::ofstream fout("POSCAR");
    if (!fout) {
        std::cerr << "❌ Cannot write POSCAR file.\n";
        std::abort();
    }

    fout << this -> to_poscar();
    fout.close();

    // 写 POTCAR
    int code0 = std::system("echo 103 | vaspkit > /dev/null");
    if (code0 != 0) {
        std::cerr << "❌ Failed to generate POTCAR.\n";
        std::abort();
    }

    // 写 INCAR
    std::ofstream incar("INCAR");
    incar << 
R"(PREC = ACCURATE
IBRION  = -1
NSW     = 0
ENCUT   = 400
EDIFF   = 1e-06
EDIFFG  = -0.02
IVDW    = 12
ISMEAR  = 0
SIGMA   = 0.05
ISIF    = 2
ISPIN   = 1
# LDIPOL  = .TRUE.
# IDIPOL  = 3
LERAL   = Auto
ISTART  = 0
ICHARG  = 2)";

    incar.close();

    // 写 KPOINTS
    std::ofstream kpoints("KPOINTS");
    kpoints << "Gamma-point only\n0\nGamma\n1 1 1\n";
    kpoints.close();

    std::cout << "📎 Expecting POTCAR in working directory\n";
    std::cout << "⚙️  Running VASP (vasp_std > out.log)...\n";

    int code = std::system("mpirun -np 16 vasp_gam > vasp.log");
    if (code != 0) {
        std::cerr << "❌ VASP execution failed.\n";
        std::abort();
    }

    // 解析 OUTCAR 中的 TOTAL-FORCE
    std::ifstream log("OUTCAR");
    if (!log) {
        std::cerr << "❌ Cannot open OUTCAR.\n";
        std::abort();
    }

    const double sc = 0.0194469;

    std::size_t nAtom = num_atoms();
    Matrix forces(nAtom, 3);

    std::string line;
    std::size_t idx = 0;
    while (std::getline(log, line)) {
        if (line.find("TOTAL-FORCE (eV/Angst)") != std::string::npos) {
            std::getline(log, line); // --------

            while (std::getline(log, line)) {
                if (line.find("--------") != std::string::npos) {
                    break;
                }

                std::istringstream iss(line);
                double x, y, z, fx, fy, fz;
                iss >> x >> y >> z >> fx >> fy >> fz;

                forces(idx, 0) = fx;
                forces(idx, 1) = fy;
                forces(idx, 2) = fz;

                ++idx;
            }

            break;
        }
    }

    if (idx != nAtom) {
        std::cerr << "❌ Parsed force number mismatch (" << idx
                  << " vs " << nAtom << ").\n";
        std::abort();
    }

    return forces;
}


Matrix Cell::calculate_gradient_from_gauss() const {
    const std::string job = "gauss";
    const std::string gjf = job + ".gjf";
    const std::string log = job + ".log";

    int charge = 0;
    int multip = 1;

    // 1) write .gjf input
    std::ofstream f(gjf);
    f << "%chk=" << job << ".chk\n"
        << "# hf/6-31G(d) force\n\n"
        << "Force calculation\n\n"
        << charge << " " << multip << "\n";

    char buff[128];

    Matrix coords_unit_A = coords / BohrPerAngstrom;

    for (size_t i = 0; i < coords_unit_A.rows(); ++i) {
        std::sprintf(buff, "  %16.8lf  %16.8lf  %16.8lf", 
                     coords_unit_A(i, 0),
                     coords_unit_A(i, 1), 
                     coords_unit_A(i, 2));

        f << symbols[i] << "  " << std::string(buff) << std::endl;
    }

    f << std::endl;

    // 2) use g16（确保设置 GAUSS_SCRDIR，某些环境未设置会导致 g16 段错误）
    std::cout << "⚙️ Running Gaussian...\n";
    const std::string scratch = "g16_scratch";
    std::system(("mkdir -p " + scratch + " > /dev/null 2>&1").c_str());
    const std::string run_cmd =
        "GAUSS_SCRDIR=./" + scratch +
        " OMP_NUM_THREADS=1 g16 " + gjf + " > " + log + " 2>&1";
    if (std::system(run_cmd.c_str()) != 0) {
        std::cerr << "❌ Gaussian failed.\n"; std::abort();
    }

    static std::size_t istep = 0;
    std::string cmdStr = "cp gauss.log   gauss" + std::to_string(istep) + ".log";
    std::system(cmdStr.data());

    istep++;


    // 3) analys .log 中的 "Forces (Hartrees/Bohr)"
    std::ifstream fin(log);
    if (!fin) { std::cerr << "❌ Cannot open " << log << "\n"; std::abort(); }


    std::size_t nAtom = num_atoms();
    Matrix grads(nAtom, 3);
    std::size_t idx = 0;

    std::string line;
    while (std::getline(fin, line)) {
        if (line.find("Forces (Hartrees/Bohr)") != std::string::npos) {
            // skip header lines
            std::getline(fin, line);
            std::getline(fin, line);

            while (std::getline(fin, line)) {
                if (line.find("----") != std::string::npos) break;

                std::istringstream iss(line);
                int atmIdx, Z;
                double fx, fy, fz;
                if (iss >> atmIdx >> Z >> fx >> fy >> fz) {
                    // gradient = –Force
                    // grads.emplace_back(-fx, -fy, -fz);
                    grads(idx, 0) = -fx;
                    grads(idx, 1) = -fy;
                    grads(idx, 2) = -fz;
                    ++idx;
                }
            }
            break;
        }
    }

    if (idx != nAtom) {
        std::cerr << "❌ Gaussian gradient count mismatch\n";
        std::abort();
    }
    return grads;
}
// ===== Gaussian能量提取 =====
double Cell::get_energy_from_gauss() const {
    const std::string log = "gauss.log";
    std::ifstream fin(log);
    if (!fin) { 
        std::cerr << "❌ Cannot open " << log << "\n"; 
        std::abort(); 
    }

    double energy = 0.0;
    std::string line;
    bool found_energy = false;
    
    while (std::getline(fin, line)) {
        // 查找SCF Done行，包含最终能量
        if (line.find("SCF Done:") != std::string::npos) {
            std::istringstream iss(line);
            std::string token;
            // 跳过 "SCF Done: E(RHF) ="
            for (int i = 0; i < 4; ++i) {
                iss >> token;
            }
            // 读取能量值
            iss >> energy;
            found_energy = true;
            break;
        }
    }
    
    if (!found_energy) {
        std::cerr << "❌ 未能在Gaussian输出中找到能量信息\n";
        std::abort();
    }
    
    return energy;  // 返回Hartree单位的能量
}

// ===== Gaussian收敛性检查 =====
bool Cell::is_converged_from_gauss(const std::string& log_file) const {
    std::ifstream fin(log_file);
    if (!fin) {
        std::cerr << "❌ Cannot open " << log_file << " for convergence check.\n";
        return false;
    }

    std::string line;
    bool in_table = false;
    int yes_count = 0;

    while (std::getline(fin, line)) {
        // 保险：Gaussian 已给出全局收敛提示
        if (line.find("Optimization completed") != std::string::npos) {
            return true;
        }

        // 找到收敛表的表头
        if (!in_table &&
            line.find("Item") != std::string::npos &&
            line.find("Converged?") != std::string::npos) {
            in_table = true;
            yes_count = 0;
            continue;
        }

        if (in_table) {
            // 表格一般在遇到空行或 Predicted change 后结束
            if (line.empty() ||
                line.find("Predicted change in Energy") != std::string::npos) {
                break;
            }

            // 统计 YES（Gaussian 会把每一行的最后一列标为 YES/NO）
            if (line.find("YES") != std::string::npos) {
                ++yes_count;
            }
        }
    }

    // 四个判据均为 YES，则视为收敛
    return (yes_count >= 4);
}

// ===== 纯Gaussian几何优化 =====
bool Cell::pure_gaussian_optimize(const std::string& method, int max_cycles) {
    const std::string gjf = "gauss_opt.gjf";
    const std::string log = "gauss_opt.log";
    
    int charge = 0;
    int multip = 1;
    
    // 1) 生成Gaussian优化输入文件
    std::ofstream f(gjf);
    if (!f) {
        std::cerr << "❌ Cannot create " << gjf << std::endl;
        return false;
    }
    
    f << "%chk=gauss_opt.chk" << std::endl;
    f << "%mem=8GB" << std::endl;
    f << "%nprocshared=8" << std::endl;
    f << "# " << method << " opt=(cartesian,loose,maxcycles=" << max_cycles << ") scf=xqc" << std::endl;
    f << std::endl;
    f << "Pure Gaussian geometry optimization" << std::endl;
    f << std::endl;
    f << charge << " " << multip << std::endl;
    
    // 输出笛卡尔坐标（转换为埃单位）
    Matrix coords_unit_A = coords / BohrPerAngstrom;
    char buff[128];
    
    for (size_t i = 0; i < coords_unit_A.rows(); ++i) {
        std::sprintf(buff, "  %16.8lf  %16.8lf  %16.8lf", 
                     coords_unit_A(i, 0),
                     coords_unit_A(i, 1), 
                     coords_unit_A(i, 2));
        f << symbols[i] << "  " << std::string(buff) << std::endl;
    }
    
    f << std::endl;
    f.close();
    
    // 2) 运行Gaussian优化
    std::cout << "🚀 启动纯Gaussian几何优化..." << std::endl;
    std::cout << "   方法: " << method << std::endl;
    std::cout << "   最大循环数: " << max_cycles << std::endl;
    std::cout << "   输入文件: " << gjf << std::endl;
    std::cout << "   输出文件: " << log << std::endl;
    
    // 设置scratch目录和环境变量以避免段错误
    const std::string scratch = "g16_scratch";
    std::system(("mkdir -p " + scratch + " > /dev/null 2>&1").c_str());
    const std::string cmd =
        "GAUSS_SCRDIR=./" + scratch +
        " OMP_NUM_THREADS=1 g16 " + gjf + " > " + log + " 2>&1";
    int result = std::system(cmd.c_str());
    
    if (result != 0) {
        std::cerr << "❌ Gaussian优化执行失败！" << std::endl;
        return false;
    }
    
    // 3) 检查优化是否收敛
    std::ifstream fin(log);
    if (!fin) {
        std::cerr << "❌ 无法打开Gaussian优化输出文件: " << log << std::endl;
        return false;
    }
    
    std::string line;
    bool optimization_completed = false;
    bool stationary_point_found = false;
    int cycles_used = 0;
    
    while (std::getline(fin, line)) {
        // 检查是否找到稳定点
        if (line.find("Stationary point found") != std::string::npos) {
            stationary_point_found = true;
            optimization_completed = true;
        }
        
        // 检查是否达到最大循环数
        if (line.find("Maximum number of steps exceeded") != std::string::npos ||
            line.find("The maximum number of steps has been exceeded") != std::string::npos) {
            std::cout << "⚠️  达到最大优化步数限制" << std::endl;
            break;
        }
        
        // 计算使用的循环数（可选）
        if (line.find("Step number") != std::string::npos && line.find("out of a maximum of") != std::string::npos) {
            std::istringstream iss(line);
            std::string word;
            while (iss >> word) {
                if (word == "number") {
                    iss >> cycles_used;
                    break;
                }
            }
        }
    }
    
    // 4) 如果优化成功，更新分子坐标
    if (stationary_point_found) {
        std::cout << "✅ Gaussian优化成功完成！" << std::endl;
        if (cycles_used > 0) {
            std::cout << "   使用优化步数: " << cycles_used << std::endl;
        }
        
        // 提取最终坐标并更新当前Cell对象
        Matrix final_coords = extract_final_coords_from_gauss_opt(log);
        if (final_coords.rows() == coords.rows()) {
            coords = final_coords;
            std::cout << "   已更新分子坐标为优化后的构型" << std::endl;
            return true;
        } else {
            std::cerr << "❌ 提取最终坐标失败：原子数不匹配" << std::endl;
        }
    } else {
        std::cout << "❌ Gaussian优化未收敛" << std::endl;
    }
    
    return false;
}

// ===== 从Gaussian优化输出中提取最终坐标 =====
Matrix Cell::extract_final_coords_from_gauss_opt(const std::string& log_file) const {
    std::ifstream fin(log_file);
    if (!fin) {
        std::cerr << "❌ Cannot open " << log_file << " for coordinate extraction" << std::endl;
        return Matrix(0, 0);
    }
    
    std::size_t nAtom = num_atoms();
    Matrix final_coords(nAtom, 3);
    
    std::string line;
    bool found_optimized_coordinates = false;
    
    // 寻找优化后的坐标
    // Gaussian通常在"Standard orientation"或"Input orientation"部分输出最终坐标
    while (std::getline(fin, line)) {
        if (line.find("Standard orientation:") != std::string::npos ||
            (line.find("Input orientation:") != std::string::npos && !found_optimized_coordinates)) {
            
            // 跳过表头
            for (int i = 0; i < 4; ++i) {
                std::getline(fin, line);
            }
            
            // 读取坐标
            for (std::size_t i = 0; i < nAtom; ++i) {
                if (!std::getline(fin, line)) break;
                
                std::istringstream iss(line);
                int center_number, atomic_number, atomic_type;
                double x, y, z;
                
                if (iss >> center_number >> atomic_number >> atomic_type >> x >> y >> z) {
                    // 将坐标从埃转换为Bohr
                    final_coords(i, 0) = x * BohrPerAngstrom;
                    final_coords(i, 1) = y * BohrPerAngstrom;
                    final_coords(i, 2) = z * BohrPerAngstrom;
                } else {
                    std::cerr << "❌ 解析坐标行失败: " << line << std::endl;
                    return Matrix(0, 0);
                }
            }
            
            found_optimized_coordinates = true;
            // 继续读取以获取最后的坐标（如果有多个orientation输出）
        }
    }
    
    if (!found_optimized_coordinates) {
        std::cerr << "❌ 未能在Gaussian输出中找到优化后的坐标" << std::endl;
        return Matrix(0, 0);
    }
    
    return final_coords;
}

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
    
    // 纯Gaussian几何优化功能
    bool pure_gaussian_optimize(const std::string& method = "hf/6-31G(d)", int max_cycles = 200);
    Matrix extract_final_coords_from_gauss_opt(const std::string& log_file = "gauss_opt.log") const;
};
