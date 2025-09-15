#include <iostream>
#include <string>
#include <algorithm>
#include "cell.hpp"
#include "geomopt.hpp"
#include "rfo.hpp"

int main(int argc, char* argv[]) {
    if (argc < 2 || argc > 3) {
        std::cerr << "Usage: opt POSCAR [vasp|gauss|both]\n";
        return 1;
    }

    Cell cell = Cell::from_poscar(argv[1]);

    // Cell cell = Cell::from_poscar("/home/lbl8/code/hartree-fock_geomopt-rev-3/example/POSCAR-bak");
    GeomOpt opt(cell);
    bool success = opt.run();

    if (cell.is_converged_from_gauss("gauss.log")) {
        std::cout << "✅ Gaussian convergence table: ALL 4 criteria = YES.\n";
    } else {
        std::cout << "⚠️ Not converged per Gaussian criteria.\n";
    }
    return success ? 0 : 2;

}