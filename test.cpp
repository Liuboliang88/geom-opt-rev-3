#include "krks.hpp"
#include "kpotset.hpp"
#include "utility.hpp"
#include "wfn_utility.hpp"
#include "pyscf-base.hpp"
#include "pyscf-int.hpp"
#include "pyscf-pbc.hpp"
#include "ftdf.hpp"
#include "gfdf.hpp"
#include "timer.hpp"
#include "control.hpp"
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <cmath>

static bool usePyscfJK = true;


static std::size_t cnt_nonzero(const MatrixZ& m, double cut = 1e-8) {
    std::size_t n = 0;
    for (std::size_t i = 0; i < m.size(); ++i) {
        if (std::abs(m(i)) > cut) n++;
    }

    return n;
}


KRKS::KRKS(const Cell& cell, const KpotSet& kpts,const BasisSet& bset, 
           const PsudoSet& pset, const std::string& psdo, const std::string& xcfn)
: cell(cell), kpts(kpts), bset(bset), pset(pset), psdo(psdo), 
  xcfn(xcfn), hybrid(is_hybrid(xcfn)), hcoeff(hybrid_coeff(xcfn))
{
    std::size_t nkpot = kpts.size();
    std::size_t naofn = bset.n_aofn();

    ovlp_S = std::vector<MatrixZ>(nkpot, MatrixZ(naofn, naofn));
    kine_T = std::vector<MatrixZ>(nkpot, MatrixZ(naofn, naofn));
    extn_V = std::vector<MatrixZ>(nkpot, MatrixZ(naofn, naofn));
    core_H = std::vector<MatrixZ>(nkpot, MatrixZ(naofn, naofn));
    sqrt_A = std::vector<MatrixZ>(nkpot, MatrixZ(naofn, naofn));
    vxcl_V = std::vector<MatrixZ>(nkpot, MatrixZ(naofn, naofn));
    coul_J = std::vector<MatrixZ>(nkpot, MatrixZ(naofn, naofn));
    xchg_K = std::vector<MatrixZ>(nkpot, MatrixZ(naofn, naofn));
    veff_V = std::vector<MatrixZ>(nkpot, MatrixZ(naofn, naofn));
    fock_F = std::vector<MatrixZ>(nkpot, MatrixZ(naofn, naofn));
    last_F = std::vector<MatrixZ>(nkpot, MatrixZ(naofn, naofn));
    wfns_C = std::vector<MatrixZ>(nkpot, MatrixZ(naofn, naofn));
    dens_D = std::vector<MatrixZ>(nkpot, MatrixZ(naofn, naofn));
    last_D = std::vector<MatrixZ>(nkpot, MatrixZ(naofn, naofn));
    erro_E = std::vector<MatrixZ>(nkpot, MatrixZ(naofn, naofn));

    eige_E = std::vector<Matrix>(nkpot, Matrix(1, naofn));
    occu_O = std::vector<Matrix>(nkpot, Matrix(1, naofn));

    diis_f = std::vector<Complex>(nkpot * naofn * naofn);
    diis_r = std::vector<Complex>(nkpot * naofn * naofn);

    std::cout << "hybrid coeff: " << hybrid_coeff(xcfn) << std::endl;
    std::cout << std::endl;

    std::cout << "number of shell funciotns: " << bset.size() << std::endl;
    std::cout << "number of basis functions: " << bset.n_aofn() << std::endl;
    std::cout << "cell mesh: " << get_cell_mesh(bset, cell.cellv).to_string() << std::endl;


    // init density fitting
    const std::vector<Vec3f64>& cellv = cell.cellv;
    assert(cellv.size() == 3);

    const double xlen = cellv[0].len();
    const double ylen = cellv[1].len();
    const double zlen = cellv[2].len();

    // interval of grid;
    const double dlen = 0.20;
    assert(xlen > dlen);
    assert(ylen > dlen);
    assert(zlen > dlen);

    std::size_t nx = std::size_t(xlen / dlen);
    std::size_t ny = std::size_t(ylen / dlen);
    std::size_t nz = std::size_t(zlen / dlen);

    if (nx % 2) nx -= 1;
    if (ny % 2) ny -= 1;
    if (nz % 2) nz -= 1;

    Vec3u64 fftmesh(nx, ny, nz);
    std::cout << "fft mesh: " << fftmesh.to_string() << std::endl;
    std::cout << std::endl;

    diis = DIIS(SCF_diis_space);
    // ftdf = FDF(bset, pset, kpts, fftmesh, cell.cellv);

    std::vector<Vec3f64> imgls = make_imgls(bset, cell.cellv);  

    Timer t0("KRKS::int1e");

    bset.pbc_int1e_ovlp_d00_kpts(kpts, imgls, ovlp_S);
    bset.pbc_int1e_kine_d00_kpts(kpts, imgls, kine_T);

    if (psdo.empty()) {
        throw std::runtime_error("only support pseudopotential calculation.");
    }

    psdo_d00_kpts_pyscf(bset, psdo, cell.cellv, kpts, extn_V);

    t0.print_time_used();

    // core Hamiltonian
    for (std::size_t ik = 0; ik < nkpot; ++ik) {
        core_H[ik] = kine_T[ik] + extn_V[ik];
    }


    /**************************************************************************************/
    /*                               ckeck Core Hamiltonian                               */
    /**************************************************************************************/

    // std::vector<MatrixZ> core_from_pyscf = std::vector<MatrixZ>(nkpot, MatrixZ(naofn, naofn));
    // core_d00_kpts_pyscf(bset, psdo, cell.cellv, kpts, core_from_pyscf);

    // const double absErr = 1e-2;
    // for (std::size_t ik = 0; ik < nkpot; ++ik) {
    //     for (std::size_t i = 0; i < naofn; ++i) {
    //     for (std::size_t j = 0; j < naofn; ++j) {
    //         double absDiff = std::abs(core_from_pyscf[ik](i,j) - core_H[ik](i,j));
    //         if (absDiff > absErr) {
    //             std::cerr << "Core Hamiltonian is not match." << std::endl;

    //             std::cout << "ik = " << ik << std::endl;
    //             std::cout << " i = " <<  i << std::endl;
    //             std::cout << " j = " <<  j << std::endl;

    //             std::cerr << core_from_pyscf[ik](i,j) << std::endl;
    //             std::cerr << core_H[ik](i,j) << std::endl;
    //             std::abort();
    //         }
    //     }}
    // }

    // std::cout << "Core Hamiltonian is checked." << std::endl;
    

    // calculate A = S^(-1/2)
    for (std::size_t ik = 0; ik < nkpot; ++ik) {
        sqrt_A[ik] = matrixz_pow(ovlp_S[ik], -0.5);
    }

    // // init density fitting
    // const std::vector<Vec3f64>& cellv = cell.cellv;
    // assert(cellv.size() == 3);

    // const double xlen = cellv[0].len();
    // const double ylen = cellv[1].len();
    // const double zlen = cellv[2].len();

    // // interval of grid;
    // const double dlen = 0.20;
    // assert(xlen > dlen);
    // assert(ylen > dlen);
    // assert(zlen > dlen);

    // std::size_t nx = std::size_t(xlen / dlen);
    // std::size_t ny = std::size_t(ylen / dlen);
    // std::size_t nz = std::size_t(zlen / dlen);

    // if (nx % 2) nx -= 1;
    // if (ny % 2) ny -= 1;
    // if (nz % 2) nz -= 1;

    // Vec3u64 fftmesh(nx, ny, nz);
    
    // std::cout << "fft mesh: " << fftmesh.to_string() << std::endl;
    // fdf = FDF(bset, pset, kpts, fftmesh, cell.cellv);

    if (! usePyscfJK) {
        gfdf = GFDF(hybrid, kpts, bset, bset, pset, cell.cellv);
        mydf = &gfdf;  // set mydf to gfdf
    }
}


// in pyscf-int.hpp

// *************************************************************************************************
//                                             int2e
// *************************************************************************************************
void coul_cart_kpts_pyscf(const BasisSet& bs, const std::vector<Vec3f64> &cellv, const KpotSet& kpots, const std::vector<MatrixZ>& dms, std::vector<MatrixZ>& mat);
void coul_sphl_kpts_pyscf(const BasisSet& bs, const std::vector<Vec3f64> &cellv, const KpotSet& kpots, const std::vector<MatrixZ>& dms, std::vector<MatrixZ>& mat);
void xchg_cart_kpts_pyscf(const BasisSet& bs, const std::vector<Vec3f64> &cellv, const KpotSet& kpots, const std::vector<MatrixZ>& dms, std::vector<MatrixZ>& mat);
void xchg_sphl_kpts_pyscf(const BasisSet& bs, const std::vector<Vec3f64> &cellv, const KpotSet& kpots, const std::vector<MatrixZ>& dms, std::vector<MatrixZ>& mat);


void coul_xchg_cart_kpts_pyscf(const BasisSet& bs, const std::vector<Vec3f64> &cellv, 
                               const KpotSet& kpots, const std::vector<MatrixZ>& dms, 
                               std::vector<MatrixZ>& coul, std::vector<MatrixZ>& xchg);

void coul_xchg_sphl_kpts_pyscf(const BasisSet& bs, const std::vector<Vec3f64> &cellv, 
                               const KpotSet& kpots, const std::vector<MatrixZ>& dms, 
                               std::vector<MatrixZ>& coul, std::vector<MatrixZ>& xchg);

// *************************************************************************************************
//                                           eval_xc
// *************************************************************************************************
void eval_xc_kpts_pyscf(const BasisSet& bs, const std::vector<Vec3f64> &cellv, 
                        const KpotSet& kpots, const std::vector<MatrixZ>& dms, 
                        const std::string& xc, std::vector<MatrixZ>& mat, double& energy);


// *************************************************************************************************
//                                encode and decode diis vector
// *************************************************************************************************                        
static void vector_matrixz_to_vector_complex(const std::vector<MatrixZ>& mats, std::vector<Complex>& list) {
    std::size_t nkpot = mats.size();
    std::size_t naofn = mats.front().rows();

    assert(list.size() == nkpot*naofn*naofn);

    for (std::size_t ik = 0; ik < nkpot; ++ik) {
        assert(mats[ik].rows() == naofn);
        assert(mats[ik].cols() == naofn);

        for (std::size_t i = 0; i < naofn; ++i) {
        for (std::size_t j = 0; j < naofn; ++j) {
            list[ik*naofn*naofn + i*naofn + j] = mats[ik](i,j);
        }}
    }
}

static void vector_complex_to_vector_matrixz(const std::vector<Complex>& list, 
                                                   std::vector<MatrixZ>& mats) {
    std::size_t nkpot = mats.size();
    std::size_t naofn = mats.front().rows();

    assert(list.size() == nkpot*naofn*naofn);

    for (std::size_t ik = 0; ik < nkpot; ++ik) {
        assert(mats[ik].rows() == naofn);
        assert(mats[ik].cols() == naofn);

        for (std::size_t i = 0; i < naofn; ++i) {
        for (std::size_t j = 0; j < naofn; ++j) {
            mats[ik](i,j) = list[ik*naofn*naofn + i*naofn + j];
        }}
    }
}

bool KRKS::scf() {
    std::size_t loop = 0;
    bool converged = false;
    
    std::size_t nkpot = kpts.size();
    std::size_t naofn = bset.n_aofn();

    Timer tGuess("KRKS::scf::init_guess");
    generate_initial_density_matrix("atom");
    tGuess.print_time_used();

    if (!SCF_chk_load.empty()) {
        load_wavefunction(dens_D, SCF_chk_load);
    }

    std::size_t nValenceElec = pset.n_valence_electron();
    std::cout << "number of electron: " << nValenceElec << std::endl;

    /******************************************************************************/
    /*                              calc Fock                                     */
    /******************************************************************************/
    eval_xc_kpts_pyscf(bset, cell.cellv, kpts, xcfn, dens_D, vxcl_V, _exc);

    if (usePyscfJK) {
        coul_kpts_pyscf(bset, cell.cellv, kpts, dens_D, coul_J);
        if (hybrid) xchg_kpts_pyscf(bset, cell.cellv, kpts, dens_D, xchg_K);
    } else {
        mydf -> pbc_int2e_coul_kpts(dens_D, coul_J);
        if (hybrid) mydf -> pbc_int2e_xchg_kpts(dens_D, xchg_K);
    }

    for (std::size_t ik = 0; ik < nkpot; ++ik) {
        veff_V[ik] = vxcl_V[ik] + coul_J[ik] - 0.5 * hcoeff * xchg_K[ik];
        fock_F[ik] = core_H[ik] + veff_V[ik];
    }

    // /******************************************************************************/
    // /*                              solve eigen                                   */
    // /******************************************************************************/
    // for (std::size_t ik = 0; ik < nkpot; ++ik) {
    //     GeneralSelfAdjointEigenSolverZ
    //     gesz(fock_F[ik], ovlp_S[ik]);

    //     wfns_C[ik] = gesz.eigen_vec();
    //     MatrixZ ev =  gesz.eigen_val();
    //     assert(ev.size() == naofn);

    //     for (std::size_t ia = 0; ia < naofn; ++ia) {
    //         eige_E[ik](ia) = ev(ia).real();
    //     }

    //     occu_O[ik] = make_occu(eige_E[ik], nValenceElec);
    //     // dens_D[ik] = make_rdm1(wfns_C[ik], occu_O[ik]);
    //     dens_D[ik] = make_rdm1(wfns_C[ik], nValenceElec);
    // }

    eCore = ::energy_core(dens_D, core_H);
    eCoul = ::energy_coul(dens_D, coul_J);
    
    // double eExclKylin = this -> energy_excl();
    eExcl = _exc;
    // be careful, there is 0.5*0.5
    if (hybrid) { eExcl -= ::energy_xchg(dens_D, xchg_K) * hcoeff; }
    
    eNucl = this -> energy_nucl();
    eTotl = eCore + eCoul + eExcl + eNucl;

    std::cout << std::endl;
    std::printf("step = %ld, kylin: eTotl = %.8lf, eCore = %.8lf, "
                "eCoul = %.8lf, eExcl = %.8lf, eNucl = %.8lf\n", 
                std::size_t(0), eTotl, eCore, eCoul, eExcl, eNucl);

    std::cout << std::endl;

    while (true) {
        if ((++loop) > SCF_max_cycle) break;

        Timer ti("SCF loop step " + std::to_string(loop));

        double eLast = eTotl;
        last_D = dens_D;
        last_F = fock_F;


        /******************************************************************************/
        /*                               calc fock                                    */
        /******************************************************************************/
        eval_xc_kpts_pyscf(bset, cell.cellv, kpts, xcfn, dens_D, vxcl_V, _exc);

        if (usePyscfJK) {
            coul_kpts_pyscf(bset, cell.cellv, kpts, dens_D, coul_J);
            if (hybrid) xchg_kpts_pyscf(bset, cell.cellv, kpts, dens_D, xchg_K);
        } else {
            mydf -> pbc_int2e_coul_kpts(dens_D, coul_J);
            if (hybrid) mydf -> pbc_int2e_xchg_kpts(dens_D, xchg_K);
        }


        // std::cout << "num of AO: " << bset.n_aofn() << std::endl;

        // std::cout << "J cnt = 1e-8: ";
        // for (std::size_t ik = 0; ik < nkpot; ++ik) {
        //     std::cout << cnt_nonzero(coul_J[ik]) << " ";
        // }
        // std::cout << std::endl;

        // std::cout << "K cnt = 1e-8: ";
        // for (std::size_t ik = 0; ik < nkpot; ++ik) {
        //     std::cout << cnt_nonzero(xchg_K[ik]) << " ";
        // }
        // std::cout << std::endl;

        for (std::size_t ik = 0; ik < nkpot; ++ik) {
            veff_V[ik] = vxcl_V[ik] + coul_J[ik] - 0.5 * hcoeff * xchg_K[ik];
            fock_F[ik] = core_H[ik] + veff_V[ik];
        }


        /******************************************************************************/
        /*                                  DIIS                                      */
        /******************************************************************************/
        if (SCF_diis_method != "NONE" && loop >= SCF_diis_start) {
            for (std::size_t ik = 0; ik < nkpot; ++ik) {
                erro_E[ik] = make_error_vector(fock_F[ik], last_D[ik],
                                               ovlp_S[ik], sqrt_A[ik]);
            }

            // push fock vector and error vector to diis space
            vector_matrixz_to_vector_complex(fock_F, diis_f);
            vector_matrixz_to_vector_complex(erro_E, diis_r);
            diis.update_diis_space(diis_f, diis_r);
            
            if (diis.diis_space_ready()) {
                std::vector<Complex> newFock = diis.extrapolate_new_vector();
                vector_complex_to_vector_matrixz(newFock, fock_F);
            } 
            
            // else if (std::abs(SCF_damp_factor) > 1e-4) {
            //     for (std::size_t ik = 0; ik < nkpot; ++ik) {
            //         fock_F[ik] = damping(fock_F[ik], last_F[ik], SCF_damp_factor);
            //     }
            // }
        }

        // /*************************************************************************************/
        // /*                                    Damping                                        */
        // /*************************************************************************************/
        // if (std::abs(SCF_damp_factor) > 1e-4 && loop > 1 && loop < SCF_diis_start) {
        //     for (std::size_t ik = 0; ik < nkpot; ++ik) {
        //         fock_F[ik] = damping(fock_F[ik], last_F[ik], SCF_damp_factor);
        //     }
        // }

        /******************************************************************************/
        /*                              level shift                                   */
        /******************************************************************************/
        if (std::abs(SCF_level_shift) > 1e-4) {
            for (std::size_t ik = 0; ik < nkpot; ++ik) {
                fock_F[ik] = level_shift(fock_F[ik], last_D[ik], 
                                         ovlp_S[ik], SCF_level_shift);
            }
        }

        /******************************************************************************/
        /*                              solve eigen                                   */
        /******************************************************************************/
        Timer t2("KRKS::scf::eigensolver");

        for (std::size_t ik = 0; ik < nkpot; ++ik) {
            GeneralSelfAdjointEigenSolverZ
            gesz(fock_F[ik], ovlp_S[ik]);

            wfns_C[ik] = gesz.eigen_vec();

            MatrixZ ev =  gesz.eigen_val();
            assert(ev.size() == naofn);

            for (std::size_t ia = 0; ia < naofn; ++ia) {
                eige_E[ik](ia) = ev(ia).real();
            }

            occu_O[ik] = make_occu(eige_E[ik], nValenceElec);
            // dens_D[ik] = make_rdm1(wfns_C[ik], occu_O[ik]);
            // dens_D[ik] = make_rdm1(wfns_C[ik], nValenceElec);
        }

        // make occupation
        if (SCF_smear_method == "NONE") { // no smearing
            for (std::size_t ik = 0; ik < nkpot; ++ik) {
                occu_O[ik] = make_occu(eige_E[ik], nValenceElec);
            }
            std::cout << "make_occu: normal." << std::endl;
            // make_occu_kpts(eige_E, nValenceElec, occu_O);
        } else if (SCF_smear_method == "GAUSS"){
            std::cout << "make_occu: gauss." << std::endl;
            occu_O = make_occu_gauss_kpts(eige_E, nValenceElec, SCF_smear_sigma);
        } else if (SCF_smear_method == "METH") {
            std::cout << "make_occu: meth." << std::endl;
            occu_O = make_occu_methf_kpts(eige_E, nValenceElec, SCF_smear_sigma);
        } else if (SCF_smear_method == "FERMI") {
            std::cout << "make_occu: fermi." << std::endl;
            occu_O = make_occu_fermi_kpts(eige_E, nValenceElec, SCF_smear_sigma);
        } else {
            std::cerr << "don't support this type of smearing." << std::endl;
            std::abort();
        }

        // make rdm
        for (std::size_t ik = 0; ik < nkpot; ++ik) {
            dens_D[ik] = make_rdm1(wfns_C[ik], occu_O[ik]);
        }

        t2.print_time_used();

        Timer tx("KRKS::scf::save_dm");
        if (!SCF_chk_save.empty() && loop % 5 == 0) {
            void save_wavefunction(const std::vector<MatrixZ>& dm, const std::string& filename);
            save_wavefunction(dens_D, SCF_chk_save);
        }
        tx.print_time_used();


        /*************************************************************************************/
        /*                    calculate energy in **new** density matrix                     */
        /*************************************************************************************/
        Timer t3("KRKS::scf::energy");

        eCore = ::energy_core(dens_D, core_H);
        eCoul = ::energy_coul(dens_D, coul_J);
        
        // double eExclKylin = this -> energy_excl();
        eExcl = _exc;
        // be careful, there is 0.5*0.5
        if (hybrid) { eExcl -= ::energy_xchg(dens_D, xchg_K) * hcoeff; }
        
        eNucl = this -> energy_nucl();
        eTotl = eCore + eCoul + eExcl + eNucl;

        t3.print_time_used();

        /*************************************************************************************/
        /*                    calculate energy in **new** density matrix                     */
        /*************************************************************************************/

        std::printf("step = %ld, kylin: eTotl = %.8lf, eCore = %.8lf, "
                    "eCoul = %.8lf, eExcl = %.8lf, eNucl = %.8lf\n", 
                    loop, eTotl, eCore, eCoul, eExcl, eNucl);

        double enDiff = eTotl - eLast;
        double dmDiff = density_matrix_diff(last_D, dens_D);
        std::printf("step = %ld, de = %.8lf, ddm = %.8lf\n", loop, enDiff, dmDiff);

        if (std::abs(SCF_level_shift) > 1e-4 && std::abs(enDiff) < SCF_eps_energy 
                                             && std::abs(dmDiff) < SCF_eps_dmat) {
            SCF_level_shift = 0.0;
            std::cout << std::endl;
            continue;
        }


        if (std::abs(enDiff) < SCF_eps_energy && std::abs(dmDiff) < SCF_eps_dmat) {
            std::printf("SCF is converged. E = %.8lf\n", eTotl);
            std::cout << std::endl;

            converged = true;
            return true;
        }

        ti.print_time_used();

        std::cout << std::endl;

    }

    std::printf("SCF didn't converge.\n");
    return false;
}



/*************************************************************************************/
/*                               kpt_band calculation                                */
/*************************************************************************************/
void KRKS::band(const KBand& kb) {
    std::size_t nband = kb.size();
    std::size_t naofn = bset.n_aofn();

    KpotSet kbs;
    kbs.rcplv = kb.rcplv;
    kbs.kpots = kb.kpots;

    std::vector<MatrixZ> ovlp_S_band(nband, MatrixZ(naofn, naofn));
    std::vector<MatrixZ> kine_T_band(nband, MatrixZ(naofn, naofn));
    std::vector<MatrixZ> extn_V_band(nband, MatrixZ(naofn, naofn));
    std::vector<MatrixZ> core_H_band(nband, MatrixZ(naofn, naofn));
    std::vector<MatrixZ> vxcl_V_band(nband, MatrixZ(naofn, naofn));
    std::vector<MatrixZ> coul_J_band(nband, MatrixZ(naofn, naofn));
    std::vector<MatrixZ> xchg_K_band(nband, MatrixZ(naofn, naofn));
    std::vector<MatrixZ> veff_V_band(nband, MatrixZ(naofn, naofn));
    std::vector<MatrixZ> fock_F_band(nband, MatrixZ(naofn, naofn));

    std::vector<Vec3f64> imgls = make_imgls(bset, cell.cellv); 

    bset.pbc_int1e_ovlp_d00_kpts(kbs, imgls, ovlp_S_band);
    bset.pbc_int1e_kine_d00_kpts(kbs, imgls, kine_T_band);

    if (psdo.empty()) {
        throw std::runtime_error("only support pseudopotential calculation.");
    }

    psdo_d00_kpts_pyscf(bset, psdo, cell.cellv, kbs, extn_V_band);

    for (std::size_t ik = 0; ik < nband; ++ik) {
        core_H_band[ik] = kine_T_band[ik] + extn_V_band[ik];
    }

    // Timer t0("KRKS::scf::eval_xc(first)");
    eval_xc_kband_pyscf(bset, cell.cellv, kpts, kb, xcfn, dens_D, vxcl_V_band, _exc);
    // t0.print_time_used();

    // Timer t1("KRKS::scf::GFDF::pbc_int2e_coul_kpts(first)");
    gfdf.pbc_int2e_coul_kband(dens_D, kb, coul_J_band);
    // t1.print_time_used();

    // if (hybrid) {
    //     // Timer t2("KRKS::scf::GFDF::pbc_int2e_xchg_kpts(first)");
    //     gfdf.pbc_int2e_xchg_kpts(dens_D, xchg_K);
    //     // t2.print_time_used();

    //     for (std::size_t ik = 0; ik < nkpot; ++ik) {
    //         xchg_K[ik] *= hcoeff;
    //     }
    // }

    for (std::size_t ik = 0; ik < nband; ++ik) {
        veff_V_band[ik] = vxcl_V_band[ik] + coul_J_band[ik] - 0.5 * xchg_K_band[ik];
        fock_F_band[ik] = core_H_band[ik] + veff_V_band[ik];
    }


    std::ofstream ofs("band_energy.txt");

    for (std::size_t ik = 0; ik < nband; ++ik) {
        GeneralSelfAdjointEigenSolverZ
        gesz(fock_F_band[ik], ovlp_S_band[ik]);

        // wfns_C[ik] = gesz.eigenVec();

        MatrixZ ev =  gesz.eigen_val();
        assert(ev.size() == naofn);

        for (std::size_t ia = 0; ia < naofn; ++ia) {
            // eige_E[ik](ia) = ev(ia).real();

            std::printf("%15.8e ", ev(ia).real());
        }

        std::printf("\n");

        // occu_O[ik] = make_occu(eige_E[ik], nValenceElec);
        // dens_D[ik] = make_rdm1(wfns_C[ik], occu_O[ik]);
        // // dens_D[ik] = make_rdm1(wfns_C[ik], nValenceElec);
    }


}


double KRKS::density_matrix_diff(const std::vector<MatrixZ>& dm1,
                                 const std::vector<MatrixZ>& dm2) const {
    assert(dm1.size() == dm2.size());
    std::size_t nkpot = dm1.size();

    double diffMax = 0.0;
    for (std::size_t k = 0; k < nkpot; ++k) {
        assert(dm1[k].rows() == dm2[k].rows());
        assert(dm1[k].cols() == dm2[k].cols());

        double diffNorm = (dm1[k] - dm2[k]).norm();
        diffMax = std::max(diffMax, diffNorm);
    }

    return diffMax;
}

// this function will generate init dm, write into last_D
void KRKS::generate_initial_density_matrix(const std::string& type) {
    std::size_t nkpot = kpts.size();
    std::size_t naofn = bset.n_aofn();

    last_D = std::vector<MatrixZ>(nkpot, MatrixZ(naofn, naofn));
    dens_D = std::vector<MatrixZ>(nkpot, MatrixZ(naofn, naofn));
    
    if (type == "core") { return ; }

    if (type == "atom") {
        guess_kpts_pyscf(bset, psdo, cell.cellv, xcfn, kpts, dens_D);
        return ;
    }
    
    if (type == "rand") {
        // TODO
    }

    std::cerr << "Don't support initial guess method: \"" << type << "\"." << std::endl;
    std::cerr << "Use one guess method in [core, atom, rand]." << std::endl;
    std::abort();
}

/******************************************************************************/
/*                            calculate energy                                */
/******************************************************************************/

double KRKS::energy_nucl() const {
    return energy_nucl_ewald_pyscf(bset, psdo, cell.cellv);
}

// double KRKS::energy_excl() const {
//     double excl = 0.0;
//     std::size_t nkpot = kpts.size();
//     std::size_t naofn = bset.n_aofn();
//     static std::vector<MatrixZ> _vxcl = std::vector<MatrixZ>(nkpot, MatrixZ(naofn, naofn));
//     eval_xc_kpts_pyscf(bset, cell.cellv, kpts, xcfn, last_D, _vxcl, excl);
//     return excl;
// }


