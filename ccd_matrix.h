#ifndef CCD_MATRIX_H
#define CCD_MATRIX_H

#include <iostream>
#include <armadillo>
#include <array>
#include <abstract_coulomb.h>
#include <coulomb_function2.h>

using std::array;
using namespace arma;
using namespace std;

class CCD_matrix
{
public:

    CCD_matrix(int sh_number, int fermi_lv, vec TBME, vec SP_energies);

    int number_of_states;
    int fermi_level;
    int n_holes;
    int n_particles;
    vec TBME;
    vec sp_energies;
    vec sp_energies_HF;
    double correlated_energy_old;
    double correlated_energy_new;

    //********  Intermediates  ****************

    mat I_1;
    mat I_2;
    mat I_3;
    mat I_4;

    mat Vab_ij;
    mat Vab_cd;
    mat Vkl_cd;
    mat Vkl_bj;
    mat Vkl_ij;
    mat Vkc_ld;
    mat Vkc_bj;

    //*********  Functions   *******************

    //vec intial_amplitudes(mat mapping);
    vec intial_amplitudes_old();
    mat intial_amplitudes_matrix();
    void one_particle_energies(vec sp_energies);
    vec one_particle_energies_new_basis(vec sp_energies);
    void CCD_update(vec amplitudes_old, vec &amplitudes_new);
    mat CCD_update_matrix(mat amplitudes_old_mat);
    void make_V_matrix();
    double CCD_solver(mat mapping);

    double CCD_energy(mat amplitudes);
    double CCD_energy_vec(vec amplitudes);

    double CCD_energy_total();
    inline int index(int p, int q, int r, int s);

    double intermediate_1(int k, int l, int i, int j, vec amplitudes);
    double intermediate_2(int j, int k, int b, int c, vec amplitudes);
    double intermediate_3(int j, int l, vec amplitudes);
    double intermediate_4(int b, int c, vec amplitudes);

    mat intermediate_1_mat(mat amplitudes);

};

#endif // CCD_H
