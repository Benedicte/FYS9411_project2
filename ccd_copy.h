#ifndef CCD_COPY_H
#define CCD_COPY_H

#include <iostream>
#include <armadillo>
#include <array>
#include <abstract_coulomb.h>
#include <coulomb_function2.h>

using std::array;
using namespace arma;
using namespace std;

class CCD_1
{
public:

    CCD_1(int sh_number, int fermi_lv, vec TBME, vec SP_energies);

    int number_of_states;
    int fermi_level;
    int n_holes;
    int n_particles;
    vec TBME;
    vec sp_energies;
    vec sp_energies_HF;
    double correlated_energy_old;
    double correlated_energy_new;

    //vec intial_amplitudes(mat mapping);
    vec intial_amplitudes_old(mat mapping);
    void one_particle_energies(vec sp_energies);
    vec one_particle_energies_new_basis(vec sp_energies, mat mapping);
    void CCD_update(mat mapping, vec amplitudes_old, vec &amplitudes_new);
    void CCD_update_matrix(vec amplitudes_old, vec &amplitudes_new);
    double CCD_solver(mat mapping);
    double CCD_energy(mat mapping, vec amplitudes);
    double CCD_energy_total(mat mapping);
    int index(int p, int q, int r, int s);

    double intermediate_1(int k, int l, int i, int j, vec amplitudes);
    double intermediate_2(int j, int k, int b, int c, vec amplitudes);
    double intermediate_3(int j, int l, vec amplitudes);
    double intermediate_4(int b, int c, vec amplitudes);


};

#endif // CCD_H
