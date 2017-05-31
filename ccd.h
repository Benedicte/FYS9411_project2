#ifndef CCD_H
#define CCD_H

#include <iostream>
#include <armadillo>
#include <array>
#include <abstract_coulomb.h>
#include <coulomb_function2.h>

using std::array;
using namespace arma;
using namespace std;

class CCD
{
public:

    CCD(int sh_number, int fermi_lv, vec TBME, vec SP_energies);

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
    vec CCD_update(mat mapping, vec amplitudes_old);
    double CCD_solver(mat mapping);
    double CCD_energy(mat mapping, vec amplitudes);
    double CCD_energy_total(mat mapping);
    double TBMEl(double p, double q, double r, double s, mat mapping);
    int index(int p, int q, int r, int s);

    //double computeTau(int a, int b, int i, int j);
    //double computeW1(int m, int n, int i, int j);
    //double computeW2(int a, int b, int e, int f);
    //double computeW3(int m, int b, int e, int j);
    //double computeF1(int m, int e);
    //double computeF2(int m, int i);
    //double computeF3(int a, int e);

};

#endif // CCD_H
