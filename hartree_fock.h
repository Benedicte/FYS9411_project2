#ifndef HARTREE_FOCK_EQUATIONS_H
#define HARTREE_FOCK_EQUATIONS_H
#include "coulomb_functions.h"
#include <armadillo>
#include <array>

using std::array;
using namespace arma;


class Hartree_fock_equations
{
public:
    Hartree_fock_equations(int sh_number, int n_states);

    int number_of_electrons;
    int number_of_states;
    int shell_number;
    mat slater_determinant;
    mat fock_matrix;
    mat density_matrix;
    vec TBME;
    vec twoBme;
    vec energies;
    vec energies_pre;
    double* coloumb_matrix_direct;
    double* coloumb_matrix_exchange;
    vec direct_elements;
    vec exchange_elements;


    void GetCoulombIntegrals(mat mapping, double hw);
    mat hartree_fock_method(mat mapping, double hw);
    vec get_energies_harmonic(mat mapping, double hw);
    double get_energy(int m, int n, double omega);
    double kroneckerDelta(int x, int y);
    double total_energy();
    double total_energy_test(mat mapping, vec energies, vec new_basis);
    vec new_basis();
    double one_particle_energies_new_basis(vec TBME_new);

    double KroneckerDelta(int i, int j);
    vec TBME_new(int shells, double hw, mat mapping);

};

#endif // HARTREE_FOCK_EQUATIONS_H
