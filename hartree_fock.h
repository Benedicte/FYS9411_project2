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
    mat slater_determinant;
    mat fock_matrix;
    mat density_matrix;
    vec TBME;
    vec energies;

    double* coloumb_matrix_direct;
    double* coloumb_matrix_exchange;
    void GetCoulombIntegrals(mat mapping, double hw);
    mat hartree_fock_method(mat mapping);
    double get_energy(int m, int n, double omega);
    double kroneckerDelta(int x, int y);
    double total_energy(mat mapping, vec energies);
    double total_energy_test(mat mapping, vec energies, vec new_basis);
    vec new_basis(mat mapping, mat fock_matrix);
    vec GetTBME();

};

#endif // HARTREE_FOCK_EQUATIONS_H
