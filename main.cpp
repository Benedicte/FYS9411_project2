#include "basis.h"
#include "hartree_fock.h"
#include "coulomb_functions.h"
#include "ccd.h"
#include "ccd_copy.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <armadillo>
#include <typeinfo>
#include <chrono>
#include <examples.h>
#include <abstract_coulomb.h>
#include <coulomb_function2.h>

using namespace std;
using namespace arma;

typedef std::chrono::high_resolution_clock Clock;

int main()
{
    int particles = 2;
    int shells = 3;
    double hw = 1;
    //int basisFunctions = shells*(shells+1);


    // ************************ HARTREE-FOCK PART    ************************
    Basis first(shells);
    mat mapping = first.map_quantum_numbers(first.number_of_states);
    vec sp_energies = first.get_energy();

    /*
    Hartree_fock_equations trial_1(shells, particles);
    trial_1.GetCoulombIntegrals(mapping, hw);
    vec energies = trial_1.hartree_fock_method(mapping);

    double ref_energy_HF = trial_1.total_energy(mapping,energies);

    cout << "ground state energy " << ref_energy_HF << endl;

    trial_1.total_energy(mapping, sp_energies);
    */

   // ************************ CCD PART 1   ************************


   // vec TBME_3 = trial_1.TBME;

    vec TBME_2 = Examples::TwoParticleDotTest(shells);


    CCD_1 trial_copy(shells, particles, TBME_2, sp_energies);


    trial_copy.CCD_solver(mapping);



    // ************************ CCD PART 2   ************************


    //CCD trial(shells, particles, sp_energies);

    //double CCD_correlated_energy = trial.CCD_solver(mapping);

    //cout << "This is the total CCD energy " << ref_energy_HF + CCD_correlated_energy << endl;
    //cout << "This is the new total CCD energy " << trial.CCD_energy_total(mapping) + trial.correlated_energy_new << endl;

    return 0;

}
