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

    Basis first(shells);
    mat mapping = first.map_quantum_numbers(first.number_of_states);
    vec sp_energies = first.get_energy();

    Hartree_fock_equations trial_1(shells, particles);
    trial_1.GetCoulombIntegrals(mapping, hw);
    vec energies = trial_1.hartree_fock_method(mapping);

    cout << "The method is done" << endl;

    double energy = trial_1.total_energy(mapping,energies);

    cout << "ground state energy" << energy << endl;

    //vec TBME = trial_1.GetTBME();

    //int TBME_length = trial_1.number_of_states*trial_1.number_of_states*trial_1.number_of_states*trial_1.number_of_states;

    //for (int i = 0; i < TBME_length; i++){
    //    if (TBME(i) != 0){
    //        cout << "not zero" << endl;
    //    }
    //}

    trial_1.total_energy(mapping, sp_energies);

    vec TBME_3 = trial_1.TBME;

    vec TBME_2 = Examples::TwoParticleDotTest();

    cout << "vec was made" << endl;

    CCD_1 trial_copy(shells, particles, TBME_3, sp_energies);

    trial_copy.CCD_solver(mapping);

    cout << "CCD_copy object was initialized" << endl;

    /*
    CCD trial(shells, particles, TBME_2, sp_energies);

    cout << "CCD object was initialized" << endl;

    trial.one_particle_energies_new_basis(sp_energies, mapping);

    cout <<"Amplitudes have been initalized" << endl;

    trial.CCD_solver(mapping);
    return 0;
    */


    //cout << "This is the total CCD energy " << total_energy + trial.correlated_energy_new << endl;
    //cout << "This is the new total CCD energy " << trial.CCD_energy_total(mapping) + trial.correlated_energy_new << endl;

    return 0;

}
