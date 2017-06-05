#include "basis.h"
#include "hartree_fock.h"
#include "coulomb_functions.h"
#include "ccd.h"
#include "ccd_copy.h"
#include "ccd_matrix.h"
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

    int number_of_states = first.number_of_states;
    int length_of_TBME = pow(number_of_states, 4);
    vec TBME = zeros<vec>(length_of_TBME);


    /*
    Hartree_fock_equations trial_1(shells, particles);
    vec energies = trial_1.hartree_fock_method(mapping, hw);

    energies.print();

    vec TBME = trial_1.new_basis();

    trial_1.total_energy();
    */

    //****************Reading TBME from file****************************

    ifstream myfile;
    myfile.open ("/Users/benedicte/Programs/fys4410/project_2/CCD/P2S2.txt");

    if (myfile.fail()) {
        cout << "Could not find file" << endl;
    }

    double value;

    for(int i = 0; i < number_of_states; i++){
        myfile >> value;
        TBME(i) = value;
    }


    myfile.close();




   // ************************ CCD PART 1   ************************


    //vec TBME = Examples::TwoParticleDotTest(shells, hw);

    //CCD_1 trial_copy(shells, particles, TBME, energies);

    //CCD_1 trial_copy(shells, particles, TBME, sp_energies);
    //trial_copy.CCD_solver(mapping);

    CCD_matrix trial_matrix(shells, particles, TBME, sp_energies);

    cout << "Success!" << endl;

    //trial_matrix.CCD_solver(mapping);


    // ************************ CCD PART 2   ************************


    //CCD trial(shells, particles, sp_energies);

    //double CCD_correlated_energy = trial.CCD_solver(mapping);


    return 0;

}
