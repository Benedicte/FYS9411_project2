#include "basis.h"
#include <iostream>

using namespace std;

Basis::Basis (int sh_number){

    shell_number = sh_number;

    for(int i = shell_number; i > 0; i--){
        number_of_states = number_of_states + i*2;
    }
}

vec Basis::get_energy(){

    double energy;
    sp_energies = zeros<vec>(number_of_states);

    for(int i = 0; i < number_of_states; i++){
      energy = (2*mapping(i,0) + abs(mapping(i,1)) + 1);//*hw;
      sp_energies(i) = energy;
    }

    for(int i = 0; i < number_of_states; i++){
        cout << "i " << i << " " << mapping(i,0) << " " << mapping(i,1) << " "<<  mapping(i,2) << " " << sp_energies(i)<<endl;
    }

    return sp_energies;
}

mat Basis::get_quantum_numbers(int shell){

/*
 *Based on the equation E = (2n - |m| + 1) where E is known from the shell number.
 *First we find all the n that can satisfy the equation, and then we find the
 *corresponding m. We then add them into a 2 dimensional array.
 *
 */
    mat quantum_numbers(shell,2);

    int counter = 0;

    for (int value_n = 0; value_n < (shell + 1)/2 ; value_n++){

        int value_m = shell - 1 - 2*value_n;

         if(value_m == 0){
            quantum_numbers(counter,0) = value_n;
            quantum_numbers(counter,1) = value_m;

            counter = counter + 1;
        } //End if

         if(value_m != 0 ){

            quantum_numbers(counter,0) = value_n;
            quantum_numbers(counter,1) = value_m;

            quantum_numbers(counter + 1,0) = value_n;
            quantum_numbers(counter + 1,1) = value_m * -1;

            counter = counter + 2;

        } // End if
    } // End for loop
    return quantum_numbers;
}

mat Basis::map_quantum_numbers(int number_of_particles){
/*
*We use the quantum numbers generated above, and
* make a mapping so that we have one mapping per state.
*/

    mapping = zeros<mat>(number_of_particles,3);
    int mapping_counter = 0;

    for(int i = 1; i <= shell_number; i++){

        mat quantum_numbers = Basis::get_quantum_numbers(i);

        for(int j = 0; j < i; j++){

            mapping(mapping_counter,0) = quantum_numbers(j,0);
            mapping(mapping_counter,1) = quantum_numbers(j,1);
            mapping(mapping_counter,2) = true;


            mapping_counter++;

            mapping(mapping_counter,0) = quantum_numbers(j,0);
            mapping(mapping_counter,1) = quantum_numbers(j,1);
            mapping(mapping_counter,2) = false;

            mapping_counter ++;

       } //end inner loop
    } //end outer loop

    return mapping;
}

