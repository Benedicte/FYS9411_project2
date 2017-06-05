#include "ccd_matrix.h"
#include "coulomb_function2.h"
#include <iomanip>

using namespace std;

CCD_matrix::CCD_matrix(int shell_number, int fermi_lv, vec two_body_matrix_elements, vec single_particle_energies)
{
    fermi_level = fermi_lv;

    for(int i = shell_number; i > 0; i--){
        number_of_states = number_of_states + i*2;
    }

    TBME = two_body_matrix_elements;
    sp_energies = single_particle_energies;

}

void CCD_matrix::one_particle_energies(vec sp_energies){

    // rotate sp energies to HF energies

    double ei,ea;
    vec spEnergies_HF = zeros<vec>(number_of_states);

    for(int i = 0; i < fermi_level; i++){
        ei = sp_energies[i];
        for(int j = 0; j < fermi_level; j++){
            //ei += calcVpqrs(i,j,i,j,SPbasis); Fix this part. make sure corrrrect
        }
        spEnergies_HF(i) = ei;
    }
    for(int a = fermi_level; a < number_of_states; a++){
        ea = sp_energies[a];
        for(int j = 0; j < fermi_level; j++){
            //ea += calcVpqrs(a,j,a,j,SPbasis); Fix this part. make sure corrrrect
        }
        spEnergies_HF(a) = ea;
    }
}

vec CCD_matrix::one_particle_energies_new_basis(vec sp_energies){

    double ei,ea;

    for(int i= 0 ; i < fermi_level; i++){
        ei = sp_energies(i);
        for(int j = 0 ; j < fermi_level; j++){
            ei += TBME(index(i,j,i,j));
        };
        sp_energies(i) = ei;
    }

    for(int a = fermi_level; a < number_of_states; a++){
        ea = sp_energies(a);
        for(int j=0; j<fermi_level; j++){
            ea += TBME(index(a,j,a,j));
        }
        sp_energies(a) = ea;
    }

    sp_energies_HF = sp_energies;

    return sp_energies;
}

vec CCD_matrix::intial_amplitudes_old(){

    cout << number_of_states << endl;
    cout << fermi_level << endl;

    vec initial_amplitudes= zeros<vec>(number_of_states*number_of_states*number_of_states*number_of_states);

    n_holes = number_of_states - fermi_level;
    n_particles = fermi_level;

    cout << size(TBME)<< endl;
    cout << size(initial_amplitudes)<< endl;

    cout << "Size of SP-Energies" << size(sp_energies) << endl;

    for(int i = 0; i < n_particles; i++){
        for(int j = 0; j < n_particles; j++){
            for(int a = n_particles; a < number_of_states; a++){
                for(int b = n_particles; b < number_of_states; b++){

                    //cout << index(a,b,i,j) << endl;

                    initial_amplitudes(index(a,b,i,j)) = TBME(index(a,b,i,j))
                            /(sp_energies(i) + sp_energies(j) - sp_energies(a) - sp_energies(b));
                }
            }
        }
    }

    return initial_amplitudes;
}

mat CCD_matrix::intial_amplitudes_matrix(){

    // ********** In Matrix-multiplication format*************//

    n_holes = number_of_states - fermi_level;
    n_particles = fermi_level;

    mat T = zeros<mat>(n_holes*n_holes,n_particles*n_particles);

    for(int i = 0; i < n_particles; i++){
        for( int j = 0; j < n_particles; j++){
            for(int a = n_particles; a < number_of_states; a++){
                for(int b = n_particles; b < number_of_states; b++){
                    T((b-n_particles)*n_holes+(a-n_particles), j*n_particles+i) = TBME(index(a,b,i,j))
                            /(sp_energies(i) + sp_energies(j) - sp_energies(a) - sp_energies(b));
                }
            }
        }
    }

    return T;
}

void CCD_matrix::CCD_update(vec amplitudes_old, vec &amplitudes_new){
    double sum = 0;
    double energy_denom = 0;

    for(int i = 0; i < n_particles; i++){
        for(int j = 0; j < n_particles; j++){
            for(int a = n_particles; a < number_of_states; a++){
                for(int b = n_particles; b < number_of_states; b++){

                    sum = TBME(index(a,b,i,j)); //First term

                    for(int c = n_particles; c < number_of_states; c++){
                        for(int d = n_particles; d < number_of_states; d++){

                            sum += 0.5*TBME(index(a,b,c,d))*amplitudes_old(index(c,d,i,j)); //Second term
                        }
                    }

                    for(int k = 0; k < n_particles; k++){
                        for(int l = 0; l < n_particles; l++){
                            sum += 0.5*amplitudes_old(index(a,b,k,l))
                                    *intermediate_1(k,l,i,j,amplitudes_old); //Third term
                        }
                    }

                    for(int k = 0; k < n_particles; k++){
                        for(int c = n_particles; c < number_of_states; c++){
                            sum += amplitudes_old(index(a,c,i,k))
                                    *intermediate_2(k,b,c,j,amplitudes_old); //fourth term

                            //sum -= amplitudes_old(index(a,c,j,k))
                              //       *intermediate_2(k,b,c,i,amplitudes_old);
                            // sum -= amplitudes_old(index(b,c,i,k))
                            //         *intermediate_2(k,a,c,j,amplitudes_old);//fifth term switch a and b
                            // sum += amplitudes_old(index(b,c,j,k))
                            //         *intermediate_2(k,a,c,i,amplitudes_old);
                        }
                    }

                    /*
                    for(int k = 0; k < n_particles; k++){
                        sum -= amplitudes_old(index(a,b,i,k))
                                *intermediate_3(k,j,amplitudes_old);//eigth term

                        sum += amplitudes_old(index(a,b,j,k))
                                *intermediate_3(k,i,amplitudes_old);
                    }

                    for(int c = n_particles; c < number_of_states; c++){
                        sum += amplitudes_old(index(a,c,i,j))
                                *intermediate_4(b,c,amplitudes_old); //tenth term

                        sum -= amplitudes_old(index(b,c,i,j))
                                *intermediate_4(a,c,amplitudes_old);

                    }
*/
                    energy_denom = sp_energies(i) + sp_energies(j) - sp_energies(a) - sp_energies(b);
                    amplitudes_new(index(a,b,i,j)) = sum/energy_denom;

                }//End of main indexes
            }
        }
    }

    return;


} //End of function

mat CCD_matrix::CCD_update_matrix(mat amplitudes_old_mat){

    make_V_matrix();
    intermediate_1_mat(amplitudes_old_mat);

    //In the CCD_update

    mat T_2 = zeros<mat>(n_holes*n_holes, n_particles*n_particles);
    mat T_3 = zeros<mat>(n_holes*n_holes, n_particles*n_particles);
    mat Tai_kc = zeros<mat>(n_holes*n_particles, n_holes*n_particles);
    mat T_5 = zeros<mat>(n_holes*n_particles, n_holes*n_particles);


    for(int i = 0; i < n_particles; i++){
        for( int j = 0; j < n_particles; j++){
            for(int a = n_particles; a < number_of_states; a++){
                for(int b = n_particles; b < number_of_states; b++){
                    for(int c = n_particles; c < number_of_states; c++){
                        for(int k = 0; k < n_particles; k++){

                            Tai_kc(i*n_holes+(a-n_particles),(c-n_particles)*n_particles+k)
                                    = amplitudes_old_mat((c-n_particles)*n_holes+(a-n_particles), k*n_particles+i)
                                    /(sp_energies(i) + sp_energies(j) - sp_energies(a) - sp_energies(b));
                            T_5(j*n_holes+(a-n_particles),(c-n_particles)*n_particles+k)
                                    = amplitudes_old_mat((c-n_particles)*n_holes+(a-n_particles), k*n_particles+j)
                                    /(sp_energies(i) + sp_energies(j) - sp_energies(a) - sp_energies(b));
                        }

                        for (int d = n_particles; d < number_of_states; d++){
                            T_2((d-n_particles)*n_holes+(c-n_particles), j*n_particles+i)
                                    = amplitudes_old_mat((d-n_particles)*n_holes+(c-n_particles), j*n_particles+i);
                        }

                        for(int k = 0; k < n_particles; k++){
                            for(int l = 0; l < n_particles; l++){
                                T_3((b-n_particles)*n_holes+(a-n_particles),l*n_particles+k)
                                        = amplitudes_old_mat((b-n_particles)*n_holes+(a-n_particles),l*n_particles+k)
                                        /(sp_energies(i) + sp_energies(j) - sp_energies(a) - sp_energies(b));
                            }
                        }
                    }
                }
            }
        }
    }

    exit(1);

    mat first_term = Vab_ij;
    mat second_term = Vab_cd*T_2;
    mat third_term = T_3*I_1;
    mat fourth_term_unaligned = Tai_kc*I_2;
    mat fourth_term = zeros<mat>(n_holes*n_holes, n_particles*n_particles);

    //******Realign*******

    for(int i = 0; i < n_particles; i++){
        for( int j = 0; j < n_particles; j++){
            for(int a = n_particles; a < number_of_states; a++){
                for(int b = n_particles; b < number_of_states; b++){
                    fourth_term_unaligned((b-n_particles)*n_holes+(a-n_particles), j*n_particles+i)
                            = fourth_term(i*n_holes+(a-n_particles),j*n_holes+(b-n_particles));
                }
            }
        }
    }

    mat amplitudes_new = first_term + 0.5*second_term + 0.5*third_term + fourth_term;


    return amplitudes_new;

}

void CCD_matrix::make_V_matrix(){

    Vab_ij = zeros<mat>(n_holes*n_holes, n_particles*n_particles);
    Vab_cd = zeros<mat>(n_holes*n_holes, n_holes*n_holes);
    Vkl_cd = zeros<mat>(n_particles*n_particles, n_holes*n_holes);
    Vkl_bj = zeros<mat>(n_particles*n_particles, n_holes*n_particles);
    Vkc_ld = zeros<mat>(n_particles*n_holes, n_particles*n_holes);
    Vkl_ij = zeros<mat>(n_particles*n_particles, n_particles*n_particles);
    Vkc_bj = zeros<mat>(n_particles*n_holes, n_particles*n_holes);

    for(int i = 0; i < n_particles; i++){
        for( int j = 0; j < n_particles; j++){
            for(int k = 0; k < n_particles; k++){
                for(int l = 0; l < n_particles; l++){
                    for(int a = n_particles; a < number_of_states; a++){
                        for(int b = n_particles; b < number_of_states; b++){
                            for(int c = n_particles; c < number_of_states; c++){
                                for (int d = n_particles; d < number_of_states; d++){
                                    Vab_ij((b-n_particles)*n_holes+(a-n_particles), j*n_particles+i) = TBME(index(a,b,i,j))
                                            /(sp_energies(i) + sp_energies(j) - sp_energies(a) - sp_energies(b));
                                    Vab_cd((b-n_particles)*n_holes+(a-n_particles),(d-n_particles)*n_holes+(c-n_particles))
                                            = TBME(index(a,b,c,d))
                                            /(sp_energies(i) + sp_energies(j) - sp_energies(a) - sp_energies(b));
                                    Vkl_cd((l*n_particles+k), (d-n_particles)*n_holes+(c-n_particles))
                                            = TBME(index(k,l,c,d));
                                    Vkl_bj((l*n_particles+k),j*n_holes+(b-n_particles))
                                            = TBME(index(k,l,b,j));
                                    Vkc_ld((c-n_particles)*n_particles+k, (d-n_particles)*n_particles+l)
                                            = TBME(index(k,c,l,d));
                                    Vkl_ij(l*n_particles+k, j*n_particles+i) = TBME(index(k,l,i,j));
                                    Vkc_bj((c-n_particles)*n_particles+k, j*n_holes+(b-n_particles))
                                            = TBME(index(k,c,b,j));
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

double CCD_matrix::CCD_energy(mat T){

    mat A = zeros<mat>(n_particles*n_particles, n_holes*n_holes);

    T = trans(T);

    for(int i = 0; i < n_particles; i++){
        for( int j = 0; j < n_particles; j++){
            for(int a = n_particles; a < number_of_states; a++){
                for(int b = n_particles; b < number_of_states; b++){
                    A(j*n_particles+i, (b-n_particles)*n_holes+(a-n_particles)) = TBME(index(i,j,a,b));
                }
            }
        }
    }

    double matrix_CCD_energy = 0.25*dot(A,T);

    return matrix_CCD_energy;
}

double CCD_matrix::CCD_energy_vec(vec amplitudes){

    double CCD_energy = 0;

    for(int i = 0; i < n_particles; i++){
        for(int j = 0; j < n_particles; j++){
            for(int a = n_particles; a < number_of_states; a++){
                for(int b = n_particles; b < number_of_states; b++){
                    CCD_energy += 0.25*TBME(index(i,j,a,b))
                            *amplitudes(index(a,b,i,j));
                }
            }
        }
    }

    return CCD_energy;
}

double CCD_matrix::CCD_solver(mat mapping){

    int CCD_counter = 0;
    int max_iterator = 1;
    double epsilon = 1.0e-6;
    double difference = 1.0;

    double CCD_energy_old = 0;
    double CCD_energy_new = 0;

    double CCD_energy_old_vec = 0;
    double CCD_energy_new_vec = 0;

    one_particle_energies_new_basis(sp_energies);

    mat amplitudes_old = intial_amplitudes_matrix();
    vec amplitudes_old_vec = intial_amplitudes_old();

    mat amplitudes_new_mat;
    vec amplitudes_new = zeros<vec>(pow(number_of_states, 4));

    cout << "initialized amplitudes" << endl;

    CCD_energy_new = CCD_energy(amplitudes_old);

    cout << "MBPT2 correction: " <<  setprecision(16) << CCD_energy_new << endl;

    while(max_iterator > CCD_counter && difference > epsilon){

        CCD_update(amplitudes_old_vec, amplitudes_new); //The test update

        cout << "far far " << endl;

        amplitudes_new_mat = CCD_update_matrix(amplitudes_old); //The matrix update

        CCD_energy_new_vec = CCD_energy_vec(amplitudes_new); //The test energy
        CCD_energy_new = CCD_energy(amplitudes_new_mat); //The matrix energy

        difference = std::abs(CCD_energy_old - CCD_energy_new);
        CCD_energy_old = CCD_energy_new;
        CCD_energy_old_vec = CCD_energy_new_vec;
        amplitudes_old = amplitudes_new;
        CCD_counter ++;

        cout << "iteration: " << CCD_counter << endl;
        cout << "CCD energy vec: " << CCD_energy_old_vec << endl;

        cout << "iteration: " << CCD_counter << endl;
        cout << "CCD energy mat: " << CCD_energy_old << endl;

    }

    if(max_iterator <= CCD_counter){
        cout << "Did not converge" << endl;
    }

    /*
    cout << "The counter reached " << CCD_counter << endl;

    double E_ref = CCD_energy_total();

    cout << "E_ref: " << E_ref << endl;

    cout << setprecision(12) << "CCD Energy " << E_ref + CCD_energy_old << endl;

    */

    return 0;
}

double CCD_matrix::CCD_energy_total(){

    int N = fermi_level;

    double Eref;
    double Eref_single = 0;
    double Eref_double = 0;

    for(int i = 0; i < N; i++) {
        Eref_single += sp_energies(i);
        for(int j = 0; j < N; j++) {
            Eref_double += 0.5*TBME(index(i,j,i,j));
            cout << 0.5*TBME(index(i,j,i,j)) << endl;
        }
    }

    cout << "Eref single: " << Eref_single << endl;
    cout << "Eref double: " << Eref_double << endl;

    Eref = Eref_single - Eref_double; //Se på dette!!

    return Eref;
}

inline
int CCD_matrix::index(int p,int q,int r,int s){

    int np = number_of_states;
    int np2 = np*np;
    int np3 = np2*np;
    int index = p*np3 + q*np2 + r*np + s;

    return index;
}

mat CCD_matrix::intermediate_1_mat(mat amplitudes){

    cout << "in the intermediate" << endl;

    mat Tcd_ij = zeros<mat>(n_holes*n_holes, n_particles*n_particles);

    //For the second intermediate

    mat Tdl_bj = zeros<mat>(n_particles*n_holes, n_holes*n_particles);
    mat Tdl_bi = zeros<mat>(n_particles*n_holes, n_holes*n_particles);


    for(int i = 0; i < n_particles; i++){
        for( int j = 0; j < n_particles; j++){
            for( int k = 0; k < n_particles; k++){
                for( int l = 0; l < n_particles; l++){
                    for(int a = n_particles; a < number_of_states; a++){
                        for(int b = n_particles; b < number_of_states; b++){
                            for (int c = n_particles; c < number_of_states; c++){
                                for (int d = n_particles; d < number_of_states; d++){
                                    Tcd_ij((d-n_particles)*n_holes+(c-n_particles),j*n_particles+i)
                                            = amplitudes((d-n_particles)*n_holes+(c-n_particles),j*n_particles+i);

                                    Tdl_bj(l*n_holes+(d-n_particles),j*n_holes+(b-n_particles))
                                            = amplitudes((b-n_particles)*n_holes+(d-n_particles),j*n_particles+l);

                                    Tdl_bi(l*n_holes+(d-n_particles),i*n_holes+(b-n_particles))
                                            = amplitudes((b-n_particles)*n_holes+(d-n_particles),i*n_particles+l);
                                }
                            }
                        }
                    }
                }
            }
        }
    }



    I_1 = Vkl_ij + 0.5*(Vkl_cd*Tcd_ij);

    I_2 = Vkc_bj + 0.5*(-1*Vkc_ld*Tdl_bj); //The minus is because of the transformation

    //I_3 = Vkc_bi + 0.5*(-1*Vkc_ld*Tdl_bi);

}

double CCD_matrix::intermediate_1(int i, int j, int k, int l, vec amplitudes){

    double value = TBME(index(k,l,i,j));

    for(int c = n_particles; c < number_of_states; c++){
        for(int d = n_particles; d < number_of_states; d++){
            value += 0.5*TBME(index(k,l,c,d))*amplitudes(index(c,d,i,j));
        }
    }

    return value;

}

double CCD_matrix::intermediate_2(int k,int b,int c,int j, vec amplitudes){

    double value = TBME(index(k,b,c,j));

    for(int l = 0; l < n_particles; l++){
        for(int d = n_particles; d < number_of_states; d++){
            value += 0.5*TBME(index(k,l,c,d))*amplitudes(index(d,b,l,j));
        }
    }

    return value;
}

double CCD_matrix::intermediate_3(int k, int j, vec amplitudes){

    double value = 0;

    for(int l = 0; l < n_particles; l++){
        for(int c = n_particles; c < number_of_states; c++){
            for(int d = n_particles; d < number_of_states; d++){

                value += 0.5*TBME(index(k,l,c,d))*amplitudes(index(c,d,j,l));
            }
        }
    }

    return value;
}

double CCD_matrix::intermediate_4(int b, int c, vec amplitudes){

    double value = 0;

    for(int k = 0; k < n_particles; k++){
        for(int l = 0; l < n_particles; l++){
            for(int d = n_particles; d < number_of_states; d++){
                value -= 0.5*TBME(index(k,l,c,d))*amplitudes(index(b,d,k,l));
            }
        }
    }

    return value;

}

