#include "ccd_copy.h"
#include "coulomb_function2.h"
#include <iomanip>

using namespace std;

CCD_1::CCD_1(int shell_number, int fermi_lv, vec two_body_matrix_elements, vec single_particle_energies)
{
    fermi_level = fermi_lv;

    for(int i = shell_number; i > 0; i--){
        number_of_states = number_of_states + i*2;
    }

    TBME = two_body_matrix_elements;
    //TBME.print(std::cout);
    sp_energies = single_particle_energies;
    //sp_energies = zeros<vec>(number_of_states*number_of_states*number_of_states*number_of_states);

}

void CCD_1::one_particle_energies(vec sp_energies){

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

vec CCD_1::one_particle_energies_new_basis(vec sp_energies, mat mapping){

    double ei,ea;

    coulomb_function2 test(1, number_of_states, fermi_level);

    /*
    cout << "Beginning of the TBME we need" << endl;

    for(int i=0; i<fermi_level; i++){
        for(int j = 0; j < fermi_level; j++){
            cout << (test.calc_TBME(i, j, i, j, mapping)-test.calc_TBME(i, j, j, i, mapping)) << " ";
        }
        cout << endl;
    }

    cout << "End of the TBME we need" << endl;
    */

    for(int i= 0 ; i < fermi_level; i++){
        ei = sp_energies(i);
        for(int j = 0 ; j < fermi_level; j++){
            //ei += calcVpqrs(i,j,i,j,SPbasis); double
            ei += TBME(index(i,j,i,j));
        }
        //sp_energies_newbasis(i) = ei;
        sp_energies(i) = ei;
    }

    for(int a = fermi_level; a < number_of_states; a++){
        //ea = sp_energies_newbasis(a);
        ea = sp_energies(a);
        for(int j=0; j<fermi_level; j++){
            //ea += calcVpqrs(a,j,a,j,SPbasis);
            ea += TBME(index(a,j,a,j));
        }
        //sp_energies_newbasis(a) = ea;
        sp_energies(a) = ea;
    }

    for(int i = 0; i < number_of_states; i++){
        cout << mapping(i,0) << " " << mapping(i,1) << " "<<  mapping(i,2) << " " << sp_energies(i)<<endl;
    }

    sp_energies_HF = sp_energies;

    return sp_energies;
}

vec CCD_1::intial_amplitudes_old(mat mapping){

    cout << number_of_states << endl;
    cout << fermi_level << endl;

    vec initial_amplitudes= zeros<vec>(number_of_states*number_of_states*number_of_states*number_of_states);

    n_holes = number_of_states - fermi_level;
    n_particles = fermi_level;

    cout << size(TBME)<< endl;
    cout << size(initial_amplitudes)<< endl;

    cout << "Size of SP-Energies" << size(sp_energies) << endl;

    cout << "HEI" << endl;

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

void CCD_1::CCD_update(mat mapping, vec amplitudes_old, vec &amplitudes_new){

    double sum = 0.0;
    double energy_denom = 0.0;

    //vec amplitudes_new = zeros<vec>(number_of_states*number_of_states*number_of_states*number_of_states);

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
                            sum += 0.5*TBME(index(k,l,i,j))
                                    *amplitudes_old(index(a,b,k,l)); //Tenth term
                        }
                    }

                    for(int k = 0; k < n_particles; k++){
                        for(int c = n_particles; c < number_of_states; c++){
                            sum += TBME(index(k,b,c,j))
                                    *amplitudes_old(index(a,c,i,k)); //Eleventh term
                            sum -= TBME(index(k,a,c,j))
                                    *amplitudes_old(index(b,c,i,k)); //Twelvth term
                            sum -= TBME(index(k,b,c,i))
                                    *amplitudes_old(index(a,c,j,k)); //Thirteenth term
                            sum += TBME(index(k,a,c,i))
                                    *amplitudes_old(index(b,c,j,k)); //Fourteenth term
                        }
                    }

                    for(int k = 0; k < n_particles; k++){
                        for(int l = 0; l < n_particles; l++){
                            for(int c = n_particles; c < number_of_states; c++){
                                for(int d = n_particles; d < number_of_states; d++){
                                    sum += 0.25*TBME(index(k,l,c,d))*amplitudes_old(index(c,d,i,j))
                                            *amplitudes_old(index(a,b,k,l)); //Third term
                                }
                            }
                        }
                    }
                    for(int k = 0; k < n_particles; k++){
                        for(int l = 0; l < n_particles; l++){
                            for(int c = n_particles; c < number_of_states; c++){
                                for(int d = n_particles; d < number_of_states; d++){
                                    sum += TBME(index(k,l,c,d))*amplitudes_old(index(a,c,i,k))
                                            *amplitudes_old(index(b,d,j,l)); //Fourth term
                                    sum -= TBME(index(k,l,c,d))*amplitudes_old(index(a,c,j,k))
                                            *amplitudes_old(index(b,d,i,l)); //Fifth term
                                }
                            }
                        }
                    }
                    for(int k = 0; k < n_particles; k++){
                        for(int l = 0; l < n_particles; l++){
                            for(int c = n_particles; c < number_of_states; c++){
                                for(int d = n_particles; d < number_of_states; d++){

                                    sum -= 0.5*TBME(index(k,l,c,d))*amplitudes_old(index(d,c,i,k))
                                            *amplitudes_old(index(a,b,l,j)); //Sixth term
                                    sum += 0.5*TBME(index(k,l,c,d))*amplitudes_old(index(d,c,j,k))
                                            *amplitudes_old(index(a,b,l,i)); //Seventh term
                                }
                            }
                        }
                    }

                    for(int k = 0; k < n_particles; k++){
                        for(int l = 0; l < n_particles; l++){
                            for(int c = n_particles; c < number_of_states; c++){
                                for(int d = n_particles; d < number_of_states; d++){

                                    sum -= 0.5*TBME(index(k,l,c,d))*amplitudes_old(index(a,c,l,k))
                                            *amplitudes_old(index(d,b,i,j)); // Eighth term
                                    sum += 0.5*TBME(index(k,l,c,d))*amplitudes_old(index(b,c,l,k))
                                            *amplitudes_old(index(d,a,i,j)); //Ninth term

                                }
                            }
                        }
                    }

                    energy_denom = sp_energies(i) + sp_energies(j) - sp_energies(a) - sp_energies(b);
                    amplitudes_new(index(a,b,i,j)) = sum/energy_denom;


                }//End of main indexes
            }
        }
    }

    return;

} //End of function

double CCD_1::CCD_energy(mat mapping, vec amplitudes){

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

double CCD_1::CCD_solver(mat mapping){

    int CCD_counter = 0;
    int max_iterator = 400;
    double epsilon = 1.0e-6;
    double difference = 1.0;

    double CCD_energy_old = 0;
    double CCD_energy_new = 0;


    vec amplitudes_old = intial_amplitudes_old(mapping);
    vec amplitudes_new(pow(number_of_states,4));

    cout << "initialized amplitudes" << endl;

    CCD_energy_new = CCD_energy(mapping, amplitudes_old);

    cout << "MBPT2 correction: " <<  setprecision(16) << CCD_energy_new << endl;

    while(max_iterator > CCD_counter && difference > epsilon){

        CCD_update(mapping, amplitudes_old,amplitudes_new);

        CCD_energy_new = CCD_energy(mapping, amplitudes_new);
        difference = std::abs(CCD_energy_old - CCD_energy_new);
        CCD_energy_old = CCD_energy_new;
        amplitudes_old = amplitudes_new;
        CCD_counter ++;

        cout << "iteration: " << CCD_counter << endl;
        cout << "CCD energy: " << CCD_energy_old << endl;
    }

    if(max_iterator <= CCD_counter){
        cout << "Did not converge" << endl;
    }
    cout << "The counter reached " << CCD_counter << endl;

    double E_ref = CCD_energy_total(mapping);

    cout << "E_ref: " << E_ref << endl;

    cout << setprecision(12) << "CCD Energy " << E_ref + CCD_energy_old << endl;

    return 0;
}

double CCD_1::CCD_energy_total(mat mapping){

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

    Eref = Eref_single + Eref_double;

    return Eref;
}

int CCD_1::index(int p,int q,int r,int s){

    int np = number_of_states;
    int np2 = np*np;
    int np3 = np2*np;
    int index = p*np3 + q*np2 + r*np + s;

    return index;
}

/*

double CCD::computeTau(int a, int b, int i, int j) {


    double val = m_t2_old[a-n_holes][b-n_holes][i][j]
            + m_t1_old[a-n_holes][i]*m_t1_old[b-n_holes][j]
            - m_t1_old[b-n_holes][i]*m_t1_old[a-n_holes][j];
    return val;

}

double CCD::computeW1(int m, int n, int i, int j) {

    int N = m_holeStates;
    int L = m_basisFunctions;

    double val = m_twoBodyElements[m][n][i][j];

    for(int e = N; e < L; e++){
        val += m_twoBodyElements[m][n][i][e]*m_t1_old[e-N][j] - m_twoBodyElements[m][n][j][e]*m_t1_old[e-N][i];
    }

    for(int e = N; e < L; e++){
        for(int f = N; f < L; f++){
            val += 0.25*m_twoBodyElements[m][n][e][f]*m_tau[e-N][f-N][i][j];
        }
    }

    return val;
}

double CCD::computeW2(int a, int b, int e, int f) {

    int N = m_holeStates;
    //int L = m_basisFunctions;

    double val = m_twoBodyElements[a][b][e][f];

    for(int m=0; m < N; m++) {
        val -= m_t1_old[b-N][m]*m_twoBodyElements[a][m][e][f] - m_t1_old[a-N][m]*m_twoBodyElements[b][m][e][f];
        for(int n=0; n < N; n++) {
            val += 0.25*m_tau[a-N][b-N][m][n]*m_twoBodyElements[m][n][e][f];
        }
    }

    return val;
}

double CCD::computeW3(int m, int b, int e, int j) {

    int N = m_holeStates;
    int L = m_basisFunctions;

    double val = m_twoBodyElements[m][b][e][j];

    for(int f = N; f < L; f++) {
        val += m_twoBodyElements[m][b][e][f]*m_t1_old[f-N][j];
    }

    for(int n = 0; n < N; n++) {
        val -= m_t1_old[b-N][n]*m_twoBodyElements[m][n][e][j];
    }

    for(int n = 0; n < N; n++) {
        for(int f = N; f < L; f++) {
            val -= (0.5*m_t2_old[f-N][b-N][j][n]+m_t1_old[f-N][j]*m_t1_old[b-N][n])*m_twoBodyElements[m][n][e][f];
        }
    }

    return val;
}

double CCD::computeF1(int m, int e) {

    int N = m_holeStates;
    int L = m_basisFunctions;

    double val = m_F[m][e];

    for(int n = 0; n < N; n++) {
        for(int f = N; f < L; f++){
            val += m_twoBodyElements[m][n][e][f]*m_t1_old[f-N][n];
        }
    }

    return val;
}

double CCD::computeF2(int m, int i) {

    int N = m_holeStates;
    int L = m_basisFunctions;

    double val = (1.0-delta(m,i))*m_F[m][i];

    for(int e = N; e < L; e++) {
        val += 0.5*m_F[m][e]*m_t1_old[e-N][i];
    }

    for(int e = N; e < L; e++) {
        for(int n = 0; n < N; n++) {
            val += m_twoBodyElements[m][n][i][e]*m_t1_old[e-N][n];
        }
    }

    for(int n = 0; n < N; n++) {
        for(int e = N; e < L; e++) {
            for(int f = N; f < L; f++) {
                val += 0.5*m_twoBodyElements[m][n][e][f]*m_tau2[e-N][f-N][i][n];
            }
        }
    }

    return val;
}

double CCD::computeF3(int a, int e) {

    int N = m_holeStates;
    int L = m_basisFunctions;

    double val = (1.0-delta(a,e))*m_F[a][e];

    for(int m = 0; m < N; m++) {
        val -= 0.5*m_F[m][e]*m_t1_old[a-N][m];
    }

    for(int m = 0; m < N; m++) {
        for(int n = 0; n < N; n++) {
            for(int f = N; f < L; f++) {
                val -= 0.5*m_twoBodyElements[m][n][e][f]*m_tau2[a-N][f-N][m][n];
            }
        }
    }

    for(int m = 0; m < N; m++) {
        for(int f = N; f < L; f++) {
            val += m_twoBodyElements[m][a][f][e]*m_t1_old[f-N][m];
        }
    }

    return val;
}
*/
