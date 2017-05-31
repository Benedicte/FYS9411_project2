#include "hartree_fock.h"

Hartree_fock_equations::Hartree_fock_equations(int sh_number, int number_of_electrons){

    int shell_number = sh_number;
    this->number_of_electrons = number_of_electrons;

    for(int i = shell_number; i > 0; i--){
        number_of_states = number_of_states + i*2;
    }
}

void Hartree_fock_equations::GetCoulombIntegrals(mat mapping, double hw){

    coloumb_matrix_direct = new double[number_of_states*number_of_states*number_of_states*number_of_states];
    coloumb_matrix_exchange = new double[number_of_states*number_of_states*number_of_states*number_of_states];

    int n_alpha, m_alpha, s_alpha;
    int n_beta, m_beta, s_beta;
    int n_gamma, m_gamma, s_gamma;
    int n_delta, m_delta, s_delta;

    int np = number_of_states;
    int np2 = np*np;
    int np3 = np2*np;

    for( int alpha = 0; alpha < number_of_states; alpha++){

        //cout << "alpha "  << alpha << endl;

        n_alpha = mapping(alpha,0);
        m_alpha = mapping(alpha,1);
        s_alpha = mapping(alpha,2);

        for( int beta = 0; beta < number_of_states; beta++){

            n_beta = mapping(beta,0);
            m_beta = mapping(beta,1);
            s_beta = mapping(beta,2);

            for (int gamma = 0; gamma < number_of_states; gamma++){

                n_gamma = mapping(gamma,0);
                m_gamma = mapping(gamma,1);
                s_gamma = mapping(gamma,2);

                for( int delta = 0; delta < number_of_states; delta++){

                    // cout << "inner "  << alpha << endl;

                    n_delta = mapping(delta,0);
                    m_delta = mapping(delta,1);
                    s_delta = mapping(delta,2);

                    if(s_gamma == s_delta && s_alpha == s_beta){

                        coloumb_matrix_direct[np3*(alpha) + np2*(beta) + np*(gamma) + delta] =
                            Coulomb_HO(hw, n_alpha, m_alpha, n_gamma, m_gamma, n_beta, m_beta, n_delta, m_delta);

                    } else{
                        coloumb_matrix_direct[np3*(alpha) + np2*(beta) + np*(gamma) + delta] = 0.0;
                    }

                    if(s_alpha == s_delta && s_gamma == s_beta){
                        coloumb_matrix_exchange[np3*(alpha) + np2*(beta) + np*(gamma) + delta] =
                                Coulomb_HO(hw, n_alpha, m_alpha, n_gamma, m_gamma, n_delta, m_delta, n_beta, m_beta);
                    } else{
                        coloumb_matrix_exchange[np3*(alpha) + np2*(beta) + np*(gamma) + delta] = 0.0;
                    }
                }
            }
        }
    }
}

mat Hartree_fock_equations::hartree_fock_method(mat mapping){

    cout << "number of particles " << number_of_states << endl;
    slater_determinant = eye<mat>(number_of_states, number_of_states);
    density_matrix = zeros<mat>(number_of_states, number_of_states);
    mat fock_matrix = zeros<mat>(number_of_states,number_of_states);

    int np = number_of_states;
    int np2 = np*np;
    int np3 = np2*np;

    double direct_elements = 0;
    double exchange_elements = 0;
    double hw = 1;

    int maxIterator = 200;
    double epsilon = 1.0e-9;
    double difference = 1.0;
    int hf_count = 1;
    double sumFockTerm = 0.0;

    int s_alpha, s_beta, s_gamma, s_delta;

    vec oldenergies(number_of_states);
    vec newenergies(number_of_states);

    for( int gamma = 0; gamma < number_of_states; gamma++){
        for( int delta = 0; delta < number_of_states; delta++){
            double sum = 0.0;
            for(int i = 0; i < number_of_electrons; i++){
                sum += slater_determinant(gamma,i)*slater_determinant(delta,i);
            }
            density_matrix(gamma, delta) = sum;
        }
    }

    while( hf_count < maxIterator && difference >= epsilon){
        cout << "############### Iteration " << hf_count <<" ###############" << endl;
        fock_matrix.zeros();

        for( int alpha = 0; alpha < number_of_states; alpha++){

            double energy_alpha = Hartree_fock_equations::get_energy(mapping(alpha,0), mapping(alpha,1), hw);

            s_alpha = mapping(alpha,2);

            fock_matrix(alpha,alpha) = energy_alpha;

            for( int beta = 0; beta < number_of_states; beta++){

                sumFockTerm = 0;
                s_beta = mapping(beta,2);

                for (int gamma = 0; gamma < number_of_states; gamma++){

                    s_gamma = mapping(gamma,2);

                    for( int delta = 0; delta < number_of_states; delta++){

                        s_delta = mapping(delta,2);
                        direct_elements = 0;
                        exchange_elements = 0;

                        if(s_gamma == s_delta && s_alpha == s_beta){
                            direct_elements = density_matrix(gamma,delta)
                                    *coloumb_matrix_direct[np3*(alpha) + np2*(beta) + np*(gamma) + delta];

                            sumFockTerm += direct_elements;
                        }

                        if(s_alpha == s_delta && s_gamma == s_beta){
                            exchange_elements = density_matrix(gamma,delta)
                                    *coloumb_matrix_exchange[np3*(alpha) + np2*(beta) + np*(gamma) + delta];
                            sumFockTerm -= exchange_elements;
                        }
                    }
                }
                fock_matrix(alpha, beta) += sumFockTerm;
                //fock_matrix(beta, alpha) += sumFockTerm;
            }
        }

        vec energies;

        eig_sym(energies, slater_determinant, fock_matrix);

        /// Setting up new density matrix

        for( int gamma = 0; gamma < number_of_states; gamma++){
            for( int delta = 0; delta < number_of_states; delta++){
                double sum = 0.0;

                for(int i = 0; i < number_of_electrons; i++){
                    sum += slater_determinant(gamma,i)*slater_determinant(delta,i);
                }
                density_matrix(gamma, delta) = sum;
            }
        }

        newenergies = energies;

        //Brute force computation of difference between previous and new sp HF energies """
        double sum = 0.0;

        for( int i = 0; i < number_of_states; i++){
            sum += (std::abs(newenergies(i)-oldenergies(i)))/number_of_states;
            difference = sum;
            oldenergies = newenergies;
        }

        hf_count ++;
        cout << "The difference " << difference << endl;
    }

    std::setprecision(7);
    cout.setf(ios::fixed);

    newenergies.raw_print(std::cout);
    energies = newenergies;


    return newenergies;
}

double Hartree_fock_equations::get_energy(int n, int m, double omega){

    double energy;

    energy = omega*(2.0*n + abs(m) + 1);

    return energy;
}

vec Hartree_fock_equations::new_basis(mat mapping, mat fock_matrix){


    int np = number_of_states;
    int np2 = np*np;
    int np3 = np2*np;

    vec new_basis = zeros<vec>(np3*number_of_states);
    double direct_elements;
    double exchange_elements;
    double anti_sym;

    //fock_matrix = fock_matrix.t();

    //for(int i = 0; i < number_of_states; i++){
    //    for(int j = 0; j < number_of_states; j++){

     //       if( fabs(fock_matrix(i,j)) < 1E-8)
                cout << 0 << " ";
     //       else
     //           cout << fock_matrix(i,j) << " ";

     //   }
       // cout << endl;
    //}

    //fock_matrix.print(std::cout);
    for(int p = 0; p < number_of_states; p++){
        for(int q = 0; q < number_of_states; q++){
            for(int r = 0; r < number_of_states; r++){
                for(int s = 0; s < number_of_states; s++){

                    double element = 0;

                    for( int alpha = 0; alpha < number_of_states; alpha++){
                        int s_alpha = mapping(alpha,2);

                        for( int beta = 0; beta < number_of_states; beta++){
                            int s_beta = mapping(beta,2);

                            for (int gamma = 0; gamma < number_of_states; gamma++){

                                int s_gamma = mapping(gamma,2);

                                for( int delta = 0; delta < number_of_states; delta++){

                                    int s_delta = mapping(delta,2);

                                    direct_elements = 0.0;
                                    exchange_elements = 0.0;

                                    //if(s_alpha == s_beta && s_gamma == s_delta){
                                    //    direct_elements =
                                    //            coloumb_matrix_direct[np3*(alpha) + np2*(beta) + np*(gamma) + delta];

                                    //}

                                    //if(s_alpha == s_delta && s_gamma == s_beta){

                                    //    exchange_elements =
                                    //            coloumb_matrix_exchange[np3*(alpha) + np2*(beta) + np*(gamma) + delta];

                                    //}

                                    //if(s_gamma == s_delta && s_alpha == s_beta){
                                        direct_elements =
                                                coloumb_matrix_direct[np3*(alpha) + np2*(beta) + np*(gamma) + delta];

                                    //}

                                    //if(s_alpha == s_delta && s_gamma == s_beta){

                                        exchange_elements =
                                                coloumb_matrix_exchange[np3*(alpha) + np2*(beta) + np*(gamma) + delta];

                                    //}

                                    anti_sym = direct_elements - exchange_elements;
                                    //fock_matrix(p, alpha)*fock_matrix(q,gamma)*fock_matrix(r,delta)*fock_matrix(s,beta)*anti_sym;
                                    if(anti_sym != 0){

                                        element +=
                                                fock_matrix(p, alpha)*fock_matrix(q,delta)*fock_matrix(r,gamma)*fock_matrix(s,beta)*anti_sym;
                                    }
                                }
                            }
                        }
                    }

                    new_basis(np3*p + np2*q + np*r + s) = element;

                }
            }
        }
    }

    return new_basis;

}

double Hartree_fock_equations::total_energy(mat mapping, vec energies){

    int fermi_level = number_of_electrons;
    double ground_state_energy;
    double single_particle_contributions = 0.0;
    double double_particle_contributions = 0.0;
    double direct_elements;
    double exchange_elements;
    double anti_sym;
    int TBME_length = number_of_states*number_of_states*number_of_states*number_of_states;
    TBME = zeros<vec>(TBME_length);

    int np = number_of_states;
    int np2 = np*np;
    int np3 = np2*np;


    for(int i = 0; i < fermi_level; i++){
        single_particle_contributions += energies(i);
    }

    for( int alpha = 0; alpha < number_of_states; alpha++){
        int s_alpha = mapping(alpha,2);

        for( int beta = 0; beta < number_of_states; beta++){
            int s_beta = mapping(beta,2);

            for (int gamma = 0; gamma < number_of_states; gamma++){

                int s_gamma = mapping(gamma,2);

                for( int delta = 0; delta < number_of_states; delta++){

                    int s_delta = mapping(delta,2);

                    direct_elements = 0.0;
                    exchange_elements = 0.0;

                    if(s_alpha == s_beta && s_gamma == s_delta){
                        direct_elements =
                                coloumb_matrix_direct[np3*(alpha) + np2*(beta) + np*(gamma) + delta];

                    }

                    if(s_alpha == s_delta && s_gamma == s_beta){

                        exchange_elements =
                                coloumb_matrix_exchange[np3*(alpha) + np2*(beta) + np*(gamma) + delta];

                    }

                    anti_sym = direct_elements - exchange_elements;

                    TBME(np3*(alpha) + np2*(beta) + np*(gamma) + delta) = anti_sym;

                   // if(anti_sym_new_basis != 0){
                   //     double_particle_contributions_HF += density_matrix(alpha,beta)*density_matrix(gamma,delta)* anti_sym_new_basis;
                   // }

                    if(anti_sym != 0){
                        double_particle_contributions += density_matrix(alpha,beta)*density_matrix(gamma,delta)*anti_sym;
                    }
                }
            }
        }
    }


    //double_particle_contributions = 19.527
    //Ground state energy 20.76692

    ground_state_energy = single_particle_contributions - 0.5 * double_particle_contributions;

    //ground_state_energy_HF = single_particle_contributions - 0.5 * double_particle_contributions_HF;

    cout << "number of particles: " << number_of_electrons << endl;

    cout << "double particle cont" << endl;
    cout << double_particle_contributions << endl;


    //cout << "double particle cont HF basis" << endl;
    //cout << double_particle_contributions_HF << endl;

    cout << "single particle contribution " << endl;


    return ground_state_energy;
}

double Hartree_fock_equations::total_energy_test(mat mapping, vec energies, vec new_basis){

    int fermi_level = number_of_electrons;
    double ground_state_energy;
    double single_particle_contributions = 0.0;
    double double_particle_contributions = 0.0;
    double direct_elements;
    double exchange_elements;
    double anti_sym;

    int np = number_of_states;
    int np2 = np*np;
    int np3 = np2*np;

    for(int i = 0; i < fermi_level; i++){
        single_particle_contributions += energies(i);
    }

    for( int alpha = 0; alpha < number_of_states; alpha++){
        int s_alpha = mapping(alpha,2);

        for( int beta = 0; beta < number_of_states; beta++){
            int s_beta = mapping(beta,2);

            for (int gamma = 0; gamma < number_of_states; gamma++){

                int s_gamma = mapping(gamma,2);

                for( int delta = 0; delta < number_of_states; delta++){

                    int s_delta = mapping(delta,2);

                    direct_elements = 0.0;
                    exchange_elements = 0.0;

                    if(s_alpha == s_beta && s_gamma == s_delta){
                        direct_elements =
                                new_basis[np3*(alpha) + np2*(beta) + np*(gamma) + delta];

                    }

                    if(s_alpha == s_delta && s_gamma == s_beta){

                        exchange_elements =
                                new_basis[np3*(alpha) + np2*(beta) + np*(gamma) + delta];

                    }

                    anti_sym = direct_elements;

                    if(anti_sym != 0){
                        double_particle_contributions += density_matrix(alpha,beta)*density_matrix(gamma,delta)*anti_sym;
                    }

                }
            }
        }
    }


    //double_particle_contributions = 19.527
    //Ground state energy 20.76692

    ground_state_energy = single_particle_contributions - 0.5 * double_particle_contributions;

    cout << "number of particles: " << number_of_electrons << endl;
    cout << "Ground state energy HF basis 2" << endl;
    cout << double_particle_contributions << endl;

    return ground_state_energy;
}

vec Hartree_fock_equations::GetTBME(){

    int TBME_length = number_of_states*number_of_states*number_of_states*number_of_states;
    vec TBME = zeros<vec>(TBME_length);

    for(int i = 0; i < TBME_length; i++){
        TBME(i) = coloumb_matrix_direct[i] - coloumb_matrix_exchange[i];
    }

    return TBME;

}
