#include "hartree_fock.h"
#include <iomanip>

Hartree_fock_equations::Hartree_fock_equations(int sh_number, int number_of_electrons){

    shell_number = sh_number;
    this->number_of_electrons = number_of_electrons;

    for(int i = shell_number; i > 0; i--){
        number_of_states = number_of_states + i*2;
    }
}

vec Hartree_fock_equations::get_energies_harmonic(mat mapping, double hw){

    energies_pre = zeros<vec>(number_of_states);

    for(int i = 0; i<number_of_states; i++){
        energies_pre(i) = get_energy(mapping(i,0), mapping(i,1), hw);
    }

    return energies_pre;

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

                      //  coloumb_matrix_direct[np3*(alpha) + np2*(beta) + np*(gamma) + delta] =
                      //          Coulomb_HO(hw, n_alpha, m_alpha, n_gamma, m_gamma, n_beta, m_beta, n_delta, m_delta);

                    } else{
                       // coloumb_matrix_direct[np3*(alpha) + np2*(beta) + np*(gamma) + delta] = 0.0;
                    }

                    if(s_alpha == s_delta && s_gamma == s_beta){
                       // coloumb_matrix_exchange[np3*(alpha) + np2*(beta) + np*(gamma) + delta] =
                       //         Coulomb_HO(hw, n_alpha, m_alpha, n_gamma, m_gamma, n_delta, m_delta, n_beta, m_beta);
                    } else{
                       // coloumb_matrix_exchange[np3*(alpha) + np2*(beta) + np*(gamma) + delta] = 0.0;
                    }
                }
            }
        }
    }
}

mat Hartree_fock_equations::hartree_fock_method(mat mapping, double hw){

    cout << "number of particles " << number_of_states << endl;
    slater_determinant = eye<mat>(number_of_states, number_of_states);
    density_matrix = zeros<mat>(number_of_states, number_of_states);
    mat fock_matrix = zeros<mat>(number_of_states,number_of_states);

    TBME_new(shell_number, hw, mapping);

    int np = number_of_states;
    int np2 = np*np;
    int np3 = np2*np;

    int maxIterator = 200;
    double epsilon = 1.0e-7;
    double difference = 1.0;
    int hf_count = 1;
    double sumFockTerm = 0.0;

    vec oldenergies(number_of_states);
    vec newenergies(number_of_states);

    while( hf_count < maxIterator && difference >= epsilon){
        cout << "############### Iteration " << hf_count <<" ###############" << endl;
        fock_matrix.zeros();

        for( int alpha = 0; alpha < number_of_states; alpha++){

            double energy_alpha = Hartree_fock_equations::get_energy(mapping(alpha,0), mapping(alpha,1), hw);


            fock_matrix(alpha,alpha) = energy_alpha;

            for( int beta = 0; beta < number_of_states; beta++){

                sumFockTerm = 0;

                for (int gamma = 0; gamma < number_of_states; gamma++){
                    for( int delta = 0; delta < number_of_states; delta++){

                        sumFockTerm += twoBme(np3*alpha + np2*gamma + np*beta + delta)
                                *density_matrix(gamma, delta);
                    }
                }
                fock_matrix(alpha, beta) += sumFockTerm;
            }
        }

        vec energies;

        eig_sym(energies, slater_determinant, fock_matrix);

        // Setting up new density matrix

        mat slater_determinant_t = slater_determinant.t();


        for( int gamma = 0; gamma < number_of_states; gamma++){
            for( int delta = 0; delta < number_of_states; delta++){
                double sum = 0.0;
                for(int i = 0; i < number_of_electrons; i++){
                    sum += slater_determinant_t(i,gamma)*slater_determinant_t(i,delta);
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

vec Hartree_fock_equations::new_basis(){

    //Definity need to use density matrices here

    int np = number_of_states;
    int np2 = np*np;
    int np3 = np2*np;

    vec TBME = zeros<vec>(number_of_states*number_of_states*number_of_states*number_of_states);

    cout << "entered new basis" << endl;

    for(int p = 0; p < number_of_states; p++){
        for(int q = 0; q < number_of_states; q++){
            for(int r = 0; r < number_of_states; r++){
                for(int s = 0; s < number_of_states; s++){
                    double element = 0.0;
                    for( int alpha = 0; alpha < number_of_states; alpha++){
                        for( int beta = 0; beta < number_of_states; beta++){
                            for (int gamma = 0; gamma < number_of_states; gamma++){
                                for( int delta = 0; delta < number_of_states; delta++){

                                    if(twoBme(np3*alpha + np2*beta + np*gamma + delta)!= 0){
                                        element += slater_determinant(alpha,p)*slater_determinant(beta,q)
                                                *slater_determinant(gamma,r)*slater_determinant(delta,s)
                                                *twoBme(np3*alpha + np2*beta + np*gamma + delta);
                                        //element += density_matrix(alpha,gamma)*density_matrix(beta,delta)
                                        //     *twoBme(np3*alpha + np2*beta + np*gamma + delta);
                                    }
                                }
                            }
                        }
                    }
                    TBME(np3*p + np2*q + np*r + s) = element;
                }
            }
        }
    }

    cout << "new basis created" << endl;

    return TBME;
}

double Hartree_fock_equations::total_energy(){

    int fermi_level = number_of_electrons;
    double ground_state_energy;
    double single_particle_contributions = 0.0;
    double single_particle_contributions_HF = 0.0;
    double double_particle_contributions = 0.0;
    double double_particle_contributions_HF = 0.0;
    double ground_state_energy_HF = 0.0;

    int np = number_of_states;
    int np2 = np*np;
    int np3 = np2*np;

    vec TBME_newbasis = new_basis();

    for(int i = 0; i < fermi_level; i++){
        single_particle_contributions += energies(i);
    }

    cout << single_particle_contributions << endl;

    single_particle_contributions_HF = one_particle_energies_new_basis(TBME_newbasis);


    for( int alpha = 0; alpha < number_of_states; alpha++){
        for( int beta = 0; beta < number_of_states; beta++){
            for (int gamma = 0; gamma < number_of_states; gamma++){
                for( int delta = 0; delta < number_of_states; delta++){

                    double_particle_contributions += density_matrix(alpha,gamma)*density_matrix(beta,delta)
                            *twoBme(np3*alpha + np2*beta + np*gamma + delta);

                    double_particle_contributions_HF += density_matrix(alpha,gamma)*density_matrix(beta,delta)
                            *TBME_newbasis(np3*alpha + np2*beta + np*gamma + delta);
                }
            }
        }
    }

    ground_state_energy = single_particle_contributions - 0.5 * double_particle_contributions;

    ground_state_energy_HF = single_particle_contributions_HF - 0.5 * double_particle_contributions_HF;

    cout << "number of particles: " << number_of_electrons << endl;

    cout << std::setprecision(24);


    cout << "single particle cont" << endl;
    cout << single_particle_contributions << endl;

    cout << "single particle cont HF basis" << endl;
    cout << single_particle_contributions_HF << endl;

    cout << "double particle cont" << endl;
    cout << double_particle_contributions << endl;

    cout << "double particle cont HF basis" << endl;
    cout << double_particle_contributions_HF << endl;

    cout << "Diff " << double_particle_contributions - double_particle_contributions_HF <<endl;

    cout << "Ground state energy" << endl;
    cout << ground_state_energy << endl;

    cout << "Ground state energy HF basis" << endl;
    cout << ground_state_energy_HF << endl;

    return ground_state_energy;
}

vec Hartree_fock_equations::TBME_new(int shells, double hw, mat mapping) {

    int basisFunctions = shells*(shells+1);

    int np3 = pow(basisFunctions,3);
    int np2 = pow(basisFunctions,2);
    int np1 = basisFunctions;

    twoBme = zeros<vec>(pow(basisFunctions,4));
    direct_elements = zeros<vec>(pow(basisFunctions,4));
    exchange_elements = zeros<vec>(pow(basisFunctions,4));

    for(int p=0; p<basisFunctions; p++) {

        int np = mapping(p,0);
        int mp = mapping(p,1);
        int sp = mapping(p,2);

        //cout << p << " " << np << " " << mp << " " << sp << endl;

        for(int q=0; q<basisFunctions; q++) {
            int nq = mapping(q,0);
            int mq = mapping(q,1);
            int sq = mapping(q,2);
            for(int r=0; r<basisFunctions; r++) {
                int nr = mapping(r,0);
                int mr = mapping(r,1);
                int sr = mapping(r,2);
                for(int s=0; s<basisFunctions; s++) {
                    int ns = mapping(s,0);
                    int ms = mapping(s,1);
                    int ss = mapping(s,2);

                    //if ( (mp + mq == mr +ms) && (sp+ sq == sr + ss) ){

                    double direct = 0;

                    //direct   = Coulomb_HO(hw, np, mp, nq, mq, nr, mr, ns, ms);

                    double exchange = 0;

                    //exchange = Coulomb_HO(hw, np, mp, nq, mq, ns, ms, nr, mr);

                    double TBMEAS = KroneckerDelta(sp, sr)*KroneckerDelta(sq, ss)*direct
                            -KroneckerDelta(sp, ss)*KroneckerDelta(sq, sr)*exchange;
                    int index = np3*p + np2*q + np1*r + s;
                    twoBme(index) = TBMEAS;

                    //cout << TBMEAS << " " << p << " " << q << " " << r << " " << s << endl;
                    //}

                }}}}
    //exit(1);
    return twoBme;

}

double Hartree_fock_equations::KroneckerDelta(int i, int j){
    if(i == j){
        return 1.0;
    }

    if(i != j){
        return 0;
    }
}

double Hartree_fock_equations::one_particle_energies_new_basis(vec TBME_new){

    double ei,ea = 0;
    double sp_energy = 0;

    int np3 = pow(number_of_states,3);
    int np2 = pow(number_of_states,2);
    int np = number_of_states;

    for(int i= 0 ; i < number_of_electrons; i++){
        ei = energies(i);
        for(int j = 0 ; j < number_of_electrons; j++){
            ei += TBME_new(np3*i + np3*j + np*i + j);
        }
        energies(i) = ei;
    }

    for(int a = number_of_electrons; a < number_of_states; a++){
        ea = energies(a);
        for(int j=0; j<number_of_electrons; j++){
            ea += TBME_new(np3*a + np2*j + np*a + j);
        }

        energies(a) = ea;
    }

    for(int i = 0; i < number_of_electrons; i++){
        sp_energy += energies(i);
    }

    return sp_energy;
}
