#include "examples.h"
#include <unordered_map>
#include <basis.h>
#include <coulomb_functions.h>
#include <armadillo>

using namespace std;
using namespace arma;


Examples::Examples()
{

}

vec Examples::TwoParticleDotTest(int shells) {



    int N = 2;
    //int shells = 4;
    int basisFunctions = shells*(shells+1);

    int np3 = pow(basisFunctions,3);
    int np2 = pow(basisFunctions,2);
    int np1 = basisFunctions;

    vec twoBme = zeros<vec>(pow(basisFunctions,4));

    Basis* basis = new Basis(shells);
    double hw = 1.0;
    mat mapping = basis->map_quantum_numbers(basisFunctions);



    for(int p=0; p<basisFunctions; p++) {

        int np = mapping(p,0);
        int mp = mapping(p,1);
        int sp = mapping(p,2);

        cout << p << " " << np << " " << mp << " " << sp << endl;

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

                        direct   = Coulomb_HO(hw, np, mp, nq, mq, nr, mr, ns, ms);

                        double exchange = 0;

                        exchange = Coulomb_HO(hw, np, mp, nq, mq, ns, ms, nr, mr);

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

double Examples::KroneckerDelta(int i, int j){
    if(i == j){
        return 1.0;
    }

    if(i != j){
        return 0;
    }
}

