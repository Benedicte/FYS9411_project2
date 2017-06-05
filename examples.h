#ifndef EXAMPLES_H
#define EXAMPLES_H
#include <armadillo>

using namespace arma;

class Examples
{
public:
    Examples();
    static vec TwoParticleDotTest(int shells, double hw);
    static double KroneckerDelta(int i, int j);

};

#endif // EXAMPLES_H
