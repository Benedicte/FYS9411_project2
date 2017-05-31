#ifndef ABSTRACT_COULOMB.H
#define ABSTRACT_COULOMB.H
#include "abstract_coulomb.h"

using namespace arma;

// This probably shouldn't be in
// the abstract. Depends on basis.
struct channelBundle{
  int chanNx;
  int chanNy;
  int chanNz;
  int chanSz;
  int chanTz;
  int chanMl;
};

class abstract_coulomb{
public:
  // double g;
  // double density;
  // double r_s;
  // int tzMax;
  // int nMax;
  // int shellMax;
  int Nspstates;
  int Nparticles;
  int Nchannels;
  int ** indexMap;
  double * spEnergy;
  channelBundle * chanValue;
  channelBundle * chanModValue;

  abstract_coulomb() {}
  // abstractSPbasis(double densityIn, int tzMaxIn, int shellMaxIn, int NparticlesIn);
  virtual void generateIndexMap() = 0;
  virtual void generateBasis() = 0;
  virtual int checkSympqrs(int p, int q, int r, int s) = 0;
  virtual int checkModSympqrs(int p, int q, int r, int s) = 0;
  virtual int checkChanSym(int p, int q, int ichan) = 0;
  virtual int checkChanModSym(int p, int q, int ichan) = 0;
  virtual void setUpTwoStateChannels() = 0;
  virtual void printBasis() = 0;
  virtual void deallocate() = 0;

  virtual double calc_TBME(int p, int q, int r, int s, mat mapping) = 0;
};

#endif /*ABSTRACTSPBASIS_HPP*/
