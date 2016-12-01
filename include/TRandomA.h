#ifndef GUARD_TRandomA_h
#define GUARD_TRandomA_h
////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Mon Aug 29 17:22:42 EDT 2016
//
////////////////////////////////////////////////////////////////////

#include <random>



class TRandomA
{
  public:
    TRandomA ();
    TRandomA (int const);
    ~TRandomA ();

    void SetSeed (int const);

    double Normal ();
    double Uniform ();

  private:
    std::random_device* fRD;
    std::mt19937* fMT;

    std::normal_distribution<double> fNormalDist;
    std::uniform_real_distribution<double> fUniformDist;


};













#endif
