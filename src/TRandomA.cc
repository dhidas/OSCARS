////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Mon Aug 29 17:22:42 EDT 2016
//
////////////////////////////////////////////////////////////////////

#include "TRandomA.h"

// Global random generator
TRandomA* gRandomA = new TRandomA();



TRandomA::TRandomA ()
{
  fRD = new std::random_device();
  fMT = new std::mt19937((*fRD)());
  fNormalDist  = std::normal_distribution<double>(0, 1);
  fUniformDist = std::uniform_real_distribution<double>(0, 1);
}



TRandomA::TRandomA (int const Seed)
{
  fMT = new std::mt19937(Seed);
  fNormalDist  = std::normal_distribution<double>(0, 1);
  fUniformDist = std::uniform_real_distribution<double>(0, 1);
}



TRandomA::~TRandomA ()
{
  delete fMT;
}



void TRandomA::SetSeed (int const Seed)
{
  delete fMT;
  fMT = new std::mt19937(Seed);

  return;
}




double TRandomA::Normal ()
{
  return fNormalDist(*fMT);
}




double TRandomA::Uniform ()
{
  return fUniformDist(*fMT);
}
