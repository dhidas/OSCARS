////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Thu Feb  2 11:55:24 EST 2017
//
////////////////////////////////////////////////////////////////////

#include "OSCARSTH.h"

#include <iostream>

#include "TOSCARSSR.h"

#include <cmath>

#include <iomanip>

OSCARSTH::OSCARSTH ()
{
  // Default constructor
}



OSCARSTH::~OSCARSTH ()
{
  // Destruction
}





double OSCARSTH::UndulatorK (double const BFieldMax, double const Period) const
{
  // Return the 'K' value for an undulator with max bfield [T], Period [m]

  return BFieldMax * Period * TOSCARSSR::Qe() / (TOSCARSSR::TwoPi() * TOSCARSSR::Me() * TOSCARSSR::C());
}


double OSCARSTH::DipoleSpectrum (double const BField, double const BeamEnergy_GeV, double const Angle, TVector2D const EnergyRange_eV) const
{
  std::cout << "BField:         " << BField << std::endl;
  std::cout << "BeamEnergy_GeV: " << BeamEnergy_GeV << std::endl;
  std::cout << "Angle:          " << Angle << std::endl;
  std::cout << "EnergyRange_eV: " << EnergyRange_eV << std::endl;
    
  double const R = BeamEnergy_GeV  * 1e9 / (BField * TOSCARSSR::C());
    
  std::cout << "R: " << R << std::endl;

  // I take out the "return" otherwise what is below it never happens
  // return R;
    
  double const Q = TOSCARSSR::Qe();
    
  std::cout << "Q: " << Q << std::endl;
    
    long double const v1 = ((TOSCARSSR::Me() * TOSCARSSR::Me()) * (TOSCARSSR::C() * TOSCARSSR::C() * TOSCARSSR::C() * TOSCARSSR::C())) / ((BeamEnergy_GeV * TOSCARSSR::Qe() * (1e-9)) * (BeamEnergy_GeV * TOSCARSSR::Qe() * (1e-9)));
    
    std::setprecision(40);
    std::cout << sizeof(long double) << std::endl;
    std::cout << "v1: " << v1 << std::endl;

    double const v2 = v1 * pow(10,28);
    
    std::cout << "v2: " << v2 << std::endl;
    
    long double const v3 = v2 * pow(10, -28);
    
    std::cout << sizeof(long double) << std::endl;
    std::cout << "v3: " << v3 << std::endl;
  
  double const v = TOSCARSSR::C() * sqrt(1 - v1);
    
    std::cout << "v: " << v << std::endl;

  double const w0 = v / R;
    
    std::cout << "w0: " << w0 << std::endl;
    
  long double const beta_sqr = 1 - v1;
    
    std::cout << std::setprecision(40) << sizeof(long double) << std::endl;
    std::cout << "beta_sqr: " << beta_sqr << std::endl;
    
  long double const beta = v / TOSCARSSR::C();
    
    std::cout << std::setprecision(40) << sizeof(long double) << std::endl;
    std::cout << "beta: " << beta << std::endl;
    
    long double const gamma = BeamEnergy_GeV / TOSCARSSR::kgToGeV( TOSCARSSR::Me());
    
    std::cout << sizeof(long double) << std::endl;
    std::cout << "gamma: " << gamma << std::endl;
    
  long double const w = (3/2) * (gamma * gamma * gamma) * ((beta * TOSCARSSR::C()) / R);
    
    std::cout << sizeof(long double) << std::endl;
    std::cout << "w: " << w << std::endl;

    
  return w;
}
