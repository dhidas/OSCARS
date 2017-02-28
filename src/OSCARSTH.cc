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

#include <math.h>

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
    
    double const v1 = ((TOSCARSSR::Me() * TOSCARSSR::Me()) * (TOSCARSSR::C() * TOSCARSSR::C() * TOSCARSSR::C() * TOSCARSSR::C())) / ((BeamEnergy_GeV) * (BeamEnergy_GeV));
    
    std::cout << std::endl;
  
  double const v = TOSCARSSR::C() * sqrt(1 - v1);
    
    std::cout << "v: " << v << std::endl;

  return v;
}
