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

#include "TOMATH.h"

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
    
    double const v2 = 1 - v1;
    
    std::cout << std::setprecision(40) << std::endl;
    std::cout << "v2: " << v2 << std::endl;
    
    double const v3 = sqrt(v2);
    
    std::cout << std::setprecision(40) << std::endl;
    std::cout << "v3: " << v3 << std::endl;
  
    double const v = TOSCARSSR::C() * v3;
    
    std::cout << std::setprecision(40) << std::endl;
    std::cout << "v: " << v << std::endl;

    
  long double const gamma = BeamEnergy_GeV / TOSCARSSR::kgToGeV( TOSCARSSR::Me());
    
    std::cout << sizeof(long double) << std::endl;
    std::cout << "gamma: " << gamma << std::endl;
    
  long double const beta_sqr = 1 - (1 / (gamma * gamma));
    
    std::cout << std::setprecision(40) << std::endl;
    std::cout << "beta_sqr: " << beta_sqr << std::endl;
    
  long double const beta = sqrt(beta_sqr);
    
    std::cout << std::setprecision(40) << sizeof(long double) << std::endl;
    std::cout << "beta: " << beta << std::endl;
    
    double const v4 = beta * TOSCARSSR::C() * 1e-8;
    
    std::cout << std::setprecision(40) << std::endl;
    std::cout << "v4: " << v4 << std::endl;
    
    double const w0 = v4 / R;
    
    std::cout << "w0: " << w0 << std::endl;
    
    long double const wc = (3/2) * (gamma * gamma * gamma) * ((beta * TOSCARSSR::C()) / R);
    
    std::cout << sizeof(long double) << std::endl;
    std::cout << "wc: " << wc << std::endl;

    double const psi = Angle;
    
    std::cout << std::setprecision(40) << std::endl;
    std::cout << "psi: " << psi << std::endl;
    
    double const xi = (1/2) * (w / wc) * sqrt(1 + (gamma*gamma) * (psi*psi)) * sqrt(1 + (gamma*gamma) * (psi*psi)) * sqrt(1 + (gamma*gamma) * (psi*psi));
    
    std::cout << std::setprecision(40) << std::endl;
    std::cout << "xi: " << xi << std::endl;
    
    double const myK = TOMATH::BesselK( 2/3, xi);
    
    std::cout << std::setprecision(40) << std::endl;
    std::cout << "myK: " << myK << std::endl;
    
  return myK;
}
