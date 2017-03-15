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


double OSCARSTH::DipoleSpectrum (double const BField, double const BeamEnergy_GeV, double const Angle, double const Energy_eV) const
{
  std::cout << "BField:         " << BField << std::endl;
  std::cout << "BeamEnergy_GeV: " << BeamEnergy_GeV << std::endl;
  std::cout << "Angle:          " << Angle << std::endl;
  std::cout << "Energy_eV:      " << Energy_eV << std::endl;
  //std::cout << "EnergyRange_eV: " << EnergyRange_eV << std::endl;
    
  double const R = BeamEnergy_GeV  * 1e9 / (BField * TOSCARSSR::C());
    
  std::cout << "R: " << R << std::endl;

  // I take out the "return" otherwise what is below it never happens
  // return R;
    
  double const Q = TOSCARSSR::Qe();
    
  std::cout << "Q: " << Q << std::endl;
   
  double const Me = 0.511e6;
    
    std::cout << "Me: " << Me << std::endl;
    
  long double const v5 = Me * Me;
    
    std::cout << std::setprecision(40) << "v5: " << v5 << std::endl;
    
    double const v6 = (BeamEnergy_GeV * 1e9) * (BeamEnergy_GeV * 1e9);
    
    std::cout << std::setprecision(40) <<"v6: " << v6 << std::endl;
    
    double const v7 = v5 / v6;
    
    std::cout << std::setprecision(40) <<"v7: " << v7 << std::endl;
    
  long double const v1 = ((TOSCARSSR::Me() * TOSCARSSR::Me()) * (TOSCARSSR::C() * TOSCARSSR::C() * TOSCARSSR::C() * TOSCARSSR::C())) / ((BeamEnergy_GeV * TOSCARSSR::Qe() * (1e-9)) * (BeamEnergy_GeV * TOSCARSSR::Qe() * (1e-9)));
    
    std::cout << std::setprecision(40) << "v1: " << v1 << std::endl;
    
    double const v2 = 1 - v7;
    
    std::cout << std::setprecision(40) << "v2: " << v2 << std::endl;
    
    double const v3 = sqrt(v2);
    
    std::cout << std::setprecision(40) << "v3: " << v3 << std::endl;
  
  double const v = TOSCARSSR::C() * v3;
    
    std::cout << std::setprecision(40) << "v: " << v << std::endl;

    
  long double const gamma = BeamEnergy_GeV / TOSCARSSR::kgToGeV( TOSCARSSR::Me());
    
    std::cout << std::setprecision(40) << std::endl;
    std::cout << "gamma: " << gamma << std::endl;
    
  long double const beta_sqr = 1 - (1 / (gamma * gamma));
    
    std::cout << std::setprecision(40) << std::endl;
    std::cout << "beta_sqr: " << beta_sqr << std::endl;
    
  long double const beta = sqrt(beta_sqr);
    
    std::cout << std::setprecision(40) << "beta: " << beta << std::endl;
    
  double const w0 = v / R;
    
    std::cout << std::setprecision(40) << std::endl;
    std::cout << "w0: " << w0 << std::endl;
    
    double const w01 = TOSCARSSR::AngularFrequencyToEv(w0);
    
    std::cout << std::setprecision(40) << "w01: " << w01 << std::endl;
    
    long double const wc1 = gamma * gamma * gamma;
    
    std::cout << std::setprecision(40) << "wc1: " << wc1 << std::endl;
    
    double const wc2 = 1.5 * wc1;
    
    std::cout << std::setprecision(40) << "wc2: " << wc2 << std::endl;

    double const wc3 = wc2 * w01;
    
    std::cout << std::setprecision(40) << "wc3: " << wc3 << std::endl;
    
  double const wc4 = TOSCARSSR::EvToAngularFrequency(wc3);
    
    std::cout << std::setprecision(40) << "wc4: " << wc4 << std::endl;
    
  double const w = TOSCARSSR::EvToAngularFrequency(2390);
    
    std::cout << std::setprecision(40) << "w: " << w << std::endl;
    
  double const psi = -(3 * Angle);
    
    std::cout << std::setprecision(10) << std::endl;
    std::cout << "psi: " << psi << std::endl;
    
  double const xi = (1./2.) * (w / wc4) * sqrt(1 + (gamma*gamma) * (psi*psi)) * sqrt(1 + (gamma*gamma) * (psi*psi)) * sqrt(1 + (gamma*gamma) * (psi*psi));
    
    std::cout << std::setprecision(40) << std::endl;
    std::cout << "xi: " << xi << std::endl;
    
    double const K2 = TOMATH::BesselK( 2/3, xi);
    
    std::cout << std::setprecision(40) << std::endl;
    std::cout << "K2: " << K2 << std::endl;
    
  double const myK2 = K2 * K2;
    
    std::cout << std::setprecision(40) << "myK2: " << myK2 << std::endl;
    
    double const K1 = TOMATH::BesselK( 1/3, xi);
    
    std::cout << std::setprecision(40) << "K1: " << K1 << std::endl;
    
  double const myK1 = K1 * K1;
    
    std::cout << std::setprecision(40) << "myK1: " << myK1 << std::endl;
    
  double const pi = TOSCARSSR::Pi();
    
    std::cout << std::endl;
    std::cout << "pi: " << pi << std::endl;
    
  double const Epsilon0 = TOSCARSSR::Epsilon0();
    
    std::cout << "Epsilon0: " << Epsilon0 << std::endl;
    
  double const c = TOSCARSSR::C();
    
    std::cout << "c: " << c << std::endl;
    
  double const d2I1 = (3./4.) * ((Q*Q*gamma*gamma) / (4. * (pi*pi*pi) * Epsilon0 * c));
    
    std::cout << std::endl;
    std::cout << "d2I1: " << d2I1 << std::endl;
    
    double const d2I2 = myK2 + (((gamma*gamma*psi*psi) / (1 + (gamma*gamma*psi*psi))) * myK1);
    
    std::cout << "d2I2: " << d2I2 << std::endl;
    
    double const d2I3 = ((w / wc4) * (w / wc4)) * ((1 + (gamma*gamma)*(psi*psi))*(1 + (gamma*gamma)*(psi*psi)));
    
    std::cout << "d2I3: " << d2I3 << std::endl;
    
  double const d2I = d2I1 * d2I2 * d2I3;
    
    std::cout << "d2I: " << d2I << std::endl;
    
    double const Hbar = TOSCARSSR::Hbar();
    
    std::cout << std::endl;
    std::cout << "Hbar: " << Hbar << std::endl;
    
    double const alpha = (Q*Q) / (4. * (pi*pi*pi) * Epsilon0 * c * Hbar);
    
    std::cout << "alpha: " << alpha << std::endl;
    
    double const I = 0.5;
    
    std::cout << "I: " << I << std::endl;
    
    double const d = (3./4.) * (alpha) * (gamma*gamma) * (I / Q) * (1.) * d2I3 * d2I2;
    
    std::cout << "d: " << d << std::endl;
    
  return d;
}
