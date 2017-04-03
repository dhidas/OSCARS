////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Thu Feb  2 11:55:24 EST 2017
//
////////////////////////////////////////////////////////////////////

#include "OSCARSTH.h"

#include <iostream>
#include <cmath>
#include <iomanip>
#include <complex>

#include "TOSCARSSR.h"
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
    
  double const w = TOSCARSSR::EvToAngularFrequency(Energy_eV);
    
    std::cout << std::setprecision(40) << "w: " << w << std::endl;
    
  double const psi = Angle;
    
    std::cout << std::setprecision(10) << std::endl;
    std::cout << "psi: " << psi << std::endl;
    
  double const xi = (1./2.) * (w / wc4) * sqrt(1 + (gamma*gamma) * (psi*psi)) * sqrt(1 + (gamma*gamma) * (psi*psi)) * sqrt(1 + (gamma*gamma) * (psi*psi));
    
    std::cout << std::setprecision(40) << std::endl;
    std::cout << "xi: " << xi << std::endl;
    
    double const K2 = TOMATH::BesselK( 2./3. , xi);
    
    std::cout << std::setprecision(40) << std::endl;
    std::cout << "K2: " << K2 << std::endl;
    
  double const myK2 = K2 * K2;
    
    std::cout << std::setprecision(40) << "myK2: " << myK2 << std::endl;
    
    double const K1 = TOMATH::BesselK( 1. /3. , xi);
    
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
    
    double const d = (3./4.) * (alpha) * (gamma*gamma) * (I / Q) * (0.001) * d2I3 * d2I2;
    
    std::cout << "d: " << d << std::endl;
    
    double const d2N = d * (1e-6);
    
    std::cout << "d2N: " << d2N << std::endl;
    
  return d2N;
}





double OSCARSTH::UndulatorFlux (double const BField, double const Period, double const NPeriods, double const BeamEnergy, double const AngleV, double const AngleH,  double const Energy_eV) const
{
  // Return the flux at a given energy and horizontal and vertical angle [photons/s/mrad^2/0.1%bw]

  // Print input fields as a check
  std::cout << "BField:       " << BField << std::endl;
  std::cout << "Period:       " << Period << std::endl;
  std::cout << "NPeriods:     " << NPeriods << std::endl;
  std::cout << "BeamEnergy:   " << BeamEnergy << std::endl;
  std::cout << "AngleV:       " << AngleV << std::endl;
  std::cout << "AngleH:       " << AngleH << std::endl;
  std::cout << "Energy_eV:    " << Energy_eV << std::endl;

  // The undulator K
  double const K = this->UndulatorK(BField, Period);
  std::cout << "K:            " << K << std::endl;

  // This is a complex number
  std::complex<double> const I(0, 1);
  std::cout << "I:            " << I << std::endl;

  return 0;
}



double J(int const n, double const x)
{
  return TOMATH::BesselJ(n, x);;
}


double OSCARSTH::UndulatorFluxKHarmonic (double const K, double const Period, double const NPeriods, double const BeamEnergyGeV, int const Harmonic) const
{
  // Return the on-axis flux for this K value and harmonic

  std::cout << "K             " << K << std::endl;
  std::cout << "Period:       " << Period << std::endl;
  std::cout << "NPeriods:     " << NPeriods << std::endl;
  std::cout << "BeamEnergyGeV:   " << BeamEnergyGeV << std::endl;
  std::cout << "Harmonic      " << Harmonic << std::endl;

    double const n = 1.;
    
    std::cout << "n: " << n << std::endl;
    
    double const nu = (n + 1.) / 2.;
    
    std::cout << "nu: " << nu << std::endl;
    
    double const JB = (n * K * K) / (4. + (2. * K * K));
    
    std::cout << "JB: " << JB << std::endl;
    
    double const E = BeamEnergyGeV * 1e9;
    
    std::cout << "E: " << E << std::endl;
    
    double const gamma = BeamEnergyGeV / TOSCARSSR::kgToGeV( TOSCARSSR::Me());
    
    std::cout << "gamma: " << gamma << std::endl;
    
    double const JBessel1 = TOMATH::BesselJ(nu, JB);
    
    std::cout << "JBessel1: " << JBessel1 << std::endl;
    
    double const z1 = n * K;
    
    std::cout << "z1: " << z1 << std::endl;
    
    double const z2 = 1. + ((K * K)/2.);
    
    std::cout << "z2: " << z2 << std::endl;
    
    double const z3 = z1 / z2;
    
    std::cout << "z3: " << z3 << std::endl;
    
    double const JBessel2 = z3 * gamma * JBessel1;
    
    std::cout << "JBessel2: " << JBessel2 << std::endl;
    
  return JBessel2;
}

double OSCARSTH::UndulatorFluxOnAxis (double const BField, double const Period, double const NPeriods, double const BeamEnergy, double const Energy_eV, int const FirstHarmonic, int const LastHarmonic) const
{
  // Return the flux at a given energy and horizontal and vertical angle [photons/s/mrad^2/0.1%bw]

  // Print input fields as a check
  std::cout << "BField:       " << BField << std::endl;
  std::cout << "Period:       " << Period << std::endl;
  std::cout << "NPeriods:     " << NPeriods << std::endl;
  std::cout << "BeamEnergy:   " << BeamEnergy << std::endl;
  std::cout << "Energy_eV:    " << Energy_eV << std::endl;
  std::cout << "FirstHarmonic " << FirstHarmonic << std::endl;
  std::cout << "LastHarmonic  " << LastHarmonic << std::endl;

  // The undulator K
  double const K = this->UndulatorK(BField, Period);
  double const K2 = K * K;

  return 0;
}








double OSCARSTH::UndulatorFluxWeak (double const K, double const Period, double const NPeriods, double const BeamEnergyGeV, int const Harmonic) const
{
  // Return the on-axis flux for this K value and harmonic

  std::cout << "K             " << K << std::endl;
  std::cout << "Period:       " << Period << std::endl;
  std::cout << "NPeriods:     " << NPeriods << std::endl;
  std::cout << "BeamEnergyGeV:   " << BeamEnergyGeV << std::endl;
  std::cout << "Harmonic      " << Harmonic << std::endl;


  double const Gamma = BeamEnergyGeV / TOSCARSSR::kgToGeV( TOSCARSSR::Me());

  return Gamma;
}





double OSCARSTH::DipoleBrightness () const
{
  double const BField = 0.4;
  double const BeamEnergy_GeV = 3;
  double const Energy_eV = 123;
  double const Angle = 0;
  double const Current = 0.500;
  double const EmittanceX = 0.9e-9;
  double const EmittanceY = 0.008e-9;
  double const SigmaE = 0.005;
  
  double const BetaX = pow(44.2e-6, 2) / EmittanceX;
  double const BetaY = pow(15.7e-6, 2) / EmittanceY;
  double const EtaX = 0;
  double const EtaY = 0;

  double const Gamma = BeamEnergy_GeV / TOSCARSSR::kgToGeV( TOSCARSSR::Me());
  double const Alpha = 1. / 137.;

  double const Omega = TOSCARSSR::EvToAngularFrequency(Energy_eV);
  double const OmegaC = 3. * Gamma * Gamma * Gamma * TOSCARSSR::C() / (2. * BeamEnergy_GeV * 1e9 * TOSCARSSR::Qe() / (TOSCARSSR::Qe() * TOSCARSSR::C() * BField));

  double const Y = Omega / OmegaC;

  double const D1F = sqrt(3) / TOSCARSSR::TwoPi() * Alpha * Gamma * 0.001 * Current / TOSCARSSR::Qe() * Y * TOMATH::BesselK_IntegralToInfty(5. / 3., Y);
  double const D2F = 3. * Alpha / (4 * TOSCARSSR::Pi2()) * Gamma * Gamma * 0.001 * Current / TOSCARSSR::Qe() * Y * Y * pow(TOMATH::BesselK(2. / 3., Y / 2.), 2);

  double const SigmaPsi = 1. / TOSCARSSR::Sqrt2Pi() * D1F / D2F;

  // Diffraction limited source size
  //double const SigmaR = (TOSCARSSR::C() * TOSCARSSR::TwoPi() / Omega) / (TOSCARSSR::FourPi() * SigmaPsi);
  double const SigmaR = (TOSCARSSR::C() / Omega) / (2. * SigmaPsi);


  double const SigmaX = sqrt(EmittanceX * BetaX + EtaX * EtaX * SigmaE * SigmaE + SigmaR * SigmaR);
  //double const SigmaY = sqrt(EmittanceY * BetaY + SigmaR * SigmaR + (EmittanceY * EmittanceY + EmittanceY * GammaY * SigmaR * SigmaR) / (SigmaPsi * SigmaPsi));


  return TOSCARSSR::AngularFrequencyToEv(OmegaC);
  return TOMATH::BesselK_IntegralToInfty(2./3., 1);
}






TVector2D OSCARSTH::UndulatorBrightness (double const BField,
                                         double const Period,
                                         int    const NPeriods,
                                         int    const N,
                                         double const BeamEnergy_GeV,
                                         double const SigmaE,
                                         double const Current,
                                         TVector2D const Beta,
                                         TVector2D const Emittance
                                         ) const
{
  // Return the on-axis theoretical brightness for a planar undulator
  if (N % 2 == 0) {
    return TVector2D(0, 0);
  }

  double const sigx = sqrt(Emittance[0] * Beta[0]);
  double const sigy = sqrt(Emittance[1] * Beta[1]);

  double const sigxp = sqrt(Emittance[0] / Beta[0]);
  double const sigyp = sqrt(Emittance[1] / Beta[1]);

  double const Gamma = BeamEnergy_GeV / TOSCARSSR::kgToGeV( TOSCARSSR::Me());

  double const K = this->UndulatorK(BField, Period);
  double const K2 = K * K;

  double const Lambda = Period / (2 * Gamma * Gamma) * (1. + K2 / 2.) / (double) N;
  double const Energy_eV = TOSCARSSR::FrequencyToEv(TOSCARSSR::C() / Lambda);

  double const Fn = K2 * N * N / pow(1. + K2 / 2., 2) * pow(
      TOMATH::BesselJ( (N - 1) / 2, N * K2 / (4 * (1. + 0.5 * K2))) - TOMATH::BesselJ( (N + 1) / 2, N * K2 / (4 * (1. + 0.5 * K2))),
      2);

  double const Qn = (1. + K2 / 2.) * Fn / (double) N;

  double const Fu = TOSCARSSR::Pi() * TOSCARSSR::Alpha() * NPeriods * 0.001 * Current / TOSCARSSR::Qe() * Qn;


  double const sigr = 1 / TOSCARSSR::FourPi() * sqrt(Lambda * Period * NPeriods);
  double const sigrp = sqrt(Lambda / (Period * NPeriods));
  double const SigmaX = sqrt(sigx * sigx + sigr * sigr);
  double const SigmaY = sqrt(sigy * sigy + sigr * sigr);
  double const SigmaXP = sqrt(sigxp * sigxp + sigrp * sigrp);
  double const SigmaYP = sqrt(sigyp * sigyp + sigrp * sigrp);

  return TVector2D(Energy_eV, Fu / (4 * TOSCARSSR::Pi2() * SigmaX * SigmaY * SigmaXP * SigmaYP) * 1.e-12);
  //return TVector2D(Energy_eV, Fu / (TOSCARSSR::TwoPi() * SigmaXP * SigmaYP) * 1e-6);
}
