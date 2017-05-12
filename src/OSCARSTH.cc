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



double OSCARSTH::UndulatorBField (double const K, double const Period) const
{
  // Return the 'BFieldMax' [T] value for an undulator with deflection parameter K, Period [m]

  return K / (Period * TOSCARSSR::Qe() / (TOSCARSSR::TwoPi() * TOSCARSSR::Me() * TOSCARSSR::C()));
}



double OSCARSTH::UndulatorPeriod (double const BFieldMax, double const K) const
{
  // Return the Period [m] value for an undulator with max bfield [T] and deflection parameter K

  return K / (BFieldMax * TOSCARSSR::Qe() / (TOSCARSSR::TwoPi() * TOSCARSSR::Me() * TOSCARSSR::C()));
}



double OSCARSTH::DipoleCriticalEnergy (double const BField, double const BeamEnergy_GeV) const
{
  // Return the critical energy in eV for dipole and electron beam
  double const Gamma = BeamEnergy_GeV / TOSCARSSR::kgToGeV( TOSCARSSR::Me());
  double const OmegaC = 3. * Gamma * Gamma * Gamma * TOSCARSSR::C() / (2. * BeamEnergy_GeV * 1e9 * TOSCARSSR::Qe() / (TOSCARSSR::Qe() * TOSCARSSR::C() * fabs(BField)));

  return TOSCARSSR::AngularFrequencyToEv(OmegaC);
}



double OSCARSTH::DipoleSpectrum (double const BField, double const BeamEnergy_GeV, double const Angle, double const Energy_eV) const
{

  double const R = BeamEnergy_GeV  * 1e9 / (BField * TOSCARSSR::C());
  double const Q = TOSCARSSR::Qe();
  double const Me = 0.511e6;
  long double const v5 = Me * Me;
  double const v6 = (BeamEnergy_GeV * 1e9) * (BeamEnergy_GeV * 1e9);
  double const v7 = v5 / v6;
  long double const v1 = ((TOSCARSSR::Me() * TOSCARSSR::Me()) * (TOSCARSSR::C() * TOSCARSSR::C() * TOSCARSSR::C() * TOSCARSSR::C())) / ((BeamEnergy_GeV * TOSCARSSR::Qe() * (1e-9)) * (BeamEnergy_GeV * TOSCARSSR::Qe() * (1e-9)));
  double const v2 = 1 - v7;
  double const v3 = sqrt(v2);
  double const v = TOSCARSSR::C() * v3;
  long double const gamma = BeamEnergy_GeV / TOSCARSSR::kgToGeV( TOSCARSSR::Me());
  long double const beta_sqr = 1 - (1 / (gamma * gamma));
  long double const beta = sqrt(beta_sqr);
  double const w0 = v / R;
  double const w01 = TOSCARSSR::AngularFrequencyToEv(w0);
  long double const wc1 = gamma * gamma * gamma;
  double const wc2 = 1.5 * wc1;
  double const wc3 = wc2 * w01;
  double const wc4 = TOSCARSSR::EvToAngularFrequency(wc3);
  double const w = TOSCARSSR::EvToAngularFrequency(Energy_eV);
  double const psi = Angle;
  double const xi = (1./2.) * (w / wc4) * sqrt(1 + (gamma*gamma) * (psi*psi)) * sqrt(1 + (gamma*gamma) * (psi*psi)) * sqrt(1 + (gamma*gamma) * (psi*psi));
  double const K2 = TOMATH::BesselK( 2./3. , xi);
  double const myK2 = K2 * K2;
  double const K1 = TOMATH::BesselK( 1. /3. , xi);
  double const myK1 = K1 * K1;
  double const pi = TOSCARSSR::Pi();
  double const Epsilon0 = TOSCARSSR::Epsilon0();
  double const c = TOSCARSSR::C();
  double const d2I1 = (3./4.) * ((Q*Q*gamma*gamma) / (4. * (pi*pi*pi) * Epsilon0 * c));
  double const d2I2 = myK2 + (((gamma*gamma*psi*psi) / (1 + (gamma*gamma*psi*psi))) * myK1);
  double const d2I3 = ((w / wc4) * (w / wc4)) * ((1 + (gamma*gamma)*(psi*psi))*(1 + (gamma*gamma)*(psi*psi)));
  double const d2I = d2I1 * d2I2 * d2I3;
  double const Hbar = TOSCARSSR::Hbar();
  double const alpha = (Q*Q) / (4. * (pi*pi*pi) * Epsilon0 * c * Hbar);
  double const I = 0.5;
  double const d = (3./4.) * (alpha) * (gamma*gamma) * (I / Q) * (0.001) * d2I3 * d2I2;
  double const d2N = d * (1e-6);

  return d2N;
}





double OSCARSTH::UndulatorFlux (double const BField, double const Period, double const NPeriods, double const BeamEnergy, double const AngleV, double const AngleH,  double const Energy_eV) const
{
  // Return the flux at a given energy and horizontal and vertical angle [photons/s/mrad^2/0.1%bw]
  return 0;

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




TVector2D OSCARSTH::UndulatorFluxOnAxisK (double const K, double const Period, double const NPeriods, int const Harmonic) const
{
  // Return the on-axis flux for this K value and harmonic

  double const BeamEnergyGeV = fParticleBeam.GetE0();

  double const n = Harmonic;
  double const nu = (n + 1.) / 2.;
  double const nu2 = (n - 1.) / 2.;
  double const JB = (n * K * K) / (4. + (2. * K * K));
  double const E = BeamEnergyGeV * 1e9;
  double const gamma = BeamEnergyGeV / TOSCARSSR::kgToGeV( TOSCARSSR::Me());
  double const JBessel1 = TOMATH::BesselJ(nu, JB);
  double const JBessel3 = TOMATH::BesselJ(nu2, JB);
  double const Q = TOSCARSSR::Qe();
  double const c = TOSCARSSR::C();
  double const Me_eV = TOSCARSSR::kgToGeV(TOSCARSSR::Me()) * 1e9;
  double const pi = TOSCARSSR::Pi();
  double const K_1 = Q * 1. * Period;
  double const K_2 = 2. * pi * Me_eV * c * c;
  double const K_0 = K_1 / K_2;
  double const z1 = n * K;
  double const z2 = 1. + ((K * K)/2.);
  double const z3 = z1 / z2;
  double const h0 = z3 * gamma * (JBessel1 - JBessel3);
  double const Epsilon0 = TOSCARSSR::Epsilon0();
  double const Hbar = TOSCARSSR::Hbar();
  double const alpha = (Q*Q) / (4. * pi * Epsilon0 * c * Hbar);
  double const I = 0.5;
  double const N = NPeriods;
  double const Nsq = N * N;
  double const h0squrd = h0 * h0;
  double const h001 = alpha * (I / Q) * Nsq * h0squrd;
  double const h002 = h001 * 0.001;
  double const h00 = h002 * 1e-6;
  double const w1 = 2. * n * c * gamma * gamma;
  double const w2 = Period * (1. + ((K * K)/ 2.));
  double const wn = w1 / w2;
  double const omega_1 = TOSCARSSR::FrequencyToEv(wn);

  return TVector2D(omega_1, h00);
}





TVector2D OSCARSTH::UndulatorFluxOnAxisB (double const BField, double const Period, double const NPeriods, int const Harmonic) const
{
  // Return the on-axis flux for this K value and harmonic

  // Undulator deflection parameter
  double const K = this->UndulatorK(BField, Period);

  return this->UndulatorFluxOnAxisK(K, Period, NPeriods, Harmonic);
}









double OSCARSTH::UndulatorFluxWeak (double const K, double const Period, double const NPeriods, double const BeamEnergyGeV, int const Harmonic) const
{
  // Return the on-axis flux for this K value and harmonic
  return 0;

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






double OSCARSTH::UndulatorEnergyAtHarmonicK (double const K,
                                             double const Period,
                                             int    const Harmonic
                                             ) const
{
  // Return the on-axis theoretical photon energy for a planar undulator

  // Even harmonics don't exist
  if (Harmonic % 2 == 0) {
    return 0;
  }

  double const Gamma = fParticleBeam.GetGamma();

  double const K2 = K * K;

  double const Lambda = Period / (2 * Gamma * Gamma) * (1. + K2 / 2.) / (double) Harmonic;
  double const Energy_eV = TOSCARSSR::FrequencyToEv(TOSCARSSR::C() / Lambda);

  return Energy_eV;
}






double OSCARSTH::UndulatorEnergyAtHarmonicB (double const BField,
                                             double const Period,
                                             int    const Harmonic
                                             ) const
{
  // Return the on-axis theoretical photon energy for a planar undulator

  // Undulator deflection parameter
  double const K = this->UndulatorK(BField, Period);

  // The value is calculated using the K value in the following function
  return this->UndulatorEnergyAtHarmonicK(K, Period, Harmonic);
}






TVector2D OSCARSTH::UndulatorBrightnessK (double const K,
                                          double const Period,
                                          int    const NPeriods,
                                          int    const N
                                          ) const
{
  // Return the on-axis theoretical brightness for a planar undulator
  if (N % 2 == 0) {
    return TVector2D(0, 0);
  }

  // Properties from beam
  double    const Gamma          = fParticleBeam.GetGamma();
  TVector2D const Beta           = fParticleBeam.GetBeta();
  TVector2D const Emittance      = fParticleBeam.GetEmittance();
  double    const Current        = fParticleBeam.GetCurrent();


  double const sigx = sqrt(Emittance[0] * Beta[0]);
  double const sigy = sqrt(Emittance[1] * Beta[1]);

  double const sigxp = sqrt(Emittance[0] / Beta[0]);
  double const sigyp = sqrt(Emittance[1] / Beta[1]);

  double const K2 = K * K;

  double const Lambda = Period / (2 * Gamma * Gamma) * (1. + K2 / 2.) / (double) N;
  double const Energy_eV = UndulatorEnergyAtHarmonicK(K, Period, N);

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
}






TVector2D OSCARSTH::UndulatorBrightnessB (double const BField,
                                          double const Period,
                                          int    const NPeriods,
                                          int    const Harmonic
                                          ) const
{
  // Return the on-axis theoretical brightness for a planar undulator

  return UndulatorBrightnessK(this->UndulatorK(BField, Period), Period, NPeriods, Harmonic);
}






void OSCARSTH::SetParticleBeam (std::string const& Beam)
{
  // Add a particle beam
  // Beam - The name of the predefined particle beam to add

  fParticleBeam.SetPredefinedBeam(Beam);
  return;
}




void OSCARSTH::SetParticleBeam (double const Energy_GeV,
                                double const Current,
                                TVector2D const& Beta,
                                TVector2D const& Emittance,
                                double const SigmaEnergyGeV)
{
  // Add a particle beam
  // Energy_GeV  - Energy of particle beam in GeV
  // Current     - Beam current in Amperes

  fParticleBeam = TParticleBeam("electron", "beam", Energy_GeV, Current);
  fParticleBeam.SetBetaEmittance(TVector3D(1, 0, 0), Beta, Emittance, TVector3D(0, 0, 0), SigmaEnergyGeV);

  return;
}




TParticleBeam& OSCARSTH::GetParticleBeam ()
{
  // Return a reference to the particle beam by a given name
  return fParticleBeam;
}




