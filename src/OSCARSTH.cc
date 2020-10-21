////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Thu Feb  2 11:55:24 EST 2017
//
////////////////////////////////////////////////////////////////////

#include "OSCARSTH.h"
#include "OSCARSTH_Cuda.h"

#include <iostream>
#include <cmath>
#include <iomanip>
#include <complex>

#include "TOSCARS.h"
#include "TOMATH.h"





OSCARSTH::OSCARSTH ()
{
  // Default constructor
  SetUseGPUGlobal(0);   // GPU off by default
  SetNThreadsGlobal(1); // Use N threads for calculations by default
}



OSCARSTH::~OSCARSTH ()
{
  // Destruction
}





double OSCARSTH::UndulatorK (double const BFieldMax, double const Period) const
{
  // Return the 'K' value for an undulator with max bfield [T], Period [m]

  return BFieldMax * Period * TOSCARS::Qe() / (TOSCARS::TwoPi() * TOSCARS::Me() * TOSCARS::C());
}



double OSCARSTH::UndulatorBField (double const K, double const Period) const
{
  // Return the 'BFieldMax' [T] value for an undulator with deflection parameter K, Period [m]

  return K / (Period * TOSCARS::Qe() / (TOSCARS::TwoPi() * TOSCARS::Me() * TOSCARS::C()));
}



double OSCARSTH::UndulatorPeriod (double const BFieldMax, double const K) const
{
  // Return the Period [m] value for an undulator with max bfield [T] and deflection parameter K

  return K / (BFieldMax * TOSCARS::Qe() / (TOSCARS::TwoPi() * TOSCARS::Me() * TOSCARS::C()));
}



double OSCARSTH::DipoleCriticalEnergy (double const BField) const
{

  // Return the critical energy in eV for dipole and electron beam
  double const BeamEnergy_GeV = fParticleBeam.GetE0();
  double const Gamma = BeamEnergy_GeV / TOSCARS::kgToGeV( fParticleBeam.GetM());
  //double const OmegaC = 3. * Gamma * Gamma * Gamma * TOSCARS::C() / (2. * BeamEnergy_GeV * 1e9 * TOSCARS::Qe() / (TOSCARS::Qe() * TOSCARS::C() * fabs(BField)));
  double const OmegaC = 3. * Gamma * Gamma * Gamma * TOSCARS::C() / (2. * BeamEnergy_GeV * 1e9 / (TOSCARS::C() * fabs(BField)));

  return TOSCARS::AngularFrequencyToEv(OmegaC);
}




void OSCARSTH::DipoleSpectrumEnergy (double const BField, 
                                     TSpectrumContainer& Spectrum,
                                     double const Angle) const
{
  // Beam energy from internal beam
  double const BeamEnergy_GeV =  fParticleBeam.GetE0();

  // Calculate flux for each point
  for (size_t i = 0; i < Spectrum.GetNPoints(); ++i) {
    Spectrum.SetFlux(i, DipoleSpectrum(BField, BeamEnergy_GeV, Angle, Spectrum.GetEnergy(i)));
  }

  return;
}




void OSCARSTH::DipoleSpectrumAngle (double const BField, 
                                    TSpectrumContainer& Spectrum,
                                    double const Energy_eV) const
{
  // Beam energy from internal beam
  double const BeamEnergy_GeV =  fParticleBeam.GetE0();

  // Calculate flux for each point
  for (size_t i = 0; i < Spectrum.GetNPoints(); ++i) {

    // Here Spectrum.GetEnergy is really getting the angle
    Spectrum.SetFlux(i, DipoleSpectrum(BField, BeamEnergy_GeV, Spectrum.GetEnergy(i), Energy_eV));
  }

  return;
}




void OSCARSTH::DipoleSpectrumEnergyAngleIntegrated (double const BField, 
                                     TSpectrumContainer& Spectrum) const
{
  // Beam energy from internal beam
  double const BeamEnergy_GeV =  fParticleBeam.GetE0();

  // Calculate flux for each point
  for (size_t i = 0; i < Spectrum.GetNPoints(); ++i) {
    Spectrum.SetFlux(i, DipoleSpectrumAngleIntegrated(BField, BeamEnergy_GeV, Spectrum.GetEnergy(i)));
  }


  return;
}




double OSCARSTH::DipoleSpectrum (double const BField, double const BeamEnergy_GeV, double const Angle, double const Energy_eV) const
{

  double const y = Energy_eV / this->DipoleCriticalEnergy(BField);
  double const Gamma = fParticleBeam.GetGamma();
  double const Current = fParticleBeam.GetCurrent();

  double const X = Gamma * Angle;
  double const Xi = y / 2 * pow(1 + X * X, 3/2);
  double const C0 = 3.0 * TOSCARS::Alpha() / (4*TOSCARS::Pi2()) * 0.001 * Gamma * Gamma * Current / fParticleBeam.GetQ() * y * y * pow(1 + X * X, 2);

  return C0 * (pow(TOMATH::BesselK( 2. / 3., Xi), 2) + X * X / (1 + X * X) * pow(TOMATH::BesselK( 1. / 3., Xi), 2)) * 1e-6;
}




double OSCARSTH::WigglerFlux (double const BField, double const Period, double const Angle, double const Energy_eV) const
{

  double const K = this->UndulatorK(BField, Period);
  double const Gamma = fParticleBeam.GetGamma();
  double const y = Energy_eV / (this->DipoleCriticalEnergy(BField) * sqrt(1 - pow(Angle / (K / Gamma), 2)));
  double const Current = fParticleBeam.GetCurrent();

  double const X = Gamma * Angle;
  double const Xi = y / 2 * pow(1 + X * X, 3/2);
  double const C0 = 3.0 * TOSCARS::Alpha() / (4*TOSCARS::Pi2()) * 0.001 * Gamma * Gamma * Current / TOSCARS::Qe() * y * y * pow(1 + X * X, 2);

  return C0 * (pow(TOMATH::BesselK( 2. / 3., Xi), 2) + X * X / (1 + X * X) * pow(TOMATH::BesselK( 1. / 3., Xi), 2)) * 1e-6;
}




double OSCARSTH::WigglerFluxAngleIntegrated (double const BField, double const Period, double const Energy_eV) const
{

  double const Gamma = fParticleBeam.GetGamma();
  double const K = this->UndulatorK(BField, Period);
  double const Current = fParticleBeam.GetCurrent();

  int const N = 1000;
  double const Range = 3 / Gamma;
  double const da = Range / N;

  double Angle = -da;
  double y;
  double X;
  double Xi;
  double C0;
  double Sum = 0;
  for (int i = 0; i != N; ++i) {
    Angle += da;
    y = Energy_eV / (this->DipoleCriticalEnergy(BField) * sqrt(1 - pow(Angle / (K / Gamma), 2)));

    X = Gamma * Angle;
    Xi = y / 2 * pow(1 + X * X, 3/2);
    C0 = 3.0 * TOSCARS::Alpha() / (4*TOSCARS::Pi2()) * 0.001 * Gamma * Gamma * Current / TOSCARS::Qe() * y * y * pow(1 + X * X, 2);

    Sum += C0 * (pow(TOMATH::BesselK( 2. / 3., Xi), 2) + X * X / (1 + X * X) * pow(TOMATH::BesselK( 1. / 3., Xi), 2)) * 1e-3;
  }

  return Sum * da * 2;
}




/*
double OSCARSTH::DipoleSpectrum_OLD (double const BField, double const BeamEnergy_GeV, double const Angle, double const Energy_eV) const
{

  double const R = BeamEnergy_GeV  * 1e9 / (BField * TOSCARS::C());
  double const Q = TOSCARS::Qe();
  double const Me = TOSCARS::Me();
  long double const v5 = Me * Me;
  double const v6 = (BeamEnergy_GeV * 1e9) * (BeamEnergy_GeV * 1e9);
  double const v7 = v5 / v6;
  //long double const v1 = ((TOSCARS::Me() * TOSCARS::Me()) * (TOSCARS::C() * TOSCARS::C() * TOSCARS::C() * TOSCARS::C())) / ((BeamEnergy_GeV * TOSCARS::Qe() * (1e-9)) * (BeamEnergy_GeV * TOSCARS::Qe() * (1e-9)));
  double const v2 = 1 - v7;
  double const v3 = sqrt(v2);
  double const v = TOSCARS::C() * v3;
  long double const gamma = BeamEnergy_GeV / TOSCARS::kgToGeV( TOSCARS::Me());
  //long double const beta_sqr = 1 - (1 / (gamma * gamma));
  //long double const beta = sqrt(beta_sqr);
  double const w0 = v / R;
  double const w01 = TOSCARS::AngularFrequencyToEv(w0);
  long double const wc1 = gamma * gamma * gamma;
  double const wc2 = 1.5 * wc1;
  double const wc3 = wc2 * w01;
  double const wc4 = TOSCARS::EvToAngularFrequency(wc3);
  double const w = TOSCARS::EvToAngularFrequency(Energy_eV);
  double const psi = Angle;
  double const xi = (1./2.) * (w / wc4) * sqrt(1 + (gamma*gamma) * (psi*psi)) * sqrt(1 + (gamma*gamma) * (psi*psi)) * sqrt(1 + (gamma*gamma) * (psi*psi));
  double const K2 = TOMATH::BesselK( 2./3. , xi);
  double const myK2 = K2 * K2;
  double const K1 = TOMATH::BesselK( 1. /3. , xi);
  double const myK1 = K1 * K1;
  double const pi = TOSCARS::Pi();
  double const Epsilon0 = TOSCARS::Epsilon0();
  double const c = TOSCARS::C();
  //double const d2I1 = (3./4.) * ((Q*Q*gamma*gamma) / (4. * (pi*pi*pi) * Epsilon0 * c));
  double const d2I2 = myK2 + (((gamma*gamma*psi*psi) / (1 + (gamma*gamma*psi*psi))) * myK1);
  double const d2I3 = ((w / wc4) * (w / wc4)) * ((1 + (gamma*gamma)*(psi*psi))*(1 + (gamma*gamma)*(psi*psi)));
  //double const d2I = d2I1 * d2I2 * d2I3;
  double const Hbar = TOSCARS::Hbar();
  double const alpha = (Q*Q) / (4. * (pi*pi*pi) * Epsilon0 * c * Hbar);
  double const I = 0.5;
  double const d = (3./4.) * (alpha) * (gamma*gamma) * (I / Q) * (0.001) * d2I3 * d2I2;
  double const d2N = d * (1e-6);

  return d2N;
}
*/




double OSCARSTH::DipoleSpectrumAngleIntegrated (double const BField, double const BeamEnergy_GeV, double const Energy_eV) const
{
  // This function calculates the dipole spectrum integrated over the vertical angle and returns the result as:
  // [0.1%bw / mrad]

  double const Gamma = fParticleBeam.GetGamma();
  double const Current = fParticleBeam.GetCurrent();
  double const y = Energy_eV / this->DipoleCriticalEnergy(BField);

  return sqrt(3.0) / TOSCARS::TwoPi() * TOSCARS::Alpha() * Gamma * 0.001 * Current / TOSCARS::Qe() * y * TOMATH::BesselK_IntegralToInfty(5./3., y) * 0.001;
}


/*
double OSCARSTH::DipoleSpectrumAngleIntegrated_OLD (double const BField, double const BeamEnergy_GeV, double const Energy_eV) const
{
  // This function calculates the dipole spectrum integrated over the vertical angle and returns the result as:
  // [0.1%bw / mrad]

  double const Gamma = BeamEnergy_GeV / TOSCARS::kgToGeV( TOSCARS::Me());
  double const OmegaC = 3. * Gamma * Gamma * Gamma * TOSCARS::C() / (2. * BeamEnergy_GeV * 1e9 * TOSCARS::Qe() / (TOSCARS::Qe() * TOSCARS::C() * fabs(BField)));
  double const Omega = TOSCARS::EvToAngularFrequency(Energy_eV);
  double const y = Omega / OmegaC;
  return sqrt(3.0) / TOSCARS::TwoPi() * TOSCARS::Alpha() * Gamma * 0.001 * y * fParticleBeam.GetCurrent() / TOSCARS::Qe() * TOMATH::BesselK_IntegralToInfty(5./3., y) * 0.001;
}
*/





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




TVector2D OSCARSTH::UndulatorFluxOnAxisK (double const K,
                                          double const Period,
                                          double const NPeriods,
                                          int const N) const
{
  // Return the on-axis flux for this K value and harmonic

  // Properties from beam
  double    const Gamma          = fParticleBeam.GetGamma();
  TVector2D const Beta           = fParticleBeam.GetTwissBeta();
  TVector2D const Emittance      = fParticleBeam.GetEmittance();
  double    const Current        = fParticleBeam.GetCurrent();

  // Check that we can do this calculation, else reject
  if (Gamma == 0 || Beta[0] == 0 || Beta[1] == 0 || Emittance[0] == 0 || Emittance[1] == 0 || Current == 0) {
    throw std::invalid_argument("Beam definition incorrect for this calculation: Check energy, current, beta, emittance");
  }

  double const sigx = sqrt(Emittance[0] * Beta[0]);
  double const sigy = sqrt(Emittance[1] * Beta[1]);

  double const sigxp = sqrt(Emittance[0] / Beta[0]);
  double const sigyp = sqrt(Emittance[1] / Beta[1]);

  double const K2 = K * K;

  double const Lambda = Period / (2 * Gamma * Gamma) * (1. + K2 / 2.) / (double) N;
  double const Energy_eV = UndulatorEnergyAtHarmonicK(K, Period, N);

  if (N % 2 != 1) {
    return TVector2D(Energy_eV, 0);
  }


  double const Fn = K2 * N * N / pow(1. + K2 / 2., 2) * pow(
      TOMATH::BesselJ( (N - 1) / 2, N * K2 / (4 * (1. + 0.5 * K2))) - TOMATH::BesselJ( (N + 1) / 2, N * K2 / (4 * (1. + 0.5 * K2))),
      2);

  double const Qn = (1. + K2 / 2.) * Fn / (double) N;

  double const Fu = TOSCARS::Pi() * TOSCARS::Alpha() * NPeriods * 0.001 * Current / TOSCARS::Qe() * Qn;


  double const sigr = 1 / TOSCARS::FourPi() * sqrt(Lambda * Period * NPeriods);
  double const sigrp = sqrt(Lambda / (Period * NPeriods));
  double const SigmaXP = sqrt(sigxp * sigxp + sigrp * sigrp);
  double const SigmaYP = sqrt(sigyp * sigyp + sigrp * sigrp);


  return TVector2D(Energy_eV, Fu / (TOSCARS::Pi2() * SigmaXP * SigmaYP) * 1.e-6);
}





TVector2D OSCARSTH::UndulatorFluxOnAxisB (double const BField,
                                          double const Period,
                                          double const NPeriods,
                                          int const N) const
{
  // Return the on-axis flux for this K value and harmonic

  // Undulator deflection parameter
  double const K = this->UndulatorK(BField, Period);

  return this->UndulatorFluxOnAxisK(K, Period, NPeriods, N);
}








TVector2D OSCARSTH::UndulatorFluxB (double const BField,
                                    double const Period,
                                    int    const NPeriods,
                                    int    const Harmonic
                                    ) const
{
  // Return the on-axis theoretical brightness for a planar undulator

  return UndulatorFluxK(this->UndulatorK(BField, Period), Period, NPeriods, Harmonic);
}








TVector2D OSCARSTH::UndulatorFluxK (double const K,
                                    double const Period,
                                    int    const NPeriods,
                                    int    const N
                                    ) const
{
  // Return the on-axis flux for a planar undulator
  if (N % 2 == 0) {
    return TVector2D(0, 0);
  }


  // Properties from beam
  double    const Gamma          = fParticleBeam.GetGamma();
  TVector2D const Beta           = fParticleBeam.GetTwissBeta();
  TVector2D const Emittance      = fParticleBeam.GetEmittance();
  double    const Current        = fParticleBeam.GetCurrent();

  // Check that we can do this calculation, else reject
  if (Gamma == 0 || Beta[0] == 0 || Beta[1] == 0 || Emittance[0] == 0 || Emittance[1] == 0 || Current == 0) {
    throw std::invalid_argument("Beam definition incorrect for this calculation: Check energy, current, beta, emittance");
  }

  double const sigx = sqrt(Emittance[0] * Beta[0]);
  double const sigy = sqrt(Emittance[1] * Beta[1]);

  double const sigxp = sqrt(Emittance[0] / Beta[0]);
  double const sigyp = sqrt(Emittance[1] / Beta[1]);

  double const K2 = K * K;

  double const Lambda = Period / (2 * Gamma * Gamma) * (1. + K2 / 2.) / (double) N;
  double const Energy_eV = UndulatorEnergyAtHarmonicK(K, Period, N);

  if (N % 2 != 1) {
    return TVector2D(Energy_eV, 0);
  }

  double const Fn = K2 * N * N / pow(1. + K2 / 2., 2) * pow(
      TOMATH::BesselJ( (N - 1) / 2, N * K2 / (4 * (1. + 0.5 * K2))) - TOMATH::BesselJ( (N + 1) / 2, N * K2 / (4 * (1. + 0.5 * K2))),
      2);

  double const Qn = (1. + K2 / 2.) * Fn / (double) N;

  double const Fu = TOSCARS::Pi() * TOSCARS::Alpha() * NPeriods * 0.001 * Current / TOSCARS::Qe() * Qn;




  return TVector2D(Energy_eV, Fu);
}






double OSCARSTH::UndulatorTotalPower (double const K,
                                      double const Period,
                                      int const NPeriods
                                      ) const
{
  // This function is taken from Kwang-Je KIM
  // Nuclear Instruments and Methods in Physics Research A246 (1986) 67-70
  // ANGULAR DISTRIBUTION OF UNDULATOR POWER FOR AN ARBITRARY DEFLECTION PARAMETER K *

  // Properties from beam
  double    const Gamma          = fParticleBeam.GetGamma();
  double    const Current        = fParticleBeam.GetCurrent();

  double const Z0 = 376.73031346177; // This is the vacuum impedence of free space

  // Check that we can do this calculation, else reject
  if (Gamma == 0 || Current == 0) {
    throw std::invalid_argument("Beam definition incorrect for this calculation: Check gamma, current");
  }

  return NPeriods / 6 * Z0 * Current * TOSCARS::TwoPi() * TOSCARS::Qe() * TOSCARS::C() / Period * Gamma*Gamma * K*K;
}



TVector2D OSCARSTH::UndulatorPowerDensity (double const        K,
                                           double const        Period,
                                           int    const        NPeriods,
                                           T3DScalarContainer& PowerDensityContainer
                                           ) const
{
  // This function is taken from Kwang-Je KIM
  // Nuclear Instruments and Methods in Physics Research A246 (1986) 67-70
  // ANGULAR DISTRIBUTION OF UNDULATOR POWER FOR AN ARBITRARY DEFLECTION PARAMETER K *

  // Return the on-axis theoretical brightness for a planar undulator


  // Properties from beam
  double    const Gamma          = fParticleBeam.GetGamma();
  double    const Current        = fParticleBeam.GetCurrent();

  // Check that we can do this calculation, else reject
  if (Gamma == 0 || Current == 0) {
    throw std::invalid_argument("Beam definition incorrect for this calculation: Check energy, current, beta, emittance");
  }

  // Get total power
  double const TotalPower = this->UndulatorTotalPower(K, Period, NPeriods);

  double const Factor = 21 * Gamma*Gamma * TotalPower / (7 * TOSCARS::Pi2());

  int const N = 1000;
  double const dalpha = TOSCARS::TwoPi() / (double) N;


  for (int ip = 0; ip != PowerDensityContainer.GetNPoints(); ++ip) {
    double const Theta = PowerDensityContainer.GetPoint(ip).GetX().GetX();
    double const Psi   = PowerDensityContainer.GetPoint(ip).GetX().GetY();

    double const D1 = 1 + pow(Gamma*Psi, 2);
    double const GammaTimesTheta = Gamma*Theta;

    double Integral = 0;
    for (int i = 0; i != N; ++i) {
      double const alpha = -TOSCARS::Pi() + i*dalpha;

      double const D = D1 + pow(GammaTimesTheta - K*cos(alpha), 2);

      Integral += (1/pow(D, 3) - (4*pow(GammaTimesTheta - K*cos(alpha), 2)/pow(D, 5))) * pow(sin(alpha), 2);
    }
    Integral *= dalpha * Factor;

    Integral *= 1e-6; // To urad^2

    PowerDensityContainer.AddToPoint(ip, Integral);
  }
  return TVector2D(0, 0);
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


  double const Gamma = BeamEnergyGeV / TOSCARS::kgToGeV( TOSCARS::Me());

  return Gamma;
}





void OSCARSTH::DipoleBrightness (double const BField,
                                 TSpectrumContainer& SpectrumContainer) const
{

  // Beam energy from internal beam
  double const BeamEnergy_GeV =  fParticleBeam.GetE0();
  double const sigma_E = fParticleBeam.GetSigmaEnergyGeV() / fParticleBeam.GetE0();

  TVector2D const Beta = fParticleBeam.GetTwissBeta();
  TVector2D const Alpha = fParticleBeam.GetTwissAlpha();
  TVector2D const Emittance = fParticleBeam.GetEmittance();


  double lam = fabs(TOSCARS::H() * TOSCARS::C() / (1000. * fParticleBeam.GetCharge()));
  double const eta_x = fParticleBeam.GetEta().GetX();
  double const beta_x = Beta.GetX();
  double const epsilon_x = Emittance.GetX();

  double const epsilon_y = Emittance.GetY();
  double const beta_y = Beta.GetY();
  double const alpha_y = Alpha.GetY();
  double const gamma_y = (1. + alpha_y*alpha_y) / beta_y;

  for (size_t i = 0; i != SpectrumContainer.GetNPoints(); ++i) {
    double const dfdt = this->DipoleSpectrumAngleIntegrated(BField, BeamEnergy_GeV, SpectrumContainer.GetEnergy(i));
    double const df2dtdp = this->DipoleSpectrum(BField, BeamEnergy_GeV, 0, SpectrumContainer.GetEnergy(i));

    double const sigma_psi = dfdt / (df2dtdp * sqrt(TOSCARS::TwoPi()));


    double const sigma_r = lam / (TOSCARS::FourPi() * sigma_psi);


    double const Sigma_x = sqrt(epsilon_x * beta_x + pow(eta_x * sigma_E, 2) + sigma_r*sigma_r);
    double const Sigma_y = sqrt(epsilon_y * beta_y + (epsilon_y*epsilon_y + epsilon_y * gamma_y * sigma_r*sigma_r) / (sigma_psi*sigma_psi));

    SpectrumContainer.SetFlux(i, df2dtdp / (TOSCARS::TwoPi() * Sigma_x * Sigma_y) / 1e6);
  }
  return;
}






void OSCARSTH::WigglerBrightnessK (double const K,
                                   double const Period,
                                   int    const NPeriods,
                                   TSpectrumContainer& SpectrumContainer) const
{
  double const BField = OSCARSTH::UndulatorBField(K, Period);
  this->WigglerBrightnessBField(BField, Period, NPeriods, SpectrumContainer);
  return;
}

void OSCARSTH::WigglerBrightnessBField (double const BField,
                                        double const Period,
                                        int    const NPeriods,
                                        TSpectrumContainer& SpectrumContainer) const
{

  // Beam energy from internal beam
  double const BeamEnergy_GeV =  fParticleBeam.GetE0();
  double const sigma_E = fParticleBeam.GetSigmaEnergyGeV() / fParticleBeam.GetE0();
  double const Gamma = fParticleBeam.GetGamma();
  double const Current = fParticleBeam.GetCurrent();

  TVector2D const Beta = fParticleBeam.GetTwissBeta();
  TVector2D const Alpha = fParticleBeam.GetTwissAlpha();
  TVector2D const Emittance = fParticleBeam.GetEmittance();

  double K = OSCARSTH::UndulatorK(BField, Period);

  double const eta_x = fParticleBeam.GetEta().GetX();
  double const beta_x = Beta.GetX();
  double const epsilon_x = Emittance.GetX();
  double const sigma_x  = sqrt(epsilon_x * beta_x + eta_x*eta_x*sigma_E*sigma_E);
  double const sigma_xp = sqrt(epsilon_x / beta_x);

  double const epsilon_y = Emittance.GetY();
  double const beta_y = Beta.GetY();
  double const sigma_y  = sqrt(epsilon_y * beta_y);
  double const sigma_yp = sqrt(epsilon_y / beta_y);

  for (size_t i = 0; i != SpectrumContainer.GetNPoints(); ++i) {
    double const dfdt1 =    2 * NPeriods * this->DipoleSpectrumAngleIntegrated(BField, BeamEnergy_GeV, SpectrumContainer.GetEnergy(i));
    //double const df2dtdp = 2 * NPeriods * this->DipoleSpectrum(BField, BeamEnergy_GeV, 0, SpectrumContainer.GetEnergy(i));
    double const dfdt =    2 * NPeriods * this->WigglerFluxAngleIntegrated(BField, Period, SpectrumContainer.GetEnergy(i));
    double const df2dtdp = 2 * NPeriods * this->WigglerFlux(BField, Period, 0, SpectrumContainer.GetEnergy(i));
    double const x0 = K * Period / (TOSCARS::TwoPi() * Gamma);

    std::cout << "Frac: " << (dfdt - dfdt1) / dfdt << std::endl;

    double const sigma_psi = dfdt / (df2dtdp * sqrt(TOSCARS::TwoPi()));


    double BSum = 0;
    for (int in = 0; in != NPeriods; ++in) {
      double zn = -Period * NPeriods / 2 + Period / 4 + in * Period;
      BSum += (1/TOSCARS::TwoPi()) * exp(-(1/2) * (x0 * x0 / (sigma_x * sigma_x + zn * zn * sigma_xp * sigma_xp)))
          / sqrt((sigma_x * sigma_x + zn * zn * sigma_xp * sigma_xp) * (epsilon_y * epsilon_y / (sigma_psi * sigma_psi) + sigma_y * sigma_y + zn * zn * sigma_yp * sigma_yp));
      zn = -Period * NPeriods / 2 + 3 * Period / 4 + in * Period;
      BSum += (1/TOSCARS::TwoPi()) * exp(-(1/2) * (x0 * x0 / (sigma_x * sigma_x + zn * zn * sigma_xp * sigma_xp)))
          / sqrt((sigma_x * sigma_x + zn * zn * sigma_xp * sigma_xp) * (epsilon_y * epsilon_y / (sigma_psi * sigma_psi) + sigma_y * sigma_y + zn * zn * sigma_yp * sigma_yp));
    }


    SpectrumContainer.SetFlux(i, Current * df2dtdp * BSum / 1e6);
  }
  return;
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
  double const Energy_eV = TOSCARS::FrequencyToEv(TOSCARS::C() / Lambda);

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


TVector2D OSCARSTH::UndulatorCoherentFluxFractionB (double const BField,
                                                    double const Period,
                                                    int    const NPeriods,
                                                    int    const Harmonic
                                                    ) const
{
  // Return the on-axis coherent flux for a planar undulator

  return UndulatorCoherentFluxFractionK(this->UndulatorK(BField, Period), Period, NPeriods, Harmonic);
}



TVector2D OSCARSTH::UndulatorCoherentFluxFractionK (double const K,
                                                    double const Period,
                                                    int    const NPeriods,
                                                    int    const N
                                                    ) const
{
  // Return the on-axis coherent flux fraction for a planar undulator
  //
  // This is the ALS/SLS method for energy spread (lambda/4pi)
  //
  if (N % 2 == 0) {
    return TVector2D(0, 0);
  }

  // Properties from beam
  double    const Gamma          = fParticleBeam.GetGamma();
  TVector2D const Beta           = fParticleBeam.GetTwissBeta();
  TVector2D const Emittance      = fParticleBeam.GetEmittance();
  double    const Current        = fParticleBeam.GetCurrent();
  double    const EnergySpreadFrac = fParticleBeam.GetSigmaEnergyGeV() / fParticleBeam.GetE0();

  // Check that we can do this calculation, else reject
  if (Gamma == 0 || Beta[0] == 0 || Beta[1] == 0 || Emittance[0] == 0 || Emittance[1] == 0 || Current == 0) {
    throw std::invalid_argument("Beam definition incorrect for this calculation: Check energy, current, beta, emittance");
  }

  double const K2 = K * K;
  double const Lambda = Period / (2 * Gamma * Gamma) * (1. + K2 / 2.) / (double) N;

  double const sigx = sqrt(Emittance[0] * Beta[0]);
  double const sigy = sqrt(Emittance[1] * Beta[1]);

  double const sigr = 1. / (4 * TOSCARS::Pi()) * sqrt(2 * Lambda * Period * NPeriods);
  double const sigrp = sqrt(Lambda / (2 * Period * NPeriods));
  double const sigw = 0.36 / (N * NPeriods);
  double const mu_lambda = sqrt(1 + pow(2. * EnergySpreadFrac / sigw, 2));

  //static std::ofstream of("delme.dat");
  //of << Lambda << "   " << N << "   " << mu_lambda << std::endl;


  double const sigxp = sqrt(Emittance[0] / Beta[0]);
  double const sigyp = sqrt(Emittance[1] / Beta[1]);

  double const Energy_eV = UndulatorEnergyAtHarmonicK(K, Period, N);


  double const SigmaX = sqrt(sigx * sigx + sigr * sigr);
  double const SigmaY = sqrt(sigy * sigy + sigr * sigr);
  double const SigmaXP = sqrt(sigxp * sigxp + mu_lambda * sigrp * sigrp);
  double const SigmaYP = sqrt(sigyp * sigyp + mu_lambda * sigrp * sigrp);


  return TVector2D(Energy_eV, Lambda * Lambda / (16 * TOSCARS::Pi2() * SigmaX * SigmaY * SigmaXP * SigmaYP));
}



TVector2D OSCARSTH::UndulatorCoherentFluxB (double const BField,
                                            double const Period,
                                            int    const NPeriods,
                                            int    const Harmonic
                                            ) const
{
  // Return the on-axis coherent flux for a planar undulator

  return UndulatorCoherentFluxK(this->UndulatorK(BField, Period), Period, NPeriods, Harmonic);
}








TVector2D OSCARSTH::UndulatorCoherentFluxK (double const K,
                                            double const Period,
                                            int    const NPeriods,
                                            int    const N
                                            ) const
{
  // Return the on-axis coherent flux for a planar undulator
  //
  // This is the ALS/SLS method for energy spread (lambda/4pi)
  //
  if (N % 2 == 0) {
    return TVector2D(0, 0);
  }

  // Properties from beam
  double    const Gamma          = fParticleBeam.GetGamma();
  TVector2D const Beta           = fParticleBeam.GetTwissBeta();
  TVector2D const Emittance      = fParticleBeam.GetEmittance();
  double    const Current        = fParticleBeam.GetCurrent();
  double    const EnergySpreadFrac = fParticleBeam.GetSigmaEnergyGeV() / fParticleBeam.GetE0();

  // Check that we can do this calculation, else reject
  if (Gamma == 0 || Beta[0] == 0 || Beta[1] == 0 || Emittance[0] == 0 || Emittance[1] == 0 || Current == 0) {
    throw std::invalid_argument("Beam definition incorrect for this calculation: Check energy, current, beta, emittance");
  }

  double const K2 = K * K;
  double const Lambda = Period / (2 * Gamma * Gamma) * (1. + K2 / 2.) / (double) N;

  double const sigx = sqrt(Emittance[0] * Beta[0]);
  double const sigy = sqrt(Emittance[1] * Beta[1]);

  double const sigr = 1. / (4 * TOSCARS::Pi()) * sqrt(2 * Lambda * Period * NPeriods);
  double const sigrp = sqrt(Lambda / (2 * Period * NPeriods));
  double const sigw = 0.36 / (N * NPeriods);
  double const mu_lambda = sqrt(1 + pow(2. * EnergySpreadFrac / sigw, 2));

  double const sigxp = sqrt(Emittance[0] / Beta[0]);
  double const sigyp = sqrt(Emittance[1] / Beta[1]);

  double const Energy_eV = UndulatorEnergyAtHarmonicK(K, Period, N);

  double const Fn = K2 * N * N / pow(1. + K2 / 2., 2) * pow(
      TOMATH::BesselJ( (N - 1) / 2, N * K2 / (4 * (1. + 0.5 * K2))) - TOMATH::BesselJ( (N + 1) / 2, N * K2 / (4 * (1. + 0.5 * K2))),
      2);

  double const Qn = (1. + K2 / 2.) * Fn / (double) N;

  double const Fu = TOSCARS::Pi() * TOSCARS::Alpha() * NPeriods * 0.001 * Current / TOSCARS::Qe() * Qn;

  double const FFnn = 1. * 1.43e14 * NPeriods * Qn * Current;

  double const SigmaX = sqrt(sigx * sigx + sigr * sigr);
  double const SigmaY = sqrt(sigy * sigy + sigr * sigr);
  double const SigmaXP = sqrt(sigxp * sigxp + mu_lambda * sigrp * sigrp);
  double const SigmaYP = sqrt(sigyp * sigyp + mu_lambda * sigrp * sigrp);


  return TVector2D(Energy_eV, Lambda * Lambda * FFnn / (16 * TOSCARS::Pi2() * SigmaX * SigmaY * SigmaXP * SigmaYP));
}






TVector2D OSCARSTH::UndulatorBrightnessK (double const K,
                                          double const Period,
                                          int    const NPeriods,
                                          int    const N
                                          ) const
{
  // Return the on-axis theoretical brightness for a planar undulator
  //
  // This is the ALS/SLS method for energy spread
  //
  if (N % 2 == 0) {
    return TVector2D(0, 0);
  }

  // Properties from beam
  double    const Gamma          = fParticleBeam.GetGamma();
  TVector2D const Beta           = fParticleBeam.GetTwissBeta();
  TVector2D const Emittance      = fParticleBeam.GetEmittance();
  double    const Current        = fParticleBeam.GetCurrent();
  double    const EnergySpreadFrac = fParticleBeam.GetSigmaEnergyGeV() / fParticleBeam.GetE0();

  // Check that we can do this calculation, else reject
  if (Gamma == 0 || Beta[0] == 0 || Beta[1] == 0 || Emittance[0] == 0 || Emittance[1] == 0 || Current == 0) {
    throw std::invalid_argument("Beam definition incorrect for this calculation: Check energy, current, beta, emittance");
  }

  double const K2 = K * K;
  double const Lambda = Period / (2 * Gamma * Gamma) * (1. + K2 / 2.) / (double) N;

  double const sigx = sqrt(Emittance[0] * Beta[0]);
  double const sigy = sqrt(Emittance[1] * Beta[1]);

  double const sigr = 1. / (4 * TOSCARS::Pi()) * sqrt(2 * Lambda * Period * NPeriods);
  double const sigrp = sqrt(Lambda / (2 * Period * NPeriods));
  double const sigw = 0.36 / (N * NPeriods);
  double const mu_lambda = sqrt(1 + pow(2. * EnergySpreadFrac / sigw, 2));

  double const sigxp = sqrt(Emittance[0] / Beta[0]);
  double const sigyp = sqrt(Emittance[1] / Beta[1]);

  double const Energy_eV = UndulatorEnergyAtHarmonicK(K, Period, N);

  double const Fn = K2 * N * N / pow(1. + K2 / 2., 2) * pow(
      TOMATH::BesselJ( (N - 1) / 2, N * K2 / (4 * (1. + 0.5 * K2))) - TOMATH::BesselJ( (N + 1) / 2, N * K2 / (4 * (1. + 0.5 * K2))),
      2);

  double const Qn = (1. + K2 / 2.) * Fn / (double) N;

  double const Fu = TOSCARS::Pi() * TOSCARS::Alpha() * NPeriods * 0.001 * Current / TOSCARS::Qe() * Qn;

  double const FFnn = 1. * 1.43e14 * NPeriods * Qn * Current;

  double const SigmaX = sqrt(sigx * sigx + sigr * sigr);
  double const SigmaY = sqrt(sigy * sigy + sigr * sigr);
  double const SigmaXP = sqrt(sigxp * sigxp + mu_lambda * sigrp * sigrp);
  double const SigmaYP = sqrt(sigyp * sigyp + mu_lambda * sigrp * sigrp);


  return TVector2D(Energy_eV, FFnn / (4 * TOSCARS::Pi2() * SigmaX * SigmaY * SigmaXP * SigmaYP) * 1e-12);
}






//TVector2D OSCARSTH::UndulatorBrightnessK (double const K,
//                                          double const Period,
//                                          int    const NPeriods,
//                                          int    const N
//                                          ) const
//{
//  // Return the on-axis theoretical brightness for a planar undulator
//  if (N % 2 == 0) {
//    return TVector2D(0, 0);
//  }
//
//  // Properties from beam
//  double    const Gamma          = fParticleBeam.GetGamma();
//  TVector2D const Beta           = fParticleBeam.GetTwissBeta();
//  TVector2D const Emittance      = fParticleBeam.GetEmittance();
//  double    const Current        = fParticleBeam.GetCurrent();
//
//  // Check that we can do this calculation, else reject
//  if (Gamma == 0 || Beta[0] == 0 || Beta[1] == 0 || Emittance[0] == 0 || Emittance[1] == 0 || Current == 0) {
//    throw std::invalid_argument("Beam definition incorrect for this calculation: Check energy, current, beta, emittance");
//  }
//
//  double const sigx = sqrt(Emittance[0] * Beta[0]);
//  double const sigy = sqrt(Emittance[1] * Beta[1]);
//
//  double const sigxp = sqrt(Emittance[0] / Beta[0]);
//  double const sigyp = sqrt(Emittance[1] / Beta[1]);
//
//  double const K2 = K * K;
//
//  double const Lambda = Period / (2 * Gamma * Gamma) * (1. + K2 / 2.) / (double) N;
//  double const Energy_eV = UndulatorEnergyAtHarmonicK(K, Period, N);
//
//  double const Fn = K2 * N * N / pow(1. + K2 / 2., 2) * pow(
//      TOMATH::BesselJ( (N - 1) / 2, N * K2 / (4 * (1. + 0.5 * K2))) - TOMATH::BesselJ( (N + 1) / 2, N * K2 / (4 * (1. + 0.5 * K2))),
//      2);
//
//  double const Qn = (1. + K2 / 2.) * Fn / (double) N;
//
//  double const Fu = TOSCARS::Pi() * TOSCARS::Alpha() * NPeriods * 0.001 * Current / TOSCARS::Qe() * Qn;
//
//
//  double const sigr = 1 / TOSCARS::FourPi() * sqrt(Lambda * Period * NPeriods);
//  double const sigrp = sqrt(Lambda / (Period * NPeriods));
//  double const SigmaX = sqrt(sigx * sigx + sigr * sigr);
//  double const SigmaY = sqrt(sigy * sigy + sigr * sigr);
//  double const SigmaXP = sqrt(sigxp * sigxp + sigrp * sigrp);
//  double const SigmaYP = sqrt(sigyp * sigyp + sigrp * sigrp);
//
//
//  return TVector2D(Energy_eV, Fu / (4 * TOSCARS::Pi2() * SigmaX * SigmaY * SigmaXP * SigmaYP) * 1.e-12);
//}






TVector2D OSCARSTH::UndulatorBrightnessB (double const BField,
                                          double const Period,
                                          int    const NPeriods,
                                          int    const Harmonic
                                          ) const
{
  // Return the on-axis theoretical brightness for a planar undulator

  return UndulatorBrightnessK(this->UndulatorK(BField, Period), Period, NPeriods, Harmonic);
}







void OSCARSTH::WigglerFluxK (double         const  K,
                             double         const  Period,
                             int            const  NPeriods,
                             TSurfacePoints const& Surface,
                             double         const  Energy_eV,
                             T3DScalarContainer&   FluxContainer) const
{
  // Extra inputs for calculation
  bool Done = false;
  size_t const iFirst = 0;
  size_t const iLast = Surface.GetNPoints() - 1;

  WigglerFluxKPoints(K,
                     Period,
                     NPeriods,
                     Surface,
                     Energy_eV,
                     FluxContainer,
                     iFirst,
                     iLast,
                     Done);

  return;
}





void OSCARSTH::WigglerFluxK (double         const  K,
                             double         const  Period,
                             int            const  NPeriods,
                             TSurfacePoints const& Surface,
                             double         const  Energy_eV,
                             T3DScalarContainer&   FluxContainer,
                             int            const  NThreads,
                             int            const  GPU) const
{
  // Calculate the wiggler flux
  //
  // THIS is the ENTRY POINT typically (or WigglerFluxB forwards to here)

  // SHOULD be given in input, fake for now
  int const Dimension = 2;

  double const BeamEnergyGeV = fParticleBeam.GetE0();

  // Number of threads to possibly use
  int const NThreadsToUse = NThreads < 1 ? fNThreadsGlobal : NThreads;
  if (NThreadsToUse <= 0) {
    throw std::out_of_range("NThreads or NThreadsGlobal must be >= 1");
  }

  // Should we use the GPU or not?
  bool const UseGPU = GPU == 0 ? false : this->GetUseGPUGlobal() && (this->CheckGPU() > 0) ? true : false;

  if (Dimension == 3) {
    for (size_t i = 0; i != Surface.GetNPoints(); ++i) {
      FluxContainer.AddPoint(Surface.GetPoint(i).GetPoint(), 0);
    }
  } else if (Dimension == 2) {
    for (size_t i = 0; i != Surface.GetNPoints(); ++i) {
      FluxContainer.AddPoint( TVector3D(Surface.GetX1(i), Surface.GetX2(i), 0), 0);
    }
  } else {
    throw std::out_of_range("wROng dimension");
  }

  // GPU will outrank NThreads...
  if (UseGPU == 0) {
    if (NThreadsToUse == 1) {
      this->WigglerFluxK(K, Period, NPeriods, Surface, Energy_eV, FluxContainer);
    } else {
      //this->CalculateFluxThreads(fParticle,
      //                           1);
    }
  } else if (UseGPU == 1) {
    //this->CalculateFluxGPU(fParticle,
  }

  return;
}





void OSCARSTH::WigglerFluxB (double         const  BField,
                             double         const  Period,
                             int            const  NPeriods,
                             TSurfacePoints const& Surface,
                             double         const  Energy_eV,
                             T3DScalarContainer&   FluxContainer,
                             int            const  NThreads,
                             int            const  GPU) const
{
  // Return the on-axis flux for this K value and harmonic

  // Undulator deflection parameter
  double const K = this->UndulatorK(BField, Period);

  this->WigglerFluxK(K, Period, NPeriods, Surface, Energy_eV, FluxContainer, NThreads, GPU);

  return;
}




void OSCARSTH::WigglerFluxKPoints (double         const  K,
                                   double         const  Period,
                                   int            const  NPeriods,
                                   TSurfacePoints const& Surface,
                                   double         const  Energy_eV,
                                   T3DScalarContainer&   FluxContainer,
                                   size_t const iFirst,
                                   size_t const iLast,
                                   bool& Done
                                  ) const
{
  // Calculate flux for the points given in surface/flux container


  double const Gamma = fParticleBeam.GetE0() / TOSCARS::kgToGeV( TOSCARS::Me());
  double const Gamma2 = Gamma*Gamma;
  double const C0 = 3. * TOSCARS::Qe() * TOSCARS::Qe() / (16. * TOSCARS::Pi3() * TOSCARS::Epsilon0() * TOSCARS::C()) * Gamma2;

  double const K2 = K * K;

  // Loop over all points in the surface container
  for (size_t i = iFirst; i <= iLast; ++i) {

    // Obs point
    TVector3D const& ObservationPoint = Surface.GetPoint(i).GetPoint();
    double const X0 = ObservationPoint.GetX();
    double const Y0 = ObservationPoint.GetY();
    double const Z0 = ObservationPoint.GetZ();
    double const ThetaX = atan2(X0, Z0);
    double const ThetaY = atan2(Y0, sqrt(Z0*Z0 + X0*X0));

    double const omega = TOSCARS::EvToAngularFrequency(Energy_eV);
    double const omega1 = TOSCARS::TwoPi() * TOSCARS::C() / Period * 2. * Gamma2 / (1. + K2 / 2. + Gamma2 * pow(ObservationPoint.GetTheta(), 2));

    double const alpha = Gamma * ThetaX / K;
    double const A = 1. + K2 / 2. + Gamma2 *(ThetaX*ThetaX + ThetaY*ThetaY);
    double const Delta = omega / omega1 * (TOSCARS::Pi() + 2. * asin(alpha) + 3. * K2 / A * alpha * sqrt(1. - alpha*alpha));

    double const omega_c0 = TOSCARS::TwoPi() * TOSCARS::C() * 2. * Gamma2 / Period;
    double const omega_c = omega_c0 * sqrt(1. - alpha*alpha);
    double const y = omega / omega_c;
    double const X = Gamma * ThetaY;
    double const Xi = 0.5 * y * pow(1. + X*X, 1.5);

    double const SinFactor = pow(sin(NPeriods * TOSCARS::Pi() * omega / omega1), 2) / pow(sin(TOSCARS::Pi() * omega / omega1), 2);

    double const AxMag2 = 4. * pow(sin(Delta/2.), 2) * pow(TOMATH::BesselK(2./3., Xi), 2) * SinFactor;
    double const AyMag2 = 4. * pow(cos(Delta/2.), 2) * X*X/(1. + X*X) * pow(TOMATH::BesselK(1./3., Xi), 2) * SinFactor;

    // Flux in terms of intensity per rad^2 ==> nphotons / mrad^2 / s / 0.1%bw
    //double const ThisFlux = C0 * y*y * pow(1. + X*X, 2) * AxMag2;
    double const ThisFlux = C0 * y*y * pow(1. + X*X, 2) * (AxMag2 + AyMag2);

    // All point to flux container
    FluxContainer.AddToPoint(i, ThisFlux);
  }

  // Set done to true before returnning
  Done = true;

  return;
}



/*
TParticleBeam& OSCARSTH::SetParticleBeam (std::string const& Beam, std::string const& Name)
{
  // Add a particle beam
  // Beam - The name of the predefined particle beam to add

  fParticleBeam.SetPredefinedBeam(Beam);
  fParticleBeam.SetName(Name);
  fParticleBeam.SetWeight(1);
  return fParticleBeam;
}




TParticleBeam& OSCARSTH::SetParticleBeam (double const Energy_GeV,
                                double const Current,
                                TVector2D const& Beta,
                                TVector2D const& Emittance,
                                double const SigmaEnergyGeV,
                                TVector2D const& Eta,
                                std::string const& Name)
{
  // Add a particle beam
  // Energy_GeV  - Energy of particle beam in GeV
  // Current     - Beam current in Amperes

  fParticleBeam = TParticleBeam("electron", Name, Energy_GeV, Current);
  fParticleBeam.SetWeight(1);
  fParticleBeam.SetEmittance(Emittance);
  fParticleBeam.SetTwissBetaAlpha(Beta, TVector2D(0, 0));
  fParticleBeam.SetSigmaEnergyGeV(SigmaEnergyGeV);
  fParticleBeam.SetEta(Eta);

  return fParticleBeam;
}
*/



TParticleBeam& OSCARSTH::AddParticleBeam (std::string const& Type,
                                          std::string const& Name,
                                          TVector3D const& X0,
                                          TVector3D const& V0,
                                          double const Energy_GeV,
                                          double const T0,
                                          double const Current,
                                          double const Weight,
                                          double const Charge,
                                          double const Mass)
{
  // Add a particle beam
  // Type        - The name of the particle type that you want to use
  // Name        - A user specified 'name' for this beam
  // X0          - Initial position in X,Y,Z
  // V0          - A vector pointing in the direction of the velocity of arbitrary magnitude
  // Energy_GeV  - Energy of particle beam in GeV
  // T0          - Time of initial conditions, specified in units of [m] (for v = c)
  // Current     - Beam current in Amperes
  // Weight      - Relative weight to give this beam when randomly sampling
  // Charge      - Charge of custom particle
  // Mass        - Mass of custom particle

  this->ClearParticleBeams();
  fParticleBeam = fParticleBeamContainer.AddNewParticleBeam(Type, Name, X0, V0, Energy_GeV, T0, Current, Weight, Charge, Mass);
  return fParticleBeam;
}




TParticleBeam& OSCARSTH::AddParticleBeam (std::string const& Beam,
                                          std::string const& Name,
                                          double const Weight)
{
  // Add a particle beam
  // Beam - The name of the predefined particle beam to add

  this->ClearParticleBeams();
  fParticleBeam = fParticleBeamContainer.AddNewParticleBeam(Beam, Name, Weight);
  return fParticleBeam;
}




void OSCARSTH::ClearParticleBeams ()
{
  // Clear the contents of the particle beam container
  fParticleBeamContainer.Clear();

  return;
}




TParticleBeam& OSCARSTH::GetParticleBeam ()
{
  // Return a reference to the particle beam by a given name
  return fParticleBeam;
}



bool OSCARSTH::CheckBeam () const
{
  // Check the beam exists and meets some minimum criteria
  if (fParticleBeam.GetE0() <= 0) {
    return false;
  }

  return true;
}




bool OSCARSTH::SetUseGPUGlobal (int const in)
{
  // Must be a 0 or a 1 at the moment.  Will return false if you tried to set the GPU but it's not available

  if (in == 0) {
    fUseGPUGlobal = 0;
    return true;
  }

  if (in != 1) {
    fUseGPUGlobal = 0;
    return false;
  }

  #ifdef CUDA
  if (OSCARSTH_Cuda_GetDeviceCount() > 0) {
    fUseGPUGlobal = 1;
    return true;
  } else {
    fUseGPUGlobal = 0;
    return false;
  }
  #endif

  fUseGPUGlobal = 0;

  return false;
}




int OSCARSTH::GetUseGPUGlobal () const
{
  return fUseGPUGlobal;
}




int OSCARSTH::CheckGPU () const
{
  #ifdef CUDA
    return OSCARSTH_Cuda_GetDeviceCount();
  #endif
  return -1;
}




void OSCARSTH::SetNThreadsGlobal (int const N)
{
  fNThreadsGlobal = N;
  return;
}




