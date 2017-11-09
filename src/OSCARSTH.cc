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

#include "TOSCARSSR.h"
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



double OSCARSTH::DipoleCriticalEnergy (double const BField) const
{

  // Return the critical energy in eV for dipole and electron beam
  double const BeamEnergy_GeV = fParticleBeam.GetE0();
  double const Gamma = BeamEnergy_GeV / TOSCARSSR::kgToGeV( TOSCARSSR::Me());
  double const OmegaC = 3. * Gamma * Gamma * Gamma * TOSCARSSR::C() / (2. * BeamEnergy_GeV * 1e9 * TOSCARSSR::Qe() / (TOSCARSSR::Qe() * TOSCARSSR::C() * fabs(BField)));

  return TOSCARSSR::AngularFrequencyToEv(OmegaC);
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

  double const R = BeamEnergy_GeV  * 1e9 / (BField * TOSCARSSR::C());
  double const Q = TOSCARSSR::Qe();
  double const Me = TOSCARSSR::Me();
  long double const v5 = Me * Me;
  double const v6 = (BeamEnergy_GeV * 1e9) * (BeamEnergy_GeV * 1e9);
  double const v7 = v5 / v6;
  //long double const v1 = ((TOSCARSSR::Me() * TOSCARSSR::Me()) * (TOSCARSSR::C() * TOSCARSSR::C() * TOSCARSSR::C() * TOSCARSSR::C())) / ((BeamEnergy_GeV * TOSCARSSR::Qe() * (1e-9)) * (BeamEnergy_GeV * TOSCARSSR::Qe() * (1e-9)));
  double const v2 = 1 - v7;
  double const v3 = sqrt(v2);
  double const v = TOSCARSSR::C() * v3;
  long double const gamma = BeamEnergy_GeV / TOSCARSSR::kgToGeV( TOSCARSSR::Me());
  //long double const beta_sqr = 1 - (1 / (gamma * gamma));
  //long double const beta = sqrt(beta_sqr);
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
  //double const d2I1 = (3./4.) * ((Q*Q*gamma*gamma) / (4. * (pi*pi*pi) * Epsilon0 * c));
  double const d2I2 = myK2 + (((gamma*gamma*psi*psi) / (1 + (gamma*gamma*psi*psi))) * myK1);
  double const d2I3 = ((w / wc4) * (w / wc4)) * ((1 + (gamma*gamma)*(psi*psi))*(1 + (gamma*gamma)*(psi*psi)));
  //double const d2I = d2I1 * d2I2 * d2I3;
  double const Hbar = TOSCARSSR::Hbar();
  double const alpha = (Q*Q) / (4. * (pi*pi*pi) * Epsilon0 * c * Hbar);
  double const I = 0.5;
  double const d = (3./4.) * (alpha) * (gamma*gamma) * (I / Q) * (0.001) * d2I3 * d2I2;
  double const d2N = d * (1e-6);

  return d2N;
}





double OSCARSTH::DipoleSpectrumAngleIntegrated (double const BField, double const BeamEnergy_GeV, double const Energy_eV) const
{
  // This function calculates the dipole spectrum integrated over the vertical angle and returns the result as:
  // [0.1%bw / mrad]

  double const Gamma = BeamEnergy_GeV / TOSCARSSR::kgToGeV( TOSCARSSR::Me());
  double const OmegaC = 3. * Gamma * Gamma * Gamma * TOSCARSSR::C() / (2. * BeamEnergy_GeV * 1e9 * TOSCARSSR::Qe() / (TOSCARSSR::Qe() * TOSCARSSR::C() * fabs(BField)));
  double const Omega = TOSCARSSR::EvToAngularFrequency(Energy_eV);
  double const y = Omega / OmegaC;
  return sqrt(3.0) / TOSCARSSR::TwoPi() * TOSCARSSR::Alpha() * Gamma * 0.001 * y * fParticleBeam.GetCurrent() / TOSCARSSR::Qe() * TOMATH::BesselK_IntegralToInfty(5./3., y) * 0.001;
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

  if (Harmonic % 2 != 1) {
    return TVector2D(0, 0);
  }

  double const BeamEnergyGeV = fParticleBeam.GetE0();

  double const n = Harmonic;
  double const nu = (n + 1.) / 2.;
  double const nu2 = (n - 1.) / 2.;
  double const JB = (n * K * K) / (4. + (2. * K * K));
  //double const E = BeamEnergyGeV * 1e9;
  double const gamma = BeamEnergyGeV / TOSCARSSR::kgToGeV( TOSCARSSR::Me());
  double const JBessel1 = TOMATH::BesselJ(nu, JB);
  double const JBessel3 = TOMATH::BesselJ(nu2, JB);
  double const Q = TOSCARSSR::Qe();
  double const c = TOSCARSSR::C();
  //double const Me_eV = TOSCARSSR::kgToGeV(TOSCARSSR::Me()) * 1e9;
  double const pi = TOSCARSSR::Pi();
  //double const K_1 = Q * 1. * Period;
  //double const K_2 = 2. * pi * Me_eV * c * c;
  //double const K_0 = K_1 / K_2;
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

  if (Harmonic % 2 != 1) {
    return TVector2D(omega_1, 0);
  }

  return TVector2D(omega_1, h00);
}





TVector2D OSCARSTH::UndulatorFluxOnAxisB (double const BField,
                                          double const Period,
                                          double const NPeriods,
                                          int const Harmonic) const
{
  // Return the on-axis flux for this K value and harmonic

  // Undulator deflection parameter
  double const K = this->UndulatorK(BField, Period);

  return this->UndulatorFluxOnAxisK(K, Period, NPeriods, Harmonic);
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
  // Return the on-axis theoretical brightness for a planar undulator
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

  double const Fn = K2 * N * N / pow(1. + K2 / 2., 2) * pow(
      TOMATH::BesselJ( (N - 1) / 2, N * K2 / (4 * (1. + 0.5 * K2))) - TOMATH::BesselJ( (N + 1) / 2, N * K2 / (4 * (1. + 0.5 * K2))),
      2);

  double const Qn = (1. + K2 / 2.) * Fn / (double) N;

  double const Fu = TOSCARSSR::Pi() * TOSCARSSR::Alpha() * NPeriods * 0.001 * Current / TOSCARSSR::Qe() * Qn;


  double const sigr = 1 / TOSCARSSR::FourPi() * sqrt(Lambda * Period * NPeriods);
  double const sigrp = sqrt(Lambda / (Period * NPeriods));
  double const SigmaXP = sqrt(sigxp * sigxp + sigrp * sigrp);
  double const SigmaYP = sqrt(sigyp * sigyp + sigrp * sigrp);


  return TVector2D(Energy_eV, Fu / (TOSCARSSR::Pi2() * SigmaXP * SigmaYP) * 1.e-6);
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





void OSCARSTH::DipoleBrightness (double const BField,
                                 TSpectrumContainer& SpectrumContainer) const
{

  // Beam energy from internal beam
  double const BeamEnergy_GeV =  fParticleBeam.GetE0();
  double const sigma_E = fParticleBeam.GetSigmaEnergyGeV() / fParticleBeam.GetE0();

  TVector2D const Beta = fParticleBeam.GetTwissBeta();
  TVector2D const Alpha = fParticleBeam.GetTwissAlpha();
  TVector2D const Emittance = fParticleBeam.GetEmittance();


  double lam = fabs(TOSCARSSR::H() * TOSCARSSR::C() / (1000. * fParticleBeam.GetCharge()));
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

    double const sigma_psi = dfdt / (df2dtdp * sqrt(TOSCARSSR::TwoPi()));


    double const sigma_r = lam / (TOSCARSSR::FourPi() * sigma_psi);


    double const Sigma_x = sqrt(epsilon_x * beta_x + pow(eta_x * sigma_E, 2) + sigma_r*sigma_r);
    double const Sigma_y = sqrt(epsilon_y * beta_y + (epsilon_y*epsilon_y + epsilon_y * gamma_y * sigma_r*sigma_r) / (sigma_psi*sigma_psi));

    SpectrumContainer.SetFlux(i, df2dtdp / (TOSCARSSR::TwoPi() * Sigma_x * Sigma_y) / 1e6);
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


  double const Gamma = fParticleBeam.GetE0() / TOSCARSSR::kgToGeV( TOSCARSSR::Me());
  double const Gamma2 = Gamma*Gamma;
  double const C0 = 3. * TOSCARSSR::Qe() * TOSCARSSR::Qe() / (16. * TOSCARSSR::Pi3() * TOSCARSSR::Epsilon0() * TOSCARSSR::C()) * Gamma2;

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

    double const omega = TOSCARSSR::EvToAngularFrequency(Energy_eV);
    double const omega1 = TOSCARSSR::TwoPi() * TOSCARSSR::C() / Period * 2. * Gamma2 / (1. + K2 / 2. + Gamma2 * pow(ObservationPoint.GetTheta(), 2));

    double const alpha = Gamma * ThetaX / K;
    double const A = 1. + K2 / 2. + Gamma2 *(ThetaX*ThetaX + ThetaY*ThetaY);
    double const Delta = omega / omega1 * (TOSCARSSR::Pi() + 2. * asin(alpha) + 3. * K2 / A * alpha * sqrt(1. - alpha*alpha));

    double const omega_c0 = TOSCARSSR::TwoPi() * TOSCARSSR::C() * 2. * Gamma2 / Period;
    double const omega_c = omega_c0 * sqrt(1. - alpha*alpha);
    double const y = omega / omega_c;
    double const X = Gamma * ThetaY;
    double const Xi = 0.5 * y * pow(1. + X*X, 1.5);

    double const SinFactor = pow(sin(NPeriods * TOSCARSSR::Pi() * omega / omega1), 2) / pow(sin(TOSCARSSR::Pi() * omega / omega1), 2);

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




