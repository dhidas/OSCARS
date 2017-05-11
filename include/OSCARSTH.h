#ifndef GUARD_OSCARSTH_h
#define GUARD_OSCARSTH_h
////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Thu Feb  2 11:55:24 EST 2017 (minus a day or so)
//
////////////////////////////////////////////////////////////////////

#include "TVector2D.h"
#include "TOMATH.h"
#include "TParticleBeamContainer.h"


class OSCARSTH
{
  public:
    OSCARSTH ();
    ~OSCARSTH ();

    double UndulatorK (double const BFieldMax,
                       double const Period) const;

    double UndulatorBField (double const K,
                            double const Period) const;

    double UndulatorPeriod (double const BField,
                            double const K) const;

    double DipoleCriticalEnergy (double const BField,
                                 double const BeamEnergy_GeV) const;

    double DipoleSpectrum (double const BField,
                           double const BeamEnergy_GeV,
                           double const Angle,
                           double const Energy_eV) const;

    double UndulatorFlux (double const BField,
                          double const Period,
                          double const NPeriods,
                          double const BeamEnergy,
                          double const AngleV,
                          double const AngleH,
                          double const Energy_eV) const;


    TVector2D UndulatorFluxOnAxisK (double const K,
                                    double const Period,
                                    double const NPeriods,
                                    int    const Harmonic) const;

    TVector2D UndulatorFluxOnAxisB (double const BField,
                                    double const Period,
                                    double const NPeriods,
                                    int    const Harmonic) const;

    double UndulatorFluxWeak (double const K,
                              double const Period,
                              double const NPeriods,
                              double const BeamEnergy,
                              int const Harmonic) const;

    double DipoleBrightness () const;

    double UndulatorEnergyAtHarmonicK (double const K,
                                       double const Period,
                                       int    const Harmonic) const;

    double UndulatorEnergyAtHarmonicB (double const BField,
                                       double const Period,
                                       int    const Harmonic) const;

    TVector2D UndulatorBrightnessK (double const BField,
                                    double const Period,
                                    int    const NPeriods,
                                    int    const N) const;

    TVector2D UndulatorBrightnessB (double const BField,
                                    double const Period,
                                    int    const NPeriods,
                                    int    const N) const;



    // Functions related to the particle beam
    void SetParticleBeam (std::string const& Beam);
    void SetParticleBeam (double const Energy_GeV,
                          double const Current,
                          TVector2D const& Beta = TVector2D(0, 0),
                          TVector2D const& Emittance = TVector2D(0, 0),
                          double const SigmaEnergyGeV = 0);
    TParticleBeam& GetParticleBeam ();

  private:
    TParticleBeam fParticleBeam;

};














#endif
