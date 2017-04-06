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


class OSCARSTH
{
  public:
    OSCARSTH ();
    ~OSCARSTH ();

    double UndulatorK (double const BFieldMax,
                       double const Period) const;

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


    TVector2D UndulatorFluxKHarmonic (double const K,
                                      double const Period,
                                      double const NPeriods,
                                      double const BeamEnergy,
                                      int    const Harmonic) const;

    double UndulatorFluxOnAxis (double const BField,
                                double const Period,
                                double const NPeriods,
                                double const BeamEnergy,
                                double const Energy_eV,
                                int const FirstHarmonic,
                                int const LastHarmonic) const;

    double UndulatorFluxWeak (double const K,
                              double const Period,
                              double const NPeriods,
                              double const BeamEnergy,
                              int const Harmonic) const;

    double DipoleBrightness () const;

    TVector2D UndulatorBrightness (double const BField,
                                   double const Period,
                                   int    const NPeriods,
                                   int    const N,
                                   double const BeamEnergy_GeV,
                                   double const SigmaE,
                                   double const Current,
                                   TVector2D const BetaX,
                                   TVector2D const Emittance) const;

};














#endif
