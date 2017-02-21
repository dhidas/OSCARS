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


class OSCARSTH
{
  public:
    OSCARSTH ();
    ~OSCARSTH ();

    double UndulatorK (double const BFieldMax, double const Period) const;
    double DipoleSpectrum (double const BField, double const BeamEnergy_GeV, double const Angle, TVector2D const EnergyRange_eV) const;
};














#endif
