#ifndef GUARD_TField3D_IdealUndulator_h
#define GUARD_TField3D_IdealUndulator_h
////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Tue Jul 19 08:09:18 EDT 2016
//
// An idealized undulator field
//
////////////////////////////////////////////////////////////////////

#include "TField.h"



class TField3D_IdealUndulator : public TField
{
  public:
    TField3D_IdealUndulator ();
    TField3D_IdealUndulator (TVector3D const&, TVector3D const&, int const, TVector3D const& Center = TVector3D(0, 0, 0), double const Phase = 0, double const Taper = 0);
    ~TField3D_IdealUndulator ();

    double    GetFx (double const, double const, double const) const;
    double    GetFy (double const, double const, double const) const;
    double    GetFz (double const, double const, double const) const;
    TVector3D GetF  (double const, double const, double const) const;
    TVector3D GetF  (TVector3D const&) const;

    void Init (TVector3D const&, TVector3D const&, int const, TVector3D const& Center = TVector3D(0, 0, 0), double const Phase = 0, double const Taper = 0);



  private:
    TVector3D fField;
    TVector3D fPeriod;
    TVector3D fPeriodUnitVector;
    double    fPeriodLength;
    int       fNPeriods;
    TVector3D fCenter;
    double    fPhase;
    double    fTaper;

    double fUndulatorLength;

};



#endif
