#ifndef GUARD_TBField3D_IdealUndulator_h
#define GUARD_TBField3D_IdealUndulator_h
////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Tue Jul 19 08:09:18 EDT 2016
//
// An idealized undulator field
//
////////////////////////////////////////////////////////////////////

#include "TBField.h"



class TBField3D_IdealUndulator : public TBField
{
  public:
    TBField3D_IdealUndulator ();
    TBField3D_IdealUndulator (TVector3D const&, TVector3D const&, int const, TVector3D const& Center = TVector3D(0, 0, 0), double const Phase = 0);
    ~TBField3D_IdealUndulator ();

    double    GetBx (double const, double const, double const) const;
    double    GetBy (double const, double const, double const) const;
    double    GetBz (double const, double const, double const) const;
    TVector3D GetB  (double const, double const, double const) const;
    TVector3D GetB  (TVector3D const&) const;

    void Init (TVector3D const&, TVector3D const&, int const, TVector3D const& Center = TVector3D(0, 0, 0), double const Phase = 0);



  private:
    TVector3D fBField;
    TVector3D fPeriod;
    TVector3D fPeriodUnitVector;
    double    fPeriodLength;
    int       fNPeriods;
    TVector3D fCenter;
    double    fPhase;

    double fUndulatorLength;

};



#endif
