#ifndef GUARD_TBField3D_Gaussian
#define GUARD_TBField3D_Gaussian
////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Wed Jun 29 08:52:22 EDT 2016
//
////////////////////////////////////////////////////////////////////

#include "TBField.h"


class TBField3D_Gaussian : public TBField
{
  public:
    TBField3D_Gaussian ();
    TBField3D_Gaussian (TVector3D const&, TVector3D const&, TVector3D const&, TVector3D const& Rotations = TVector3D(0, 0, 0));
    ~TBField3D_Gaussian ();

    double    GetBx (double const, double const, double const) const;
    double    GetBy (double const, double const, double const) const;
    double    GetBz (double const, double const, double const) const;
    TVector3D GetB  (double const, double const, double const) const;
    TVector3D GetB  (TVector3D const&) const;

    bool IsWithinRange (double const, double const, double const) const;

  private:
    TVector3D fBField;
    TVector3D fCenter;
    TVector3D fSigma;
    TVector3D fRotated;

    bool fIgnoreAxisX;
    bool fIgnoreAxisY;
    bool fIgnoreAxisZ;

};




















#endif
