#ifndef GUARD_TField3D_Gaussian
#define GUARD_TField3D_Gaussian
////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Wed Jun 29 08:52:22 EDT 2016
//
////////////////////////////////////////////////////////////////////

#include "TField.h"


class TField3D_Gaussian : public TField
{
  public:
    TField3D_Gaussian ();
    TField3D_Gaussian (TVector3D const&, TVector3D const&, TVector3D const&, TVector3D const& Rotations = TVector3D(0, 0, 0));
    ~TField3D_Gaussian ();

    double    GetFx (double const, double const, double const) const;
    double    GetFy (double const, double const, double const) const;
    double    GetFz (double const, double const, double const) const;
    TVector3D GetF  (double const, double const, double const) const;
    TVector3D GetF  (TVector3D const&) const;

    bool IsWithinRange (double const, double const, double const) const;

  private:
    TVector3D fField;
    TVector3D fCenter;
    TVector3D fSigma;
    TVector3D fRotated;

    bool fIgnoreAxisX;
    bool fIgnoreAxisY;
    bool fIgnoreAxisZ;

};




















#endif
