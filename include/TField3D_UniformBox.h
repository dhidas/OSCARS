#ifndef GUARD_TField3D_UniformBox_h
#define GUARD_TField3D_UniformBox_h
////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Thu Jun 30 08:09:53 EDT 2016
//
// UPDATE: Comments
//
////////////////////////////////////////////////////////////////////

#include "TField.h"
#include "TVector3D.h"

class TField3D_UniformBox : public TField
{
  public:
    TField3D_UniformBox (double const, double const, double const);
    TField3D_UniformBox (TVector3D const&, TVector3D const& Width = TVector3D(0, 0, 0), TVector3D const& Center = TVector3D(0, 0, 0), TVector3D const& Rotations = TVector3D(0, 0, 0));
    ~TField3D_UniformBox ();

    double    GetFx (double const, double const, double const) const;
    double    GetFy (double const, double const, double const) const;
    double    GetFz (double const, double const, double const) const;
    TVector3D GetF  (double const, double const, double const) const;
    TVector3D GetF  (TVector3D const&) const;


  private:
    TVector3D fField;
    TVector3D fWidth;
    TVector3D fRotated;
    TVector3D fCenter;

    bool fIgnoreAxisX;
    bool fIgnoreAxisY;
    bool fIgnoreAxisZ;
};






#endif

