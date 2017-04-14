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
    TField3D_UniformBox (double const Fx, double const Fy, double const Fz);
    TField3D_UniformBox (TVector3D const& Field,
                         TVector3D const& Width = TVector3D(0, 0, 0),
                         TVector3D const& Center = TVector3D(0, 0, 0),
                         TVector3D const& Rotations = TVector3D(0, 0, 0));

    ~TField3D_UniformBox ();

    double    GetFx (double const X, double const Y, double const Z) const;
    double    GetFy (double const X, double const Y, double const Z) const;
    double    GetFz (double const X, double const Y, double const Z) const;
    TVector3D GetF  (double const X, double const Y, double const Z) const;
    TVector3D GetF  (TVector3D const& X) const;

    TVector3D GetField () const;
    TVector3D GetWidth () const;
    TVector3D GetRotated () const;
    TVector3D GetCenter () const;

    void Print (std::ostream& os) const;

  private:
    TVector3D fField;
    TVector3D fWidth;
    TVector3D fRotated;
    TVector3D fCenter;

    bool fIgnoreAxisX;
    bool fIgnoreAxisY;
    bool fIgnoreAxisZ;
};





inline std::ostream& operator << (std::ostream& os, TField3D_UniformBox const& o)
{
  // For easy printing
  os << "TField3D_UniformBox " << "\n"
     << "Field               " << o.GetField() << "\n"
     << "Width               " << o.GetWidth() << "\n"
     << "Rotations           " << o.GetRotated() << "\n"
     << "Center              " << o.GetCenter() << "\n";

  return os;
}





#endif

