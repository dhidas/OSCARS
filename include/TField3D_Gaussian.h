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
    TField3D_Gaussian (TVector3D const& PeakField,
                       TVector3D const& Center,
                       TVector3D const& Sigma,
                       TVector3D const& Rotations);
    ~TField3D_Gaussian ();

    double    GetFx (double const X, double const Y, double const Z) const;
    double    GetFy (double const X, double const Y, double const Z) const;
    double    GetFz (double const X, double const Y, double const Z) const;
    TVector3D GetF  (double const X, double const Y, double const Z) const;
    TVector3D GetF  (TVector3D const& X) const;

    bool IsWithinRange (double const X, double const Y, double const Z) const;

    TVector3D const& GetPeakField () const;
    TVector3D const& GetCenter () const;
    TVector3D const& GetSigma () const;
    TVector3D const& GetRotations () const;

    void Print (std::ostream& os) const;



  private:
    TVector3D fPeakField;
    TVector3D fCenter;
    TVector3D fSigma;
    TVector3D fRotated;

    bool fIgnoreAxisX;
    bool fIgnoreAxisY;
    bool fIgnoreAxisZ;

};


inline std::ostream& operator << (std::ostream& os, TField3D_Gaussian const& o)
{
  // For easy printing
  os << "Gaussian  " << "\n"
     << "Peak      " << o.GetPeakField() << "\n"
     << "Center    " << o.GetCenter() << "\n"
     << "Sigma     " << o.GetSigma() << "\n"
     << "Rotations " << o.GetRotations() << "\n";

  return os;
}























#endif
