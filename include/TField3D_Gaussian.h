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
    TField3D_Gaussian (std::string const Name = "");

    TField3D_Gaussian (TVector3D const& PeakField,
                       TVector3D const& Center,
                       TVector3D const& Sigma,
                       TVector3D const& Rotations,
                       double    const  Frequency = 0,
                       double    const  FrequencyPhase = 0,
                       double    const  TimeOffset = 0,
                       std::string const& Name = "");

    ~TField3D_Gaussian ();

    TVector3D GetF  (double const X, double const Y, double const Z, double const T = 0) const;
    TVector3D GetF  (TVector3D const& X, double const T = 0) const;

    bool IsWithinRange (double const X, double const Y, double const Z) const;

    TVector3D const& GetPeakField () const;
    TVector3D const& GetCenter () const;
    TVector3D const& GetSigma () const;
    TVector3D const& GetRotations () const;

    double GetFrequency () const;
    double GetFrequencyPhase () const;
    double GetTimeOffset () const;

    void Print (std::ostream& os) const;



  private:
    TVector3D fPeakField;
    TVector3D fCenter;
    TVector3D fSigma;
    TVector3D fRotated;

    double fFrequency;
    double fFrequencyPhase;
    double fTimeOffset;

    bool fIgnoreAxisX;
    bool fIgnoreAxisY;
    bool fIgnoreAxisZ;

};


inline std::ostream& operator << (std::ostream& os, TField3D_Gaussian const& o)
{
  // For easy printing
  os << "Gaussian       " << "\n"
     << "Name           " << o.GetName() << "\n"
     << "Peak           " << o.GetPeakField() << "\n"
     << "Center         " << o.GetCenter() << "\n"
     << "Sigma          " << o.GetSigma() << "\n"
     << "Rotations      " << o.GetRotations() << "\n"
     << "Frequency      " << o.GetFrequency() << "\n"
     << "FrequencyPhase " << o.GetFrequencyPhase() << "\n"
     << "TimeOffset     " << o.GetTimeOffset() << "\n";

  return os;
}























#endif
