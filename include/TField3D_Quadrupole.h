#ifndef GUARD_TField3D_Quadrupole_h
#define GUARD_TField3D_Quadrupole_h
////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Fri May 19 14:41:31 CEST 2017
//
////////////////////////////////////////////////////////////////////


#include "TField.h"

class TField3D_Quadrupole : public TField
{
  public:
    TField3D_Quadrupole (std::string const& Name = "");

    TField3D_Quadrupole (double const K,
                         double const Width,
                         TVector3D const& Rotations,
                         TVector3D const& Translation,
                         double    const  Frequency = 0,
                         double    const  FrequencyPhase = 0,
                         double    const  TimeOffset = 0,
                         std::string const& Name = "");

    ~TField3D_Quadrupole ();


    TVector3D GetF  (double const, double const, double const, double const T = 0) const;
    TVector3D GetF  (TVector3D const&, double const T = 0) const;

    void      Print (std::ostream&) const;

    double GetK () const;
    double GetWidth () const;
    TVector3D const& GetRotations () const;
    TVector3D const& GetTranslation () const;

    double GetFrequency () const;
    double GetFrequencyPhase () const;
    double GetTimeOffset () const;


  private:
    double fK;
    double fWidth;
    TVector3D fRotations;
    TVector3D fTranslation;

    double fFrequency;
    double fFrequencyPhase;
    double fTimeOffset;
};



inline std::ostream& operator << (std::ostream& os, TField3D_Quadrupole const& o)
{
  // For easy printing
  os << "TField3D_Quadrupole " << "\n"
     << "Name                " << o.GetName() << "\n"
     << "K                   " << o.GetK() << "\n"
     << "Width               " << o.GetWidth() << "\n"
     << "Rotations           " << o.GetRotations() << "\n"
     << "Translation         " << o.GetTranslation() << "\n"
     << "Frequency           " << o.GetFrequency() << "\n"
     << "FrequencyPhase      " << o.GetFrequencyPhase() << "\n"
     << "TimeOffset          " << o.GetTimeOffset() << "\n";


  return os;
}





#endif
