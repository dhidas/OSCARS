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
    TField3D_UniformBox (double      const  Fx,
                         double      const  Fy,
                         double      const  Fz,
                         double      const  Frequency = 0,
                         double      const  FrequencyPhase = 0,
                         double      const  TimeOffset = 0,
                         std::string const& Name = "");

    TField3D_UniformBox (TVector3D   const& Field,
                         TVector3D   const& Width = TVector3D(0, 0, 0),
                         TVector3D   const& Center = TVector3D(0, 0, 0),
                         TVector3D   const& Rotations = TVector3D(0, 0, 0),
                         double      const  Frequency = 0,
                         double      const  FrequencyPhase = 0,
                         double      const  TimeOffset = 0,
                         std::string const& Name = "");

    ~TField3D_UniformBox ();

    TVector3D GetF  (double const X, double const Y, double const Z, double const T = 0) const;
    TVector3D GetF  (TVector3D const& X, double const T = 0) const;

    TVector3D GetField () const;
    TVector3D GetWidth () const;
    TVector3D GetRotated () const;
    TVector3D GetCenter () const;

    double    GetFrequency () const;
    double    GetFrequencyPhase () const;
    double    GetTimeOffset () const;

    void Print (std::ostream& os) const;

  private:
    TVector3D fField;
    TVector3D fWidth;
    TVector3D fRotated;
    TVector3D fCenter;
    double    fFrequency;
    double    fFrequencyPhase;
    double    fTimeOffset;


    bool fIgnoreAxisX;
    bool fIgnoreAxisY;
    bool fIgnoreAxisZ;
};





inline std::ostream& operator << (std::ostream& os, TField3D_UniformBox const& o)
{
  // For easy printing
  os << "TField3D_UniformBox " << "\n"
     << "Name                " << o.GetName() << "\n"
     << "Field               " << o.GetField() << "\n"
     << "Width               " << o.GetWidth() << "\n"
     << "Rotations           " << o.GetRotated() << "\n"
     << "Center              " << o.GetCenter() << "\n"
     << "Frequency           " << o.GetFrequency() << "\n"
     << "FrequencyPhase      " << o.GetFrequency() << "\n"
     << "TimeOffset          " << o.GetTimeOffset() << "\n";

  return os;
}





#endif

