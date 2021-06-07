#ifndef GUARD_TField3D_IdealEPU_h
#define GUARD_TField3D_IdealEPU_h
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



class TField3D_IdealEPU : public TField
{
  public:
    TField3D_IdealEPU (std::string const& Name);

    TField3D_IdealEPU (TVector3D   const& FieldA,
                       TVector3D   const& FieldB,
                       TVector3D   const& Period,
                       int         const  NPeriods,
                       TVector3D   const& Center    = TVector3D(0, 0, 0),
                       double      const  Phase     = 0,
                       double      const  Taper     = 0,
                       double      const  Frequency = 0,
                       double      const  FrequencyPhase = 0,
                       double      const  TimeOffset = 0,
                       std::string const& Name      = "");

    ~TField3D_IdealEPU ();

    TVector3D GetF  (double const X, double const Y, double const Z, double const T = 0) const;
    TVector3D GetF  (TVector3D const& X, double const T = 0) const;

    void Init (TVector3D   const& FieldA,
               TVector3D   const& FieldB,
               TVector3D   const& Period,
               int         const  NPeriods,
               TVector3D   const& Center    = TVector3D(0, 0, 0),
               double      const  Phase     = 0,
               double      const  Taper     = 0,
               double      const  Frequency = 0,
               double      const  FrequencyPhase = 0,
               double      const  TimeOffset = 0,
               std::string const& Name      = "");

    TVector3D GetFieldA () const;
    TVector3D GetFieldB () const;
    TVector3D GetPeriod () const;
    int       GetNPeriods () const;
    TVector3D GetCenter () const;
    double    GetPhase () const;
    double    GetTaper () const;

    double    GetFrequency () const;
    double    GetFrequencyPhase () const;
    double    GetTimeOffset () const;

    void Print (std::ostream& os) const;


  private:
    TVector3D fFieldA;
    TVector3D fFieldB;
    TVector3D fPeriod;
    TVector3D fPeriodUnitVector;
    double    fPeriodLength;
    int       fNPeriods;
    TVector3D fCenter;
    double    fPhase;
    double    fTaper;

    double fFrequency;
    double fFrequencyPhase;
    double fTimeOffset;

    double fUndulatorLength;

};


inline std::ostream& operator << (std::ostream& os, TField3D_IdealEPU const& o)
{
  // For easy printing
  os << "TField3D_IdealEPU " << "\n"
     << "Name                    " << o.GetName() << "\n"
     << "FieldA                  " << o.GetFieldA() << "\n"
     << "FieldB                  " << o.GetFieldB() << "\n"
     << "Period                  " << o.GetPeriod() << "  (" << o.GetPeriod().Mag() << " [m])\n"
     << "NPeriods                " << o.GetNPeriods() << "\n"
     << "Center                  " << o.GetCenter() << "\n"
     << "Phase                   " << o.GetPhase() << "\n"
     << "Taper                   " << o.GetTaper() << "\n"
     << "Frequency               " << o.GetFrequency() << "\n"
     << "FrequencyPhase          " << o.GetFrequencyPhase() << "\n"
     << "TimeOffset              " << o.GetTimeOffset() << "\n";

  return os;
}






#endif
