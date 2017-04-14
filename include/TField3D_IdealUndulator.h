#ifndef GUARD_TField3D_IdealUndulator_h
#define GUARD_TField3D_IdealUndulator_h
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



class TField3D_IdealUndulator : public TField
{
  public:
    TField3D_IdealUndulator ();
    TField3D_IdealUndulator (TVector3D const& Field,
                             TVector3D const& Period,
                             int       const  NPeriods,
                             TVector3D const& Center    = TVector3D(0, 0, 0),
                             double    const  Phase     = 0,
                             double    const  Taper     = 0);

    ~TField3D_IdealUndulator ();

    double    GetFx (double const X, double const Y, double const Z) const;
    double    GetFy (double const X, double const Y, double const Z) const;
    double    GetFz (double const X, double const Y, double const Z) const;
    TVector3D GetF  (double const X, double const Y, double const Z) const;
    TVector3D GetF  (TVector3D const& X) const;

    void Init (TVector3D const& Field,
               TVector3D const& Period,
               int       const  NPeriods,
               TVector3D const& Center    = TVector3D(0, 0, 0),
               double    const  Phase     = 0,
               double    const Taper      = 0);

    TVector3D GetField () const;
    TVector3D GetPeriod () const;
    int       GetNPeriods () const;
    TVector3D GetCenter () const;
    double    GetPhase () const;
    double    GetTaper () const;

    void Print (std::ostream& os) const;


  private:
    TVector3D fField;
    TVector3D fPeriod;
    TVector3D fPeriodUnitVector;
    double    fPeriodLength;
    int       fNPeriods;
    TVector3D fCenter;
    double    fPhase;
    double    fTaper;

    double fUndulatorLength;

};


inline std::ostream& operator << (std::ostream& os, TField3D_IdealUndulator const& o)
{
  // For easy printing
  os << "TField3D_IdealUndulator " << "\n"
     << "Field                   " << o.GetField() << "\n"
     << "Period                  " << o.GetPeriod() << "  (" << o.GetPeriod().Mag() << " [m])\n"
     << "NPeriods                " << o.GetNPeriods() << "\n"
     << "Center                  " << o.GetCenter() << "\n"
     << "Phase                   " << o.GetPhase() << "\n"
     << "Taper                   " << o.GetTaper() << "\n";

  return os;
}






#endif
