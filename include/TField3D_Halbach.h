#ifndef GUARD_TField3D_Halbach_h
#define GUARD_TField3D_Halbach_h
////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Tue Nov 21 13:52:28 EST 2017
//
////////////////////////////////////////////////////////////////////

#include "TField.h"

class TField3D_Halbach : public TField
{
  public:
    TField3D_Halbach ();
    TField3D_Halbach (double const Field,
                      double const Period,
                      int    const NPeriods,
                      double const Gap,
                      double const MagnetHeight,
                      double const MagnetWidth,
                      int    const NPerPeriod,
                      std::string const& Name = "");
    ~TField3D_Halbach ();

    TVector3D GetF  (double const, double const, double const, double const T = 0) const;
    TVector3D GetF  (TVector3D const&, double const T = 0) const;

    void Print (std::ostream&) const;

    double GetField () const;

  private:
    double fField;
    double fPeriod;
    int    fNPeriods;
    double fGap;
    double fMagnetHeight;
    double fMagnetWidth;
    int    fNPerPeriod;
};



inline std::ostream& operator << (std::ostream& os, TField3D_Halbach const& o)
{
  // For easy printing
  os << "TField3D_Halbach" << "\n"
     << "Name:           " << o.GetName() << "\n"
     << "Field:          " << o.GetField() << "\n";

  return os;
}










#endif
