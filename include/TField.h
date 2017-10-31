#ifndef GUARD_TField_h
#define GUARD_TField_h
////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Tue Sep 20 07:47:09 EDT 2016
//
////////////////////////////////////////////////////////////////////

#include "TVector3D.h"

#include <string>
#include <stdexcept>

class TField
{
  // This class is designed to be a base class for a field object.

  public:
    virtual TVector3D GetF  (double const, double const, double const, double const T = 0) const = 0;
    virtual TVector3D GetF  (TVector3D const&, double const T = 0) const = 0;

    void SetName (std::string const& Name)
    {
      fName = Name;
      return;
    }
    std::string const& GetName () const
    {
      return fName;
    }

    // Scale factor
    void SetScaleFactor (double const& s)
    {
      if (fScaleFactorMinimum > fScaleFactorMaximum || (s >= fScaleFactorMinimum && s <= fScaleFactorMaximum)) {
        fScaleFactor = s;
      } else {
        throw std::out_of_range("TField::SetScaleFactor input is incorrect");
      }
      return;
    }
    double GetScaleFactor () const
    {
      return fScaleFactor;
    }
    void SetScaleFactorMinimumMaximum(double const& Min = 1, double const& Max = 1)
    {
      fScaleFactorMinimum = Min;
      fScaleFactorMaximum = Max;
      if (fScaleFactorMinimum == fScaleFactorMaximum) {
        this->SetScaleFactor(fScaleFactorMaximum);
      }
      return;
    }

    virtual void Print (std::ostream&) const = 0;

    virtual ~TField () {};

  private:
    std::string fName;

    // Scale factor for field magnitude.  This must be implemented in the following way:
    //  If min > max, assume the limits are +/- infinity for fScaleFactor.  fScaleFactor
    //  should always be initialized to 1.
    double fScaleFactor;
    double fScaleFactorMinimum;
    double fScaleFactorMaximum;

};





#endif
