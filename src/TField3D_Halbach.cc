////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Tue Nov 21 13:52:28 EST 2017
//
////////////////////////////////////////////////////////////////////


#include "TField3D_Halbach.h"
#include "math.h"

TField3D_Halbach::TField3D_Halbach ()
{
  // Default constructor
}



TField3D_Halbach::TField3D_Halbach (double const Field,
                                    double const Period,
                                    int    const NPeriods,
                                    double const Gap,
                                    double const MagnetHeight,
                                    double const MagnetWidth,
                                    int    const NPerPeriod,
                                    std::string const& Name)
{
  fField = Field;
  fPeriod = Period;
  fNPeriods = NPeriods;
  fGap = Gap;
  fMagnetHeight = MagnetHeight;
  fMagnetWidth = MagnetWidth;
  fNPerPeriod = NPerPeriod;
  this->SetName(Name);
}




TField3D_Halbach::~TField3D_Halbach ()
{
  // Destruction!
}




TVector3D TField3D_Halbach::GetF (double const X, double const Y, double const Z, double const T) const
{
  return this->GetF(TVector3D(X, Y , Z), T);
}




TVector3D TField3D_Halbach::GetF (TVector3D const& X, double const T) const
{
  // Field for Halback undulator!

  double fLength = fPeriod * fNPeriods / 2.0;

  double fEpsilon =((double) fNPerPeriod) * fMagnetWidth / fPeriod;
  double fBY = -2 * fField * cos(2 * M_PI * X.GetZ() / fPeriod) * sin(fEpsilon * M_PI / fNPerPeriod) * fNPerPeriod / M_PI * exp(-M_PI * fGap / fPeriod) * (1 - exp(-2 * M_PI * fMagnetHeight / fPeriod));
  
  if (X.GetZ() < fLength + fPeriod && X.GetZ() > -fLength - fPeriod) {
    double a = 1.0;
    if (X.GetZ() > fLength || X.GetZ() < -fLength) {
      if (X.GetZ() > fLength + fPeriod / 2.0 || X.GetZ() < -fLength - fPeriod / 2.0) {
        a = 0.25;
      } else {
        a = 0.75;
      }
    }
    return TVector3D(0.0, a * fBY, 0.0); 
  } else {
    return TVector3D(0.0, 0.0, 0.0);  
  }
}


double TField3D_Halbach::GetField () const
{
  // Return the remnant field
  return fField;
}


double TField3D_Halbach::GetPeriod () const
{
  // Return the period
  return fPeriod;
}


int TField3D_Halbach::GetNPeriods () const
{
  // Return the number of periods
  return fNPeriods;
}


double TField3D_Halbach::GetGap () const
{
  // Return the gap length
  return fGap;
}


double TField3D_Halbach::GetMagnetHeight () const
{
  // Return the magnet height
  return fMagnetHeight;
}


double TField3D_Halbach::GetMagnetWidth () const
{
  // Return the magnet width
  return fMagnetWidth;
}


int TField3D_Halbach::GetNPerPeriod () const
{
  // Return the number of magnets per period
  return fNPerPeriod;
}


void TField3D_Halbach::Print (std::ostream& os) const
{
  os << *this << std::endl;
  return;
}
