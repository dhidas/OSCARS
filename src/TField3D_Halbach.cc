////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Tue Nov 21 13:52:28 EST 2017
//
////////////////////////////////////////////////////////////////////


#include "TField3D_Halbach.h"


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

  // assume s is z

  return TVector3D(0, 0, 0);
}





void TField3D_Halbach::Print (std::ostream& os) const
{
  os << *this << std::endl;
  return;
}



double TField3D_Halbach::GetField () const
{
  return fField;
}

