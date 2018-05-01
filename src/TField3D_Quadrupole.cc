#include "TField3D_Quadrupole.h"

#include "TOSCARSSR.h"

#include <cmath>

TField3D_Quadrupole::TField3D_Quadrupole (std::string const& Name)
{
  // Default constructor

  // Set the name and default scale factors
  this->SetName(Name);
  this->SetScaleFactorMinimumMaximum();
}




TField3D_Quadrupole::TField3D_Quadrupole (double const K,
                                          double const Width,
                                          TVector3D const& Rotations,
                                          TVector3D const& Translation,
                                          double    const  Frequency,
                                          double    const  FrequencyPhase,
                                          double    const  TimeOffset,
                                          std::string const& Name) {
  // Constructor

  // Set the name and default scale factors
  this->SetName(Name);
  this->SetScaleFactorMinimumMaximum();

  fK = K;
  fWidth = Width;
  fRotations = Rotations;
  fTranslation = Translation;

  fFrequency = Frequency;
  fFrequencyPhase = FrequencyPhase;
  fTimeOffset = TimeOffset;
}




TField3D_Quadrupole::~TField3D_Quadrupole ()
{
  // Default constructor
}




TVector3D TField3D_Quadrupole::GetF (double const X, double const Y, double const Z, double const T) const
{
  return this->GetF(TVector3D(X, Y, Z));
}




TVector3D TField3D_Quadrupole::GetF (TVector3D const& X, double const T) const
{
  // Get the magnetic field at a point in space.

  TVector3D P = X;
  P.RotateSelfXYZ(fRotations);

  P -= fTranslation;

  if (fabs(P.GetZ()) > fWidth) {
    return TVector3D(0, 0, 0);
  }

  TVector3D Ret(fK * P.GetY(), fK * P.GetX(), 0);
  Ret.RotateSelfXYZ(fRotations);
  if (fFrequency == 0) {
    // It has no time dependence
    return Ret;
  }
  return Ret * cos(TOSCARSSR::TwoPi() * fFrequency * (T + fTimeOffset) + fFrequencyPhase);
}









double TField3D_Quadrupole::GetFrequency () const
{
  // Return the frequency
  return fFrequency;
}




double TField3D_Quadrupole::GetFrequencyPhase () const
{
  // Return the frequency
  return fFrequencyPhase;
}




double TField3D_Quadrupole::GetTimeOffset () const
{
  // Return the time offset for the frequency
  return fTimeOffset;
}




void TField3D_Quadrupole::Print (std::ostream& os) const
{
  os << *this << std::endl;
  return;
}




double TField3D_Quadrupole::GetK () const
{
  return fK;
}




double TField3D_Quadrupole::GetWidth () const
{
  return fWidth;
}




TVector3D const& TField3D_Quadrupole::GetRotations () const
{
  return fRotations;
}




TVector3D const& TField3D_Quadrupole::GetTranslation () const
{
  return fTranslation;
}
