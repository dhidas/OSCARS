#include "TField3D_Quadrupole.h"

#include <cmath>

TField3D_Quadrupole::TField3D_Quadrupole (std::string const& Name)
{
  // Default constructor

  // Set Name
  this->SetName(Name);
}




TField3D_Quadrupole::TField3D_Quadrupole (double const K,
                                          double const Width,
                                          TVector3D const& Rotations,
                                          TVector3D const& Translation,
                                          std::string const& Name = "") {
  // Constructor

  // Set Name
  this->SetName(Name);

  fK = K;
  fWidth = Width;
  fRotations = Rotations;
  fTranslation = Translation;
}




TField3D_Quadrupole::~TField3D_Quadrupole ()
{
  // Default constructor
}




double TField3D_Quadrupole::GetFx (double const X, double const Y, double const Z) const
{
  return this->GetF(TVector3D(X, Y, Z)).GetX();
}




double TField3D_Quadrupole::GetFy (double const X, double const Y, double const Z) const
{
  return this->GetF(TVector3D(X, Y, Z)).GetY();
}




double TField3D_Quadrupole::GetFz (double const X, double const Y, double const Z) const
{
  return this->GetF(TVector3D(X, Y, Z)).GetZ();
}




TVector3D TField3D_Quadrupole::GetF (double const X, double const Y, double const Z) const
{
  return this->GetF(TVector3D(X, Y, Z));
}




TVector3D TField3D_Quadrupole::GetF (TVector3D const& X) const
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
  return Ret;
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
