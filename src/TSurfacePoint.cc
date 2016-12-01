#include "TSurfacePoint.h"


TSurfacePoint::TSurfacePoint ()
{
  // Default constructor
}




TSurfacePoint::TSurfacePoint (TVector3D const& x, TVector3D const& n)
{
  // Constructor when you give me vectors
  fX = x;
  fN = n.UnitVector();
}




TSurfacePoint::TSurfacePoint (double const& x, double const& y, double const& z, double const& nx, double const& ny, double const& nz)
{
  // Constructor when you give me vectors
  fX.SetXYZ(x, y, z);
  fN = TVector3D(nx, ny, nz).UnitVector();
}



TSurfacePoint::~TSurfacePoint ()
{
  // Destruction! Nothing to do here
}



TVector3D const& TSurfacePoint::GetPoint () const
{
  return fX;
}



double TSurfacePoint::GetX () const
{
  return fX.GetX();
}



double TSurfacePoint::GetY () const
{
  return fX.GetY();
}



double TSurfacePoint::GetZ () const
{
  return fX.GetZ();
}



void TSurfacePoint::SetXYZ (TVector3D const& x)
{
  // Set the xyz of the location
  fX = x;
  return;
}



void TSurfacePoint::SetXYZ (double const& x, double const& y, double const& z)
{
  // Set the xyz of the location
  fX.SetXYZ(x, y, z);
  return;
}





TVector3D const& TSurfacePoint::GetNormal () const
{
  return fN;
}



double TSurfacePoint::GetNormalX () const
{
  return fN.GetX();
}



double TSurfacePoint::GetNormalY () const
{
  return fN.GetY();
}



double TSurfacePoint::GetNormalZ () const
{
  return fN.GetZ();
}







void TSurfacePoint::SetNormalXYZ (TVector3D const& n)
{
  // Set the xyz of the location
  fN = n.UnitVector();
  return;
}




void TSurfacePoint::SetNormalXYZ (double const& x, double const& y, double const& z)
{
  // Set the xyz of the location
  fN = TVector3D(x, y, z).UnitVector();
  return;
}




