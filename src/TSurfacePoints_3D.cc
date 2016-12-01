#include "TSurfacePoints_3D.h"

TSurfacePoints_3D::TSurfacePoints_3D ()
{
  // Constructor
}




TSurfacePoints_3D::~TSurfacePoints_3D ()
{
  // Self destruction
}




TSurfacePoint const TSurfacePoints_3D::GetPoint (size_t const i) const
{
  return fPoints[i];
}




size_t TSurfacePoints_3D::GetNPoints () const
{
  return fPoints.size();
}




double TSurfacePoints_3D::GetX1 (size_t const i) const
{
  return 0;
}




double TSurfacePoints_3D::GetX2 (size_t const i) const
{
  return 0;
}




void TSurfacePoints_3D::AddPoint(TSurfacePoint const& P)
{
  fPoints.push_back(P);
  return;
}




void TSurfacePoints_3D::AddPoint(TVector3D const& X, TVector3D const& N)
{
  fPoints.push_back(TSurfacePoint(X, N));
  return;
}




void TSurfacePoints_3D::AddPoint(double const& X, double const& Y, double const& Z, double const& NX, double const& NY, double const& NZ)
{
  fPoints.push_back(TSurfacePoint(X, Y, Z, NX, NY, NZ));
  return;
}
