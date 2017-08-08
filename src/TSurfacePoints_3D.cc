#include "TSurfacePoints_3D.h"

TSurfacePoints_3D::TSurfacePoints_3D ()
{
  // Constructor
  fHasNormal = true;
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




bool TSurfacePoints_3D::HasNormal () const
{
  return fHasNormal;
}



void TSurfacePoints_3D::AddPoint(TVector3D const& X)
{
  fHasNormal = false;
  fPoints.push_back(TSurfacePoint(X, TVector3D(0, 0, 0)));
  return;
}




void TSurfacePoints_3D::AddPoint(TVector3D const& X, TVector3D const& N)
{
  fPoints.push_back(TSurfacePoint(X, N));
  return;
}
