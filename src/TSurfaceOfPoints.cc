#include "TSurfaceOfPoints.h"


TSurfaceOfPoints::TSurfaceOfPoints ()
{
  // Default constructor
}




TSurfaceOfPoints::~TSurfaceOfPoints ()
{
  // Destruction
}



TSurfacePoint const& TSurfaceOfPoints::GetPoint (size_t const i) const
{
  return fPoints[i];
}



size_t TSurfaceOfPoints::GetNPoints () const
{
  return fPoints.size();
}



void TSurfaceOfPoints::AddPoint (TSurfacePoint const& XN)
{
  fPoints.push_back(XN);
  return;
}



void TSurfaceOfPoints::AddPoint (TVector3D const& x, TVector3D const& n)
{
  fPoints.push_back( TSurfacePoint(x, n) );
  return;
}



void TSurfaceOfPoints::AddPoint (double const& x, double const& y, double const& z, double const& nx, double const& ny, double const& nz)
{
  fPoints.push_back( TSurfacePoint(x, y, z, nx, ny, nz) );
  return;
}
