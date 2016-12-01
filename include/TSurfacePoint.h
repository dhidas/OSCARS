#ifndef GUARD_TSurfacePoint_h
#define GUARD_TSurfacePoint_h
////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Mon May  2 15:50:30 EDT 2016
//
// This is a class for a surface point.  It should consist of a 3D
// point in space and a unit vector normal to the surface at that
// point, pointing "outward" if you will.
//
////////////////////////////////////////////////////////////////////

#include "TVector3D.h"


class TSurfacePoint
{
  public:
    TSurfacePoint ();
    TSurfacePoint (TVector3D const&, TVector3D const&);
    TSurfacePoint (double const&, double const&, double const&, double const&, double const&, double const&);
    ~TSurfacePoint ();

    TVector3D const& GetPoint () const;
    double GetX () const;
    double GetY () const;
    double GetZ () const;
    void SetXYZ (TVector3D const&);
    void SetXYZ (double const&, double const&, double const&);

    TVector3D const& GetNormal () const;
    double GetNormalX () const;
    double GetNormalY () const;
    double GetNormalZ () const;
    void SetNormalXYZ (TVector3D const&);
    void SetNormalXYZ (double const&, double const&, double const&);

  private:
    TVector3D fX;
    TVector3D fN;

};








inline std::ostream& operator << (std::ostream& os, TSurfacePoint const& o)
{
  os << "(" << o.GetPoint() << ", " << o.GetNormal() <<  ")";
  return os;
}







#endif
