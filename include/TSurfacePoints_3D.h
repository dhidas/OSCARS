#ifndef GUARD_TSurfacePoints_3D_h
#define GUARD_TSurfacePoints_3D_h
////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Tue Jul 12 17:09:31 EDT 2016
//
// Surface of arbitrary points
//
////////////////////////////////////////////////////////////////////


#include "TSurfacePoints.h"

#include <vector>

class TSurfacePoints_3D : public TSurfacePoints
{
  public:
    TSurfacePoints_3D ();
    ~TSurfacePoints_3D ();

    TSurfacePoint const GetPoint (size_t const) const;
    size_t GetNPoints () const;

    // Don't have meaning for this object
    double GetX1 (size_t const) const;
    double GetX2 (size_t const) const;

    void AddPoint (TSurfacePoint const&);
    void AddPoint (TVector3D const&, TVector3D const&);
    void AddPoint (double const&, double const&, double const&, double const&, double const&, double const&);

  private:
    std::vector<TSurfacePoint> fPoints;

};

#endif
