#ifndef GUARD_TSurfaceOfPoints_h
#define GUARD_TSurfaceOfPoints_h
////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Mon May  2 17:14:01 EDT 2016
//
// This class is actually intended to hold a list of points and
// their normals (as TSurfacePoint objects)
//
////////////////////////////////////////////////////////////////////

#include "TSurfacePoints.h"

#include <vector>

class TSurfaceOfPoints
{
  public:
    TSurfaceOfPoints ();
    ~TSurfaceOfPoints ();

    TSurfacePoint const& GetPoint (size_t const) const;
    size_t GetNPoints () const;


    void AddPoint (TSurfacePoint const&);
    void AddPoint (TVector3D const&, TVector3D const&);
    void AddPoint (double const&, double const&, double const&, double const&, double const&, double const&);


  private:
    std::vector<TSurfacePoint> fPoints;

};

#endif
