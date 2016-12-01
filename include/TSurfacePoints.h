#ifndef GUARD_TSurfacePoints_h
#define GUARD_TSurfacePoints_h
////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Mon May  2 17:14:01 EDT 2016
//
// This is a pure base class.  The name may be misleading because
// one can hide a function behind a derivative class, which in the
// end looks like discritized points.  This is the purpose of this
// class, to handle both
//
////////////////////////////////////////////////////////////////////


#include "TSurfacePoint.h"

#include <vector>

class TSurfacePoints
{
  public:
    virtual TSurfacePoint const GetPoint (size_t const) const = 0;
    virtual size_t GetNPoints () const = 0;

    virtual double GetX1 (size_t const) const = 0;
    virtual double GetX2 (size_t const) const = 0;

};

#endif
