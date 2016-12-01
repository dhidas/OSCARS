#ifndef GUARD_TSurfaceElement_h
#define GUARD_TSurfaceElement_h
////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Fri May 13 14:52:05 EDT 2016
//
// This is a base class for any surface element.
//
////////////////////////////////////////////////////////////////////

#include "TSurfacePoint.h"

class TSurfaceElement
{
  public:
    virtual TSurfacePoint const& GetSurfacePoint () const = 0;
    virtual TVector3D     const& GetPoint () const = 0;
    virtual TVector3D     const& GetNormal () const = 0;
    virtual double        GetArea () const = 0;

  protected:
    TSurfacePoint fSurfacePoint;
    double fArea;

};







#endif
