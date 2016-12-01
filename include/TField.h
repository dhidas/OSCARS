#ifndef GUARD_TField_h
#define GUARD_TField_h
////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Tue Sep 20 07:47:09 EDT 2016
//
////////////////////////////////////////////////////////////////////

#include "TVector3D.h"

class TField
{
  // This class is designed to be a base class for a field object.

  public:
    virtual double    GetFx (double const, double const, double const) const = 0;
    virtual double    GetFy (double const, double const, double const) const = 0;
    virtual double    GetFz (double const, double const, double const) const = 0;
    virtual TVector3D GetF  (double const, double const, double const) const = 0;
    virtual TVector3D GetF  (TVector3D const&) const = 0;

    virtual ~TField () {};


};





#endif
