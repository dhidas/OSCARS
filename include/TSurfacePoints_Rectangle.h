#ifndef GUARD_TSurfacePoints_Rectangle_h
#define GUARD_TSurfacePoints_Rectangle_h
////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Mon May  2 17:55:52 EDT 2016
//
// Surface of a simple rectangle
//
////////////////////////////////////////////////////////////////////

#include "TSurfacePoints.h"

#include <string>
#include <stdexcept>

class TSurfacePoints_Rectangle : public TSurfacePoints
{
  public:
    TSurfacePoints_Rectangle ();
    TSurfacePoints_Rectangle (int const, int const, TVector3D const&, TVector3D const&, TVector3D const&, int const);
    TSurfacePoints_Rectangle (std::string const&, int const, int const, double const, double const, TVector3D const&, TVector3D const&, int const);
    ~TSurfacePoints_Rectangle ();

    void Init (int const, int const, TVector3D const&, TVector3D const&, TVector3D const&, int const);
    void Init (std::string const&, int const, int const, double const, double const, TVector3D const&, TVector3D const&, int const);
    TSurfacePoint const GetPoint (size_t const) const;
    TVector3D GetXYZ (size_t const) const;
    size_t GetNPoints () const;

    double GetX1 (size_t const) const;
    double GetX2 (size_t const) const;

    double GetElementArea () const;

  private:
    int fNX1;
    int fNX2;


    int fNormal;

    double fX1StepSize;
    double fX2StepSize;

    size_t fNPoints;

    TVector3D fNormalVector;
    TVector3D fX1Vector;
    TVector3D fX2Vector;
    TVector3D fStartVector;


};




#endif
