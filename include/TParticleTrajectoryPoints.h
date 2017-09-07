#ifndef GUARD_TParticleTrajectoryPoints_h
#define GUARD_TParticleTrajectoryPoints_h
////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Wed May 18 18:08:30 EDT 2016
//
////////////////////////////////////////////////////////////////////

#include <vector>

#include "TParticleTrajectoryPoint.h"
#include "TVector3D.h"

class TParticleTrajectoryPoints
{
  public:
    TParticleTrajectoryPoints ();
    TParticleTrajectoryPoints (double const);
    ~TParticleTrajectoryPoints ();

    TVector3D const& GetX      (size_t const) const;
    TVector3D const& GetB      (size_t const) const;
    TVector3D        GetV      (size_t const) const;
    TVector3D const& GetAoverC (size_t const) const;
    TVector3D        GetA      (size_t const) const;

    double GetDeltaT () const;
    void   SetDeltaT (double const);
    size_t GetNPoints () const;

    void AddPoint (TParticleTrajectoryPoint const& P, double const T = 0);
    void AddPoint (TVector3D const&, TVector3D const&, TVector3D const&, double const T = 0);
    void AddPoint (double const, double const, double const, double const, double const, double const, double const, double const, double const, double const T = 0);

    void ReverseArrays ();

    void WriteToFile       (std::string const&) const;
    void WriteToFileBinary (std::string const&) const;

    void ReadFromFile       (std::string const&);
    void ReadFromFileBinary (std::string const&);

    void Clear ();



  private:
    std::vector<TParticleTrajectoryPoint> fP;  // Trajectory points (x, beta, a/c)
    //std::vector<double> fT;       // Time.  Not used, but could be in the future


    // For equidistant time steps use single DeltaT
    double fDeltaT;


};





















#endif
