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
#include <mutex>

#include "TParticleTrajectoryPoint.h"
#include "TVector3D.h"

class TParticleTrajectoryPoints
{
  public:
    TParticleTrajectoryPoints ();
    TParticleTrajectoryPoints (const TParticleTrajectoryPoints& TPTP);
    TParticleTrajectoryPoints (double const);
    ~TParticleTrajectoryPoints ();

    TParticleTrajectoryPoint const& GetPoint  (size_t const) const;

    TVector3D const& GetX      (size_t const) const;
    TVector3D const& GetB      (size_t const) const;
    void             SetB      (size_t const, TVector3D const&);
    TVector3D        GetV      (size_t const) const;
    TVector3D const& GetAoverC (size_t const) const;
    void             SetAoverC (size_t const, TVector3D const&);
    TVector3D        GetA      (size_t const) const;
    double           GetT      (size_t const) const;

    double           GetTStart () const;
    double           GetTStop  () const;

    double GetDeltaT () const;
    void   SetDeltaT (double const);
    size_t GetNPoints () const;

    std::vector<TParticleTrajectoryPoint> const& GetTrajectory() const;
    std::vector<double> const& GetTimePoints () const;

    void AddPoint (TParticleTrajectoryPoint const& P, double const T = 0);
    void AddPoint (TVector3D const&, TVector3D const&, TVector3D const&, double const T = 0);
    void AddPoint (double const, double const, double const, double const, double const, double const, double const, double const, double const, double const T = 0);

    void Reserve (size_t const);
    void ReverseArrays ();

    void WriteToFile       (std::string const& FileName, std::string const& FormatIn = "default") const;
    void WriteToFileBinary (std::string const& FileName, std::string const& FormatIn = "default") const;

    void ReadFromFileFormat (std::string const& FileName, std::string const& FormatIn = "default");
    void ReadFromFileBinary (std::string const& FileName, std::string const& FormatIn = "");

    void ConstructBetaAtPoints ();
    void ConstructAoverCAtPoints ();

    void Lock ();
    void UnLock ();

    void Clear ();

    TParticleTrajectoryPoints& operator=( const TParticleTrajectoryPoints& other ) {
      fLock_mutex = new std::mutex();
      return *this;
    }

    TParticleTrajectoryPoints& operator=( TParticleTrajectoryPoints& rhs ) {
      fLock_mutex = new std::mutex();
      return *this;
    };



  private:
    std::vector<TParticleTrajectoryPoint> fP;  // Trajectory points (x, beta, a/c)
    std::vector<double> fT;                    // Time in [s]


    // For equidistant time steps use single DeltaT
    double fDeltaT;

    // Mutex lock.  If locked don't write or read.  Lock only to be used for writing.
    // If not locked, assume trajectory is complete and will not change
    std::mutex* fLock_mutex;

};





















#endif
