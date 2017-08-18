#ifndef GUARD_TParticleBeamContainer_h
#define GUARD_TParticleBeamContainer_h
////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Thu Jun  9 14:46:05 EDT 2016
//
// This class is meant to contain several different particle beams.
// It should be the main class used by the module
//
////////////////////////////////////////////////////////////////////

#include <map>

#include "TVector3D.h"
#include "TParticleA.h"
#include "TParticleBeam.h"


class TParticleBeamContainer
{
  public:
    TParticleBeamContainer ();
    ~TParticleBeamContainer ();

    TParticleBeam& AddNewParticleBeam (std::string const& Type,
                                       std::string const& Name,
                                       TVector3D const& X0,
                                       TVector3D const& D0,
                                       double const E0,
                                       double const T0,
                                       double const Current,
                                       double const Weight = 1,
                                       double const Charge = 0,
                                       double const Mass = 0);

    TParticleBeam& AddNewParticleBeam (std::string const& Beam,
                                       std::string const& Name,
                                       double const Weight = 1);

    TParticleA GetNewParticle ();
    TParticleBeam& GetParticleBeam (size_t const);
    TParticleBeam& GetParticleBeam (std::string const&);
    TParticleBeam& GetRandomBeam ();
    size_t GetRandomBeamIndexByWeight () const;
    size_t GetNParticleBeams () const;
    void Clear ();

    void SetEmittance (std::string const& Beam,
                       TVector2D const& Emittance);

    void SetTwissParameters (std::string const& Beam,
                             TVector2D const& Beta,
                             TVector2D const& Alpha,
                             TVector2D const& Gamma,
                             TVector3D const& Lattice_Reference,
                             bool const HasReferencePoint);

  private:

    std::vector<double>        fParticleBeamWeightSums;
    std::vector<TParticleBeam> fParticleBeams;

    // A map between string and beam index
    std::map<std::string, size_t> fParticleBeamMap;

};


inline std::ostream& operator << (std::ostream& os, TParticleBeamContainer& o)
{
  // For easy printing
  os << "TParticleBeamContainer has " << o.GetNParticleBeams() << " beams" << std::endl;

  size_t const N = o.GetNParticleBeams();

  for (size_t i = 0; i != N; ++i) {
    TParticleBeam B = o.GetParticleBeam(i);

    os << B << std::endl;
  }

  return os;
}









#endif
