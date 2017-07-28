#ifndef GUARD_TDriftVolumeContainer
#define GUARD_TDriftVolumeContainer
////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Thu Jul 27 12:10:05 EDT 2017
//
// This class is meant to be a container class for all drift volumes
// a given simulation.
//
////////////////////////////////////////////////////////////////////

#include <vector>

#include "TDriftVolume.h"

class TDriftVolumeContainer
{
  public:
    TDriftVolumeContainer ();
    ~TDriftVolumeContainer ();

    bool IsInside (TVector3D const& X);

    void AddDriftVolume (TDriftVolume*);
    TDriftVolume const& GetDriftVolume (size_t const i) const;

    size_t GetNDriftVolumes () const;

    void RemoveDriftVolume (std::string const& Name);
    void Clear ();


  private:
    std::vector<TDriftVolume*> fDriftVolumes;
};


inline std::ostream& operator << (std::ostream& os, TDriftVolumeContainer const& o)
{
  // For easy printing
  //
  size_t const N = o.GetNDriftVolumes();

  os << "TDriftVolumeContainer has " << N << " DriftVolumes" << std::endl;


  for (size_t i = 0; i != N; ++i) {
    o.GetDriftVolume(i).Print(os);
  }

  return os;
}




























#endif

