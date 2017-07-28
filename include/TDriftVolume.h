#ifndef GUARD_TDriftVolume_h
#define GUARD_TDriftVolume_h
////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Thu Jul 27 10:55:02 EDT 2017
//
// Basic class for drift volumes.  All drift volumes should inherit
// from this class
//
////////////////////////////////////////////////////////////////////

#include <string>

#include "TVector3D.h"

class TDriftVolume
{
  public:
    virtual ~TDriftVolume()
    {
    }

    virtual bool IsInside (TVector3D const& X) const = 0;
    virtual void Print (std::ostream&) const = 0;

    void SetName (std::string const& Name)
    {
      fName = Name;
      return;
    }
    std::string const& GetName () const
    {
      return fName;
    }

    void SetRecordTrajectory (bool const RT)
    {
      fRecordTrajectory = RT;
      return;
    }
    bool RecordTrajectory () const
    {
      return fRecordTrajectory;
    }

  private:
    std::string fName;
    bool fRecordTrajectory;

};







#endif
