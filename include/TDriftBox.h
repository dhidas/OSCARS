#ifndef GUARD_TDriftBox_h
#define GUARD_TDriftBox_h
////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Thu Jul 27 10:55:02 EDT 2017
//
// Class for a drift box
//
////////////////////////////////////////////////////////////////////

#include "TDriftVolume.h"

class TDriftBox : public TDriftVolume
{
  public:
    TDriftBox (TVector3D   const& Width,
               TVector3D   const& Center = TVector3D(0, 0, 0),
               TVector3D   const& Rotations = TVector3D(0, 0, 0),
               std::string const& Name = "",
               bool        const  RecordTrajectory = false);

    ~TDriftBox ();

    bool IsInside (TVector3D const& X) const;
    void Print (std::ostream& os) const;

    TVector3D GetWidth () const;
    TVector3D GetCenter () const;
    TVector3D GetRotated () const;

  private:
    TVector3D fWidth;
    TVector3D fCenter;
    TVector3D fRotated;

    bool fIgnoreAxisX;
    bool fIgnoreAxisY;
    bool fIgnoreAxisZ;

};


inline std::ostream& operator << (std::ostream& os, TDriftBox const& o)
{
  // For easy printing
  os << "TDriftBox           " << "\n"
     << "Name                " << o.GetName() << "\n"
     << "Width               " << o.GetWidth() << "\n"
     << "Rotations           " << o.GetRotated() << "\n"
     << "Center              " << o.GetCenter() << "\n"
     << "RecordTrajectory    " << o.RecordTrajectory() << "\n";

  return os;
}










#endif

