#ifndef GUARD_TVector3DC_h
#define GUARD_TVector3DC_h
////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Wed May  4 08:43:33 EDT 2016
//
// This class is a basic 3-dimensional complex vector with some simple
// operators defined to make the math more transparent.  This is
// largly modeled after the TVector3 class of ROOT.
//
////////////////////////////////////////////////////////////////////


#include <complex>

#include "TVector3D.h"

class TVector3DC
{
  public:
    TVector3DC ();
    TVector3DC (TVector3D const&);
    TVector3DC (std::complex<double> const&, std::complex<double> const&, std::complex<double> const&);
    ~TVector3DC ();

    std::complex<double> GetX () const;
    std::complex<double> GetY () const;
    std::complex<double> GetZ () const;

    void SetX (std::complex<double> const&);
    void SetY (std::complex<double> const&);
    void SetZ (std::complex<double> const&);
    void SetXYZ (std::complex<double> const&, std::complex<double> const&, std::complex<double> const&);

    double Perp2 (TVector3DC const&) const;
    std::complex<double> Dot(TVector3DC const&) const;
    TVector3DC Cross (TVector3DC const&) const;
    TVector3DC UnitVector () const;
    TVector3DC CC () const;
    double Mag2 () const;
    double Mag () const;
    std::complex<double> MagC2 () const;
    std::complex<double> MagC () const;


    // Operators
    TVector3DC  operator  + (TVector3DC const&) const;
    TVector3DC  operator  - (TVector3DC const&) const;
    //TVector3DC  operator  * (double const&) const;
    TVector3DC  operator  / (double const&) const;
    TVector3DC  operator  - ();
    TVector3DC& operator  = (TVector3DC const&);
    TVector3DC& operator += (TVector3DC const&);
    TVector3DC& operator -= (TVector3DC const&);
    TVector3DC& operator *= (double const&);
    TVector3DC& operator /= (double const&);
    TVector3DC& operator *= (std::complex<double> const&);
    TVector3DC& operator /= (std::complex<double> const&);
    bool       operator == (TVector3DC const&) const;
    bool       operator != (TVector3DC const&) const;

  private:
    std::complex<double> fX;
    std::complex<double> fY;
    std::complex<double> fZ;

};



inline std::ostream& operator << (std::ostream& os, TVector3DC const& o)
{
  os << "(" << o.GetX() << ", " << o.GetY() << ", " << o.GetZ() << ")";
  return os;
}





inline TVector3DC operator * (double const& V, TVector3DC const& R)
{
  // Multiply vector by some scalar
  return TVector3DC(R.GetX() * V, R.GetY() * V, R.GetZ() * V);
}

inline TVector3DC operator * (TVector3DC const& L, double const& V)
{
  // Multiply vector by some scalar
  return TVector3DC(L.GetX() * V, L.GetY() * V, L.GetZ() * V);
}



inline TVector3DC operator * (std::complex<double> const& V, TVector3DC const& R)
{
  // Multiply vector by some scalar
  return TVector3DC(R.GetX() * V, R.GetY() * V, R.GetZ() * V);
}

inline TVector3DC operator * (TVector3DC const& L, std::complex<double> const& V)
{
  // Multiply vector by some scalar
  return TVector3DC(L.GetX() * V, L.GetY() * V, L.GetZ() * V);
}


inline TVector3DC operator + (TVector3DC const& L, TVector3D& R)
{
  // Multiply vector with and without complex
  return TVector3DC(L.GetX() + R.GetX(), L.GetY() + R.GetY(), L.GetZ() + R.GetZ());
}


inline TVector3DC operator + (TVector3D const& L, TVector3DC& R)
{
  // Multiply vector with and without complex
  return TVector3DC(L.GetX() + R.GetX(), L.GetY() + R.GetY(), L.GetZ() + R.GetZ());
}













#endif
