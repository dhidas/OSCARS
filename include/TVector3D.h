#ifndef GUARD_TVector3D_h
#define GUARD_TVector3D_h
////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Fri Mar 25 10:24:13 EDT 2016
//
// This class is a basic 3-dimensional vector with some simple
// operators defined to make the math more transparent.  This is
// largly modeled after the TVector3 class of ROOT.
//
////////////////////////////////////////////////////////////////////


#include <iostream>


class TVector3D
{
  public:
    TVector3D ();
    TVector3D (double const, double const, double const);
    ~TVector3D ();

    inline double GetX () const;
    inline double GetY () const;
    inline double GetZ () const;
    inline TVector3D Orthogonal () const;

    void SetX (double const);
    void SetY (double const);
    void SetZ (double const);
    void SetXYZ (double const, double const, double const);
    void SetXYZ (TVector3D const&);

    double Mag () const;
    double Mag2 () const;
    double Dot (TVector3D const&) const;
    double Perp2 (TVector3D const&) const;
    TVector3D Cross (TVector3D const&) const;
    TVector3D UnitVector () const;

    void RotateSelfX (double const);
    void RotateSelfY (double const);
    void RotateSelfZ (double const);
    void RotateSelfXYZ (TVector3D const&);
    void RotateSelf (double const, TVector3D const&);


    // Operators
    bool       operator  < (TVector3D const&) const;
    TVector3D  operator  + (TVector3D const&) const;
    TVector3D  operator  - (TVector3D const&) const;
    TVector3D  operator  / (double const) const;
    TVector3D  operator  - ();
    TVector3D& operator  = (TVector3D const&);
    TVector3D& operator += (TVector3D const&);
    TVector3D& operator -= (TVector3D const&);
    TVector3D& operator *= (double const);
    TVector3D& operator /= (double const);
    bool       operator == (TVector3D const&) const;
    bool       operator != (TVector3D const&) const;
    double     operator [] (int const) const;
    double&    operator [] (int const);



  private:
    double fX;
    double fY;
    double fZ;

};





inline std::ostream& operator << (std::ostream& os, TVector3D const& o)
{
  // For easy printing
  os << "(" << o.GetX() << ", " << o.GetY() << ", " << o.GetZ() << ")";
  return os;
}








inline double TVector3D::GetX () const
{
  // Return X-component of vector
  return fX;
}




inline double TVector3D::GetY () const
{
  // Return Y-component of vector
  return fY;
}




inline double TVector3D::GetZ () const
{
  // Return Z-component of vector
  return fZ;
}




inline TVector3D TVector3D::Orthogonal() const
{
  // Return a vector which is orthogonal to this vector
  double xx = fX < 0.0 ? -fX : fX;
  double yy = fY < 0.0 ? -fY : fY;
  double zz = fZ < 0.0 ? -fZ : fZ;
  if (xx < yy) {
    return xx < zz ? TVector3D(0, fZ, -fY) : TVector3D(fY, -fX, 0);
  } else {
    return yy < zz ? TVector3D(-fZ, 0, fX) : TVector3D(fY, -fX, 0);
  }
}




inline TVector3D operator * (double const V, TVector3D const& R)
{
  // Multiply vector by some scalar
  return TVector3D(R.GetX() * V, R.GetY() * V, R.GetZ() * V);
}




inline TVector3D operator * (TVector3D const& L, double const V)
{
  // Multiply vector by some scalar
  return TVector3D(L.GetX() * V, L.GetY() * V, L.GetZ() * V);
}








#endif
