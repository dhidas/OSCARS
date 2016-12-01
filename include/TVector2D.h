#ifndef GUARD_TVector2D_h
#define GUARD_TVector2D_h
////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Thu Jul 28 11:45:52 EDT 2016
//
// A 2D Vector class.
//
////////////////////////////////////////////////////////////////////


#include <iostream>


class TVector2D
{
  public:
    TVector2D ();
    TVector2D (double const, double const);
    ~TVector2D ();

    inline double GetX () const;
    inline double GetY () const;
    inline TVector2D Orthogonal () const;

    void SetX (double const);
    void SetY (double const);
    void SetXY (double const, double const);

    double Mag () const;
    double Mag2 () const;
    double Dot (TVector2D const&) const;
    double Perp2 (TVector2D const&) const;
    TVector2D Cross (TVector2D const&) const;
    TVector2D UnitVector () const;

    void RotateSelf (double const);


    // Operators
    bool       operator  < (TVector2D const&) const;
    TVector2D  operator  + (TVector2D const&) const;
    TVector2D  operator  - (TVector2D const&) const;
    TVector2D  operator  / (double const) const;
    TVector2D  operator  - ();
    TVector2D& operator  = (TVector2D const&);
    TVector2D& operator += (TVector2D const&);
    TVector2D& operator -= (TVector2D const&);
    TVector2D& operator *= (double const);
    TVector2D& operator /= (double const);
    bool       operator == (TVector2D const&) const;
    bool       operator != (TVector2D const&) const;
    double     operator [] (int const) const;
    double&    operator [] (int const);



  private:
    double fX;
    double fY;
    double fZ;

};





inline std::ostream& operator << (std::ostream& os, TVector2D const& o)
{
  // For easy printing
  os << "(" << o.GetX() << ", " << o.GetY() << ")";
  return os;
}








inline double TVector2D::GetX () const
{
  // Return X-component of vector
  return fX;
}




inline double TVector2D::GetY () const
{
  // Return Y-component of vector
  return fY;
}










inline TVector2D operator * (double const V, TVector2D const& R)
{
  // Multiply vector by some scalar
  return TVector2D(R.GetX() * V, R.GetY() * V);
}




inline TVector2D operator * (TVector2D const& L, double const V)
{
  // Multiply vector by some scalar
  return TVector2D(L.GetX() * V, L.GetY() * V);
}








#endif
