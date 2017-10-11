#ifndef GUARD_TTriangle3D_h
#define GUARD_TTriangle3D_h
////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Fri Sep 22 17:40:11 EDT 2017
//
// Class to represent a 3D triangle
//
////////////////////////////////////////////////////////////////////

#include "TVector3D.h"


class TTriangle3D
{
  public:
    TTriangle3D (TVector3D const& A,
                 TVector3D const& B,
                 TVector3D const& C,
                 TVector3D const& N);

    TTriangle3D (double const Ax, double const Ay, double const Az,
                 double const Bx, double const By, double const Bz,
                 double const Cx, double const Cy, double const Cz,
                 double const Nx, double const Ny, double const Nz);

    ~TTriangle3D ();


    void Set (TVector3D const& A,
              TVector3D const& B,
              TVector3D const& C,
              TVector3D const& N);

    void Set (double const Ax, double const Ay, double const Az,
              double const Bx, double const By, double const Bz,
              double const Cx, double const Cy, double const Cz,
              double const Nx, double const Ny, double const Nz);

    void Translate (TVector3D const& T);
    void RotateSelfXYZ    (TVector3D const& R);

    TVector3D GetCenter () const;
    TVector3D GetNormal () const;

    TVector3D  operator [] (int const) const;
    TVector3D& operator [] (int const);
    TTriangle3D& operator += (TVector3D const&);


  private:
    TVector3D fA;
    TVector3D fB;
    TVector3D fC;
    TVector3D fN;
};




#endif
