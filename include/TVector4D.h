#ifndef GUARD_TVector4D_h
#define GUARD_TVector4D_h
////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Tue Jun 28 08:28:51 EDT 2016
//
// This class provides functionality for a 4-vector, or otherwise
// called a lorentz vector.  Operations are
//
////////////////////////////////////////////////////////////////////

#include "TVector3D.h"


class TVector4D
{
  public:
    TVector4D ();
    TVector4D (double const, double const, double const, double const);
    TVector4D (TVector3D const&, double const);
    ~TVector4D ();

    inline double GetX () const;
    inline double GetY () const;
    inline double GetZ () const;
    inline double GetT () const;

    inline double GetPx () const;
    inline double GetPy () const;
    inline double GetPz () const;
    inline double GetE  () const;

    inline TVector3D GetVector () const;


  private:
    TVector3D fP;
    double    fE;

};



inline double TVector4D::GetX () const
{
  return fP.GetX();
}



inline double TVector4D::GetY () const
{
  return fP.GetY();
}



inline double TVector4D::GetZ () const
{
  return fP.GetZ();
}



inline double TVector4D::GetT () const
{
  return fE;
}





inline double TVector4D::GetPx () const
{
  return fP.GetX();
}



inline double TVector4D::GetPy () const
{
  return fP.GetY();
}



inline double TVector4D::GetPz () const
{
  return fP.GetZ();
}



inline double TVector4D::GetE () const
{
  return fE;
}







#endif
