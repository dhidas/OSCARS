#ifndef GUARD_T3DScalarContainer_h
#define GUARD_T3DScalarContainer_h
////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Tue Jun 21 17:28:24 EDT 2016
//
////////////////////////////////////////////////////////////////////

#include "TVector3D.h"

#include <fstream>
#include <vector>
#include <string>
#include <stdexcept>

class T3DScalar
{
  public:
    T3DScalar (TVector3D const& X, double const V)
    {
      fX = X;
      fV = V;
    }

    ~T3DScalar ()
    {
    }

    TVector3D const& GetX () const
    {
      return fX;
    }

    double GetV () const
    {
      return fV;
    }

    void SetV (double const V)
    {
      fV = V;
      return;
    }

  private:
    TVector3D fX;
    double    fV;
};





class T3DScalarContainer
{
  public:
    T3DScalarContainer ();
    ~T3DScalarContainer ();

    void AddPoint (TVector3D const&, double const);
    void AddToPoint (size_t const, double const);

    void Clear ();
    void AverageFromFilesText (std::vector<std::string> const&, int const Dimension);
    void AverageFromFilesBinary (std::vector<std::string> const&, int const Dimension);

    size_t GetNPoints () const;

    T3DScalar const& GetPoint (size_t const) const;

    void WriteToFileText (std::string const&, int const);
    void WriteToFileBinary (std::string const&, int const);

  private:
    std::vector<T3DScalar> fValues;
    std::vector<double> fCompensation;
};











#endif
