#ifndef GUARD_TFieldContainer
#define GUARD_TFieldContainer
////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Thu Sep 22 08:19:53 EDT 2016
//
// This class is meant to be a container class for all fields in
// a given simulation.  It will sum all B contributions and return
// the sum Fx, Fy, Fz, or zero where there is no defined field
//
////////////////////////////////////////////////////////////////////

#include <vector>

#include "TField.h"
#include "TVector2D.h"
#include "TVector3D.h"

class TFieldContainer
{
  public:
    TFieldContainer ();
    TFieldContainer (TField*);
    ~TFieldContainer ();

    void AddField (TField*);

    double    GetFx (double const, double const, double const) const;
    double    GetFy (double const, double const, double const) const;
    double    GetFz (double const, double const, double const) const;
    TVector3D GetF  (double const, double const, double const) const;
    TVector3D GetF  (TVector3D const&) const;

    size_t GetNFields () const;

    void      Clear ();

    void WriteToFile (std::string const& OutFileName, std::string const& OutFormat, TVector2D const& XLim, int const NX, TVector2D const& YLim, int const NY, TVector2D const& ZLim, int const NZ, std::string const Comment = "");

  private:
    std::vector<TField*> fFields;
};




























#endif
