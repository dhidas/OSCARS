#ifndef GUARD_TTriangle3DContainer_h
#define GUARD_TTriangle3DContainer_h
////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Fri Sep 22 19:12:02 EDT 2017
//
////////////////////////////////////////////////////////////////////

#include "TTriangle3D.h"

#include <vector>
#include <iostream>
#include <fstream>

class TTriangle3DContainer
{
  public:
    TTriangle3DContainer ();
    ~TTriangle3DContainer ();

    void Add (TTriangle3D const&);

    size_t GetNPoints () const;
    TTriangle3D GetPoint (size_t const i) const;
    void Clear ();

    void RotateSelfXYZ (TVector3D const&);
    void TranslateSelf (TVector3D const&);

    void ReadSTLFile (std::string const& FileName,
                      double const Scale);

    void WriteSTLFile (std::string const& FileName);



  private:
    std::vector<TTriangle3D> fT;

    TVector3D fBBox[2];
    double fScale;
};




















#endif
