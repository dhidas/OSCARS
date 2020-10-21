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

#include "TVector2D.h"

#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <ostream>

class TTriangle3DContainer
{
  public:
    TTriangle3DContainer ();
    ~TTriangle3DContainer ();

    void Add (TTriangle3D const&);

    size_t GetNPoints () const;
    TTriangle3D GetPoint (size_t const i) const;
    void Clear ();
    void ClearValues ();

    void AddToPoint (size_t const, double const);

    void TranslateSelf (TVector3D const&);

    void ReadSTLFile (std::string const& FileName,
                      double const Scale = 1,
                      TVector3D const& Rotations = TVector3D(0, 0, 0),
                      TVector3D const& Translation = TVector3D(0, 0, 0),
                      std::string const& Name = "");

    void WriteSTLFile (std::string const& FileName);

    std::string GetFileName () const;
    double      GetScale () const;
    TVector3D   GetRotations () const;
    TVector3D   GetTranslation () const;
    std::string GetName () const;




  private:
    void RotateSelfXYZ (TVector3D const&);

    std::vector<TTriangle3D> fT;  // List of triangles

    TVector3D fBBox[4];           // Bounding box
    std::string fFileName;        // Name of input file (if any)
    double fScale;                // Scale factor on geometry
    TVector3D   fRotations;       // Rotations used
    TVector3D   fTranslation;     // Translation used
    std::string fName;            // Object name
};








inline std::ostream& operator << (std::ostream& os, TTriangle3DContainer& o)
{
  // For easy printing
  os << "Name:             " << o.GetName() << "\n"
     << "FilenName:        " << o.GetFileName() << "\n"
     << "Scale:            " << o.GetScale() << "\n"
     << "Rotations:        " << o.GetRotations() << "\n"
     << "Translation:      " << o.GetTranslation() << "\n"
     << "NTriangles:       " << o.GetNPoints() << std::endl;


  return os;
}















#endif
