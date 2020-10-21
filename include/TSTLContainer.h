#ifndef GUARD_TSTLContainer_h
#define GUARD_TSTLContainer_h
////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Thu Sep 20 17:32:12 EDT 2018
//
// A container for all imported STL files
//
////////////////////////////////////////////////////////////////////

#include <string>
#include <vector>
#include <map>


#include "TTriangle3DContainer.h"


class TSTLContainer
{
  public:
    TSTLContainer ();
    ~TSTLContainer ();

    void AddFile (std::string const& InFileName,
                  double      const  Scale = 1,
                  TVector3D   const& Rotations = TVector3D(0, 0, 0),
                  TVector3D   const& Translation = TVector3D(0, 0, 0),
                  std::string const& Name = "");

    TTriangle3DContainer const& GetTTriangle3DContainer (size_t const) const;
    void WriteSTLFile (std::string const& FileName) const;

    size_t GetNSTL () const;
    size_t GetNPoints () const;

    TTriangle3D GetPoint (size_t const i) const;

    void AddToPoint (size_t const i, double const Value);

    void ClearValues ();
    void Clear ();
    void Remove (std::string const&);
    

  private:
    std::vector<TTriangle3DContainer> fSTL;
    std::vector<size_t> fNPoints;
    std::map<std::string, size_t> fSTLMap;

};
























inline std::ostream& operator << (std::ostream& os, TSTLContainer const& o)
{
  // For easy printing
  os << "TSTLContainer has " << o.GetNSTL() << " STL Objects" << std::endl;

  size_t const N = o.GetNSTL();

  for (size_t i = 0; i != N; ++i) {
    TTriangle3DContainer T = o.GetTTriangle3DContainer(i);

    os << T << std::endl;
  }

  return os;
}










#endif
