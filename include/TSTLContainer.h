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
    

  private:
    std::vector<TTriangle3DContainer> fSTL;
    std::map<std::string, size_t> fSTLMap;

};































#endif
