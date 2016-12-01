#ifndef GUARD_TBField3D_h
#define GUARD_TBField3D_h
////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Fri Mar 18 17:19:58 EDT 2016
//
////////////////////////////////////////////////////////////////////

#include "TBField.h"

#include <string>
#include <vector>

class TBField3D : public TBField
{

  public:
    TBField3D ();
    ~TBField3D ();

    double GetBx (double const, double const, double const) const;
    double GetBy (double const, double const, double const) const;
    double GetBz (double const, double const, double const) const;

    void ReadFile (std::string const&, bool const, bool const, bool const, bool const, bool const, bool const);

    // Regularize

    //bool operator < (std::vector<double>&, std::vector<double>&);

  private:
    bool fHasX;
    bool fHasY;
    bool fHasZ;
    bool fHasBX;
    bool fHasBY;
    bool fHasBZ;
    std::vector< std::vector<double> > fData;

};



















#endif

