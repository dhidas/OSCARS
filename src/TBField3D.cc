////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Fri Mar 18 17:19:58 EDT 2016
//
////////////////////////////////////////////////////////////////////

#include "TBField3D.h"

#include <fstream>
#include <sstream>

TBField3D::TBField3D ()
{
}




TBField3D::~TBField3D ()
{
}




double TBField3D::GetBx(double const X, double const Y, double const Z) const
{
  return 0.0;
}




double TBField3D::GetBy(double const X, double const Y, double const Z) const
{
  return 0.0;
}




double TBField3D::GetBz(double const X, double const Y, double const Z) const
{
  return 0.0;
}




void TBField3D::ReadFile (std::string const& InFileName, bool const HasX, bool const HasY, bool const HasZ, bool const HasBX, bool const HasBY, bool const HasBZ)
{
  fHasX  = HasX;
  fHasY  = HasY;
  fHasZ  = HasZ;
  fHasBX = HasBX;
  fHasBY = HasBY;
  fHasBZ = HasBZ;

  int const DIM = HasX + HasY + HasZ + HasBX + HasBY + HasBZ;


  std::vector<double> D(DIM, 0);


  std::ifstream fi(InFileName);
  if (!fi) {
    throw;
  }

  std::istringstream S;
  for (std::string Line; std::getline(fi, Line); ) {
    S.str("");
    S.str(Line);

    for (size_t i = 0; i != DIM; ++i) {
      S >> D[i];
    }

    fData.push_back(D);

  }

  return;
}











//bool TBField3D::operator < (std::vector<double>& A, std::vector<double>& B)
//{
//  return true;
//}
