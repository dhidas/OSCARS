#include "TSTLContainer.h"


TSTLContainer::TSTLContainer ()
{
  // Constructor
}



TSTLContainer::~TSTLContainer ()
{
  // Destruction!!!
}



void TSTLContainer::AddFile (std::string const& InFileName,
                             double      const  Scale,
                             TVector3D   const& Rotations,
                             TVector3D   const& Translation,
                             std::string const& Name)
{
  std::string const MyName = Name != "" ? Name : "_stl" + std::to_string(fSTL.size());

  if (fSTLMap.count(MyName) != 0) {
    std::cerr << "fSTLMap.count(Name) != 0" << std::endl;
    throw std::invalid_argument("stl with this name already exists");
  }

  fSTL.push_back( TTriangle3DContainer() );
  fSTL.back().ReadSTLFile(InFileName, Scale, Rotations, Translation, MyName);

  fSTLMap[MyName] = fSTL.size() - 1;

  fNPoints.push_back(fSTL.back().GetNPoints());

  return;
}




TTriangle3DContainer const& TSTLContainer::GetTTriangle3DContainer (size_t const i) const
{
  if (i >= fSTL.size()) {
    throw std::out_of_range("TSTLContainer::GetTTriangle3DContainer index out of range");
  }

  return fSTL[i];
}





void TSTLContainer::WriteSTLFile (std::string const& FileName) const
{
  // Write an STL file and add all points to self

  // Open the output file
  std::ofstream fi(FileName.c_str(), std::ios::binary);
  if (!fi.is_open()) {
    throw std::ofstream::failure("output file cannot be opened");
  }

  // Get total number of triangles
  uint32_t NTrianglesTotal = 0;
  for (size_t itt = 0; itt != fSTL.size(); ++itt) {
    NTrianglesTotal += this->GetTTriangle3DContainer(itt).GetNPoints();
  }

  char H[80];
  sprintf(H, "%s", "OSCARS OSCARS OSCARS OSCARS OSCARS   ");
  fi.write(H, 80 * sizeof(char));

  fi.write((char*) &NTrianglesTotal, sizeof(uint32_t));

  float N[3];
  float A[3];
  float B[3];
  float C[3];
  short S = 0;

  for (size_t itt = 0; itt != fSTL.size(); ++itt) {
    TTriangle3DContainer TT = this->GetTTriangle3DContainer(itt);
    size_t const NTriangles = TT.GetNPoints();

    double const Scale = TT.GetScale();
    for (size_t i = 0; i < NTriangles; ++i) {
      TTriangle3D T = TT.GetPoint(i);

      for (int j = 0; j != 3; ++j) {
        A[j] = T[0][j] * Scale;
        B[j] = T[1][j] * Scale;
        C[j] = T[2][j] * Scale;
        N[j] = T[3][j];
      }

      fi.write((char*)  N, 3 * sizeof(float));
      fi.write((char*)  A, 3 * sizeof(float));
      fi.write((char*)  B, 3 * sizeof(float));
      fi.write((char*)  C, 3 * sizeof(float));
      fi.write((char*) &S ,    sizeof(short));
    }
  }

  fi.close();

  return;
}




size_t TSTLContainer::GetNSTL () const
{
  return fSTL.size();
}





size_t TSTLContainer::GetNPoints () const
{
  // Return total number of points (triangles) from all STL objects
  size_t N = 0;
  for (std::vector<size_t>::const_iterator it = fNPoints.begin(); it != fNPoints.end(); ++it){
    N += *it;
  }
  return N;
}





TTriangle3D TSTLContainer::GetPoint (size_t const i) const
{
  size_t Sum = 0;
  for (size_t is = 0; is != fNPoints.size(); ++is) {
    if (i < Sum + fNPoints[is]) {
      return fSTL[is].GetPoint(i - Sum);
    }
    Sum += fNPoints[is];
  }

  throw std::length_error("TSTLContainer::GetPoint request is out of range");
}





void TSTLContainer::AddToPoint (size_t const i, double const Value)
{
  size_t Sum = 0;
  for (size_t is = 0; is != fNPoints.size(); ++is) {
    if (i < Sum + fNPoints[is]) {
      fSTL[is].AddToPoint(i - Sum, Value);
      return;
    }
    Sum += fNPoints[is];
  }

  throw std::length_error("TSTLContainer::AddToPoint request is out of range");
  return;
}





void TSTLContainer::ClearValues ()
{
  // Clear all values from TTriangle3DContainer objects in vector
  for (std::vector<TTriangle3DContainer>::iterator it = fSTL.begin(); it != fSTL.end(); ++it){
    it->ClearValues();
  }

  return;
}




void TSTLContainer::Clear ()
{
  fSTL.clear();
  fNPoints.clear();
  fSTLMap.clear();

  return;
}




void TSTLContainer::Remove (std::string const& Name)
{
  // Remove field that matches the input name

  while (fSTLMap.count(Name) != 0) {
    size_t const i = fSTLMap[Name];

    fSTL.erase(fSTL.begin() + i);
    fNPoints.erase(fNPoints.begin() + i);
    fSTLMap.erase(Name);
  }

  return;
}
