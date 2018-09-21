#include "TTriangle3DContainer.h"

#include <cstdint>

TTriangle3DContainer::TTriangle3DContainer ()
{
  fScale = 1;
}




TTriangle3DContainer::~TTriangle3DContainer ()
{
}




void TTriangle3DContainer::Add (TTriangle3D const& T)
{
  fT.push_back(T);
  return;
}



size_t TTriangle3DContainer::GetNPoints () const
{
  return fT.size();
}




TTriangle3D TTriangle3DContainer::GetPoint (size_t const i) const
{
  return fT[i];
}


void TTriangle3DContainer::Clear ()
{
  fT.clear();
  fScale = 1;
  return;
}



void TTriangle3DContainer::ClearValues ()
{
  // Set all values to zero
  for (std::vector<TTriangle3D>::iterator it = fT.begin(); it != fT.end(); ++it) {
    it->SetValue(0);
  }
  return;
}




void TTriangle3DContainer::AddToPoint (size_t const i, double const V)
{
  // Compensated sum for adding to points

  // Check that the point is within range
  if (i >= fT.size()) {
    throw std::length_error("TTriangle3DContainer::AddtoPoint index out of range");
  }

  fT[i].AddToValue(V);

  return;
}








void TTriangle3DContainer::RotateSelfXYZ (TVector3D const& R)
{
  for (std::vector<TTriangle3D>::iterator it = fT.begin(); it != fT.end(); ++it) {
    it->RotateSelfXYZ(R);
  }

  fRotations = R;
  return;
}




void TTriangle3DContainer::TranslateSelf (TVector3D const& T)
{
  for (std::vector<TTriangle3D>::iterator it = fT.begin(); it != fT.end(); ++it) {
    *it += T;
  }

  fTranslation += T;
  return;
}



void TTriangle3DContainer::ReadSTLFile (std::string const& FileName,
                                        double const Scale,
                                        TVector3D const& Rotations,
                                        TVector3D const& Translation,
                                        std::string const& Name)
{
  // Read an STL file and add all points to self

  // Open the input file
  std::ifstream fi(FileName.c_str(), std::ios::binary);
  if (!fi.is_open()) {
    throw std::ifstream::failure("input file cannot be opened");
  }

  char H[80];
  fi.read(H, 80 * sizeof(char));

  uint32_t NTriangles;
  fi.read((char*) &NTriangles, sizeof(uint32_t));

  float N[3];
  float A[3];
  float B[3];
  float C[3];
  short S;

  for (int i = 0; i < NTriangles; ++i) {
    fi.read((char*)  N, 3 * sizeof(float));
    fi.read((char*)  A, 3 * sizeof(float));
    fi.read((char*)  B, 3 * sizeof(float));
    fi.read((char*)  C, 3 * sizeof(float));
    fi.read((char*) &S ,    sizeof(short));

    fT.push_back(TTriangle3D(A[0] * Scale, A[1] * Scale, A[2] * Scale,
                             B[0] * Scale, B[1] * Scale, B[2] * Scale,
                             C[0] * Scale, C[1] * Scale, C[2] * Scale,
                             N[0], N[1], N[2]));


  }

  fi.close();

  fFileName = FileName;
  fScale = Scale;
  fRotations = Rotations;
  fTranslation = Translation;
  fName = Name;

  this->RotateSelfXYZ(fRotations);
  this->TranslateSelf(Translation);

  return;
}





void TTriangle3DContainer::WriteSTLFile (std::string const& FileName)
{
  // Write an STL file and add all points to self

  // Open the output file
  std::ofstream fi(FileName.c_str(), std::ios::binary);
  if (!fi.is_open()) {
    throw std::ofstream::failure("output file cannot be opened");
  }

  char H[80];
  sprintf(H, "%s", "OSCARS OSCARS OSCARS OSCARS OSCARS   ");
  fi.write(H, 80 * sizeof(char));

  uint32_t NTriangles = (uint32_t) fT.size();;
  fi.write((char*) &NTriangles, sizeof(uint32_t));

  float N[3];
  float A[3];
  float B[3];
  float C[3];
  short S = 0;

  for (int i = 0; i < NTriangles; ++i) {

    for (int j = 0; j != 3; ++j) {
      A[j] = fT[i][0][j] * fScale;
      B[j] = fT[i][1][j] * fScale;
      C[j] = fT[i][2][j] * fScale;
      N[j] = fT[i][3][j];
    }

    fi.write((char*)  N, 3 * sizeof(float));
    fi.write((char*)  A, 3 * sizeof(float));
    fi.write((char*)  B, 3 * sizeof(float));
    fi.write((char*)  C, 3 * sizeof(float));
    fi.write((char*) &S ,    sizeof(short));
  }

  fi.close();

  return;
}








std::string TTriangle3DContainer::GetFileName () const
{
  return fFileName;
}



double TTriangle3DContainer::GetScale () const
{
  return fScale;
}



TVector3D TTriangle3DContainer::GetRotations () const
{
  return fRotations;
}



TVector3D TTriangle3DContainer::GetTranslation () const
{
  return fTranslation;
}



std::string TTriangle3DContainer::GetName () const
{
  return fName;
}
