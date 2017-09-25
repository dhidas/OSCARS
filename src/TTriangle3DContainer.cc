#include "TTriangle3DContainer.h"


TTriangle3DContainer::TTriangle3DContainer ()
{
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
  return;
}



void TTriangle3DContainer::RotateSelfXYZ (TVector3D const& R)
{
  for (std::vector<TTriangle3D>::iterator it = fT.begin(); it != fT.end(); ++it) {
    it->RotateSelfXYZ(R);
  }
  return;
}




void TTriangle3DContainer::TranslateSelf (TVector3D const& T)
{
  for (std::vector<TTriangle3D>::iterator it = fT.begin(); it != fT.end(); ++it) {
    *it += T;
  }
  return;
}



void TTriangle3DContainer::ReadSTLFile (std::string const& FileName)
{
  // Read an STL file and add all points to self

  // Open the input file
  std::ifstream fi(FileName.c_str(), std::ios::binary);
  if (!fi.is_open()) {
    throw;
  }

  char H[80];
  fi.read(H, 80 * sizeof(char));

  unsigned int NTriangles;
  fi.read((char*) &NTriangles, sizeof(int));

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

    fT.push_back(TTriangle3D(A[0]/1000., A[1]/1000., A[2]/1000.,
                             B[0]/1000., B[1]/1000., B[2]/1000.,
                             C[0]/1000., C[1]/1000., C[2]/1000.,
                             N[0], N[1], N[2]));


  }

  return;
}
