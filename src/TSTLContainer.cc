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

  return;
}




TTriangle3DContainer const& TSTLContainer::GetTTriangle3DContainer (size_t const i) const
{
  if (i >= fSTL.size()) {
    throw std::out_of_range("TSTLContainer::GetTTriangle3DContainer index out of range");
  }

  return fSTL[i];
}





size_t TSTLContainer::GetNSTL () const
{
  return fSTL.size();
}

