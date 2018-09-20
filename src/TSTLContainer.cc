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
  fSTL.back().ReadSTLFile(InFileName, Scale);
  fSTL.back().RotateSelfXYZ(Rotations);
  fSTL.back().TranslateSelf(Translation);

  fSTLMap[MyName] = fSTL.size() - 1;

  return;
}
