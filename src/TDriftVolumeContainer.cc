////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Thu Jul 27 12:16:05 EDT 2017
//
////////////////////////////////////////////////////////////////////

#include "TDriftVolumeContainer.h"

#include <iostream>
#include <algorithm>



TDriftVolumeContainer::TDriftVolumeContainer ()
{
  // Default constructor
}



TDriftVolumeContainer::~TDriftVolumeContainer ()
{
  // Default destructor.  I own everything you have passed me.  Make no mistake there!!
  this->Clear();
}



bool TDriftVolumeContainer::IsInside (TVector3D const& X)
{
  for (std::vector<TDriftVolume*>::iterator it = fDriftVolumes.begin(); it != fDriftVolumes.end(); ++it) {
    if ((*it)->IsInside(X)) {
      if (it != fDriftVolumes.begin()) {
        // If this is not first let's put it first to save time
        std::iter_swap(fDriftVolumes.begin(), it);
      }
      return true;
    }
  }
  return false;

}





void TDriftVolumeContainer::AddDriftVolume (TDriftVolume* DV)
{
  // Construct me with a field
  fDriftVolumes.push_back(DV);
}




void TDriftVolumeContainer::RemoveDriftVolume (std::string const& Name)
{
  // Remove all fields that match the input name

  size_t i = 0;
  while (i < fDriftVolumes.size()) {
    if (fDriftVolumes[i]->GetName() == Name) {

      // Delete TDriftVolume
      delete fDriftVolumes[i];

      // Remove pointer from vector
      fDriftVolumes.erase( fDriftVolumes.begin() + i );
    } else {
      ++i;
    }
  }
  return;
}







TDriftVolume const& TDriftVolumeContainer::GetDriftVolume (size_t const i) const
{
  // Return reference to drift volume (const)
  return *fDriftVolumes[i];
}




size_t TDriftVolumeContainer::GetNDriftVolumes () const
{
  // Return the number of drift volumes input
  return fDriftVolumes.size();
}




void TDriftVolumeContainer::Clear ()
{
  for (std::vector<TDriftVolume*>::iterator it = fDriftVolumes.begin(); it != fDriftVolumes.end(); ++it) {
    if (*it != 0x0) {
      delete *it;
    }
  }

  fDriftVolumes.clear();

  return;
}





