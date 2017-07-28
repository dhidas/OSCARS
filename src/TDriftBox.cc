#include "TDriftBox.h"


TDriftBox::TDriftBox (TVector3D   const& Width,
                      TVector3D   const& Center,
                      TVector3D   const& Rotations,
                      std::string const& Name,
                      bool        const  RecordTrajectory) {

  fWidth = Width;
  fCenter = Center;
  fRotated = Rotations;
  this->SetName(Name);
  this->SetRecordTrajectory(RecordTrajectory);

  fIgnoreAxisX = false;
  fIgnoreAxisY = false;
  fIgnoreAxisZ = false;

  if (fWidth.GetX() <= 0) {
    fIgnoreAxisX = true;
  }
  if (fWidth.GetY() <= 0) {
    fIgnoreAxisY = true;
  }
  if (fWidth.GetZ() <= 0) {
    fIgnoreAxisZ = true;
  }

}

TDriftBox::~TDriftBox ()
{
}




bool TDriftBox::IsInside (TVector3D const& X) const
{
  // Translate back into box frame
  TVector3D XInBoxCoordinates = X;
  XInBoxCoordinates.RotateSelfXYZ(fRotated);

  // Position in the box frame with respect to the center
  TVector3D const RX = XInBoxCoordinates - fCenter;

  if ((!fIgnoreAxisX && fabs(RX.GetX()) > fabs(fWidth.GetX() / 2.)) || (!fIgnoreAxisY && fabs(RX.GetY()) > fabs(fWidth.GetY() / 2.)) || (!fIgnoreAxisZ && fabs(RX.GetZ()) > fabs(fWidth.GetZ() / 2.))) {
    return false;
  }

  return true;
}





TVector3D TDriftBox::GetWidth () const
{
  // Return the width
  return fWidth;
}




TVector3D TDriftBox::GetRotated () const
{
  // Return the rotations
  return fRotated;
}




TVector3D TDriftBox::GetCenter () const
{
  // Return the Center postion
  return fCenter;
}




void TDriftBox::Print (std::ostream& os) const
{
  os << *this << std::endl;
  return;
}
