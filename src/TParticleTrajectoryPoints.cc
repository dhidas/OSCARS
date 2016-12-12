////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Wed May 18 18:08:30 EDT 2016
//
////////////////////////////////////////////////////////////////////


#include "TParticleTrajectoryPoints.h"

#include "TOSCARSSR.h"

#include <algorithm>
#include <fstream>
#include <sstream>


TParticleTrajectoryPoints::TParticleTrajectoryPoints ()
{
  // Default constructor

  // Default DeltaT for this mode to zero
  fDeltaT = 0;
}



TParticleTrajectoryPoints::TParticleTrajectoryPoints (double const dt)
{
  // Default constructor

  fDeltaT = dt;
}



TParticleTrajectoryPoints::~TParticleTrajectoryPoints ()
{
  // I own pointers in my own vectors/arrays

  this->Clear();

}




TVector3D const& TParticleTrajectoryPoints::GetX (size_t const i) const
{
  return fX[i];
}





TVector3D const& TParticleTrajectoryPoints::GetB (size_t const i) const
{
  return fB[i];
}



TVector3D TParticleTrajectoryPoints::GetV (size_t const i) const
{
  return this->GetB(i) * TOSCARSSR::C();
}





TVector3D const& TParticleTrajectoryPoints::GetAoverC (size_t const i) const
{
  return fAoverC[i];
}





TVector3D TParticleTrajectoryPoints::GetA (size_t const i) const
{
  return this->GetAoverC(i) * TOSCARSSR::C();
}




double TParticleTrajectoryPoints::GetDeltaT () const
{
  return fDeltaT;
}




void TParticleTrajectoryPoints::SetDeltaT (double const DT)
{
  fDeltaT = DT;
  return;
}




size_t TParticleTrajectoryPoints::GetNPoints () const
{
  return fX.size();
}



void TParticleTrajectoryPoints::AddPoint (TVector3D const& X, TVector3D const& B, TVector3D const& AoverC, double const T)
{
  fX.push_back( TVector3D(X) );
  fB.push_back( TVector3D(B) );
  fAoverC.push_back( TVector3D(AoverC) );

  if (fDeltaT != 0) {
    fT.push_back(T);
  }

  return;
}



void TParticleTrajectoryPoints::AddPoint (double const X1, double const X2, double const X3, double const B1, double const B2, double const B3, double const AoverC1, double const AoverC2, double const AoverC3, double const T)
{
  this->AddPoint(TVector3D(X1, X2, X3), TVector3D(B1, B2, B3), TVector3D(AoverC1, AoverC2, AoverC3), T);

  return;
}




void TParticleTrajectoryPoints::ReverseArrays ()
{
  std::reverse(fX.begin(), fX.end());
  std::reverse(fB.begin(), fB.end());
  std::reverse(fAoverC.begin(), fAoverC.end());

  return;
}




void TParticleTrajectoryPoints::WriteToFile (std::string const& FileName) const
{
  // Write the trajectory data to a file in text format

  // Open file for writing
  std::ofstream f(FileName.c_str());
  if (!f.is_open()) {
    throw;
  }

  // Write header
  f << "# T X Y Z BX BY BZ" << std::endl;

  // Set the ofstream in scientific output mode
  f << std::scientific;
  f.precision(35);

  // Loop over all points and print to file
  for (size_t i = 0; i != fB.size(); ++i) {
    f << fDeltaT * (double) i
      << " "
      << fX[i].GetX()
      << " "
      << fX[i].GetY()
      << " "
      << fX[i].GetZ()
      << " "
      << fB[i].GetX()
      << " "
      << fB[i].GetY()
      << " "
      << fB[i].GetZ()
      << std::endl;
  }

  // Close file
  f.close();

  return;
}




void TParticleTrajectoryPoints::WriteToFileBinary (std::string const& FileName) const
{
  // Write the trajectory data to a file in text format

  // Open file for writing
  std::ofstream f(FileName.c_str(), std::ios::binary);
  if (!f.is_open()) {
    throw;
  }

  // Loop over all points and print to file
  for (size_t i = 0; i != fB.size(); ++i) {
    double value;

    value = fDeltaT;
    f.write((char*) &value, sizeof(double));
    value = fX[i].GetX();
    f.write((char*) &value, sizeof(double));
    value = fX[i].GetY();
    f.write((char*) &value, sizeof(double));
    value = fX[i].GetZ();
    f.write((char*) &value, sizeof(double));
    value = fB[i].GetX();
    f.write((char*) &value, sizeof(double));
    value = fB[i].GetY();
    f.write((char*) &value, sizeof(double));
    value = fB[i].GetZ();
    f.write((char*) &value, sizeof(double));
  }

  // Close file
  f.close();

  return;
}




void TParticleTrajectoryPoints::ReadFromFile (std::string const& FileName)
{
  // Read a text file of the trajectory

  // Open the input file for reading
  std::ifstream f(FileName.c_str());
  if (!f.is_open()) {
    throw;
  }

  // Variables to read from file
  double t;
  double x, y, z;
  double bx, by, bz;

  size_t n = 0;

  // Loop over all lines in file, skip anything that starts with somethign which is not
  // a floating point number
  std::istringstream LineStream;
  for (std::string Line; std::getline(f, Line); ++n) {
    if (f.eof()) {
      break;
    }

    LineStream.clear();
    LineStream.str(Line);

    LineStream >> t;
    if (LineStream.fail()) {
      // Presumably this is a comment line, although it could be
      // a data formatting error, but I'll ignore that for now.
      continue;
    }

    LineStream >> x >> y >> z >> bx >> by >> bz;

    if (!LineStream.fail()) { // check for error
      fX.push_back( TVector3D(x, y, z) );
      fB.push_back( TVector3D(bx, by, bz) );
      fAoverC.push_back( TVector3D(0, 0, 0) );
    } else {
      throw;
    }

  }

  return;
}




void TParticleTrajectoryPoints::ReadFromFileBinary (std::string const& FileName)
{
  // Read a text file of the trajectory

  // Open the input file for reading
  std::ifstream f(FileName, std::ios::binary);
  if (!f.is_open()) {
    throw;
  }

  // Variables to read from file
  double t;
  double x, y, z;
  double bx, by, bz;

  while (!f.eof()) {

    f.read( (char*)  &t, sizeof(double));
    f.read( (char*)  &x, sizeof(double));
    f.read( (char*)  &y, sizeof(double));
    f.read( (char*)  &z, sizeof(double));
    f.read( (char*) &bx, sizeof(double));
    f.read( (char*) &by, sizeof(double));
    f.read( (char*) &bz, sizeof(double));


    if (f.eof()) {
      break;
    } else {
      fX.push_back( TVector3D(x, y, z) );
      fB.push_back( TVector3D(bx, by, bz) );
      fAoverC.push_back( TVector3D(0, 0, 0) );
    }
  }

  return;
}




void TParticleTrajectoryPoints::Clear ()
{
  // Clear all vectors and free all memory that this object owns

  fX.clear();
  fB.clear();
  fAoverC.clear();

  return;
}

