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

  // Create mutex locking structure
  fLock_mutex = new std::mutex();
}




TParticleTrajectoryPoints::TParticleTrajectoryPoints (const TParticleTrajectoryPoints& TPTP)
{
  // Copy constructor

  // Create mutex locking structure
  fLock_mutex = new std::mutex();
}



TParticleTrajectoryPoints::TParticleTrajectoryPoints (double const dt)
{
  // Default constructor

  fDeltaT = dt;

  // Create mutex locking structure
  fLock_mutex = new std::mutex();
}



TParticleTrajectoryPoints::~TParticleTrajectoryPoints ()
{
  // I own pointers in my own vectors/arrays

  this->Clear();

  // Delete locking if it exists
  if (fLock_mutex != 0x0) {
    delete fLock_mutex;
    //fLock_mutex=NULL;
  }

}




TParticleTrajectoryPoint const& TParticleTrajectoryPoints::GetPoint (size_t const i) const
{
  return fP[i];
}





TVector3D const& TParticleTrajectoryPoints::GetX (size_t const i) const
{
  return fP[i].GetX();
}





TVector3D const& TParticleTrajectoryPoints::GetB (size_t const i) const
{
  return fP[i].GetB();
}



TVector3D TParticleTrajectoryPoints::GetV (size_t const i) const
{
  return this->GetB(i) * TOSCARSSR::C();
}





TVector3D const& TParticleTrajectoryPoints::GetAoverC (size_t const i) const
{
  return fP[i].GetAoverC();
}





TVector3D TParticleTrajectoryPoints::GetA (size_t const i) const
{
  return this->GetAoverC(i) * TOSCARSSR::C();
}




double TParticleTrajectoryPoints::GetT (size_t const i) const
{
  return fT[i];
}




double TParticleTrajectoryPoints::GetTStart () const
{
  return fT.front();
}




double TParticleTrajectoryPoints::GetTStop () const
{
  return fT.back();
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
  return fP.size();
}




std::vector<TParticleTrajectoryPoint> const& TParticleTrajectoryPoints::GetTrajectory() const
{
  return fP;
}




std::vector<double> const& TParticleTrajectoryPoints::GetTimePoints () const
{
  return fT;
}




void TParticleTrajectoryPoints::AddPoint (TParticleTrajectoryPoint const& P, double const T)
{
  fP.push_back(P);
  fT.push_back(T);

  return;
}




void TParticleTrajectoryPoints::AddPoint (TVector3D const& X, TVector3D const& B, TVector3D const& AoverC, double const T)
{
  fP.push_back( TParticleTrajectoryPoint(X, B, AoverC) );
  fT.push_back(T);

  return;
}




void TParticleTrajectoryPoints::AddPoint (double const X1, double const X2, double const X3, double const B1, double const B2, double const B3, double const AoverC1, double const AoverC2, double const AoverC3, double const T)
{
  this->AddPoint( TVector3D(X1, X2, X3), TVector3D(B1, B2, B3), TVector3D(AoverC1, AoverC2, AoverC3), T );

  return;
}




void TParticleTrajectoryPoints::Reserve (size_t const n)
{
  fP.reserve(n);
  fT.reserve(n);
  return;
}




void TParticleTrajectoryPoints::ReverseArrays ()
{
  std::reverse(fP.begin(), fP.end());
  std::reverse(fT.begin(), fT.end());

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
  for (size_t i = 0; i != fP.size(); ++i) {
    f << fT[i]
      << " "
      << fP[i].GetX().GetX()
      << " "
      << fP[i].GetX().GetY()
      << " "
      << fP[i].GetX().GetZ()
      << " "
      << fP[i].GetB().GetX()
      << " "
      << fP[i].GetB().GetY()
      << " "
      << fP[i].GetB().GetZ()
      << " "
      << fP[i].GetAoverC().GetX()
      << " "
      << fP[i].GetAoverC().GetY()
      << " "
      << fP[i].GetAoverC().GetZ()
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

  // For writing data
  float v = 0;

  // Loop over all points and print to file
  for (size_t i = 0; i != fP.size(); ++i) {

    v = (float) fT[i];
    f.write((char*) &v, sizeof(float));
    v = (float) fP[i].GetX().GetX();
    f.write((char*) &v, sizeof(float));
    v = (float) fP[i].GetX().GetY();
    f.write((char*) &v, sizeof(float));
    v = (float) fP[i].GetX().GetZ();
    f.write((char*) &v, sizeof(float));
    v = (float) fP[i].GetB().GetX();
    f.write((char*) &v, sizeof(float));
    v = (float) fP[i].GetB().GetY();
    f.write((char*) &v, sizeof(float));
    v = (float) fP[i].GetB().GetZ();
    f.write((char*) &v, sizeof(float));
    v = (float) fP[i].GetAoverC().GetX();
    f.write((char*) &v, sizeof(float));
    v = (float) fP[i].GetAoverC().GetY();
    f.write((char*) &v, sizeof(float));
    v = (float) fP[i].GetAoverC().GetZ();
    f.write((char*) &v, sizeof(float));
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
  double aocx, aocy, aocz;

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

    LineStream >> x >> y >> z >> bx >> by >> bz >> aocx >> aocy >> aocz;

    if (!LineStream.fail()) { // check for error
      this->AddPoint( TVector3D(x, y, z), TVector3D(bx, by, bz), TVector3D(aocx, aocy, aocz), t );
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
  float   t  = 0;
  float   x  = 0;
  float   y  = 0;
  float   z  = 0;
  float   bx = 0;
  float   by = 0;
  float   bz = 0;
  float aocx = 0;
  float aocy = 0;
  float aocz = 0;

  while (!f.eof()) {

    // Read data
    f.read( (char*)    &t, sizeof(float));
    f.read( (char*)    &x, sizeof(float));
    f.read( (char*)    &y, sizeof(float));
    f.read( (char*)    &z, sizeof(float));
    f.read( (char*)   &bx, sizeof(float));
    f.read( (char*)   &by, sizeof(float));
    f.read( (char*)   &bz, sizeof(float));
    f.read( (char*) &aocx, sizeof(float));
    f.read( (char*) &aocy, sizeof(float));
    f.read( (char*) &aocz, sizeof(float));

    // If end of file, stop, otherwise add to trajectory
    if (f.eof()) {
      break;
    } else {
      this->AddPoint( TVector3D(x, y, z), TVector3D(bx, by, bz), TVector3D(aocx, aocy, aocz), t );
    }
  }

  return;
}




void TParticleTrajectoryPoints::Lock ()
{
  // Lock mutex for trajectory writing, test, etc
  fLock_mutex->lock();
  return;
}




void TParticleTrajectoryPoints::UnLock ()
{
  // unlock mutex
  fLock_mutex->unlock();
  return;
}




void TParticleTrajectoryPoints::Clear ()
{
  // Clear all vectors and free all memory that this object owns

  fP.clear();
  fT.clear();

  return;
}

