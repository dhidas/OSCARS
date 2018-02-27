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




void TParticleTrajectoryPoints::WriteToFile (std::string const& FileName, std::string const& FormatIn) const
{
  // Write the trajectory data to a file in text format

  // Open file for writing
  std::ofstream f(FileName.c_str());
  if (!f.is_open()) {
    throw std::ifstream::failure("cannot open output file");
  }

  std::string Format = FormatIn;
  std::transform(Format.begin(), Format.end(), Format.begin(), ::toupper);

  // Write in default format
  if (Format == "DEFAULT") {
    Format = "T X Y Z BX BY BZ AX AY AZ";

    // Write header
    f << "# " << Format << std::endl;

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
        << fP[i].GetAoverC().GetX() * TOSCARSSR::C()
        << " "
        << fP[i].GetAoverC().GetY() * TOSCARSSR::C()
        << " "
        << fP[i].GetAoverC().GetZ() * TOSCARSSR::C()
        << std::endl;
    }
  } else {

    // String stream
    std::istringstream FormatStream(Format);
    std::vector<std::string> FormatWords;
    for (std::string ss; FormatStream >> ss; ) {
      FormatWords.push_back(ss.c_str());
    }

    if (FormatWords.size() < 1) {
      throw std::length_error("Format must contain at least one element");
    }

    // Write header
    f << "# " << FormatIn << std::endl;

    // Set the ofstream in scientific output mode
    f << std::scientific;
    f.precision(35);

    // Loop over all points and print to file
    for (size_t i = 0; i != fP.size(); ++i) {

      for (std::vector<std::string>::iterator it = FormatWords.begin(); it != FormatWords.end(); ++it) {
        
        // Which component to write
        if (*it == "T") {
          f << fT[i];
        } else if (*it == "X") {
          f << fP[i].GetX().GetX();
        } else if (*it == "Y") {
          f << fP[i].GetX().GetY();
        } else if (*it == "Z") {
          f << fP[i].GetX().GetZ();
        } else if (*it == "BX") {
          f << fP[i].GetB().GetX();
        } else if (*it == "BY") {
          f << fP[i].GetB().GetY();
        } else if (*it == "BZ") {
          f << fP[i].GetB().GetZ();
        } else if (*it == "AX") {
          f << fP[i].GetAoverC().GetX() * TOSCARSSR::C();
        } else if (*it == "AY") {
          f << fP[i].GetAoverC().GetY() * TOSCARSSR::C();
        } else if (*it == "AZ") {
          f << fP[i].GetAoverC().GetZ() * TOSCARSSR::C();
        } else {
          throw std::invalid_argument("format specifier not recognized");
        }

        if (it < FormatWords.end() - 1) {
          f << " ";
        }
      }
      f << std::endl;
    }
  }

  // Close file
  f.close();

  return;
}




void TParticleTrajectoryPoints::WriteToFileBinary (std::string const& FileName, std::string const & FormatIn) const
{
  // Write the trajectory data to a file in text format

  // Open file for writing
  std::ofstream f(FileName.c_str(), std::ios::binary);
  if (!f.is_open()) {
    throw std::ifstream::failure("cannot open output file");
  }


  // Transform input format
  std::string Format = FormatIn;
  std::transform(Format.begin(), Format.end(), Format.begin(), ::toupper);


  // Write in default format
  if (Format == "DEFAULT") {
    Format = "T X Y Z BX BY BZ AX AY AZ";

    // Length of format string to be stored in binary file as char(s)
    int const FormatLength = (int) Format.size();

    // Write binary header
    f.write((char*) &FormatLength, sizeof(int));
    f.write((char*) Format.c_str(), FormatLength * sizeof(char));

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
      v = (float) fP[i].GetAoverC().GetX() * TOSCARSSR::C();
      f.write((char*) &v, sizeof(float));
      v = (float) fP[i].GetAoverC().GetY() * TOSCARSSR::C();
      f.write((char*) &v, sizeof(float));
      v = (float) fP[i].GetAoverC().GetZ() * TOSCARSSR::C();
      f.write((char*) &v, sizeof(float));
    }
  } else {

    // String stream
    std::istringstream FormatStream(Format);
    std::vector<std::string> FormatWords;
    for (std::string ss; FormatStream >> ss; ) {
      FormatWords.push_back(ss.c_str());
    }

    if (FormatWords.size() < 1) {
      throw std::length_error("Format must contain at least one element");
    }

    // Length of format string to be stored in binary file as char(s)
    int const FormatLength = (int) Format.size();

    // Write binary header
    f.write((char*) &FormatLength, sizeof(int));
    f.write((char*) Format.c_str(), FormatLength * sizeof(char));

    // For writing data
    float v = 0;


    // Loop over all points and print to file
    for (size_t i = 0; i != fP.size(); ++i) {

      for (std::vector<std::string>::iterator it = FormatWords.begin(); it != FormatWords.end(); ++it) {
        
        // Which component to write
        if (*it == "T") {
          v = (float) fT[i];
        } else if (*it == "X") {
          v = (float) fP[i].GetX().GetX();
        } else if (*it == "Y") {
          v = (float) fP[i].GetX().GetY();
        } else if (*it == "Z") {
          v = (float) fP[i].GetX().GetZ();
        } else if (*it == "BX") {
          v = (float) fP[i].GetB().GetX();
        } else if (*it == "BY") {
          v = (float) fP[i].GetB().GetY();
        } else if (*it == "BZ") {
          v = (float) fP[i].GetB().GetZ();
        } else if (*it == "AX") {
          v = (float) fP[i].GetAoverC().GetX() * TOSCARSSR::C();
        } else if (*it == "AY") {
          v = (float) fP[i].GetAoverC().GetY() * TOSCARSSR::C();
        } else if (*it == "AZ") {
          v = (float) fP[i].GetAoverC().GetZ() * TOSCARSSR::C();
        } else {
          throw std::invalid_argument("format specifier not recognized");
        }
        f.write((char*) &v, sizeof(float));

      }
    }
  }

  // Close file
  f.close();

  return;
}




void TParticleTrajectoryPoints::ReadFromFile (std::string const& FileName)
{
  // Read a text file of the trajectory
  // UPDATE: Read Text File - Add format
  throw;

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

  // Read format length as first word
  int FormatLength = 0;
  f.read((char*) &FormatLength, sizeof(int));

  // Check format length
  if (FormatLength < 1) {
    throw;
  }

  // Read input format
  char* FormatIn = new char[FormatLength+1];

  FormatIn[FormatLength] = '\0';
  f.read((char*) FormatIn, FormatLength * sizeof(char));

  // Transform input format
  std::string Format = FormatIn;
  std::transform(Format.begin(), Format.end(), Format.begin(), ::toupper);

  delete FormatIn;

  // Ordered vector of inputs
  std::istringstream FormatStream(Format);
  std::vector<std::string> FormatWords;
  for (std::string ss; FormatStream >> ss; ) {
    FormatWords.push_back(ss.c_str());
  }

  if (FormatWords.size() < 1) {
    throw std::length_error("Format must contain at least one element");
  }


  // For data reading
  float* v = new float[FormatWords.size()];

  // Variables to read from file
  float   t  = 0;
  float   x  = 0;
  float   y  = 0;
  float   z  = 0;
  float   bx = 0;
  float   by = 0;
  float   bz = 0;
  float   ax = 0;
  float   ay = 0;
  float   az = 0;

  int   it  = -1;
  int   ix  = -1;
  int   iy  = -1;
  int   iz  = -1;
  int   ibx = -1;
  int   iby = -1;
  int   ibz = -1;
  int   iax = -1;
  int   iay = -1;
  int   iaz = -1;

  for (size_t i = 0; i != FormatWords.size(); ++i) {
    std::string const& Word = FormatWords[i];

    // Which component index is which variable
    if (Word == "T") {
      it = i;
    } else if (Word == "X") {
      ix = i;
    } else if (Word == "Y") {
      iy = i;
    } else if (Word == "Z") {
      iz = i;
    } else if (Word == "BX") {
      ibx = i;
    } else if (Word == "BY") {
      iby = i;
    } else if (Word == "BZ") {
      ibz = i;
    } else if (Word == "AX") {
      iax = i;
    } else if (Word == "AY") {
      iay = i;
    } else if (Word == "AZ") {
      iaz = i;
    } else {
      throw std::invalid_argument("format specifier not recognized");
    }
  }

  while (!f.eof()) {

    for (size_t i = 0; i != FormatWords.size(); ++i) {
      f.read( (char*)    &v[i], sizeof(float));
    }

    if (it != -1) { t = v[it]; }
    if (ix != -1) { x = v[ix]; }
    if (iy != -1) { y = v[iy]; }
    if (iz != -1) { z = v[iz]; }
    if (ibx != -1) { bx = v[ibx]; }
    if (iby != -1) { by = v[iby]; }
    if (ibz != -1) { bz = v[ibz]; }
    if (iax != -1) { ax = v[iax]; }
    if (iay != -1) { ay = v[iay]; }
    if (iaz != -1) { az = v[iaz]; }

    // If end of file, stop, otherwise add to trajectory
    if (f.eof()) {
      break;
    } else {
      this->AddPoint( TVector3D(x, y, z), TVector3D(bx, by, bz), TVector3D(ax, ay, az) * TOSCARSSR::C(), t );
    }
  }

  // Delete data vector
  delete v;

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

