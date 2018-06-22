////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Wed May 18 18:08:30 EDT 2016
//
////////////////////////////////////////////////////////////////////


#include "TParticleTrajectoryPoints.h"

#include "TOSCARSSR.h"
#include "TOMATH.h"

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



void TParticleTrajectoryPoints::SetB (size_t const i, TVector3D const& B)
{
  fP[i].SetB(B);
}



TVector3D TParticleTrajectoryPoints::GetV (size_t const i) const
{
  return this->GetB(i) * TOSCARSSR::C();
}





TVector3D const& TParticleTrajectoryPoints::GetAoverC (size_t const i) const
{
  return fP[i].GetAoverC();
}





void TParticleTrajectoryPoints::SetAoverC (size_t const i, TVector3D const& AoverC)
{
  fP[i].SetAoverC(AoverC);
  return;
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
  if (Format == "") {
    Format = "T X Y Z BX BY BZ BPX BPY BPZ";

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
        << fP[i].GetAoverC().GetX()
        << " "
        << fP[i].GetAoverC().GetY()
        << " "
        << fP[i].GetAoverC().GetZ()
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
        } else if (*it == "BPX") {
          f << fP[i].GetAoverC().GetX();
        } else if (*it == "BPY") {
          f << fP[i].GetAoverC().GetY();
        } else if (*it == "BPZ") {
          f << fP[i].GetAoverC().GetZ();
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
  if (Format == "") {
    Format = "T X Y Z BX BY BZ BPX BPY BPZ";

    // Length of format string to be stored in binary file as char(s)
    int const FormatLength = (int) Format.size();

    // Write binary header
    f.write((char*) &FormatLength, sizeof(int));
    f.write((char*) Format.c_str(), FormatLength * sizeof(char));

    // For writing data
    double v = 0;

    // Loop over all points and print to file
    for (size_t i = 0; i != fP.size(); ++i) {
      v = (double) fT[i];
      f.write((char*) &v, sizeof(double));
      v = (double) fP[i].GetX().GetX();
      f.write((char*) &v, sizeof(double));
      v = (double) fP[i].GetX().GetY();
      f.write((char*) &v, sizeof(double));
      v = (double) fP[i].GetX().GetZ();
      f.write((char*) &v, sizeof(double));
      v = (double) fP[i].GetB().GetX();
      f.write((char*) &v, sizeof(double));
      v = (double) fP[i].GetB().GetY();
      f.write((char*) &v, sizeof(double));
      v = (double) fP[i].GetB().GetZ();
      f.write((char*) &v, sizeof(double));
      v = (double) fP[i].GetAoverC().GetX();
      f.write((char*) &v, sizeof(double));
      v = (double) fP[i].GetAoverC().GetY();
      f.write((char*) &v, sizeof(double));
      v = (double) fP[i].GetAoverC().GetZ();
      f.write((char*) &v, sizeof(double));
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
    double v = 0;


    // Loop over all points and print to file
    for (size_t i = 0; i != fP.size(); ++i) {

      for (std::vector<std::string>::iterator it = FormatWords.begin(); it != FormatWords.end(); ++it) {
        
        // Which component to write
        if (*it == "T") {
          v = (double) fT[i];
        } else if (*it == "X") {
          v = (double) fP[i].GetX().GetX();
        } else if (*it == "Y") {
          v = (double) fP[i].GetX().GetY();
        } else if (*it == "Z") {
          v = (double) fP[i].GetX().GetZ();
        } else if (*it == "BX") {
          v = (double) fP[i].GetB().GetX();
        } else if (*it == "BY") {
          v = (double) fP[i].GetB().GetY();
        } else if (*it == "BZ") {
          v = (double) fP[i].GetB().GetZ();
        } else if (*it == "BPX") {
          v = (double) fP[i].GetAoverC().GetX();
        } else if (*it == "BPY") {
          v = (double) fP[i].GetAoverC().GetY();
        } else if (*it == "BPZ") {
          v = (double) fP[i].GetAoverC().GetZ();
        } else {
          throw std::invalid_argument("format specifier not recognized");
        }
        f.write((char*) &v, sizeof(double));

      }
    }
  }

  // Close file
  f.close();

  return;
}




void TParticleTrajectoryPoints::ReadFromFileFormat (std::string const& FileName, std::string const& FormatIn)
{
  // Read a text file of the trajectory
  // UPDATE: THIS IS TO REDO FOR SPECIFIED FORMATS
  // Open the input file for reading
  std::ifstream f(FileName.c_str());
  if (!f.is_open()) {
    throw std::ifstream::failure("cannot open input file specified");
  }

  std::string HeaderLine;
  std::getline(f, HeaderLine);
  std::vector<std::string> HeaderVector;
  std::istringstream HeaderStream(HeaderLine);
  for (std::string a; HeaderStream >> a; ) {
    HeaderVector.push_back(a);
  }

  std::string FormatFound = "";

  if (HeaderVector.size() != 0 && std::string(HeaderVector[0], 0, 1) != "#") {
    f.seekg(0, f.beg);
  } else {
    if (HeaderVector[0].size() != 1) {
      HeaderVector[0] = std::string(HeaderVector[0], 1, HeaderVector[0].size()-1);
    } else {
      HeaderVector.erase(HeaderVector.begin());
    }
    bool const HasT   =  std::find(HeaderVector.begin(), HeaderVector.end(), "T") != HeaderVector.end();
    bool const HasX   =  std::find(HeaderVector.begin(), HeaderVector.end(), "X") != HeaderVector.end();
    bool const HasY   =  std::find(HeaderVector.begin(), HeaderVector.end(), "Y") != HeaderVector.end();
    bool const HasZ   =  std::find(HeaderVector.begin(), HeaderVector.end(), "Z") != HeaderVector.end();

    // Does it have the minimum to be a format specified
    if (HasT && (HasX || HasY || HasZ)) {
      for (size_t i = 0; i != HeaderVector.size(); ++i) {
        FormatFound += " " + HeaderVector[i];
      }
    } else {
      // I guess it was just a comment
    }
  }

  // Transform input format
  std::string Format = FormatIn;
  std::transform(Format.begin(), Format.end(), Format.begin(), ::toupper);
  if (FormatIn == "" && FormatFound != "") {
    Format = FormatFound;
  }
  std::transform(Format.begin(), Format.end(), Format.begin(), ::toupper);

  // Check and set format
  if (Format == "") {
    Format = "T X Y Z BX BY BZ BPX BPY BPZ";
  }

  // Figure out how many columns in the data input
  std::istringstream FormatStringStream(Format);
  std::vector<std::string> FormatVector;
  for (std::string one; FormatStringStream >> one; ) {
    FormatVector.push_back(one);
  }

  // Number of columns
  size_t const NColumns = FormatVector.size();


  // Check which components exist
  bool const HasT   =  std::find(FormatVector.begin(), FormatVector.end(), "T") != FormatVector.end();
  bool const HasX   =  std::find(FormatVector.begin(), FormatVector.end(), "X") != FormatVector.end();
  bool const HasY   =  std::find(FormatVector.begin(), FormatVector.end(), "Y") != FormatVector.end();
  bool const HasZ   =  std::find(FormatVector.begin(), FormatVector.end(), "Z") != FormatVector.end();
  bool const HasBX  =  std::find(FormatVector.begin(), FormatVector.end(), "BX") != FormatVector.end();
  bool const HasBY  =  std::find(FormatVector.begin(), FormatVector.end(), "BY") != FormatVector.end();
  bool const HasBZ  =  std::find(FormatVector.begin(), FormatVector.end(), "BZ") != FormatVector.end();
  bool const HasBPX =  std::find(FormatVector.begin(), FormatVector.end(), "BPX") != FormatVector.end();
  bool const HasBPY =  std::find(FormatVector.begin(), FormatVector.end(), "BPY") != FormatVector.end();
  bool const HasBPZ =  std::find(FormatVector.begin(), FormatVector.end(), "BPZ") != FormatVector.end();

  bool const HasPosition  = (HasX || HasY || HasZ);
  bool const HasBeta      = (HasBX || HasBY || HasBZ);
  bool const HasBetaPrime = (HasBPX || HasBPY || HasBPZ);


  // If you do not have a time of position you are indeed done for.
  if (!(HasT && HasPosition)) {
    throw std::length_error("trajectory must contain time and position in 1d at a minimum");
  }

  // A place for every column (although some are not used)
  double Input;

  std::string IgnoredRead;

  // Out variables which will be put in trajectory
  double T;
  TVector3D X(0, 0, 0);
  TVector3D B(0, 0, 0);
  TVector3D BP(0, 0, 0);

  // Loop over all lines in file
  for (std::string line; std::getline(f, line); ) {

    X.SetXYZ(0, 0, 0);
    B.SetXYZ(0, 0, 0);
    BP.SetXYZ(0, 0, 0);

    // Parse one line
    std::istringstream sline(line);
    for (size_t i = 0; i != NColumns; ++i) {
      if (FormatVector[i] == "*") {
        sline >> IgnoredRead;
      } else {
        sline >> Input;

        // If there is a fail the format specifier does not match the input
        if (sline.fail()) {
          throw std::length_error("When parsing the input trajectory file the format specified does not match the data in length");
        }

        if (FormatVector[i] == "T") {
          T = Input;
        } else if (FormatVector[i] == "X") {
          X.SetX(Input);
        } else if (FormatVector[i] == "Y") {
          X.SetY(Input);
        } else if (FormatVector[i] == "Z") {
          X.SetZ(Input);
        } else if (FormatVector[i] == "BX") {
          B.SetX(Input);
        } else if (FormatVector[i] == "BY") {
          B.SetY(Input);
        } else if (FormatVector[i] == "BZ") {
          B.SetZ(Input);
        } else if (FormatVector[i] == "BPX") {
          BP.SetX(Input);
        } else if (FormatVector[i] == "BPY") {
          BP.SetY(Input);
        } else if (FormatVector[i] == "BPZ") {
          BP.SetZ(Input);
        }
      }
    }

    if (B.Mag() >= 1.0) {
      throw std::out_of_range("Magnitude of beta >= 1 in input file.  Sorry, this is beyond einstein.  Line:\n" + line);
    }
    // We have a line, a point, add it!!
    this->AddPoint(X, B, BP, T);
  }

  // Close file
  f.close();

  // Trajectory data is all read.  If no beta, deduce from X, T
  if (!HasBeta) {
    this->ConstructBetaAtPoints();
  }

  // If no BP deduce from B, T
  if (!HasBetaPrime) {
    this->ConstructAoverCAtPoints();
  }

  return;
}




void TParticleTrajectoryPoints::ReadFromFileBinary (std::string const& FileName, std::string const& FormatIn)
{
  // Read a text file of the trajectory

  // Open the input file for reading
  std::ifstream f(FileName, std::ios::binary);
  if (!f.is_open()) {
    throw std::ifstream::failure("cannot open input file: " + FileName);
  }

  // For reading format from file
  char* FormatRead = 0x0;

  if (FormatIn == "") {
    // Read format length as first word
    int FormatReadLength = 0;
    f.read((char*) &FormatReadLength, sizeof(int));

    // Check format length
    if (FormatReadLength < 1) {
      throw std::ifstream::failure("Trying to read binary file format from file header failed.  Incorrect format.");
    }

    // Read input format
    FormatRead = new char[FormatReadLength+1];

    FormatRead[FormatReadLength] = '\0';
    f.read((char*) FormatRead, FormatReadLength * sizeof(char));
  }

  // Transform input format
  std::string Format = FormatRead == 0x0 ? FormatIn : FormatRead;
  std::transform(Format.begin(), Format.end(), Format.begin(), ::toupper);

  // If we created it let's delete it
  if (FormatRead != 0x0) {
    delete FormatRead;
  }

  // Ordered vector of inputs
  std::istringstream FormatStream(Format);
  std::vector<std::string> FormatVector;
  for (std::string ss; FormatStream >> ss; ) {
    FormatVector.push_back(ss.c_str());
  }

  if (FormatVector.size() < 1) {
    throw std::length_error("Format must contain at least one element");
  }


  // Number of columns
  size_t const NColumns = FormatVector.size();


  // Check which components exist
  bool const HasT   =  std::find(FormatVector.begin(), FormatVector.end(), "T") != FormatVector.end();
  bool const HasX   =  std::find(FormatVector.begin(), FormatVector.end(), "X") != FormatVector.end();
  bool const HasY   =  std::find(FormatVector.begin(), FormatVector.end(), "Y") != FormatVector.end();
  bool const HasZ   =  std::find(FormatVector.begin(), FormatVector.end(), "Z") != FormatVector.end();
  bool const HasBX  =  std::find(FormatVector.begin(), FormatVector.end(), "BX") != FormatVector.end();
  bool const HasBY  =  std::find(FormatVector.begin(), FormatVector.end(), "BY") != FormatVector.end();
  bool const HasBZ  =  std::find(FormatVector.begin(), FormatVector.end(), "BZ") != FormatVector.end();
  bool const HasBPX =  std::find(FormatVector.begin(), FormatVector.end(), "BPX") != FormatVector.end();
  bool const HasBPY =  std::find(FormatVector.begin(), FormatVector.end(), "BPY") != FormatVector.end();
  bool const HasBPZ =  std::find(FormatVector.begin(), FormatVector.end(), "BPZ") != FormatVector.end();

  bool const HasPosition  = (HasX || HasY || HasZ);
  bool const HasBeta      = (HasBX || HasBY || HasBZ);
  bool const HasBetaPrime = (HasBPX || HasBPY || HasBPZ);


  // If you do not have a time of position you are indeed done for.
  if (!(HasT && HasPosition)) {
    throw std::length_error("trajectory must contain time and position in 1d at a minimum");
  }

  // For reading the binary input
  double Input;

  // Out variables which will be put in trajectory
  double T;
  TVector3D X(0, 0, 0);
  TVector3D B(0, 0, 0);
  TVector3D BP(0, 0, 0);


  size_t ipoint = 0;
  while (!f.eof()) {

    // Setting this point to zero
    X.SetXYZ(0, 0, 0);
    B.SetXYZ(0, 0, 0);
    BP.SetXYZ(0, 0, 0);

    for (size_t i = 0; i != FormatVector.size(); ++i) {
      f.read( (char*)    &Input, sizeof(double));

      if (FormatVector[i] == "T") {
        T = Input;
      } else if (FormatVector[i] == "X") {
        X.SetX(Input);
      } else if (FormatVector[i] == "Y") {
        X.SetY(Input);
      } else if (FormatVector[i] == "Z") {
        X.SetZ(Input);
      } else if (FormatVector[i] == "BX") {
        B.SetX(Input);
      } else if (FormatVector[i] == "BY") {
        B.SetY(Input);
      } else if (FormatVector[i] == "BZ") {
        B.SetZ(Input);
      } else if (FormatVector[i] == "BPX") {
        BP.SetX(Input);
      } else if (FormatVector[i] == "BPY") {
        BP.SetY(Input);
      } else if (FormatVector[i] == "BPZ") {
        BP.SetZ(Input);
      }

    }


    // If end of file, stop, otherwise add to trajectory
    if (f.eof()) {
      break;
    } else {
      if (B.Mag() >= 1.0) {
        throw std::out_of_range("Magnitude of beta >= 1 in input file.  Sorry, this is beyond einstein.  ipoint: " + std::to_string(ipoint));
      }
      // We have a line, a point, add it!!
      this->AddPoint(X, B, BP, T);
    }

    ++ipoint;
  }

  // Trajectory data is all read.  If no beta, deduce from X, T
  if (!HasBeta) {
    this->ConstructBetaAtPoints();
  }

  // If no BP deduce from B, T
  if (!HasBetaPrime) {
    this->ConstructAoverCAtPoints();
  }

  // Close file nicely
  f.close();

  return;
}




void TParticleTrajectoryPoints::ConstructBetaAtPoints ()
{
  // This functino will construct Beta from the positions and time.
  // It should only be used when Beta information is missing

  std::vector<double> Times;
  std::vector<TVector3D> Pos;
  for (size_t i = 0; i != this->GetNPoints(); ++i) {
    Times.push_back(this->GetT(i));
    Pos.push_back(this->GetX(i));
  }
  TOMATH::TSpline1D3<TVector3D> S(Times, Pos);
  for (size_t i = 0; i != this->GetNPoints(); ++i) {
    this->SetB(i, S.GetDerivative(i) / TOSCARSSR::C());
  }

  return;
}




void TParticleTrajectoryPoints::ConstructAoverCAtPoints ()
{
  // This function will construct AoverC from the Beta and time.
  // It should only be used when AoverC information is missing

  std::vector<double> Times;
  std::vector<TVector3D> Betas;
  for (size_t i = 0; i != this->GetNPoints(); ++i) {
    Times.push_back(this->GetT(i));
    Betas.push_back(this->GetB(i));
  }
  TOMATH::TSpline1D3<TVector3D> S(Times, Betas);
  for (size_t i = 0; i != this->GetNPoints(); ++i) {
    this->SetAoverC(i, S.GetDerivative(i));
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

