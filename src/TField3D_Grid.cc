////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Tue Sep 20 07:55:30 EDT 2016
//
////////////////////////////////////////////////////////////////////

#include "TField3D_Grid.h"

#include <fstream>
#include <sstream>
#include <algorithm>
#include <array>

TField3D_Grid::TField3D_Grid ()
{
  fRotated.SetXYZ(0, 0, 0);
  fTranslation.SetXYZ(0, 0, 0);
}




TField3D_Grid::TField3D_Grid (std::string const& InFileName, std::string const& FileFormat, TVector3D const& Rotations, TVector3D const& Translation, std::vector<double> const& Scaling, char const CommentChar)
{
  // I will accept lower-case
  std::string format = FileFormat;
  std::transform(format.begin(), format.end(), format.begin(), ::toupper);

  // Which file format are you looking at?
  if (format == "OSCARS") {
    this->ReadFile(InFileName, Rotations, Translation, Scaling);
  } else if (std::string(format.begin(), format.begin() + 8) == "OSCARS1D") {
    this->ReadFile_OSCARS1D(InFileName, FileFormat, Rotations, Translation, Scaling, CommentChar);
  } else if (format == "SPECTRA") {
    this->ReadFile_SPECTRA(InFileName, Rotations, Translation, CommentChar);
  } else if (format == "SRW") {
    this->ReadFile_SRW(InFileName, Rotations, Translation, CommentChar);
  } else {
    std::cerr << "TField3D_Grid::TField3D_Grid format error format: " << FileFormat << std::endl;
    throw std::invalid_argument("incorrect format given");
  }
}




TField3D_Grid::TField3D_Grid (std::vector<std::pair<double, std::string> > Mapping, std::string const& FileFormat, double const Parameter, TVector3D const& Rotations, TVector3D const& Translation, std::vector<double> const& Scaling, char const CommentChar)
{
  // This one is for interpolated fields from a mapping vector and parameter value.
  // It is meant for interpolating between different undulator gaps, but it is generalized
  // to interpolate any fields

  // I will accept lower-case
  std::string format = FileFormat;
  std::transform(format.begin(), format.end(), format.begin(), ::toupper);

  // Which file format are you looking at?
  if (format == "OSCARS") {
    this->InterpolateFromFiles(Mapping, Parameter, Rotations, Translation, Scaling);
  } else if (std::string(format.begin(), format.begin() + 8) == "OSCARS1D") {
    //this->ReadFile_OSCARS1D(InFileName, FileFormat, Rotations, Translation, Scaling, CommentChar);
  } else if (format == "SPECTRA") {
    //this->ReadFile_SPECTRA(InFileName, Rotations, Translation, CommentChar);
  } else if (format == "SRW") {
    //this->ReadFile_SRW(InFileName, Rotations, Translation, CommentChar);
  } else {
    std::cerr << "TField3D_Grid::TField3D_Grid format error format: " << FileFormat << std::endl;
    throw std::invalid_argument("incorrect format given");
  }
}




TField3D_Grid::~TField3D_Grid ()
{
  // Destruction is my goal
}




double TField3D_Grid::GetFx (double const X, double const Y, double const Z) const
{
  return this->GetF(TVector3D(X, Y, Z)).GetX();
}




double TField3D_Grid::GetFy (double const X, double const Y, double const Z) const
{
  return this->GetF(TVector3D(X, Y, Z)).GetY();
}




double TField3D_Grid::GetFz (double const X, double const Y, double const Z) const
{
  return this->GetF(TVector3D(X, Y, Z)).GetZ();
}





TVector3D TField3D_Grid::GetF (double const X, double const Y, double const Z) const
{
  return this->GetF(TVector3D(X, Y, Z));
}



size_t TField3D_Grid::GetIndex (size_t const ix, size_t const iy, size_t const iz) const
{
  return ix * fNY * fNZ + iy * fNZ + iz;
}



TVector3D TField3D_Grid::GetF (TVector3D const& XIN) const
{
  // Get the field at a point in space.  Must rotate point into coordinate system, then translate it.

  // Rotate and Translate
  TVector3D X(XIN);
  X.RotateSelfXYZ(fRotated);
  X -= fTranslation;

  // If outside the range, return a zero
  if (fNX > 1 && (X.GetX() <= fXStart || X.GetX() >= fXStop)) {
    return TVector3D(0, 0, 0);
  }
  if (fNY > 1 && (X.GetY() <= fYStart || X.GetY() >= fYStop)) {
    return TVector3D(0, 0, 0);
  }
  if (fNZ > 1 && (X.GetZ() <= fZStart || X.GetZ() >= fZStop)) {
    return TVector3D(0, 0, 0);
  }
      
  // Get index in each dimension relative to start and stop
  size_t const nx = fNX > 1 ? (X.GetX() - fXStart) / fXStep      : 0;
  double const dx = fNX > 1 ? (X.GetX() - fXStart) - nx * fXStep : 0;

  size_t const ny = fNY > 1 ? (X.GetY() - fYStart) / fYStep      : 0;
  double const dy = fNY > 1 ? (X.GetY() - fYStart) - ny * fYStep : 0;

  size_t const nz = fNZ > 1 ? (X.GetZ() - fZStart) / fZStep      : 0;
  double const dz = fNZ > 1 ? (X.GetZ() - fZStart) - nz * fZStep : 0;


  switch (fDIMX) {
    case kDIMX_XYZ:
      {
        // First move in X to find the 4 points of square
        size_t const i000 = GetIndex(nx + 0, ny + 0, nz + 0);
        size_t const i100 = GetIndex(nx + 1, ny + 0, nz + 0);
        TVector3D const v00 = fData[i000] + dx * (fData[i100] - fData[i000]) / fXStep;

        size_t const i010 = GetIndex(nx + 0, ny + 1, nz + 0);
        size_t const i110 = GetIndex(nx + 1, ny + 1, nz + 0);
        TVector3D const v10 = fData[i010] + dx * (fData[i110] - fData[i010]) / fXStep;

        size_t const i001 = GetIndex(nx + 0, ny + 0, nz + 1);
        size_t const i101 = GetIndex(nx + 1, ny + 0, nz + 1);
        TVector3D const v01 = fData[i001] + dx * (fData[i101] - fData[i001]) / fXStep;

        size_t const i011 = GetIndex(nx + 0, ny + 1, nz + 1);
        size_t const i111 = GetIndex(nx + 1, ny + 1, nz + 1);
        TVector3D const v11 = fData[i011] + dx * (fData[i111] - fData[i011]) / fXStep;

        // Step in Y to find 2 points
        TVector3D const v0 = v00 + dy * (v10 - v00) / fYStep;
        TVector3D const v1 = v01 + dy * (v11 - v01) / fYStep;

        // Step in Z to find point
        return v0 + dz * (v1 - v0) / fZStep;
      }
      break;
    case kDIMX_X:
      {
        size_t const i0 = nx + 0;
        size_t const i1 = nx + 1;
        return fData[i0] + dx * (fData[i1] - fData[i0]) / fXStep;
      }
    case kDIMX_Y:
      {
        size_t const i0 = ny + 0;
        size_t const i1 = ny + 1;
        return fData[i0] + dy * (fData[i1] - fData[i0]) / fYStep;
      }
    case kDIMX_Z:
      {
        size_t const i0 = nz + 0;
        size_t const i1 = nz + 1;
        return fData[i0] + dz * (fData[i1] - fData[i0]) / fZStep;
      }
    case kDIMX_XY:
      {
        size_t const i00 = GetIndex(nx + 0, ny + 0, 0);
        size_t const i10 = GetIndex(nx + 1, ny + 0, 0);
        TVector3D const v0 = fData[i00] + dx * (fData[i10] - fData[i00]) / fXStep;

        size_t const i01 = GetIndex(nx + 0, ny + 1, 0);
        size_t const i11 = GetIndex(nx + 1, ny + 1, 0);
        TVector3D const v1 = fData[i01] + dx * (fData[i11] - fData[i01]) / fXStep;

        return v0 + dy * (v1 - v0) / fYStep;
      }
    case kDIMX_XZ:
      {
        size_t const i00 = GetIndex(nx + 0, 0, nz + 0);
        size_t const i10 = GetIndex(nx + 1, 0, nz + 0);
        TVector3D const v0 = fData[i00] + dx * (fData[i10] - fData[i00]) / fXStep;

        size_t const i01 = GetIndex(nx + 0, 0, nz + 1);
        size_t const i11 = GetIndex(nx + 1, 0, nz + 1);
        TVector3D const v1 = fData[i01] + dx * (fData[i11] - fData[i01]) / fXStep;

        return v0 + dz * (v1 - v0) / fZStep;
      }
    case kDIMX_YZ:
      {
        size_t const i00 = GetIndex(0, ny + 0, nz + 0);
        size_t const i10 = GetIndex(0, ny + 1, nz + 0);
        TVector3D const v0 = fData[i00] + dy * (fData[i10] - fData[i00]) / fYStep;

        size_t const i01 = GetIndex(0, ny + 0, nz + 1);
        size_t const i11 = GetIndex(0, ny + 1, nz + 1);
        TVector3D const v1 = fData[i01] + dy * (fData[i11] - fData[i01]) / fYStep;

        return v0 + dz * (v1 - v0) / fZStep;
      }
      break;
    default:
      throw std::out_of_range("unknown dimension");
  }

  throw std::out_of_range("unknown dimension");
}










double TField3D_Grid::GetHeaderValue (std::string const& L) const
{
  // Get the value after comment character

  // Easy reading
  std::istringstream S;
  S.str(L);

  // For the two values of interest
  double Value;

  S >> Value;

  // Check read state of input
  if (S.bad()) {
    std::cerr << "ERROR: S is bad" << std::endl;
    throw std::ifstream::failure("cannot read header value");
  }

  return Value;
}







double TField3D_Grid::GetHeaderValueSRW (std::string const& L, const char CommentChar) const
{
  // Get the value after comment character

  // Easy reading
  std::istringstream S;
  S.str(L);

  // For the two values of interest
  char CC;
  double Value;

  S.get(CC);

  if (CC != CommentChar) {
    std::cerr << "ERROR: bad format in header" << std::endl;
    throw std::ifstream::failure("something is wrong with the comment character, it was not seen");
  }

  S >> Value;

  // Check read state of input
  if (S.bad()) {
    std::cerr << "ERROR: S is bad" << std::endl;
    throw std::ifstream::failure("cannot read header value SRW format");
  }

  return Value;
}







void TField3D_Grid::ReadFile (std::string const& InFileName, TVector3D const& Rotations, TVector3D const& Translation, std::vector<double> const& Scaling, char const CommentChar)
{
  // Read file with the best format in the entire world, OSCARSv1.0

  // Open the input file and throw exception if not open
  std::ifstream fi(InFileName);
  if (!fi) {
    std::cerr << "ERROR: cannot open file" << std::endl;
    throw std::ifstream::failure("cannot open file for reading");
  }

  // For reading lines of file
  std::istringstream S;
  std::string L;

  // Initial line for comment
  std::getline(fi, L);


  // Initial X
  std::getline(fi, L);
  double const XStartIN = GetHeaderValue(L);

  // Step X
  std::getline(fi, L);
  double const XStepIN = GetHeaderValue(L);

  // Number of points X
  std::getline(fi, L);
  int const NX = (int) GetHeaderValue(L);


  // Initial Y
  std::getline(fi, L);
  double const YStartIN = GetHeaderValue(L);

  // Step Y
  std::getline(fi, L);
  double const YStepIN = GetHeaderValue(L);

  // Number of points Y
  std::getline(fi, L);
  int const NY = (int) GetHeaderValue(L);


  // Initial Z
  std::getline(fi, L);
  double const ZStartIN = GetHeaderValue(L);

  // Step Z
  std::getline(fi, L);
  double const ZStepIN = GetHeaderValue(L);

  // Number of points Z
  std::getline(fi, L);
  int const NZ = (int) GetHeaderValue(L);

  // If we're doing any scaling, scale spatial dimensions and fields.  Start with stepsize change
  double const XStep = Scaling.size() > 0 ? XStepIN * Scaling[0] : XStepIN;
  double const YStep = Scaling.size() > 1 ? YStepIN * Scaling[1] : YStepIN;
  double const ZStep = Scaling.size() > 2 ? ZStepIN * Scaling[2] : ZStepIN;

  // Get field scaling if it exists
  double const FxScaling = Scaling.size() > 3 ? Scaling[3] : 1;
  double const FyScaling = Scaling.size() > 4 ? Scaling[4] : 1;
  double const FzScaling = Scaling.size() > 5 ? Scaling[5] : 1;

  // Calculate new start point
  double const MiddleX = XStartIN + XStepIN * (NX - 1) / 2.;
  double const XStart  = MiddleX - XStep * (NX - 1) / 2.;
  double const MiddleY = YStartIN + YStepIN * (NY - 1) / 2.;
  double const YStart  = MiddleY - YStep * (NY -1) / 2.;
  double const MiddleZ = ZStartIN + ZStepIN * (NZ - 1) / 2.;
  double const ZStart  = MiddleZ - ZStep * (NZ - 1) / 2.;


  // Check Number of points is > 0 for all
  if (NX < 1 || NY < 1 || NY < 1) {
    std::cerr << "ERROR: invalid npoints" << std::endl;
    throw std::out_of_range("invalid number of points in at least one dimension");
  }

  // Save position data to object variables
  fNX = NX;
  fNY = NY;
  fNZ = NZ;
  fXStart = XStart;
  fYStart = YStart;
  fZStart = ZStart;
  fXStep  = XStep;
  fYStep  = YStep;
  fZStep  = ZStep;
  fXStop  = fXStart + (fNX - 1) * fXStep;
  fYStop  = fYStart + (fNY - 1) * fYStep;
  fZStop  = fZStart + (fNZ - 1) * fZStep;

  fHasX = NX > 1 ? true : false;
  fHasY = NY > 1 ? true : false;
  fHasZ = NZ > 1 ? true : false;

  if (fHasX && fHasY && fHasZ) {
    fDIMX = kDIMX_XYZ;
  } else if (fHasX && fHasY) {
    fDIMX = kDIMX_XY;
  } else if (fHasX && fHasZ) {
    fDIMX = kDIMX_XZ;
  } else if (fHasY && fHasZ) {
    fDIMX = kDIMX_YZ;
  } else if (fHasX) {
    fDIMX = kDIMX_X;
  } else if (fHasY) {
    fDIMX = kDIMX_Y;
  } else if (fHasZ) {
    fDIMX = kDIMX_Z;
  } else {
    std::cerr << "ERROR: error in file header format" << std::endl;
    throw std::out_of_range("invalid dimensions");
  }

  fXDIM = 0;
  if (fHasX) {
    ++fXDIM;
  }
  if (fHasY) {
    ++fXDIM;
  }
  if (fHasZ) {
    ++fXDIM;
  }

  // Reserve correct number of points in vector (slightly faster)
  fData.reserve(fNX * fNY * fNZ);

  // Temp variables for field
  double fx;
  double fy;
  double fz;


  // Loop over all points
  for (int ix = 0; ix != NX; ++ix) {
    for (int iy = 0; iy != NY; ++iy) {
      for (int iz = 0; iz != NZ; ++iz) {

        // Grab a line from input file
        std::getline(fi, L);

        // Check we did not hit an EOF
        if (fi.eof()) {
          std::cerr << "ERROR: bad input file" << std::endl;
          throw std::ifstream::failure("error reading file.  Check format");
        }

        // Read data
        S.clear();
        S.str(L);
        S >> fx >> fy >> fz;

        // Scale values?
        if (FxScaling != 1) {
          fx *= FxScaling;
        }
        if (FyScaling != 1) {
          fy *= FyScaling;
        }
        if (FzScaling != 1) {
          fz *= FzScaling;
        }

        // Check the stream did not hit an EOF
        if (S.fail()) {
          std::cerr << "ERRROR: input stream bad" << std::endl;
          throw std::ifstream::failure("error reading file.  Check format");
        }

        // Push data to storage
        TVector3D F(fx, fy, fz);
        F.RotateSelfXYZ(Rotations);
        fData.push_back(F);
      }
    }
  }

  // Close file
  fi.close();

  // Store Rotations and Translation
  fRotated = Rotations;
  fTranslation = Translation;

  return;
}























void TField3D_Grid::ReadFile_OSCARS1D (std::string const& InFileName, std::string const& InFormat, TVector3D const& Rotations, TVector3D const& Translation, std::vector<double> const& Scaling, char const CommentChar)
{
  // Read file with OSCARS1D format

  // Open the input file and throw exception if not open
  std::ifstream fi(InFileName);
  if (!fi) {
    std::cerr << "ERROR: cannot open file" << std::endl;
    throw std::ifstream::failure("cannot open file for reading");
  }

  // And this is for which order they come in
  std::vector<int> Order(4, -1);

  // Make it a stream and set it to the format string minus the OSCARS1D
  std::string const FormatString(InFormat.begin() + 8, InFormat.end());
  std::istringstream s;
  s.str(FormatString);

  // String for identifier
  std::string c;

  // Counts
  int index = 0;
  int XDIM = 0;
  int FDIM = 0;

  // Which fields we have
  bool HasFx = false;
  bool HasFy = false;
  bool HasFz = false;

  // Which axis is used
  std::string Axis = "";

  // Look at format string
  while (s >> c) {

    if (index > 3) {
      std::cerr << "ERROR: spatial or B-field dimensions are too large(index>3)" << std::endl;
      throw std::out_of_range("spatial or B-field dimensions are too large(index>3)");
    }

    // Check if it is XYZBxByBz and in which order
    if (c == "X" || c == "Y" || c == "Z") {
      Axis = c;
      ++XDIM;
      Order[index] = 0;
      ++index;
    } else if (c == "Bx" || c == "Ex" || c == "Fx") {
      HasFx = true;
      ++FDIM;
      Order[index] = 1;
      ++index;
    } else if (c == "By" || c == "Ey" || c == "Fy") {
      HasFy = true;
      ++FDIM;
      Order[index] = 2;
      ++index;
    } else if (c == "Bz" || c == "Ez" || c == "Fz") {
      HasFz = true;
      ++FDIM;
      Order[index] = 3;
      ++index;
    } else {
      std::cerr << "ERROR: Incorrect format" << std::endl;
      throw std::invalid_argument("only excepts X Y Z Bx By Bz");
    }
  }

  // How many cols to ouout
  int InputCount = XDIM + FDIM;

  // At the moment only support 1D irregular grid
  if (XDIM != 1) {
    std::cerr << "ERROR: spatial or B-field dimensions are too large(>3)" << std::endl;
    throw std::out_of_range("spatial or B-field dimensions are too large");
  }


  // For reading lines of file
  std::istringstream S;
  std::string L;

  // Initial line for comment
  std::getline(fi, L);

  // Vector for the data inputs
  std::vector<std::array<double, 4> > InputData;

  // Loop over all lines in file
  for ( ; std::getline(fi, L); ) {

    // Look for a blank line or comment line and skip if found.  You should never use tab btw.
    size_t FirstChar = L.find_first_not_of(" \t");
    if (FirstChar == std::string::npos || L[FirstChar] == CommentChar) {
      continue;
    }

    // Clear and set the istringstream for this line
    S.clear();
    S.str(L);

    // Vector for this line of input
    std::vector<double> Value(4);

    // Read this line of data
    for (int i = 0; i < InputCount; ++i) {
      S >> Value[Order[i]];
      if (S.fail()) {
        throw std::length_error("something is incorrect with data format or iformat string");
      }
    }

    // Scale the input as requested
    for (size_t iscale = 0; iscale != Scaling.size() && iscale < 4; ++iscale) {
      Value[Order[iscale]] *= Scaling[iscale];
    }

    // Add to data array.
    std::array<double, 4> a = { {Value[0], Value[1], Value[2], Value[3]} };
    InputData.push_back(a);
  }

  // Close file
  fi.close();

  // Sort the field
  std::sort(InputData.begin(), InputData.end(), this->CompareField1D);


  // Now we need to make a regular 3D grid...

  // Check to see there are at least 2 data points
  if (InputData.size() < 2) {
    std::cerr << "ERROR: not enough data points" << std::endl;
    throw std::length_error("not enough data points");
  }

  // Grab the first and last Z
  double const First = InputData[0][0];
  double const Last  = InputData[InputData.size() - 1][0];


  // UPDATE: VAR
  int const NPointsPerMeter = 10000;

  // Get the number of points and the step size
  size_t const NPoints  = (Last - First) * NPointsPerMeter;
  double const StepSize = (Last - First) / (double) (NPoints - 1);

  // Set all to default values
  fNX = 1;
  fNY = 1;
  fNZ = 1;
  fXStart = 0;
  fYStart = 0;
  fZStart = 0;
  fXStep = 0;
  fYStep = 0;
  fZStep = 0;
  fXStop = 0;
  fYStop = 0;
  fZStop = 0;

  // Set correct one for limits based on which axis was given
  if (Axis == "X") {
    fNX = NPoints;
    fXStart = First;
    fXStop = Last;
    fXStep = StepSize;
  } else if  (Axis == "Y") {
    fNY = NPoints;
    fYStart = First;
    fYStop = Last;
    fYStep = StepSize;
  } else if  (Axis == "Z") {
    fNZ = NPoints;
    fZStart = First;
    fZStop = Last;
    fZStep = StepSize;
  }


  // Get dimensions correct
  fHasX = fNX > 1 ? true : false;
  fHasY = fNY > 1 ? true : false;
  fHasZ = fNZ > 1 ? true : false;

  if (fHasX && fHasY && fHasZ) {
    fDIMX = kDIMX_XYZ;
  } else if (fHasX && fHasY) {
    fDIMX = kDIMX_XY;
  } else if (fHasX && fHasZ) {
    fDIMX = kDIMX_XZ;
  } else if (fHasY && fHasZ) {
    fDIMX = kDIMX_YZ;
  } else if (fHasX) {
    fDIMX = kDIMX_X;
  } else if (fHasY) {
    fDIMX = kDIMX_Y;
  } else if (fHasZ) {
    fDIMX = kDIMX_Z;
  } else {
    std::cerr << "ERROR: error in file header format" << std::endl;
    throw std::out_of_range("invalid dimensions");
  }

  fXDIM = 0;
  if (fHasX) {
    ++fXDIM;
  }
  if (fHasY) {
    ++fXDIM;
  }
  if (fHasZ) {
    ++fXDIM;
  }



  // Clear the internal data and reserve the number of points
  fData.clear();
  fData.reserve(NPoints);


  // Variables to hold the slope between two real points and new By (linear interpolated)
  double SlopeX;
  double SlopeY;
  double SlopeZ;
  double NewBx;
  double NewBy;
  double NewBz;

  // I only initialize to avoid compile-time warnings.  These are for finding the
  // adj bins for linear interpolation
  size_t MinBin    = 0;
  size_t AfterBin  = 0;
  size_t BeforeBin = 0;

  // Variable to hold new calculated Z position
  double ThisZ;

  // For each desired point find the bin before and after the desired Z position
  for (size_t i = 0; i != NPoints; ++i) {
    ThisZ = i * StepSize + First;
    for (size_t j = MinBin + 1; j != InputData.size(); ++j) {
      if (InputData[j][0] > ThisZ) {
        AfterBin  = j;
        BeforeBin = j - 1;
        MinBin = j - 1;
        break;
      }
    }


    // Calculate the By at desired position based on linear interpolation
    SlopeX = !HasFx ? 0 : (InputData[AfterBin][1] - InputData[BeforeBin][1]) /  (InputData[AfterBin][0] - InputData[BeforeBin][0]);
    SlopeY = !HasFy ? 0 : (InputData[AfterBin][2] - InputData[BeforeBin][2]) /  (InputData[AfterBin][0] - InputData[BeforeBin][0]);
    SlopeZ = !HasFz ? 0 : (InputData[AfterBin][3] - InputData[BeforeBin][3]) /  (InputData[AfterBin][0] - InputData[BeforeBin][0]);
    NewBx = !HasFx ? 0 : InputData[BeforeBin][1] + (ThisZ - InputData[BeforeBin][0]) * SlopeX;
    NewBy = !HasFy ? 0 : InputData[BeforeBin][2] + (ThisZ - InputData[BeforeBin][0]) * SlopeY;
    NewBz = !HasFz ? 0 : InputData[BeforeBin][3] + (ThisZ - InputData[BeforeBin][0]) * SlopeZ;


    // Append the new BxByBz to the output vector
    TVector3D F(NewBx, NewBy, NewBz);
    F.RotateSelfXYZ(Rotations);
    fData.push_back(F);
  }

  // Clear array data
  InputData.clear();

  // Store Rotations and Translation
  fRotated = Rotations;
  fTranslation = Translation;

  return;
}























void TField3D_Grid::ReadFile_SRW (std::string const& InFileName, TVector3D const& Rotations, TVector3D const& Translation, char const CommentChar)
{
  // Read file with SRW field input format

  // Open the input file and throw exception if not open
  std::ifstream fi(InFileName);
  if (!fi) {
    std::cerr << "ERROR: cannot open file" << std::endl;
    throw std::ifstream::failure("cannot open file for reading SRW format");
  }

  // For reading lines of file
  std::istringstream S;
  std::string L;

  // Initial line for comment
  std::getline(fi, L);


  // Initial X
  std::getline(fi, L);
  double const XStart = GetHeaderValueSRW(L, CommentChar);

  // Step X
  std::getline(fi, L);
  double const XStep = GetHeaderValueSRW(L, CommentChar);

  // Number of points X
  std::getline(fi, L);
  int const NX = (int) GetHeaderValueSRW(L, CommentChar);


  // Initial Y
  std::getline(fi, L);
  double const YStart = GetHeaderValueSRW(L, CommentChar);

  // Step Y
  std::getline(fi, L);
  double const YStep = GetHeaderValueSRW(L, CommentChar);

  // Number of points Y
  std::getline(fi, L);
  int const NY = (int) GetHeaderValueSRW(L, CommentChar);


  // Initial Z
  std::getline(fi, L);
  double const ZStart = GetHeaderValueSRW(L, CommentChar);

  // Step Z
  std::getline(fi, L);
  double const ZStep = GetHeaderValueSRW(L, CommentChar);

  // Number of points Z
  std::getline(fi, L);
  int const NZ = (int) GetHeaderValueSRW(L, CommentChar);


  // Check Number of points is > 0 for all
  if (NX < 1 || NY < 1 || NY < 1) {
    std::cerr << "ERROR: invalid npoints" << std::endl;
    throw std::out_of_range("invalid dimensions");
  }

  // Save position data to object variables
  fNX = NX;
  fNY = NY;
  fNZ = NZ;
  fXStart = XStart;
  fYStart = YStart;
  fZStart = ZStart;
  fXStep  = XStep;
  fYStep  = YStep;
  fZStep  = ZStep;
  fXStop  = fXStart + (fNX - 1) * fXStep;
  fYStop  = fYStart + (fNY - 1) * fYStep;
  fZStop  = fZStart + (fNZ - 1) * fZStep;

  fHasX = NX > 1 ? true : false;
  fHasY = NY > 1 ? true : false;
  fHasZ = NZ > 1 ? true : false;

  if (fHasX && fHasY && fHasZ) {
    fDIMX = kDIMX_XYZ;
  } else if (fHasX && fHasY) {
    fDIMX = kDIMX_XY;
  } else if (fHasX && fHasZ) {
    fDIMX = kDIMX_XZ;
  } else if (fHasY && fHasZ) {
    fDIMX = kDIMX_YZ;
  } else if (fHasX) {
    fDIMX = kDIMX_X;
  } else if (fHasY) {
    fDIMX = kDIMX_Y;
  } else if (fHasZ) {
    fDIMX = kDIMX_Z;
  } else {
    std::cerr << "ERROR: error in file header format" << std::endl;
    throw std::out_of_range("invalid dimensions");
  }

  fXDIM = 0;
  if (fHasX) {
    ++fXDIM;
  }
  if (fHasY) {
    ++fXDIM;
  }
  if (fHasZ) {
    ++fXDIM;
  }

  // Resize the vector because we will need to insert values non-sequentially
  // ie can't use push_back
  fData.resize(fNX * fNY * fNZ);

  // Temp variables for field
  double fx;
  double fy;
  double fz;


  // Loop over all points
  for (int iz = 0; iz != NZ; ++iz) {
    for (int iy = 0; iy != NY; ++iy) {
      for (int ix = 0; ix != NX; ++ix) {

        // Grab a line from input file
        std::getline(fi, L);

        // Read data
        S.clear();
        S.str(L);
        S >> fx >> fy >> fz;


        // Check the stream did not hit an EOF
        if (S.fail() || fi.fail()) {
          std::cerr << "ERRROR: input stream bad" << std::endl;
          throw std::ifstream::failure("input file stream failure");
        }

        // Push data to storage
        TVector3D F(fx, fy, fz);
        F.RotateSelfXYZ(Rotations);

        size_t const Index = this->GetIndex(ix, iy, iz);
        if (Index >= fData.size()) {
          throw;
        }
        fData[Index] = F;
      }
    }
  }

  // Close file
  fi.close();

  // Store Rotations and Translation
  fRotated = Rotations;
  fTranslation = Translation;

  return;
}








void TField3D_Grid::ReadFile_SPECTRA (std::string const& InFileName, TVector3D const& Rotations, TVector3D const& Translation, char const CommentChar)
{
  // Read file with SPECTRA field input format

  // Open the input file and throw exception if not open
  std::ifstream fi(InFileName);
  if (!fi) {
    std::cerr << "ERROR: cannot open file" << std::endl;
    throw std::ifstream::failure("cannot open file");
  }

  // For reading lines of file
  std::istringstream S;
  std::string L;

  // Initial line for comment
  std::getline(fi, L);

  // Now header information
  std::getline(fi, L);
  S.str(L);

  // Grab parameters and correct for [mm] -> [m] conversion.
  S >> fXStep >> fYStep >> fZStep >> fNX >> fNY >> fNZ;
  fXStep /= 1000.;
  fYStep /= 1000.;
  fZStep /= 1000.;

  if (S.bad()) {
    throw std::ifstream::failure("file stream failure");
  }


  // Check Number of points is > 0 for all
  if (fNX < 1 || fNY < 1 || fNY < 1) {
    std::cerr << "ERROR: invalid npoints" << std::endl;
    throw std::out_of_range("invalid number of points in at least one dimension");
  }

  // Save position data to object variables
  fXStart = -(fXStep * (fNX - 1)) / 2;
  fYStart = -(fYStep * (fNY - 1)) / 2;
  fZStart = -(fZStep * (fNZ - 1)) / 2;
  fXStop  = fXStart + (fNX - 1) * fXStep;
  fYStop  = fYStart + (fNY - 1) * fYStep;
  fZStop  = fZStart + (fNZ - 1) * fZStep;

  fHasX = fNX > 1 ? true : false;
  fHasY = fNY > 1 ? true : false;
  fHasZ = fNZ > 1 ? true : false;

  if (fHasX && fHasY && fHasZ) {
    fDIMX = kDIMX_XYZ;
  } else if (fHasX && fHasY) {
    fDIMX = kDIMX_XY;
  } else if (fHasX && fHasZ) {
    fDIMX = kDIMX_XZ;
  } else if (fHasY && fHasZ) {
    fDIMX = kDIMX_YZ;
  } else if (fHasX) {
    fDIMX = kDIMX_X;
  } else if (fHasY) {
    fDIMX = kDIMX_Y;
  } else if (fHasZ) {
    fDIMX = kDIMX_Z;
  } else {
    std::cerr << "ERROR: error in file header format" << std::endl;
    throw std::out_of_range("invalid dimensions");
  }

  fXDIM = 0;
  if (fHasX) {
    ++fXDIM;
  }
  if (fHasY) {
    ++fXDIM;
  }
  if (fHasZ) {
    ++fXDIM;
  }

  // Reserve correct number of points in vector (slightly faster)
  fData.reserve(fNX * fNY * fNZ);

  // Temp variables for field
  double fx;
  double fy;
  double fz;

  // Loop over all points
  for (int ix = 0; ix != fNX; ++ix) {
    for (int iy = 0; iy != fNY; ++iy) {
      for (int iz = 0; iz != fNZ; ++iz) {

        // Grab a line from input file
        std::getline(fi, L);

        // Check we did not hit an EOF
        if (fi.eof()) {
          std::cerr << "ERROR: bad input file" << std::endl;
          throw std::ifstream::failure("file stream failure");
        }

        // Read data
        S.clear();
        S.str("");
        S.str(L);
        S >> fx >> fy >> fz;

        // Check the stream did not hit an EOF
        if (S.fail()) {
          std::cerr << "ERRROR: input stream bad" << std::endl;
          throw std::ifstream::failure("file stream failure");
        }

        // Push data to storage
        TVector3D F(fx, fy, fz);
        F.RotateSelfXYZ(Rotations);
        fData.push_back(F);
      }
    }
  }

  // Close file
  fi.close();

  // Store Rotations and Translation
  fRotated = Rotations;
  fTranslation = Translation;

  return;
}
















void TField3D_Grid::InterpolateFromFiles (std::vector<std::pair<double, std::string> > const& Mapping, double const Parameter, TVector3D const& Rotations, TVector3D const& Translation, std::vector<double> const& Scaling, char const CommentChar)
{
  // Get interpolated field based on input files

  // First sort the input mapping vector
  std::vector<std::pair<double, std::string> > MyMapping = Mapping;
  std::sort(MyMapping.begin(), MyMapping.end(), this->CompareMappingElements);

  // Open all files and push files and parameters to vectors
  std::vector<std::ifstream*> InFiles;
  std::vector<double>  Parameters;
  for (std::vector<std::pair<double, std::string> >::iterator it = MyMapping.begin(); it != MyMapping.end(); ++it) {
    Parameters.push_back(it->first);
    InFiles.push_back(new std::ifstream(it->second));

    if (!InFiles.back()->is_open()) {
      std::cerr << "ERROR: cannot open file" << std::endl;
      //throw std::ifstream::failure("cannot open file for reading");
    }
  }

  // For reading lines of file
  std::istringstream S;
  std::string L;

  // Header values
  std::vector<double> HeaderValues;
  for (size_t ih = 0; ih != 10; ++ih) {
    std::getline(*(InFiles[0]), L);
    HeaderValues.push_back(GetHeaderValue(L));

    for (size_t i = 1; i < InFiles.size(); ++i) {
      std::getline(*(InFiles[i]), L);
      if (HeaderValues[ih] != GetHeaderValue(L)) {
        throw;
      }
    }
  }


  // Initial X
  double const XStartIN = HeaderValues[1];

  // Step X
  double const XStepIN = HeaderValues[2];

  // Number of points X
  int const NX = (int) HeaderValues[3];


  // Initial Y
  double const YStartIN = HeaderValues[4];

  // Step Y
  double const YStepIN = HeaderValues[5];

  // Number of points Y
  int const NY = (int) HeaderValues[6];


  // Initial Z
  double const ZStartIN = HeaderValues[7];

  // Step Z
  double const ZStepIN = HeaderValues[8];

  // Number of points Z
  int const NZ = (int) HeaderValues[9];



  // If we're doing any scaling, scale spatial dimensions and fields.  Start with stepsize change
  double const XStep = Scaling.size() > 0 ? XStepIN * Scaling[0] : XStepIN;
  double const YStep = Scaling.size() > 1 ? YStepIN * Scaling[1] : YStepIN;
  double const ZStep = Scaling.size() > 2 ? ZStepIN * Scaling[2] : ZStepIN;

  // Get field scaling if it exists
  double const FxScaling = Scaling.size() > 3 ? Scaling[3] : 1;
  double const FyScaling = Scaling.size() > 4 ? Scaling[4] : 1;
  double const FzScaling = Scaling.size() > 5 ? Scaling[5] : 1;

  // Calculate new start point
  double const MiddleX = XStartIN + XStepIN * (NX - 1) / 2.;
  double const XStart  = MiddleX - XStep * (NX - 1) / 2.;
  double const MiddleY = YStartIN + YStepIN * (NY - 1) / 2.;
  double const YStart  = MiddleY - YStep * (NY -1) / 2.;
  double const MiddleZ = ZStartIN + ZStepIN * (NZ - 1) / 2.;
  double const ZStart  = MiddleZ - ZStep * (NZ - 1) / 2.;


  // Check Number of points is > 0 for all
  if (NX < 1 || NY < 1 || NY < 1) {
    std::cerr << "ERROR: invalid npoints" << std::endl;
    throw std::out_of_range("invalid number of points in at least one dimension");
  }

  // Save position data to object variables
  fNX = NX;
  fNY = NY;
  fNZ = NZ;
  fXStart = XStart;
  fYStart = YStart;
  fZStart = ZStart;
  fXStep  = XStep;
  fYStep  = YStep;
  fZStep  = ZStep;
  fXStop  = fXStart + (fNX - 1) * fXStep;
  fYStop  = fYStart + (fNY - 1) * fYStep;
  fZStop  = fZStart + (fNZ - 1) * fZStep;

  fHasX = NX > 1 ? true : false;
  fHasY = NY > 1 ? true : false;
  fHasZ = NZ > 1 ? true : false;

  if (fHasX && fHasY && fHasZ) {
    fDIMX = kDIMX_XYZ;
  } else if (fHasX && fHasY) {
    fDIMX = kDIMX_XY;
  } else if (fHasX && fHasZ) {
    fDIMX = kDIMX_XZ;
  } else if (fHasY && fHasZ) {
    fDIMX = kDIMX_YZ;
  } else if (fHasX) {
    fDIMX = kDIMX_X;
  } else if (fHasY) {
    fDIMX = kDIMX_Y;
  } else if (fHasZ) {
    fDIMX = kDIMX_Z;
  } else {
    std::cerr << "ERROR: error in file header format" << std::endl;
    throw std::out_of_range("invalid dimensions");
  }

  fXDIM = 0;
  if (fHasX) {
    ++fXDIM;
  }
  if (fHasY) {
    ++fXDIM;
  }
  if (fHasZ) {
    ++fXDIM;
  }

  // Reserve correct number of points in vector (slightly faster)
  fData.reserve(fNX * fNY * fNZ);




  // Temp variables for field
  double fx;
  double fy;
  double fz;

  // Fields at each parameter value
  std::vector<TVector3D> F(InFiles.size());

  // Loop over all points
  for (int ix = 0; ix != NX; ++ix) {
    for (int iy = 0; iy != NY; ++iy) {
      for (int iz = 0; iz != NZ; ++iz) {

        // Loop over all files
        for (size_t ifile = 0; ifile != InFiles.size(); ++ifile) {
          // Grab a line from input file
          std::getline(*(InFiles[ifile]), L);

          // Check we did not hit an EOF
          if (InFiles[ifile]->eof()) {
            std::cerr << "ERROR: bad input file" << std::endl;
            throw std::ifstream::failure("error reading file.  Check format");
          }

          // Read data
          S.clear();
          S.str(L);
          S >> fx >> fy >> fz;

          // Scale values?
          if (FxScaling != 1) {
            fx *= FxScaling;
          }
          if (FyScaling != 1) {
            fy *= FyScaling;
          }
          if (FzScaling != 1) {
            fz *= FzScaling;
          }

          // Check the stream did not hit an EOF
          if (S.fail()) {
            std::cerr << "ERRROR: input stream bad" << std::endl;
            throw std::ifstream::failure("error reading file.  Check format");
          }

          F[ifile].SetXYZ(fx, fy, fz);
        }

        // Interpolate and push data to memory
        TVector3D FInterpolated = this->InterpolateFields(Parameters, F, Parameter);
        FInterpolated.RotateSelfXYZ(Rotations);
        fData.push_back(FInterpolated);

      }
    }
  }


  // Close all files
  for (std::vector<std::ifstream*>::iterator it = InFiles.begin(); it != InFiles.end(); ++it) {
    (*it)->close();
    delete *it;
  }
  InFiles.clear();

  // Store Rotations and Translation
  fRotated = Rotations;
  fTranslation = Translation;

  return;
}
















TVector3D TField3D_Grid::InterpolateFields (std::vector<double> const& Parameters, std::vector<TVector3D> const& Fields, double const Parameter)
{
  // Interpolate to find field at given Parameter value

  // Check size of parameters vector.  Must be at least 2
  if (Parameters.size() < 2) {
    throw;
  }

  // Before and After index of Parameter
  int    Before;
  int    After;
  bool   HasAfter = false;

  // For each desired point find the bin before and after the desired Z position
  for (size_t i = 0; i != Parameters.size(); ++i) {
    double const V = Parameters[i];

    // If outside of the range throw for now
    if (i == 0 && Parameter < V) {
      throw;
    }

    if (V < Parameter) {
      Before= i;
      
    } else if (V >= Parameter) {
      HasAfter = true;
      After = i;
      break;
    }
  }

  // Make sure that an after point was found, otherwise outside of range
  if (!HasAfter) {
    throw;
  }

  // How far in and slope
  double const Difference = Parameter - Parameters[Before];
  TVector3D const Slope = (Fields[After] - Fields[Before]) / (Parameters[After] - Parameters[Before]);

  return Fields[Before] + Slope * Difference;
}
















bool TField3D_Grid::CompareField1D (std::array<double, 4> const& A, std::array<double, 4> const& B)
{
  // This function is used for sorting the field in 'position' cood.  It is a comparison function
  return A[0] < B[0];
}

bool TField3D_Grid::CompareMappingElements (std::pair<double, std::string> const& A, std::pair<double, std::string> const& B)
{
  // This function is used for sorting the field in 'position' cood.  It is a comparison function
  return A.first < B.first;
}

