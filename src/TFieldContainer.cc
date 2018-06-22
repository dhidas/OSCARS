////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Thu Sep 22 08:19:53 EDT 2016
//
////////////////////////////////////////////////////////////////////

#include "TFieldContainer.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <stdexcept>



TFieldContainer::TFieldContainer ()
{
  // Default constructor
}



TFieldContainer::TFieldContainer (TField* F)
{
  // Construct me with a field
  this->AddField(F);
}



TFieldContainer::~TFieldContainer ()
{
  // Default destructor.  I own everything you have passed me.  Make no mistake there!!
  this->Clear();
}



void TFieldContainer::AddField (TField* F)
{
  // Construct me with a field
  fFields.push_back(F);
}




void TFieldContainer::RemoveField (std::string const& Name)
{
  // Remove all fields that match the input name

  size_t i = 0;
  while (i < fFields.size()) {
    if (fFields[i]->GetName() == Name) {

      // Delete TField
      delete fFields[i];

      // Remove pointer from vector
      fFields.erase( fFields.begin() + i );
    } else {
      ++i;
    }
  }
  return;
}




TVector3D TFieldContainer::GetF (double const X, double const Y, double const Z, double const T, std::string const& Name) const
{
  TVector3D Sum(0, 0, 0);

  return this->GetF(TVector3D(X, Y, Z), T, Name);
}




TVector3D TFieldContainer::GetF (TVector3D const& X, double const T, std::string const& Name) const
{
  TVector3D Sum(0, 0, 0);

  // Is this named or not?
  if (Name == "") {
    // Loop over all Fields for summing fields
    for (std::vector<TField*>::const_iterator it = fFields.begin(); it != fFields.end(); ++it) {
      Sum += (*it)->GetF(X, T);
    }
  } else {
    // Loop over Fields for summing fields, select any with this exact name
    for (std::vector<TField*>::const_iterator it = fFields.begin(); it != fFields.end(); ++it) {
      if ((*it)->GetName() == Name) {
        Sum += (*it)->GetF(X, T);
      }
    }
  }

  return Sum;
}




TField const& TFieldContainer::GetField (size_t const i) const
{
  // Return reference to field (const)
  return *fFields[i];
}




size_t TFieldContainer::GetNFields () const
{
  // Return the number of fields input
  return fFields.size();
}




void TFieldContainer::Clear ()
{
  for (std::vector<TField*>::iterator it = fFields.begin(); it != fFields.end(); ++it) {
    if (*it != 0x0) {
      delete *it;
    }
  }

  fFields.clear();

  return;
}





void TFieldContainer::WriteToFile (std::string const& OutFileName,
                                   std::string const& OutFormat,
                                   TVector2D const& XLim,
                                   int const NX,
                                   TVector2D const& YLim,
                                   int const NY,
                                   TVector2D const& ZLim,
                                   int const NZ,
                                   std::string const Comment)
{
  // Write the magnetic field in a given range to an output file of the chosen format

  // Field for writing
  TVector3D B;

  // Position
  TVector3D X;

  // Open file for output
  std::ofstream of(OutFileName.c_str());
  if (!of.is_open()) {
    throw std::ifstream::failure("TFieldContainer::WriteToFile failed to open file for writing: " + OutFileName);
  }

  std::string CommentNoCRLF = Comment;
  std::replace(CommentNoCRLF.begin(), CommentNoCRLF.end(), '\n', ' ');
  std::replace(CommentNoCRLF.begin(), CommentNoCRLF.end(), '\r', ' ');

  // OSCARS format is text by default
  if (OutFormat == "OSCARS") {

    // If no comment print default comment
    if (CommentNoCRLF == "") {
      of << "# " << OutFormat << std::endl;
    } else {
      of << "# " << CommentNoCRLF << std::endl;
    }

    int const MyNX = NX == 0 ? 1 : NX;
    int const MyNY = NY == 0 ? 1 : NY;
    int const MyNZ = NZ == 0 ? 1 : NZ;

    double const XStep = MyNX == 1 ? 0 : (XLim[1] - XLim[0]) / (NX - 1);
    double const YStep = MyNY == 1 ? 0 : (YLim[1] - YLim[0]) / (NY - 1);
    double const ZStep = MyNZ == 1 ? 0 : (ZLim[1] - ZLim[0]) / (NZ - 1);

    int const Width = 15;
    of << std::setw(Width) << std::left << XLim[0]     << "  # X Start position"      << std::endl;
    of << std::setw(Width) << std::left << XStep       << "  # Step size in X"        << std::endl;
    of << std::setw(Width) << std::left << MyNX        << "  # Number of points in X" << std::endl;
    of << std::setw(Width) << std::left << YLim[0]     << "  # Y Start position"      << std::endl;
    of << std::setw(Width) << std::left << YStep       << "  # Step size in Y"        << std::endl;
    of << std::setw(Width) << std::left << MyNY        << "  # Number of points in Y" << std::endl;
    of << std::setw(Width) << std::left << ZLim[0]     << "  # Z Start position"      << std::endl;
    of << std::setw(Width) << std::left << ZStep       << "  # Step size in Z"        << std::endl;
    of << std::setw(Width) << std::left << MyNZ        << "  # Number of points in Z" << std::endl;

    // Set output format
    of << std::scientific;

    // Loop over all points and output
    for (int i = 0; i < MyNX; ++i) {
      for (int j = 0; j < MyNY; ++j) {
        for (int k = 0; k < MyNZ; ++k) {

          // Set current position
          X.SetXYZ(XLim[0] + XStep * i, YLim[0] + YStep * j, ZLim[0] + ZStep * k);

          // Get B Field
          B = this->GetF(X);

          // Print field to file
          of << B.GetX() << " " << B.GetY() << " " << B.GetZ() << std::endl;
        }
      }
    }
  } else if (std::string(OutFormat.begin(), OutFormat.begin() + 8) == "OSCARS1D") {

    // Determine output format

    // And this is for which order they come in
    std::vector<int> Order(6, -1);

    // Make it a stream and set it to the format string minus the OSCARS1D
    std::string const FormatString(OutFormat.begin() + 9, OutFormat.end());
    std::istringstream s;
    s.str(FormatString);

    // String for identifier
    std::string c;

    // Counts
    int index = 0;
    int XDIM = 0;
    int BDIM = 0;
    int N = 0;


    // Starting point and step
    TVector3D StartPoint(0, 0, 0);
    TVector3D Step(0, 0, 0);

    // Look at all input
    while (s >> c) {

      if (index > 3) {
        std::cerr << "ERROR: spatial or B-field dimensions are too large(index>3)" << std::endl;
        throw std::out_of_range("spatial or B-field dimensions are too large(index>3)");
      }

      // Check if it is XYZBxByBz and in which order
      if (c == "X") {
        ++XDIM;
        StartPoint.SetXYZ(XLim[0], 0, 0);
        Step.SetXYZ( (XLim[1] - XLim[0]) / (NX - 1), 0, 0);
        Order[index] = 0;
        N = NX;
        ++index;
      } else if (c == "Y") {
        ++XDIM;
        StartPoint.SetXYZ(0, YLim[0], 0);
        Step.SetXYZ(0, (YLim[1] - YLim[0]) / (NY - 1), 0);
        Order[index] = 1;
        N = NY;
        ++index;
      } else if (c == "Z") {
        ++XDIM;
        StartPoint.SetXYZ(0, 0, ZLim[0]);
        Step.SetXYZ(0, 0, (ZLim[1] - ZLim[0]) / (NZ - 1));
        Order[index] = 2;
        N = NZ;
        ++index;
      } else if (c == "Bx" || c == "Ex" || c == "Fx") {
        ++BDIM;
        Order[index] = 3;
        ++index;
      } else if (c == "By" || c == "Ey" || c == "Fy") {
        ++BDIM;
        Order[index] = 4;
        ++index;
      } else if (c == "Bz" || c == "Ez" || c == "Fz") {
        ++BDIM;
        Order[index] = 5;
        ++index;
      } else {
        std::cerr << "ERROR: Incorrect format" << std::endl;
        throw std::invalid_argument("only excepts X Y Z Fx Fy Fz");
      }
    }


    // How many cols to ouout
    int OutputCount = XDIM + BDIM;

    // At the moment only support 1D irregular grid
    if (XDIM > 1) {
      std::cerr << "ERROR: spatial or B-field dimensions are too large(>3)" << std::endl;
      throw std::out_of_range("spatial or B-field dimensions are too large");
    }

    // Check for correct input configuration
    if (StartPoint.Mag() == 0 || Step.Mag() == 0 || N <= 0) {
      std::cerr << "ERROR: limits or npoints not correctly defined" << std::endl;
      throw std::out_of_range("limits or npoints not correctly defined");
    }


    // Vector of outputs
    std::vector<double> Outputs(6, 0);

    // If no comment print default comment
    if (CommentNoCRLF == "") {
      of << "# " << OutFormat << std::endl;
    } else {
      of << "# " << CommentNoCRLF << std::endl;
    }

    // Set output format
    of << std::scientific;

    // Loop over all points and output
    for (int i = 0; i < N; ++i) {

      // Set current position
      X = StartPoint + Step * (double) i;

      // Get B Field
      B = this->GetF(X);

      Outputs[0] = X.GetX();
      Outputs[1] = X.GetY();
      Outputs[2] = X.GetZ();
      Outputs[3] = B.GetX();
      Outputs[4] = B.GetY();
      Outputs[5] = B.GetZ();

      for (int i = 0; i != 6; ++i) {
        if (Order[i] != -1) {
          of << Outputs[Order[i]];
          if (i != OutputCount -1) {
            of << " ";
          }
        }
      }
      of << std::endl;

    }

  } else if (OutFormat == "SRW") {

    if (CommentNoCRLF == "") {
      of << "# SRW format file generated by OSCARS" << std::endl;
    } else {
      of << "# " << CommentNoCRLF << std::endl;
    }

    int const MyNX = NX == 0 ? 1 : NX;
    int const MyNY = NY == 0 ? 1 : NY;
    int const MyNZ = NZ == 0 ? 1 : NZ;

    double const XStep = MyNX == 1 ? 0 : (XLim[1] - XLim[0]) / (NX - 1);
    double const YStep = MyNY == 1 ? 0 : (YLim[1] - YLim[0]) / (NY - 1);
    double const ZStep = MyNZ == 1 ? 0 : (ZLim[1] - ZLim[0]) / (NZ - 1);


    of << "#\t" << XLim.GetX() << "\t# X Start position" << std::endl;
    of << "#\t" << XStep << std::endl;
    of << "#\t" << MyNX << std::endl;
    of << "#\t" << YLim.GetX() << "  #\tY Start position" << std::endl;
    of << "#\t" << YStep << std::endl;
    of << "#\t" << MyNY << std::endl;
    of << "#\t" << ZLim.GetX() << "  #\tZ Start position" << std::endl;
    of << "#\t" << ZStep << std::endl;
    of << "#\t" << MyNZ << std::endl;

    // Set output format
    of << std::scientific;

    // Loop over all points and output
    for (int k = 0; k < MyNZ; ++k) {
      for (int j = 0; j < MyNY; ++j) {
        for (int i = 0; i < MyNX; ++i) {

          // Set current position
          X.SetXYZ(XLim[0] + XStep * i, YLim[0] + YStep * j, ZLim[0] + ZStep * k);

          // Get B Field
          B = this->GetF(X);

          // Print field to file
          of << B.GetX() << "\t" << B.GetY() << "\t" << B.GetZ() << std::endl;
        }
      }
    }
  } else if (OutFormat == "SPECTRA") {

    if (CommentNoCRLF == "") {
      of << "# SPECTRA format file generated by OSCARS" << std::endl;
    } else {
      of << "# " << CommentNoCRLF << std::endl;
    }

    int const MyNX = NX == 0 ? 1 : NX;
    int const MyNY = NY == 0 ? 1 : NY;
    int const MyNZ = NZ == 0 ? 1 : NZ;

    double const XStep = MyNX == 1 ? 0 : (XLim[1] - XLim[0]) / (NX - 1);
    double const YStep = MyNY == 1 ? 0 : (YLim[1] - YLim[0]) / (NY - 1);
    double const ZStep = MyNZ == 1 ? 0 : (ZLim[1] - ZLim[0]) / (NZ - 1);


    // Header information
    of << XStep*1000. << " " << YStep*1000. << " " << ZStep*1000. << " " << MyNX << " " << MyNY << " " << MyNZ << std::endl;

    // Set output format
    of << std::scientific;

    // Loop over all points and output
    for (int i = 0; i < MyNX; ++i) {
      for (int j = 0; j < MyNY; ++j) {
        for (int k = 0; k < MyNZ; ++k) {

          // Set current position
          X.SetXYZ(XLim[0] + XStep * i, YLim[0] + YStep * j, ZLim[0] + ZStep * k);

          // Get B Field
          B = this->GetF(X);

          // Print field to file
          of << B.GetX() << " " << B.GetY() << " " << B.GetZ() << std::endl;
        }
      }
    }
  }


  // Close output file
  of.close();

  return;
}




void TFieldContainer::WriteToFileBinary (std::string const& OutFileName,
                                         std::string const& OutFormat,
                                         TVector2D const& XLim,
                                         int const NX,
                                         TVector2D const& YLim,
                                         int const NY,
                                         TVector2D const& ZLim,
                                         int const NZ,
                                         std::string const Comment,
                                         int const Version)
{
  // Zero and negative default to the most current version
  if (Version == 1 || Version <= 0) {
    this->WriteToFileBinary_v1(OutFileName, OutFormat, XLim, NX, YLim, NY, ZLim, NZ, Comment);
  } else {
    throw std::invalid_argument("version number for output is unknown");
  }
  
  return;
}


void TFieldContainer::WriteToFileBinary_v1 (std::string const& OutFileName,
                                            std::string const& OutFormat,
                                            TVector2D const& XLim,
                                            int const NX,
                                            TVector2D const& YLim,
                                            int const NY,
                                            TVector2D const& ZLim,
                                            int const NZ,
                                            std::string const Comment)
{
  // Write the magnetic field in a given range to an output file of the chosen format

  // File version number
  int const FileVersionNumber = 1;

  // Field for writing
  TVector3D B;

  // Position
  TVector3D X;

  // Open file for output
  std::ofstream of(OutFileName.c_str(), std::ios::binary);
  if (!of.is_open()) {
    throw std::ofstream::failure("cannot open file for writing binary format");
  }

  std::string CommentNoCRLF = Comment;
  std::replace(CommentNoCRLF.begin(), CommentNoCRLF.end(), '\n', ' ');
  std::replace(CommentNoCRLF.begin(), CommentNoCRLF.end(), '\r', ' ');

  if (CommentNoCRLF == "") {
    CommentNoCRLF = "OSCARS";
  }

  // Header 0: Number of characters in comment, then comment
  int const NCommentChars = (int) CommentNoCRLF.size();
  of.write((char*) &NCommentChars, sizeof(int));
  of.write(CommentNoCRLF.c_str() , CommentNoCRLF.size());


  // Header 1: Version number
  of.write((char*) &FileVersionNumber, sizeof(int));

  // Header 2: Formatting
  int const NFormatChars = (int) OutFormat.size();
  of.write((char*) &NFormatChars, sizeof(int));
  of.write(OutFormat.c_str() , OutFormat.size() * sizeof(char));


  // OSCARS format is text by default
  if (OutFormat == "OSCARS") {
    int const MyNX = NX == 0 ? 1 : NX;
    int const MyNY = NY == 0 ? 1 : NY;
    int const MyNZ = NZ == 0 ? 1 : NZ;

    float const XStep = MyNX == 1 ? 0 : (XLim[1] - XLim[0]) / (NX - 1);
    float const YStep = MyNY == 1 ? 0 : (YLim[1] - YLim[0]) / (NY - 1);
    float const ZStep = MyNZ == 1 ? 0 : (ZLim[1] - ZLim[0]) / (NZ - 1);

    float const XStart = XLim[0];
    float const YStart = YLim[0];
    float const ZStart = ZLim[0];

    // Header 3: start, step, n for x, y, z
    of.write((char*) &XStart,  sizeof(float));
    of.write((char*) &XStep,   sizeof(float));
    of.write((char*) &MyNX,    sizeof(int));
    of.write((char*) &YStart,  sizeof(float));
    of.write((char*) &YStep,   sizeof(float));
    of.write((char*) &MyNY,    sizeof(int));
    of.write((char*) &ZStart,  sizeof(float));
    of.write((char*) &ZStep,   sizeof(float));
    of.write((char*) &MyNZ,    sizeof(int));

    // Loop over all points and output
    for (int i = 0; i < MyNX; ++i) {
      for (int j = 0; j < MyNY; ++j) {
        for (int k = 0; k < MyNZ; ++k) {

          // Set current position
          X.SetXYZ(XLim[0] + XStep * i, YLim[0] + YStep * j, ZLim[0] + ZStep * k);

          // Get B Field
          B = this->GetF(X);

          float const Bx = B.GetX();
          float const By = B.GetY();
          float const Bz = B.GetZ();

          // Print field to file
          of.write((char*) &Bx, sizeof(float));
          of.write((char*) &By, sizeof(float));
          of.write((char*) &Bz, sizeof(float));
        }
      }
    }
  } else if (std::string(OutFormat.begin(), OutFormat.begin() + 8) == "OSCARS1D") {

    // Determine output format

    // And this is for which order they come in
    std::vector<int> Order(6, -1);

    // Make it a stream and set it to the format string minus the OSCARS1D
    std::string const FormatString(OutFormat.begin() + 9, OutFormat.end());
    std::istringstream s;
    s.str(FormatString);

    // String for identifier
    std::string c;

    // Counts
    int index = 0;
    int XDIM = 0;
    int BDIM = 0;
    int N = 0;


    // Starting point and step
    TVector3D StartPoint(0, 0, 0);
    TVector3D Step(0, 0, 0);

    // Look at all input
    while (s >> c) {

      if (index > 3) {
        std::cerr << "ERROR: spatial or field dimensions are too large(index>3)" << std::endl;
        throw std::out_of_range("spatial or field dimensions are too large(index>3)");
      }

      // Check if it is XYZBxByBz and in which order
      if (c == "X") {
        ++XDIM;
        StartPoint.SetXYZ(XLim[0], 0, 0);
        Step.SetXYZ( (XLim[1] - XLim[0]) / (NX - 1), 0, 0);
        Order[index] = 0;
        N = NX;
        ++index;
      } else if (c == "Y") {
        ++XDIM;
        StartPoint.SetXYZ(0, YLim[0], 0);
        Step.SetXYZ(0, (YLim[1] - YLim[0]) / (NY - 1), 0);
        Order[index] = 1;
        N = NY;
        ++index;
      } else if (c == "Z") {
        ++XDIM;
        StartPoint.SetXYZ(0, 0, ZLim[0]);
        Step.SetXYZ(0, 0, (ZLim[1] - ZLim[0]) / (NZ - 1));
        Order[index] = 2;
        N = NZ;
        ++index;
      } else if (c == "Bx" || c == "Ex" || c == "Fx") {
        ++BDIM;
        Order[index] = 3;
        ++index;
      } else if (c == "By" || c == "Ey" || c == "Fy") {
        ++BDIM;
        Order[index] = 4;
        ++index;
      } else if (c == "Bz" || c == "Ez" || c == "Fz") {
        ++BDIM;
        Order[index] = 5;
        ++index;
      } else {
        std::cerr << "ERROR: Incorrect format" << std::endl;
        throw std::invalid_argument("only excepts X Y Z Fx Fy Fz");
      }
    }

    // At the moment only support 1D irregular grid
    if (XDIM > 1) {
      std::cerr << "ERROR: spatial or B-field dimensions are too large(>3)" << std::endl;
      throw std::out_of_range("spatial or B-field dimensions are too large");
    }

    // Check for correct input configuration
    if (StartPoint.Mag() == 0 || Step.Mag() == 0 || N <= 0) {
      std::cerr << "ERROR: limits or npoints not correctly defined" << std::endl;
      throw std::out_of_range("limits or npoints not correctly defined");
    }

    // Vector of outputs
    std::vector<float> Outputs(6, 0);

    // Loop over all points and output
    for (int i = 0; i < N; ++i) {

      // Set current position
      X = StartPoint + Step * (float) i;

      // Get B Field
      B = this->GetF(X);

      Outputs[0] = X.GetX();
      Outputs[1] = X.GetY();
      Outputs[2] = X.GetZ();
      Outputs[3] = B.GetX();
      Outputs[4] = B.GetY();
      Outputs[5] = B.GetZ();

      for (int i = 0; i != 6; ++i) {
        if (Order[i] != -1) {
          of.write((char*) &(Outputs[Order[i]]), sizeof(float));
        }
      }
    }
  } else {
    // Unknown format
    throw std::invalid_argument("Unknown output format specified");
  }

  // Close output file
  of.close();

  return;
}

