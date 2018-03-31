////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Thu May 19 16:05:08 EDT 2016
//
////////////////////////////////////////////////////////////////////


#include "TSpectrumContainer.h"

#include "TOSCARSSR.h"

#include <iostream>
#include <fstream>
#include <stdexcept>



TSpectrumContainer::TSpectrumContainer ()
{
  // I guess you'll build it yourself by using AddPoint()
}




TSpectrumContainer::TSpectrumContainer (std::vector<double> const& V)
{
  // Constructor for an arbitrary list of points
  // V - vector of energy points in [eV]
  this->Init(V);
}


TSpectrumContainer::TSpectrumContainer (size_t const N, double const EFirst, double const ELast)
{
  // Constructor for evenly spaced points in a given energy range.
  // N - Number of points
  // EFirst - First energy point [eV]
  // ELast  - Last energy point [eV]

  this->Init(N, EFirst, ELast);
}




TSpectrumContainer::~TSpectrumContainer ()
{
  // Destruction!
}




void TSpectrumContainer::Init (size_t const N, double const EFirst, double const ELast)
{
  // If you call this with N==1 it will only use EFirst

  // N - Number of points
  // EFirst - First energy point [eV]
  // ELast  - Last energy point [eV]

  // Clear existing data and resize member to the correct size for input
  fSpectrumPoints.clear();
  fSpectrumPoints.resize(N, std::make_pair(0.0, 0.0));
  fCompensation.resize(N, 0);


  // If you have zero elements I don't see the point of this
  if (N < 1) {
    throw std::length_error("no points specified");
  }

  // If only one point just set it to the 'First' energy
  if (N == 1) {
    fSpectrumPoints[0].first = EFirst;
    return;
  }

  // Set the energy in each element
  for (size_t i = 0; i != fSpectrumPoints.size(); ++i) {
    fSpectrumPoints[i].first = EFirst + (ELast - EFirst) / (N - 1) * (double) (i);
  }

  fNotConverged.clear();
  fNotConverged.resize(fSpectrumPoints.size() / (8 * sizeof(int)) + 1, 0);

  return;
}




void TSpectrumContainer::Init (std::vector<double> const& V)
{
  // Initialize spectrum with an arbitrary vector of points
  // V - vector of energy points in [eV]

  // Clear existing data and reserve the correct amount for input
  fSpectrumPoints.clear();
  fSpectrumPoints.reserve(V.size());
  fCompensation.resize(V.size(), 0);

  // Add each input from V to the internal vector
  for (size_t i = 0; i != V.size(); ++i) {
    fSpectrumPoints.push_back( std::make_pair(V[i], 0.0) );
  }

  fNotConverged.clear();
  fNotConverged.resize(fSpectrumPoints.size() / (8 * sizeof(int)) + 1, 0);

  return;
}



void TSpectrumContainer::SetFlux (size_t const i, double const Flux)
{
  // Set the flux for a given index

  // Simple check
  if (i >= fSpectrumPoints.size()) {
    throw std::out_of_range("index beyond fSpectrum points range");
  }

  fSpectrumPoints[i].second = Flux;
  return;
}





void TSpectrumContainer::SetPoint (size_t const i, double const Energy, double const Flux)
{
  // Set the energy and flux for a given index

  // I can't decide if I want to be nice and resize or just throw...
  if (i >= fSpectrumPoints.size()) {
    throw std::out_of_range("index beyond fSpectrum points range");
  }

  fSpectrumPoints[i].first  = Energy;
  fSpectrumPoints[i].second = Flux;

  return;
}




size_t TSpectrumContainer::AddPoint (double const Energy, double const Flux)
{
  // Add an energy point to the end of the vector.
  fSpectrumPoints.push_back( std::make_pair(Energy, Flux) );
  fCompensation.push_back(0);

  if (fSpectrumPoints.size() > fNotConverged.size() * 8 * sizeof(int)) {
    fNotConverged.push_back(0);
  }

  return fSpectrumPoints.size();
}




void TSpectrumContainer::AddToFlux (size_t const i, double const Flux)
{
  // Add to the flux, using compensated summation

  // Simple check
  if (i >= fSpectrumPoints.size()) {
    throw std::out_of_range("index beyond fSpectrum points range");
  }

  double Sum = fSpectrumPoints[i].second;
  double y = Flux - fCompensation[i];
  double t = Sum + y;
  fCompensation[i] = (t - Sum) - y;
  fSpectrumPoints[i].second = t;

  return;
}




void TSpectrumContainer::SetNotConverged (size_t const i)
{
  // Set the converged bit for this point

  size_t const VectorIndex = i / (8 * sizeof(int));
  if (VectorIndex >= fNotConverged.size()) {
    throw std::length_error("not enough elements in fNotConverged");
  }

  int const Bit = (0x1 << (i % (8 * sizeof(int))));

  fNotConverged[VectorIndex] |= Bit;

  return;
}





void TSpectrumContainer::Scale (double const ScaleFactor)
{
  // Scale all flux by the input factor
  //
  for (size_t i = 0; i != fSpectrumPoints.size(); ++i) {
    fSpectrumPoints[i].second *= ScaleFactor;
  }

  return;
}




bool TSpectrumContainer::AllConverged () const
{
  for (std::vector<int>::const_iterator it = fNotConverged.begin(); it != fNotConverged.end(); ++it) {
    if (*it != 0x0) {
      return false;
    }
  }
  return true;
}






double TSpectrumContainer::GetFlux (size_t const i) const
{
  // Get flux at a given index
  return fSpectrumPoints[i].second;
}






double TSpectrumContainer::GetEnergy (size_t const i) const
{
  // Get energy at a given index
  return fSpectrumPoints[i].first;
}






double TSpectrumContainer::GetAngularFrequency (size_t const i) const
{
  // Get the angular frequency of this index (from energy)
  // UPDATE: consider removing this function
  return TOSCARSSR::EvToAngularFrequency(fSpectrumPoints[i].first);
}





size_t TSpectrumContainer::GetNPoints () const
{
  // Return the number of points in this spectrum
  return fSpectrumPoints.size();
}




void TSpectrumContainer::WriteToFileText (std::string const FileName, std::string const Header) const
{
  // Write this spectrum to a file in text format.
  // FileName - File name to write to
  // Header   - Header to print in file

  // Open output file
  std::ofstream f(FileName.c_str());

  // Check if file is open
  if (!f.is_open()) {
    throw std::ifstream::failure("cannot open file for writing");
  }

  // If the header is specified, write one!
  if (Header != "") {
    f << Header << std::endl;
  }

  // Set in scientific mode for printing
  // UPDATE: Could change this to accept c-stype formatting
  f << std::scientific;

  // Loop over spectrum and print to file
  for (std::vector<std::pair<double, double> >::const_iterator it = fSpectrumPoints.begin(); it != fSpectrumPoints.end(); ++it) {
    f << it->first << " " << it->second << std::endl;
  }

  // Close file
  f.close();

  return;
}






void TSpectrumContainer::WriteToFileBinary (std::string const FileName, std::string const Header) const
{
  // UPDATE: Actually make this a binary output
  // Write this spectrum to a file in text format.
  // FileName - File name to write to
  // Header   - Header to print in file

  // Open output file
  std::ofstream f(FileName.c_str(), std::ios::binary);

  // Check if file is open
  // UPDATE: try a more robust check
  if (!f.is_open()) {
    throw std::ifstream::failure("cannot open file for binary write");
  }

  // I suspect there is only a need for float precision for saved data
  float a = 0;
  float b = 0;

  // Loop over spectrum and print to file
  for (std::vector<std::pair<double, double> >::const_iterator it = fSpectrumPoints.begin(); it != fSpectrumPoints.end(); ++it) {
    a = (float) it->first;
    b = (float) it->second;
    f.write((char*) &a, sizeof(float));
    f.write((char*) &b, sizeof(float));
  }

  // Close file
  f.close();

  return;
}




void TSpectrumContainer::Clear ()
{
  // Clear contents
  fSpectrumPoints.clear();
  fCompensation.clear();
  fNotConverged.clear();

  return;
}





void TSpectrumContainer::AverageFromFilesText (std::vector<std::string> const& Files, std::vector<double> const& Weights)
{
  // Average the inpput from all text files.  Text files must be the same points in energy
  // arranged in exactly the same format.  The energy points are taken from only the first file
  // in the vector

  // Clear my contents
  this->Clear();

  // Check that we have at least one file!
  if (Files.size() < 1) {
    throw std::length_error("no files specified");
  }

  // Check length of weights with number of files
  if (Weights.size() != 0 && Files.size() != Weights.size()) {
    throw std::length_error("Incorrect size for weights given the number of files");
  }

  // If no weights specified
  double const W = 1. / (double) Files.size();

  // Open all files
  std::vector<std::ifstream> f(Files.size());
  for (size_t i = 0; i != Files.size(); ++i) {

    // Open file
    f[i].open(Files[i].c_str());

    // Check that file is actually open
    if (!f[i].is_open()) {
      throw std::invalid_argument("Cannot open one or more files of input");
    }
  }

  // Vector of weights for weighting
  std::vector<double> ActualWeights = Weights.size() != 0 ? Weights : std::vector<double>(Files.size(), W);

  // Variables used for writing to file
  double X = 0;
  double V = 0;

  // Are we done reading yet
  bool NotDone = true;

  // Keep track of which point we are on
  size_t ip = 0;

  // Loop over all points until done
  while (NotDone) {

    // For each point loop over files and average
    for (size_t i = 0; i != f.size(); ++i) {

      // Read data from current file
      f[i] >> X >> V;

      // If we hit an eof we are done.
      if (f[i].fail()) {

        // Change the done state to stop reading files
        NotDone = false;

        // This must be file index 0
        if (i != 0) {
          throw std::length_error("files are not the same length");
        }

        break;
      }

      // Add point to self
      if (i == 0) {
        this->AddPoint(X, V * ActualWeights[i]);
      } else {
        this->AddToFlux(ip, V * ActualWeights[i]);
      }
    }

    // Increment point counter
    ++ip;
  }

  // Close all files
  for (size_t i = 0; i != Files.size(); ++i) {
    f[i].close();
  }

  return;
}





void TSpectrumContainer::AverageFromFilesBinary (std::vector<std::string> const& Files, std::vector<double> const& Weights)
{
  // Average the inpput from all text files.  Text files must be the same points in energy
  // arranged in exactly the same format.  The energy points are taken from only the first file
  // in the vector

  // Clear my contents
  this->Clear();

  // Check that we have at least one file!
  if (Files.size() < 1) {
    throw std::length_error("no files specified");
  }

  // Check length of weights with number of files
  if (Weights.size() != 0 && Files.size() != Weights.size()) {
    throw std::length_error("Incorrect size for weights given the number of files");
  }

  // Open all files
  std::vector<std::ifstream> f(Files.size());
  for (size_t i = 0; i != Files.size(); ++i) {

    // Open file
    f[i].open(Files[i].c_str(), std::ios::binary);

    // Check that file is actually open
    if (!f[i].is_open()) {
      throw std::invalid_argument("Cannot open one or more files of input");
    }
  }

  // If no weights specified
  double const W = 1. / (double) Files.size();

  // Vector of weights for weighting
  std::vector<double> ActualWeights = Weights.size() != 0 ? Weights : std::vector<double>(Files.size(), W);

  // Variables used for writing to file
  double X = 0;
  double V = 0;

  // Are we done reading yet
  bool NotDone = true;

  // Keep track of which point we are on
  size_t ip = 0;

  // For reading data
  float a = 0;
  float b = 0;

  // Loop over all points until done
  while (NotDone) {

    // For each point loop over files and average
    for (size_t i = 0; i != f.size(); ++i) {

      // Read data from current file
      f[i].read( (char*)  &a, sizeof(float));
      f[i].read( (char*)  &b, sizeof(float));

      X = (double) a;
      V = (double) b;

      // If we hit an eof we are done.
      if (f[i].fail()) {

        // Change the done state to stop reading files
        NotDone = false;

        // This must be file index 0
        if (i != 0) {
          throw std::length_error("files are not the same length");
        }

        break;
      }

      // Add point to self
      if (i == 0) {
        this->AddPoint(X, V * ActualWeights[i]);
      } else {
        this->AddToFlux(ip, V * ActualWeights[i]);
      }
    }

    // Increment point counter
    ++ip;
  }

  // Close all files
  for (size_t i = 0; i != Files.size(); ++i) {
    f[i].close();
  }

  return;
}


void TSpectrumContainer::AverageFromSpectra (std::vector<TSpectrumContainer> const& Spectra, std::vector<double> const& Weights)
{

  // Clear my contents
  this->Clear();

  if (Weights.size() != 0 && Spectra.size() != Weights.size()) {
    throw std::length_error("Incorrect size for weights given the spectra");
  }

  double const W = 1. / (double) Spectra.size();

  // Check all npoints in spectra
  size_t const NPoints = Spectra[0].GetNPoints();
  for (std::vector<TSpectrumContainer>::const_iterator it = Spectra.begin(); it != Spectra.end(); ++it) {
    if (it->GetNPoints() != NPoints) {
      throw std::length_error("Incorrect size in one of the spectra");
    }
  }

  for (size_t is = 0; is != Spectra.size(); ++is) {

    double const ThisWeight = Weights.size() != 0 ? Weights[is] : W;

    for (size_t ip = 0; ip != NPoints; ++ip) {

      double const X = Spectra[is].GetEnergy(ip);
      double const Y = Spectra[is].GetFlux(ip) * ThisWeight;

      // Add point to self
      if (is == 0) {
        this->AddPoint(X, Y);
      } else {
        this->AddToFlux(ip, Y);
      }
    }
  }



  return;
}








