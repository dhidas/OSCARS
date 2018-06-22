#include "T3DScalarContainer.h"



T3DScalarContainer::T3DScalarContainer ()
{
  // Default constructor
}




T3DScalarContainer::~T3DScalarContainer ()
{
  // Destruction!!!
}




void T3DScalarContainer::AddPoint (TVector3D const& X, double const V)
{
  fValues.push_back( T3DScalar(X, V) );
  fCompensation.push_back(0);

  if (fValues.size() > fNotConverged.size() * 8 * sizeof(int)) {
    fNotConverged.push_back(0);
  }
  return;
}




void T3DScalarContainer::AddToPoint (size_t const i, double const V)
{
  // Compensated sum for adding to points

  // Check that the point is within range
  if (i >= fValues.size()) {
    throw std::length_error("T3DScalarContainer::AddtoPoint index out of range");
  }

  double Sum = fValues[i].GetV();
  double y = V - fCompensation[i];
  double t = Sum + y;
  fCompensation[i] = (t - Sum) - y;
  fValues[i].SetV(t);

  return;
}




void T3DScalarContainer::SetNotConverged (size_t const i)
{
  // Set the converged bit for this point

  size_t const VectorIndex = i / (8 * sizeof(int));
  if (VectorIndex >= fNotConverged.size()) {
    throw std::out_of_range("T3DScalarContainer::SetNotConverged index out of range");
  }

  int const Bit = (0x1 << (i % (8 * sizeof(int))));

  fNotConverged[VectorIndex] |= Bit;

  return;
}



bool T3DScalarContainer::AllConverged () const
{
  for (std::vector<int>::const_iterator it = fNotConverged.begin(); it != fNotConverged.end(); ++it) {
    if (*it != 0x0) {
      return false;
    }
  }
  return true;
}



void T3DScalarContainer::Clear ()
{
  // Clear all contents from the container

  fValues.clear();
  fCompensation.clear();
  fNotConverged.clear();

  return;
}





void T3DScalarContainer::AverageFromFilesText (std::vector<std::string> const& Files, int const Dimension)
{
  // Average the inpput from all text files.  Text files must be the same points in space
  // arranged in exactly the same format.  The spatial points are taken from only the first file
  // in the vector

  // Clear my contents
  this->Clear();

  // Check that we have at least one file!
  if (Files.size() < 1) {
    throw std::length_error("no files specified");
  }

  // Double number of files for averaging
  double const N = (double) Files.size();

  // Open all files
  std::vector<std::ifstream> f(Files.size());
  for (size_t i = 0; i != Files.size(); ++i) {
    f[i].open(Files[i].c_str());
  }

  // Variables used for writing to file
  double X, Y, Z, V;

  // Are we done reading yet
  bool NotDone = true;

  // Keep track of which point we are on
  size_t ip = 0;

  // Must be 2 or 3 D at the moment
  if (Dimension == 2) {

    // Loop over all points until done
    while (NotDone) {

      // For each point loop over files and average
      for (size_t i = 0; i != f.size(); ++i) {

        // Read data from current file
        f[i] >> X >> Y >> V;

        // If we hit an eof we are done.
        if (f[i].fail()) {

          // Change the done state to stop reading files
          NotDone = false;

          // This must be file index 0
          if (i != 0) {
            throw std::length_error("Inconsistent file format found in file " + Files[i]);
          }

          break;
        }

        // Add point to self
        if (i == 0) {
          this->AddPoint(TVector3D(X, Y, 0), V/N);
        } else {
          this->AddToPoint(ip, V/N);
        }
      }

      // Increment point counter
      ++ip;
    }
  } else if (Dimension == 3) {

    // Loop over all points until done
    while (NotDone) {

    // Loop over all points until done
      for (size_t i = 0; i != f.size(); ++i) {

        // Read data from current file
        f[i] >> X >> Y >> Z >> V;

        // If we hit an eof we are done.
        if (f[i].fail()) {

          // Change the done state to stop reading files
          NotDone = false;

          // This must be file index 0
          if (i != 0) {
            throw std::length_error("Inconsistent file format found in file " + Files[i]);
          }

          break;
        }

        // Add point to self
        if (i == 0) {
          this->AddPoint(TVector3D(X, Y, Z), V/N);
        } else {
          this->AddToPoint(ip, V/N);
        }
      }

      // Increment point counter
      ++ip;
    }
  } else {
    throw std::invalid_argument("'dim' must be 2 or 3");
  }


  // Close all files
  for (size_t i = 0; i != Files.size(); ++i) {
    f[i].close();
  }

  fNotConverged.clear();
  fNotConverged.resize(fValues.size() / (8 * sizeof(int)) + 1, 0);

  return;
}







void T3DScalarContainer::AverageFromFilesBinary (std::vector<std::string> const& Files, int const Dimension)
{
  // Average the inpput from all binary files.  Binary files must be the same points in space
  // arranged in exactly the same format.  The spatial points are taken from only the first file
  // in the vector

  // Clear my contents
  this->Clear();

  // Check that we have at least one file!
  if (Files.size() < 1) {
    throw std::length_error("T3DScalarContainer::AverageFromFilesBinary no files given");
  }

  // Double number of files for averaging
  double const N = (double) Files.size();

  // Open all files in a vector, check they are open,
  std::vector<std::ifstream> f(Files.size());
  for (size_t i = 0; i != Files.size(); ++i) {
    f[i].open(Files[i].c_str(), std::ios::binary);
    if (!f[i].is_open()) {
      throw std::ifstream::failure("Cannot open input file: " + Files[i]);
    }
  }

  // Variables used for reading from file
  double X;
  double Y;
  double Z;
  double V;

  // Are we done reading yet
  bool NotDone = true;

  // Keep track of which point we are on
  size_t ip = 0;

  // For reading data
  float a = 0;
  float b = 0;
  float c = 0;
  float d = 0;

  // Must be 2 or 3 D at the moment
  if (Dimension == 2) {

    // Loop over all points until done
    while (NotDone) {

      // For each point loop over files and average
      for (size_t i = 0; i != f.size(); ++i) {

        // Read data from current file
        f[i].read( (char*)  &a, sizeof(float));
        f[i].read( (char*)  &b, sizeof(float));
        f[i].read( (char*)  &d, sizeof(float));

        X = (double) a;
        Y = (double) b;
        V = (double) d;

        // If we hit an eof we are done.
        if (f[i].fail()) {

          // Change the done state to stop reading files
          NotDone = false;

          // This must be file index 0
          if (i != 0) {
            throw std::logic_error("T3DScalarContainer::AverageFromFilesBinary index here must be 0.  Please report this.");
          }

          break;
        }

        // Add point to self
        if (i == 0) {
          this->AddPoint(TVector3D(X, Y, 0), V/N);
        } else {
          this->AddToPoint(ip, V/N);
        }
      }

      // Increment point counter
      ++ip;
    }
  } else if (Dimension == 3) {

    // Loop over all points until done
    while (NotDone) {

    // Loop over all points until done
      for (size_t i = 0; i != f.size(); ++i) {

        // Read data from current file
        f[i].read( (char*)  &a, sizeof(float));
        f[i].read( (char*)  &b, sizeof(float));
        f[i].read( (char*)  &c, sizeof(float));
        f[i].read( (char*)  &d, sizeof(float));

        X = (double) a;
        Y = (double) b;
        Z = (double) c;
        V = (double) d;

        // If we hit an eof we are done.
        if (f[i].fail()) {

          // Change the done state to stop reading files
          NotDone = false;

          // This must be file index 0
          if (i != 0) {
            throw std::logic_error("T3DScalarContainer::AverageFromFilesBinary index here must be 0..  Please report this.");
          }

          break;
        }

        // Add point to self
        if (i == 0) {
          this->AddPoint(TVector3D(X, Y, Z), V/N);
        } else {
          this->AddToPoint(ip, V/N);
        }
      }

      // Increment point counter
      ++ip;
    }
  } else {
    throw std::out_of_range("T3DScalarContainer::AverageFromFilesBinary I do not know how to calculate in this dimension");
  }


  // Close all files
  for (size_t i = 0; i != Files.size(); ++i) {
    f[i].close();
  }

  fNotConverged.clear();
  fNotConverged.resize(fValues.size() / (8 * sizeof(int)) + 1, 0);

  return;
}





void T3DScalarContainer::WeightAll (double const Weight)
{
  for (size_t i = 0; i != fValues.size(); ++i) {
    fValues[i].SetV( fValues[i].GetV() * Weight );
    fCompensation[i] = 0;
  }
}







void T3DScalarContainer::WriteToFileText (std::string const& OutFileName,
                                          int const Dimension)
{
  // Write to file in text format
  // If writing to a file, open it and set to scientific output
  std::ofstream of(OutFileName.c_str());
  if (!of.is_open()) {
    throw std::ofstream::failure("cannot open output file");
  }
  of << std::scientific;

  for (size_t i = 0; i != this->GetNPoints(); ++i) {
    TVector3D const& Obs = this->GetPoint(i).GetX();

    if (Dimension == 2) {
      of << Obs.GetX() << " " << Obs.GetY() << " " << this->GetPoint(i).GetV() << "\n";
    } else if (Dimension == 3) {
      of << Obs.GetX() << " " << Obs.GetY() << " " << Obs.GetZ() << " " << this->GetPoint(i).GetV() << "\n";
    } else {
      throw std::out_of_range("incorrect dimensions");
    }

  }

  of.close();

  return;
}


void T3DScalarContainer::WriteToFileBinary (std::string const& OutFileName,
                                            int const Dimension)
{
  // Write the data in simple binary format.  I would like to use machine independent types for this eventually

  // Open file for writing
  std::ofstream of(OutFileName.c_str(), std::ios::binary);
  if (!of.is_open()) {
    throw std::ifstream::failure("cannot open out file named: " + OutFileName);
  }

  // Variables for writing
  float X = 0;
  float Y = 0;
  float Z = 0;
  float V = 0;

  if (Dimension == 2) {
    for (size_t i = 0; i != this->GetNPoints(); ++i) {
      TVector3D const& Obs = this->GetPoint(i).GetX();
      X = (float) Obs.GetX();
      Y = (float) Obs.GetY();
      V = (float) this->GetPoint(i).GetV();

      of.write((char*) &X, sizeof(float));
      of.write((char*) &Y, sizeof(float));
      of.write((char*) &V, sizeof(float));
    }
  } else if (Dimension == 3) {
    for (size_t i = 0; i != this->GetNPoints(); ++i) {
      TVector3D const& Obs = this->GetPoint(i).GetX();
      X = (float) Obs.GetX();
      Y = (float) Obs.GetY();
      Z = (float) Obs.GetZ();
      V = (float) this->GetPoint(i).GetV();

      of.write((char*) &X, sizeof(float));
      of.write((char*) &Y, sizeof(float));
      of.write((char*) &Z, sizeof(float));
      of.write((char*) &V, sizeof(float));
    }
  } else {
    throw std::out_of_range("incorrect dimensions");
  }

  // Close file
  of.close();

  return;
}


size_t T3DScalarContainer::GetNPoints () const
{
  return fValues.size();
}




T3DScalar const& T3DScalarContainer::GetPoint (size_t const i) const
{
  if (i >= fValues.size()) {
    throw std::out_of_range("T3DScalarContainer::GetPoint index is out of range");
  }

  return fValues[i];
}







