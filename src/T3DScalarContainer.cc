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




void T3DScalarContainer::Clear ()
{
  // Clear all contents from the container

  fValues.clear();
  fCompensation.clear();

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
    throw;
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
        if (f[i].eof()) {

          // Change the done state to stop reading files
          NotDone = false;

          // This must be file index 0
          if (i != 0) {
            throw;
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
        if (f[i].eof()) {

          // Change the done state to stop reading files
          NotDone = false;

          // This must be file index 0
          if (i != 0) {
            throw;
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
    throw;
  }


  // Close all files
  for (size_t i = 0; i != Files.size(); ++i) {
    f[i].close();
  }

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
    throw;
  }

  // Double number of files for averaging
  double const N = (double) Files.size();

  // Open all files in a vector, check they are open,
  std::vector<std::ifstream> f(Files.size());
  for (size_t i = 0; i != Files.size(); ++i) {
    f[i].open(Files[i].c_str(), std::ios::binary);
    if (!f[i].is_open()) {
      throw;
    }
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
        f[i].read( (char*)  &X, sizeof(double));
        f[i].read( (char*)  &Y, sizeof(double));
        f[i].read( (char*)  &V, sizeof(double));

        // If we hit an eof we are done.
        if (f[i].eof()) {

          // Change the done state to stop reading files
          NotDone = false;

          // This must be file index 0
          if (i != 0) {
            throw;
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
        f[i].read( (char*)  &X, sizeof(double));
        f[i].read( (char*)  &Y, sizeof(double));
        f[i].read( (char*)  &Z, sizeof(double));
        f[i].read( (char*)  &V, sizeof(double));

        // If we hit an eof we are done.
        if (f[i].eof()) {

          // Change the done state to stop reading files
          NotDone = false;

          // This must be file index 0
          if (i != 0) {
            throw;
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
    throw;
  }


  // Close all files
  for (size_t i = 0; i != Files.size(); ++i) {
    f[i].close();
  }

  return;
}







void T3DScalarContainer::WriteToFileText (std::string const& OutFileName, int const Dimension)
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


void T3DScalarContainer::WriteToFileBinary (std::string const& OutFileName, int const Dimension)
{
  // Write the data in simple binary format.  I would like to use machine independent types for this eventually

  // Open file for writing
  std::ofstream of(OutFileName.c_str(), std::ios::binary);
  if (!of.is_open()) {
    throw;
  }

  // Variables for writing
  double X;
  double Y;
  double Z;
  double V;

  if (Dimension == 2) {
    for (size_t i = 0; i != this->GetNPoints(); ++i) {
      TVector3D const& Obs = this->GetPoint(i).GetX();
      X = Obs.GetX();
      Y = Obs.GetY();
      V = this->GetPoint(i).GetV();

      of.write((char*) &X, sizeof(double));
      of.write((char*) &Y, sizeof(double));
      of.write((char*) &V, sizeof(double));
    }
  } else if (Dimension == 3) {
    for (size_t i = 0; i != this->GetNPoints(); ++i) {
      TVector3D const& Obs = this->GetPoint(i).GetX();
      X = Obs.GetX();
      Y = Obs.GetY();
      Z = Obs.GetZ();
      V = this->GetPoint(i).GetV();

      of.write((char*) &X, sizeof(double));
      of.write((char*) &Y, sizeof(double));
      of.write((char*) &Z, sizeof(double));
      of.write((char*) &V, sizeof(double));
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
    throw;
  }

  return fValues[i];
}







