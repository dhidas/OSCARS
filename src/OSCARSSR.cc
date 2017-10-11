////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Tue May 24 11:33:19 EDT 2016
//
////////////////////////////////////////////////////////////////////

#include "OSCARSSR.h"

#include <cmath>
#include <complex>
#include <fstream>
#include <sstream>
#include <thread>
#include <chrono>
#include <algorithm>
#include <fstream>

#include "TVector3DC.h"
#include "TField3D_Grid.h"
#include "TField3D_Gaussian.h"
#include "TSpectrumContainer.h"
#include "TSurfacePoints_Rectangle.h"



// External global random generator
extern TRandomA* gRandomA;




OSCARSSR::OSCARSSR ()
{
  // Default constructor
  //fParticle.SetParticleType("");
  fCTStart = 0;
  fCTStop  = 0;
  fNPointsTrajectory = 0;
  fNPointsPerMeter = 10000;

  // Set derivs function default to E&B (to avoid anything nasty)
  SetDerivativesFunction();

  // Set Global compute settings
  SetUseGPUGlobal(0);   // GPU off by default
  SetNThreadsGlobal(2); // Use N threads for calculations by default
}



OSCARSSR::~OSCARSSR ()
{
  // Destructor
}




void OSCARSSR::AddMagneticField (std::string const FileName,
                                 std::string const Format,
                                 TVector3D const& Rotations,
                                 TVector3D const& Translation,
                                 std::vector<double> const& Scaling,
                                 std::string const& Name)
{
  // Add a magnetic field from a file to the field container

  // Format string all upper-case (just in case you like to type L.C.).
  std::string FormatUpperCase = Format;
  std::transform(FormatUpperCase.begin(), FormatUpperCase.end(), FormatUpperCase.begin(), ::toupper);

  // Check that the format name is correct
  if ( (FormatUpperCase == "OSCARS" || FormatUpperCase == "SRW" || FormatUpperCase == "SPECTRA") || FormatUpperCase == "BINARY" ||
       (FormatUpperCase.size() > 8 && std::string(FormatUpperCase.begin(), FormatUpperCase.begin() + 8) == std::string("OSCARS1D"))) {

    this->fBFieldContainer.AddField( new TField3D_Grid(FileName, Format, Rotations, Translation, Scaling, Name) );

  } else {
    throw std::invalid_argument("Incorrect format in format string");
  }

  // Set the derivs function accordingly
  this->SetDerivativesFunction();

  // Clear any previous fparticle trajectory data
  this->ClearTrajectory();

  return;
}



void OSCARSSR::AddMagneticFieldInterpolated (std::vector<std::pair<double, std::string> > const& Mapping,
                                             std::string const Format,
                                             double const Parameter,
                                             TVector3D const& Rotations,
                                             TVector3D const& Translation,
                                             std::vector<double> const& Scaling,
                                             std::string const& Name)
{
  // Add a magnetic field from a file to the field container

  // Format string all upper-case (just in case you like to type L.C.).
  std::string FormatUpperCase = Format;
  std::transform(FormatUpperCase.begin(), FormatUpperCase.end(), FormatUpperCase.begin(), ::toupper);


  // Check that the format name is correct
  if ( (FormatUpperCase == "OSCARS"  || FormatUpperCase == "SRW" || FormatUpperCase == "SPECTRA") ||
     (FormatUpperCase.size() > 8 && std::string(FormatUpperCase.begin(), FormatUpperCase.begin() + 8) == std::string("OSCARS1D")) ) {

    this->fBFieldContainer.AddField( new TField3D_Grid(Mapping, Format, Parameter, Rotations, Translation, Scaling, Name) );

  } else {
    throw std::invalid_argument("Incorrect format in format string");
  }

  // Set the derivs function accordingly
  this->SetDerivativesFunction();

  // Clear any previous fparticle trajectory data
  this->ClearTrajectory();

  return;
}




void OSCARSSR::AddMagneticField (TField* Field)
{
  // Add a magnetic field from a file to the field container

  this->fBFieldContainer.AddField(Field);

  // Set the derivs function accordingly
  this->SetDerivativesFunction();

  // Clear any previous fparticle trajectory data
  this->ClearTrajectory();

  return;
}



void OSCARSSR::RemoveMagneticField (std::string const& Name)
{
  // Remove all fields with the given name

  this->fBFieldContainer.RemoveField(Name);

  return;
}




void OSCARSSR::ClearMagneticFields ()
{
  // Add a magnetic field from a file to the field container

  this->fBFieldContainer.Clear();

  // Set the derivs function accordingly
  this->SetDerivativesFunction();

  return;
}




TVector3D OSCARSSR::GetB (double const X, double const Y, double const Z) const
{
  // Return summed Bx from container
  return this->fBFieldContainer.GetF(X, Y, Z);
}




TVector3D OSCARSSR::GetB (TVector3D const& X) const
{
  // Return summed Bx from container
  return this->fBFieldContainer.GetF(X);
}










void OSCARSSR::AddElectricField (std::string const FileName,
                                 std::string const Format,
                                 TVector3D const& Rotations,
                                 TVector3D const& Translation,
                                 std::vector<double> const& Scaling,
                                 std::string const& Name)
{
  // Add a electric field from a file to the field container
  this->fEFieldContainer.AddField( new TField3D_Grid(FileName, Format, Rotations, Translation, Scaling, Name) );

  // Set the derivs function accordingly
  this->SetDerivativesFunction();

  // Clear any previous fparticle trajectory data
  this->ClearTrajectory();

  return;
}




void OSCARSSR::AddElectricFieldInterpolated (std::vector<std::pair<double, std::string> > const& Mapping,
                                             std::string const Format,
                                             double const Parameter,
                                             TVector3D const& Rotations,
                                             TVector3D const& Translation,
                                             std::vector<double> const& Scaling,
                                             std::string const& Name)
{
  // Add an electric field from a file to the field container

  // Format string all upper-case (just in case you like to type L.C.).
  std::string FormatUpperCase = Format;
  std::transform(FormatUpperCase.begin(), FormatUpperCase.end(), FormatUpperCase.begin(), ::toupper);

  // Check that the format name is correct
  if ( (FormatUpperCase == "OSCARS"  || FormatUpperCase == "SRW" || FormatUpperCase == "SPECTRA") ||
     (FormatUpperCase.size() > 8 && std::string(FormatUpperCase.begin(), FormatUpperCase.begin() + 8) == std::string("OSCARS1D")) ) {

    this->fEFieldContainer.AddField( new TField3D_Grid(Mapping, Format, Parameter, Rotations, Translation, Scaling, Name) );

  } else {
    throw std::invalid_argument("Incorrect format in format string");
  }

  // Set the derivs function accordingly
  this->SetDerivativesFunction();

  // Clear any previous fparticle trajectory data
  this->ClearTrajectory();

  return;
}




void OSCARSSR::AddElectricField (TField* F)
{
  // Add a electric field from a file to the field container
  this->fEFieldContainer.AddField(F);

  // Set the derivs function accordingly
  this->SetDerivativesFunction();

  // Clear any previous fparticle trajectory data
  this->ClearTrajectory();

  return;
}



void OSCARSSR::RemoveElectricField (std::string const& Name)
{
  // Remove all fields with the given name

  this->fEFieldContainer.RemoveField(Name);

  return;
}




void OSCARSSR::ClearElectricFields ()
{
  this->fEFieldContainer.Clear();

  // Set the derivs function accordingly
  this->SetDerivativesFunction();

  return;
}




void OSCARSSR::WriteField (std::string const& BorE,
                           std::string const& OutFileName,
                           std::string const& OutFormat,
                           TVector2D const& XLim,
                           int const NX,
                           TVector2D const& YLim,
                           int const NY,
                           TVector2D const& ZLim,
                           int const NZ,
                           std::string const& Comment)
{
  if (BorE == "B") {
    fBFieldContainer.WriteToFile(OutFileName, OutFormat, XLim, NX, YLim, NY, ZLim, NZ, Comment);
  } else if (BorE == "E") {
    fEFieldContainer.WriteToFile(OutFileName, OutFormat, XLim, NX, YLim, NY, ZLim, NZ, Comment);
  } else {
    throw std::ifstream::failure("Write failure because not B or E");
  }

  return;
}




void OSCARSSR::WriteFieldBinary (std::string const& BorE,
                                 std::string const& OutFileName,
                                 std::string const& OutFormat,
                                 TVector2D const& XLim,
                                 int const NX,
                                 TVector2D const& YLim,
                                 int const NY,
                                 TVector2D const& ZLim,
                                 int const NZ,
                                 std::string const& Comment,
                                 int const Version)
{
  if (BorE == "B") {
    fBFieldContainer.WriteToFileBinary(OutFileName, OutFormat, XLim, NX, YLim, NY, ZLim, NZ, Comment, Version);
  } else if (BorE == "E") {
    fEFieldContainer.WriteToFileBinary(OutFileName, OutFormat, XLim, NX, YLim, NY, ZLim, NZ, Comment, Version);
  } else {
    throw std::ifstream::failure("Write failure because not B or E");
  }

  return;
}




TVector3D OSCARSSR::GetE (double const X, double const Y, double const Z) const
{
  // Return summed E from container
  return this->fEFieldContainer.GetF(TVector3D(X, Y, Z));
}




TVector3D OSCARSSR::GetE (TVector3D const& X) const
{
  // Return summed E from container
  return this->fEFieldContainer.GetF(X);
}
























TParticleBeam& OSCARSSR::AddParticleBeam (std::string const& Type,
                                          std::string const& Name,
                                          TVector3D const& X0,
                                          TVector3D const& V0,
                                          double const Energy_GeV,
                                          double const T0,
                                          double const Current,
                                          double const Weight,
                                          double const Charge,
                                          double const Mass)
{
  // Add a particle beam
  // Type        - The name of the particle type that you want to use
  // Name        - A user specified 'name' for this beam
  // X0          - Initial position in X,Y,Z
  // V0          - A vector pointing in the direction of the velocity of arbitrary magnitude
  // Energy_GeV  - Energy of particle beam in GeV
  // T0          - Time of initial conditions, specified in units of [m] (for v = c)
  // Current     - Beam current in Amperes
  // Weight      - Relative weight to give this beam when randomly sampling
  // Charge      - Charge of custom particle
  // Mass        - Mass of custom particle

  return fParticleBeamContainer.AddNewParticleBeam(Type, Name, X0, V0, Energy_GeV, T0, Current, Weight, Charge, Mass);
}




TParticleBeam& OSCARSSR::AddParticleBeam (std::string const& Beam,
                                          std::string const& Name,
                                          double const Weight)
{
  // Add a particle beam
  // Beam - The name of the predefined particle beam to add

  return fParticleBeamContainer.AddNewParticleBeam(Beam, Name, Weight);
}




TParticleBeam& OSCARSSR::GetParticleBeam (std::string const& Name)
{
  // Return a reference to the particle beam by a given name
  return fParticleBeamContainer.GetParticleBeam(Name);
}




size_t OSCARSSR::GetNParticleBeams () const
{
  // Return the number of particle beams defined
  return fParticleBeamContainer.GetNParticleBeams();
}




TParticleA OSCARSSR::GetNewParticle ()
{
  // Get a new particle.  Randomly sampled according to input beam parameters and beam weights
  return fParticleBeamContainer.GetNewParticle();
}




TParticleA const& OSCARSSR::GetCurrentParticle () const
{
  // Get a new particle.  Randomly sampled according to input beam parameters and beam weights
  return fParticle;
}




void OSCARSSR::SetNewParticle ()
{
  // Get a new particle.  Randomly sampled according to input beam parameters and beam weights.
  // Set this new particle as *the* particle in OSCARSSR fParticle
  fParticle = fParticleBeamContainer.GetNewParticle();

  fParticle.ResetTrajectoryData();
  return;
}




void OSCARSSR::SetNewParticle (std::string const& BeamName, std::string const& IdealOrRandom)
{
  // Get a new particle.  Randomly sampled according to input beam parameters and beam weights.
  // Set this new particle as *the* particle in OSCARSSR fParticle

  if (BeamName == "") {
    fParticle = fParticleBeamContainer.GetRandomBeam().GetNewParticle(IdealOrRandom);
  } else {
    fParticle = fParticleBeamContainer.GetParticleBeam(BeamName).GetNewParticle(IdealOrRandom);
  }

  fParticle.ResetTrajectoryData();

  return;
}




void OSCARSSR::ClearParticleBeams ()
{
  // Clear the contents of the particle beam container
  fParticleBeamContainer.Clear();

  return;
}




void OSCARSSR::SetEmittance (std::string const& Beam,
                             TVector2D const& Emittance)
{
  fParticleBeamContainer.SetEmittance(Beam, Emittance);
  return;
}




void OSCARSSR::SetTwissParameters (std::string const& Beam,
                                   TVector2D const& Beta,
                                   TVector2D const& Alpha,
                                   TVector2D const& Gamma,
                                   TVector3D const& Lattice_Reference,
                                   bool const HasReferencePoint)
{
  fParticleBeamContainer.SetTwissParameters(Beam, Beta, Alpha, Gamma, Lattice_Reference, HasReferencePoint);
  return;
}




void OSCARSSR::AddDriftVolume (TDriftVolume* DV)
{
  // Add a magnetic field from a file to the field container

  this->fDriftVolumeContainer.AddDriftVolume(DV);

  return;
}



void OSCARSSR::RemoveDriftVolume (std::string const& Name)
{
  // Remove all drift volumes with a given name

  fDriftVolumeContainer.RemoveDriftVolume(Name);

  return;
}




void OSCARSSR::ClearDriftVolumes ()
{
  // Add a magnetic field from a file to the field container

  this->fDriftVolumeContainer.Clear();

  return;
}




void OSCARSSR::SetNPointsTrajectory (size_t const N)
{
  // Set this number of points for any trajectory calculations
  fNPointsTrajectory = N;
  return;
}




void OSCARSSR::SetNPointsPerMeterTrajectory (size_t const N)
{
  // Set this number of points for any trajectory calculations
  fNPointsPerMeter = N;

  // Set NPoints trajectory if CTStart and CTStop are set
  if (fCTStart != fCTStop) {
    this->SetCTStartStop(fCTStart, fCTStop);
  }


  return;
}




void OSCARSSR::SetCTStartStop (double const Start, double const Stop)
{
  // Set the start and stop time in units of m (where v = c)

  // Use number of points per meter to scare the trajectory
  fNPointsTrajectory = fNPointsPerMeter * (Stop - Start);

  fCTStart = Start;
  fCTStop  = Stop;
  return;
}




size_t OSCARSSR::GetNPointsTrajectory () const
{
  // Return the number of points being used for trajectory calculations
  return fNPointsTrajectory;
}




double OSCARSSR::GetCTStart () const
{
  // Return the start time in units of m (where v = c)
  return fCTStart;
}




double OSCARSSR::GetCTStop () const
{
  // Return the stop time in units of m (where v = c)
  return fCTStop;
}




bool OSCARSSR::SetUseGPUGlobal (int const in)
{
  // Must be a 0 or a 1 at the moment.  Will return false if you tried to set the GPU but it's not available

  if (in == 0) {
    fUseGPUGlobal = 0;
    return true;
  }

  if (in != 1) {
    fUseGPUGlobal = 0;
    return false;
  }

  #ifdef CUDA
  if (this->CheckGPU() > 0) {
    fUseGPUGlobal = 1;
    return true;
  } else {
    fUseGPUGlobal = 0;
    return false;
  }
  #endif

  fUseGPUGlobal = 0;

  return false;
}




int OSCARSSR::GetUseGPUGlobal () const
{
  return fUseGPUGlobal;
}




int OSCARSSR::CheckGPU () const
{
  #ifdef CUDA
    static int const N = OSCARSSR_Cuda_GetDeviceCount();
    return N;
  #endif
  return -1;
}




std::string OSCARSSR::GetGPUInfo(int const i) const
{
  #ifdef CUDA
    return OSCARSSR_Cuda_GetDeviceProperties(i);
  #endif
  return "INFO: GPU not enabled in this version";
}




void OSCARSSR::SetNThreadsGlobal (int const N)
{
  fNThreadsGlobal = N;
  return;
}




int OSCARSSR::GetNThreadsGlobal () const
{
  return fNThreadsGlobal;
}




void OSCARSSR::SetSeed (int const Seed) const
{
  gRandomA->SetSeed(Seed);
  return;
}




double OSCARSSR::GetRandomNormal () const
{
  return gRandomA->Normal();
}




double OSCARSSR::GetRandomUniform () const
{
  return gRandomA->Uniform();
}




void OSCARSSR::CorrectTrajectory ()
{
  // Correct the ideal trajectory so that the position and direction at
  // a specified point is as close as possible to the specified values

  TVector3D Direction(0, 0, 1);
  TVector3D ExitPoint(0, 0, 2);

  TVector3D V0 = Direction.Perp();
  TVector3D V1 = Direction.Cross(V0);

  TParticleTrajectoryPoints& ParticleTrajectory = fParticle.GetTrajectory();

  // Closest approach to point and direction
  double MinimumDistance2;
  size_t MinimumDistanceIndex = 0;
  for (size_t i = 0; i != ParticleTrajectory.GetNPoints(); ++i) {
    double const Distance2 = (ParticleTrajectory.GetX(i) - ExitPoint).Mag2();

    if (i == 0) {
      MinimumDistance2 = Distance2;
    }

    if (Distance2 < MinimumDistance2) {
      MinimumDistance2 = Distance2;
      MinimumDistanceIndex = i;
    }
  }
  std::cout << "MinimumDistance " << sqrt(MinimumDistance2) << std::endl;

  this->AddMagneticField(new TField3D_Gaussian(TVector3D(0, 0.0001, 0),
                                               TVector3D(0, 0, -1.8),
                                               TVector3D(0, 0, 0.010),
                                               TVector3D(0, 0, 0),
                                               "_autocorrector_entry"));

  this->AddMagneticField(new TField3D_Gaussian(TVector3D(0, 0.0001, 0),
                                               TVector3D(0, 0, +1.8),
                                               TVector3D(0, 0, 0.010),
                                               TVector3D(0, 0, 0),
                                               "_autocorrector_exit"));

  return;
}


void OSCARSSR::CalculateTrajectory ()
{
  // Function to calculate the particle trajectory of the member particle fParticle

  // Check that particle has been set yet.  If fType is "" it has not been set yet
  if (fParticle.GetType() == "") {
    this->SetNewParticle();
  }

  this->CalculateTrajectory(fParticle);

  return;
}





void OSCARSSR::CalculateTrajectory (TParticleA& P)
{
  // Function to calculate the particle trajectory given initial conditions.
  // This function uses the internal Trajectory member to store results

  // Check that CTStart is not after T0 of particle
  if (this->GetCTStart() > P.GetT0()) {
    std::cerr << "GetCTStart() P.GetT0(): " << this->GetCTStart() << " " << P.GetT0() << std::endl;
    throw std::out_of_range("start time is greater than T0");
  }

  // Check that CTStart and CTStop are not the same (probably not defined if this is the case)
  if (this->GetCTStart() >= this->GetCTStop()) {
    throw std::out_of_range("start time is greater than stop time.  check that set_ctstartstop is set");
  }

  // Check that particle has been set yet.  If fType is "" it has not been set yet
  if (P.GetType() == "") {
    throw std::out_of_range("particle not initialized.  make sure you have a particle or beam defined");
  }

  // Clear any current trajectory
  P.ResetTrajectoryData();

  // Calculate the total DeltaT in seconds
  double const DeltaT = ((this->GetCTStop() - this->GetCTStart()) / TOSCARSSR::C() / (fNPointsTrajectory - 1));


  // The number of points in the forward and backward direction
  // UPDATE: Check the details of these numbers
  size_t NPointsForward  = 1 + (this->GetCTStop() - P.GetT0()) / TOSCARSSR::C() / DeltaT;
  size_t NPointsBackward = (P.GetT0() - this->GetCTStart()) / TOSCARSSR::C() / DeltaT;


  // Number of dimensions of the array to be sent for RK calculation
  int const N = 6;

  // Arrays to be sent for RK calculation
  double x[N];
  double dxdt[N];

  // Initial conditions for the forward propogation
  x[0] = P.GetX0().GetX();
  x[1] = P.GetB0().GetX() * TOSCARSSR::C();
  x[2] = P.GetX0().GetY();
  x[3] = P.GetB0().GetY() * TOSCARSSR::C();
  x[4] = P.GetX0().GetZ();
  x[5] = P.GetB0().GetZ() * TOSCARSSR::C();


  // Grap the particle trajectory object
  TParticleTrajectoryPoints& ParticleTrajectory = P.GetTrajectory();
  ParticleTrajectory.Reserve(NPointsForward + NPointsBackward);

  // Set delta T for the trajectory
  ParticleTrajectory.SetDeltaT(DeltaT);

  // Loop over points in the forward direction
  for (size_t i = 0; i != NPointsForward; ++i) {

    // This time
    double t = P.GetT0() / TOSCARSSR::C() + DeltaT * i;


    // UPDATE: dhidas dhidas dhidas
    if (fDriftVolumeContainer.IsInside(TVector3D(x[0], x[2], x[4]))) {
      x[0] += DeltaT * x[1];
      x[2] += DeltaT * x[3];
      x[4] += DeltaT * x[5];
    } else {
      // Add this point to the trajectory
      ParticleTrajectory.AddPoint(x[0], x[2], x[4], x[1] / TOSCARSSR::C(), x[3] / TOSCARSSR::C(), x[5] / TOSCARSSR::C(), dxdt[1] / TOSCARSSR::C(), dxdt[3] / TOSCARSSR::C(), dxdt[5] / TOSCARSSR::C(), t);

      // Propogate
      (this->*fDerivativesFunction)(t, x, dxdt, P);
      RK4(x, dxdt, N, t, DeltaT, x, P);
    }
  }

  // Reverse trajectory elements for backward propogation
  ParticleTrajectory.ReverseArrays();

  // Set initial conditions for propogating backwards
  x[0] =  P.GetX0().GetX();
  x[1] =  P.GetB0().GetX() * TOSCARSSR::C();
  x[2] =  P.GetX0().GetY();
  x[3] =  P.GetB0().GetY() * TOSCARSSR::C();
  x[4] =  P.GetX0().GetZ();
  x[5] =  P.GetB0().GetZ() * TOSCARSSR::C();

  // Reverse time
  double const DeltaTReversed = -DeltaT;

  // Loop over all points "before" the initial point
  for (size_t i = 0; i != NPointsBackward; ++i) {

    // This time
    double t = P.GetT0() / TOSCARSSR::C() + DeltaTReversed * (i + 1);

    if (fDriftVolumeContainer.IsInside(TVector3D(x[0], x[2], x[4]))) {
      x[0] += DeltaTReversed * x[1];
      x[2] += DeltaTReversed * x[3];
      x[4] += DeltaTReversed * x[5];
    } else {
      // Propogate backward in time!
      (this->*fDerivativesFunction)(t, x, dxdt, P);
      RK4(x, dxdt, N, t, DeltaTReversed, x, P);

      // Add the point to the trajectory
      ParticleTrajectory.AddPoint(x[0], x[2], x[4], x[1] / TOSCARSSR::C(), x[3] / TOSCARSSR::C(), x[5] / TOSCARSSR::C(), dxdt[1] / TOSCARSSR::C(), dxdt[3] / TOSCARSSR::C(), dxdt[5] / TOSCARSSR::C(), t);
    }
  }

  // Re-Reverse the trajectory to be in the proper time order
  ParticleTrajectory.ReverseArrays();

  P.SetupTrajectoryInterpolated();

  return;
}




TParticleTrajectoryPoints const& OSCARSSR::GetTrajectory ()
{
  // Get the trajectory for *the* current particle in fParticle

  return fParticle.GetTrajectory();
}




void OSCARSSR::ClearTrajectory ()
{
  fParticle.ResetTrajectoryData();
  return;
}





void OSCARSSR::SetDerivativesFunction ()
{
  // Set the derivatives function for RK4 depending on what fields exist
  // I am using a function pointer so that we avoid computations that are not needed (save CPU cycles)

  if (fBFieldContainer.GetNFields() == 0 && fEFieldContainer.GetNFields() > 0) {
    fDerivativesFunction = &OSCARSSR::DerivativesE;
  } else if (fBFieldContainer.GetNFields() > 0 && fEFieldContainer.GetNFields() == 0) {
    fDerivativesFunction = &OSCARSSR::DerivativesB;
  } else {
    fDerivativesFunction = &OSCARSSR::DerivativesEB;
  }

  return;
}









void OSCARSSR::DerivativesE (double t, double x[], double dxdt[], TParticleA const& P)
{
  // This is a second order differential equation.  It does not account for the loss in energy due to
  // radiation.  Although 't' is not used it would be easy to implement a time dependent field

  // The values correspond to:
  // x[0] - x
  // x[1] - Vx
  // x[2] - y
  // x[3] - Vy
  // x[4] - z
  // x[5] - Vz

  // BField at this point
  TVector3D const E = this->GetE(x[0], x[2], x[4]);

  double const QoverMTimesSqrtOneMinusBetaSquared = P.GetQ() / P.GetM() * sqrt(1. - (x[1]*x[1] + x[3]*x[3] + x[5]*x[5]) / (TOSCARSSR::C() * TOSCARSSR::C()));
  double const BetaDotE = (x[1] * E.GetX() + x[3] * E.GetY() + x[5] * E.GetZ()) / TOSCARSSR::C();

  dxdt[0] = x[1];
  dxdt[1] = QoverMTimesSqrtOneMinusBetaSquared * (E.GetX() - x[1] * BetaDotE / TOSCARSSR::C());
  dxdt[2] = x[3];                                          
  dxdt[3] = QoverMTimesSqrtOneMinusBetaSquared * (E.GetY() - x[3] * BetaDotE / TOSCARSSR::C());
  dxdt[4] = x[5];                                          
  dxdt[5] = QoverMTimesSqrtOneMinusBetaSquared * (E.GetZ() - x[5] * BetaDotE / TOSCARSSR::C());

  return;
}



void OSCARSSR::DerivativesB (double t, double x[], double dxdt[], TParticleA const& P)
{
  // This is a second order differential equation.  It does not account for the loss in energy due to
  // radiation.  Although 't' is not used it would be easy to implement a time dependent field

  // The values correspond to:
  // x[0] - x
  // x[1] - Vx
  // x[2] - y
  // x[3] - Vy
  // x[4] - z
  // x[5] - Vz

  // BField at this point
  TVector3D const B = this->GetB(x[0], x[2], x[4]);

  double const QoverMGamma = P.GetQoverMGamma();

  dxdt[0] = x[1];
  dxdt[1] = QoverMGamma * (-x[5] * B.GetY() + x[3] * B.GetZ());
  dxdt[2] = x[3];                                                                                          
  dxdt[3] = QoverMGamma * ( x[5] * B.GetX() - x[1] * B.GetZ());
  dxdt[4] = x[5];                                                                                          
  dxdt[5] = QoverMGamma * ( x[1] * B.GetY() - x[3] * B.GetX());

  return;
}



void OSCARSSR::DerivativesEB (double t, double x[], double dxdt[], TParticleA const& P)
{
  // This is a second order differential equation.  It does not account for the loss in energy due to
  // radiation.  Although 't' is not used it would be easy to implement a time dependent field

  // The values correspond to:
  // x[0] - x
  // x[1] - Vx
  // x[2] - y
  // x[3] - Vy
  // x[4] - z
  // x[5] - Vz

  // BField at this point
  TVector3D const B = this->GetB(x[0], x[2], x[4]);
  TVector3D const E = this->GetE(x[0], x[2], x[4]);

  double const QoverMTimesSqrtOneMinusBetaSquared = P.GetQ() / P.GetM() * sqrt(1. - (x[1]*x[1] + x[3]*x[3] + x[5]*x[5]) / (TOSCARSSR::C() * TOSCARSSR::C()));
  double const BetaDotE = (x[1] * E.GetX() + x[3] * E.GetY() + x[5] * E.GetZ()) / TOSCARSSR::C();

  dxdt[0] = x[1];
  dxdt[1] = QoverMTimesSqrtOneMinusBetaSquared * (E.GetX() - x[5] * B.GetY() + x[3] * B.GetZ() - x[1] * BetaDotE / TOSCARSSR::C());
  dxdt[2] = x[3];                                                                          
  dxdt[3] = QoverMTimesSqrtOneMinusBetaSquared * (E.GetY() + x[5] * B.GetX() - x[1] * B.GetZ() - x[3] * BetaDotE / TOSCARSSR::C());
  dxdt[4] = x[5];                                                                          
  dxdt[5] = QoverMTimesSqrtOneMinusBetaSquared * (E.GetZ() + x[1] * B.GetY() - x[3] * B.GetX() - x[5] * BetaDotE / TOSCARSSR::C());


  return;
}



void OSCARSSR::Derivatives (double t, double x[], double dxdt[], TParticleA const& P)
{
  // This is a second order differential equation.  It does not account for the loss in energy due to
  // radiation.  Although 't' is not used it would be easy to implement a time dependent field

  // The values correspond to:
  // x[0] - x
  // x[1] - Vx
  // x[2] - y
  // x[3] - Vy
  // x[4] - z
  // x[5] - Vz

  // BField at this point
  TVector3D const B = this->GetB(x[0], x[2], x[4]);

  double const QoverMGamma = P.GetQoverMGamma();

  dxdt[0] = x[1];
  dxdt[1] = QoverMGamma * (-x[5] * B.GetY() + x[3] * B.GetZ());
  dxdt[2] = x[3];                                                                                          
  dxdt[3] = QoverMGamma * ( x[5] * B.GetX() - x[1] * B.GetZ());
  dxdt[4] = x[5];                                                                                          
  dxdt[5] = QoverMGamma * ( x[1] * B.GetY() - x[3] * B.GetX());

  return;
}



void OSCARSSR::RK4 (double y[], double dydx[], int n, double x, double h, double yout[], TParticleA const& P)
{
  // Runge-Kutta 4th order method propogation
  // UPDATE: syntax

  int i;
  double xh, hh, h6;
  double *dym, *dyt, *yt;

  dym = new double[n];
  dyt = new double[n];
  yt  = new double[n];

  hh = h * 0.5;
  h6 = h / 6.0;
  xh = x + hh;

  for (i = 0; i < n; ++i) {
    yt[i] = y[i] + hh * dydx[i];
  }

  (this->*fDerivativesFunction)(xh, yt, dyt, P);
  //this->DerivativesB(xh, yt, dyt, P);

  for (i = 0; i < n; ++i) {
    yt[i] = y[i] + hh * dyt[i];
  }

  (this->*fDerivativesFunction)(xh, yt, dym, P);
  //this->DerivativesB(xh, yt, dym, P);

  for (i = 0; i < n; ++i) {
    yt[i] = y[i] + h * dym[i];
    dym[i] += dyt[i];
  }

  (this->*fDerivativesFunction)(x + h, yt, dyt, P);
  //this->DerivativesB(x + h, yt, dyt, P);

  for (i = 0; i < n; ++i) {
    yout[i] = y[i] + h6 * (dydx[i] + dyt[i] + 2.0 * dym[i]);
  }
  
  delete [] dym;
  delete [] dyt;
  delete [] yt;

  return;
}





void OSCARSSR::CalculateSpectrum (TParticleA& Particle,
                                  TVector3D const& ObservationPoint,
                                  TSpectrumContainer& Spectrum,
                                  std::string const& Polarization,
                                  double const Angle,
                                  TVector3D const& HorizontalDirection,
                                  TVector3D const& PropogationDirection,
                                  double const Precision,
                                  int    const MaxLevel,
                                  int    const MaxLevelExtended,
                                  double const Weight,
                                  int    const ReturnQuantity)
{
  // Calculates the single particle spectrum at a given observation point
  // in units of [photons / second / 0.001% BW / mm^2]
  // Save this in the spectrum container.
  //
  // Particle - the Particle.. with a Trajectory structure hopefully
  // ObservationPoint - Observation Point
  // Spectrum - Spectrum container

  // Calculate trajectory if it doesn't exist
  if (Particle.GetTrajectory().GetNPoints() == 0) {
    this->CalculateTrajectory(Particle);
  }


  // Extra inpts for calculation
  bool Done = false;

  this->CalculateSpectrumPoints(Particle,
                                ObservationPoint,
                                Spectrum,
                                0,
                                1,
                                Done,
                                Polarization,
                                Angle,
                                HorizontalDirection,
                                PropogationDirection,
                                Precision,
                                MaxLevel,
                                MaxLevelExtended,
                                Weight,
                                ReturnQuantity);

  return;
}





void OSCARSSR::CalculateSpectrum (TVector3D const& ObservationPoint,
                                  TSpectrumContainer& Spectrum,
                                  std::string const& Polarization,
                                  double const Angle,
                                  TVector3D const& HorizontalDirection,
                                  TVector3D const& PropogationDirection,
                                  int const NParticles,
                                  int const NThreads,
                                  int const GPU,
                                  int const NGPU,
                                  std::vector<int> VGPU,
                                  double const Precision,
                                  int    const MaxLevel,
                                  int    const MaxLevelExtended,
                                  int    const ReturnQuantity)
{
  // Calculate the spectrum at an observaton point.
  // THIS is the ENTRY POINT typically

  // Check that particle has been set yet.  If fType is "" it has not been set yet
  if (fParticle.GetType() == "") {
    try {
      this->SetNewParticle();
    } catch (std::exception e) {
      throw std::out_of_range("no beam defined");
    }
  }

  // Number of threads to possibly use
  int const NThreadsToUse = NThreads < 1 ? fNThreadsGlobal : NThreads;
  if (NThreadsToUse <= 0) {
    throw std::out_of_range("NThreads or NThreadsGlobal must be >= 1");
  }

  // Should we use the GPU or not?
  int const NGPUAvailable = this->CheckGPU();
  bool const UseGPU = GPU == 0 ? false : this->GetUseGPUGlobal() && (NGPUAvailable > 0) ? true : (GPU == 1) && (NGPUAvailable > 0);

  // If use GPU, let's check the GPU vector, or construct one if we must
  std::vector<int> GPUVector;
  for (std::vector<int>::const_iterator it = VGPU.begin(); it != VGPU.end(); ++it) {
    GPUVector.push_back(*it);
  }
  if (GPUVector.size() == 0) {
    for (int i = 0; i < NGPU; ++i) {
      GPUVector.push_back(i);
    }
  }
  if (NGPU != -1 && NGPU < (int) GPUVector.size()) {
    GPUVector.resize(NGPU);
  }


  // Check polarization
  if (Polarization == "all" ||
      Polarization == "linear-horizontal" ||
      Polarization == "linear-vertical"   ||
      Polarization == "linear"            ||
      Polarization == "circular-left"     ||
      Polarization == "circular-right") {
    // Do nothing
  } else {
    throw std::invalid_argument("Polarization requested not recognized");
  }


  // Which cpmpute method will we use, gpu, multi-thread, or single-thread
  if (UseGPU) {
    // Send to GPU function
    this->CalculateSpectrumGPU(fParticle,
                               ObservationPoint,
                               Spectrum,
                               Polarization,
                               Angle,
                               HorizontalDirection,
                               PropogationDirection,
                               NParticles,
                               GPUVector,
                               Precision,
                               MaxLevel,
                               MaxLevelExtended,
                               ReturnQuantity);
  } else {
    if (NParticles == 0) {
      if (NThreadsToUse == 1) {
        this->CalculateSpectrum(fParticle,
                                ObservationPoint,
                                Spectrum,
                                Polarization,
                                Angle,
                                HorizontalDirection,
                                PropogationDirection,
                                Precision,
                                MaxLevel,
                                MaxLevelExtended,
                                1,
                                ReturnQuantity);
      } else {
        this->CalculateSpectrumThreads(fParticle,
                                       ObservationPoint,
                                       Spectrum,
                                       NThreadsToUse,
                                       Polarization,
                                       Angle,
                                       HorizontalDirection,
                                       PropogationDirection,
                                       Precision,
                                       MaxLevel,
                                       MaxLevelExtended,
                                       1,
                                       ReturnQuantity);
      }
    } else {
      // Weight this by the number of particles
      double const Weight = 1.0 / (double) NParticles;

      // Loop over particles
      for (int i = 0; i != NParticles; ++i) {

        // Set a new random particle
        this->SetNewParticle();
        this->CalculateTrajectory();

        if (NThreadsToUse == 1) {
          this->CalculateSpectrum(fParticle,
                                  ObservationPoint,
                                  Spectrum,
                                  Polarization,
                                  Angle,
                                  HorizontalDirection,
                                  PropogationDirection,
                                  Precision,
                                  MaxLevel,
                                  MaxLevelExtended,
                                  Weight,
                                  ReturnQuantity);
        } else {
          this->CalculateSpectrumThreads(fParticle,
                                         ObservationPoint,
                                         Spectrum,
                                         NThreadsToUse,
                                         Polarization,
                                         Angle,
                                         HorizontalDirection,
                                         PropogationDirection,
                                         Precision,
                                         MaxLevel,
                                         MaxLevelExtended,
                                         Weight,
                                         ReturnQuantity);
        }
      }
    }
  }


  return;
}





void OSCARSSR::CalculateSpectrumPoints (TParticleA& Particle,
                                        TVector3D const& ObservationPoint,
                                        TSpectrumContainer& Spectrum,
                                        size_t const iThread,
                                        size_t const NThreads,
                                        bool& Done,
                                        std::string const& Polarization,
                                        double const Angle,
                                        TVector3D const& HorizontalDirection,
                                        TVector3D const& PropogationDirection,
                                        double const Precision,
                                        int    const MaxLevel,
                                        int    const MaxLevelExtended,
                                        double const Weight,
                                        int    const ReturnQuantity)
{
  // Calculates the single particle spectrum at a given observation point
  // in units of [photons / second / 0.001% BW / mm^2]
  // Save this in the spectrum container.
  //
  // Particle - the Particle.. with a Trajectory structure hopefully
  // ObservationPoint - Observation Point
  // Spectrum - Spectrum container

  // Check that particle has been set yet.  If fType is "" it has not been set yet
  if (Particle.GetType() == "") {
    throw std::out_of_range("no particle defined");
  }

  // Check you are not requesting a level above the maximum
  if (MaxLevel > TParticleA::kMaxTrajectoryLevel) {
    std::cerr << "WARNING: MaxLevel > TParticleA::kMaxTrajectoryLevel.  Setting MaxLevel to TParticleA::kMaxTrajectoryLevel" << std::endl;
  }

  // Set the level to stop at if requested, but not above the hard limit
  int const LevelStopMemory = MaxLevel >= -1  && MaxLevel <= TParticleA::kMaxTrajectoryLevel ? MaxLevel : TParticleA::kMaxTrajectoryLevel;
  int const LevelStopWithExtended = MaxLevelExtended > LevelStopMemory ? MaxLevelExtended : LevelStopMemory;

  // Number of points in spectrum
  int const NSpectrumPoints = Spectrum.GetNPoints();

  // Check input spectrum range numbers
  if (iThread >= NThreads) {
    throw std::out_of_range("spectrum thread range is incorrect.  Please report this error.");
  }

  // Constant C0 for calculation
  double const C0 = Particle.GetQ() / (TOSCARSSR::FourPi() * TOSCARSSR::C() * TOSCARSSR::Epsilon0() * TOSCARSSR::Sqrt2Pi());

  // Constant for flux calculation at the end
  double const C2 = TOSCARSSR::FourPi() * Particle.GetCurrent() / (TOSCARSSR::H() * fabs(Particle.GetQ()) * TOSCARSSR::Mu0() * TOSCARSSR::C()) * 1e-6 * 0.001;

  // Imaginary "i" and complxe 1+0i
  std::complex<double> const I(0, 1);
  std::complex<double> const One(1, 0);

  // Photon vertical direction and positive and negative helicity
  TVector3D  const VerticalDirection = PropogationDirection.Cross(HorizontalDirection).UnitVector();
  TVector3DC const Positive = 1. / sqrt(2) * (TVector3DC(HorizontalDirection) + VerticalDirection * I );
  TVector3DC const Negative = 1. / sqrt(2) * (TVector3DC(HorizontalDirection) - VerticalDirection * I );

  TVector3DC PolarizationVector(0, 0, 0);
  if (Polarization == "all") {
    // Do Nothing
  } else if (Polarization == "linear-horizontal") {
    PolarizationVector = HorizontalDirection;
  } else if (Polarization == "linear-vertical") {
    PolarizationVector = VerticalDirection;
  } else if (Polarization == "linear") {
    TVector3D PolarizationAngle = HorizontalDirection;
    PolarizationAngle.RotateSelf(Angle, PropogationDirection);
    PolarizationVector = PolarizationAngle;
  } else if (Polarization == "circular-left") {
    PolarizationVector = Positive;
  } else if (Polarization == "circular-right") {
    PolarizationVector = Negative;
  } else {
    // Throw invalid argument if polarization is not recognized
    throw std::invalid_argument("Polarization requested not recognized");
  }


  // Extended trajectory (not using memory for storage of arrays
  TParticleTrajectoryInterpolatedPoints TE;

  // Alternative outputs
  double Result_Precision = -1;
  int    Result_Level     = -1;


  // Loop over all points in the spectrum container
  for (size_t i = iThread; i < NSpectrumPoints; i += NThreads) {

    // Angular frequency
    double const Omega = Spectrum.GetAngularFrequency(i);

    // Electric field summation in frequency space
    TVector3DC SumE(0, 0, 0);

    double ThisMag = -1;
    double LastMag = -1;
    double ThisPhase = -1;
    double LastPhase = -1;
    double MaxDPhase = 0;
    int    LastLevel = 0;

    for (int iLevel = 0; iLevel <= LevelStopWithExtended; ++iLevel) {
      LastLevel = iLevel;

      // Grab the Trajectory (using memory arrays) if below level threshold, else set NULL
      TParticleTrajectoryPoints const& TM = Particle.GetTrajectoryLevel(iLevel <= LevelStopMemory ? iLevel : 0);
      if (iLevel > LevelStopMemory) {
        TE = Particle.GetTrajectoryExtendedLevel(iLevel);
      }


      // Number of points in the trajectory
      size_t const NTPoints = iLevel <= LevelStopMemory ? TM.GetNPoints() : TE.GetNPoints();


      MaxDPhase = 0;
      // Loop over trajectory points
      for (int iT = 0; iT != NTPoints; ++iT) {
        TParticleTrajectoryPoint const& PP = (iLevel <= LevelStopMemory ? TM.GetPoint(iT) : TE.GetTrajectoryPoint(iT));

        // Get position, Beta, and Acceleration (over c)
        TVector3D const& X = PP.GetX();
        TVector3D const& B = PP.GetB();
        TVector3D const& AoverC = PP.GetAoverC();
        double    const  Time = iLevel <= LevelStopMemory ? TM.GetT(iT) : TE.GetT(iT);

        // Define R and unit vector in direction of R, and D (distance to observer)
        TVector3D const R = ObservationPoint - X;
        TVector3D const N = R.UnitVector();
        double const D = R.Mag();

        // Exponent in transformed field
        ThisPhase = -Omega * (Time + D / TOSCARSSR::C());
        double const PhaseTestValue = fabs(ThisPhase - LastPhase);
        if (iT != 0 && PhaseTestValue > MaxDPhase) {
          MaxDPhase = PhaseTestValue;
        }
        LastPhase = ThisPhase;
        std::complex<double> Exponent(0, ThisPhase);

        TVector3DC const ThisEw = ( ( (1 - (B).Mag2()) * (N - B) ) / ( D * D * (pow(1 - N.Dot(B), 2)) )
            + ( N.Cross( (N - B).Cross(AoverC) ) ) / ( D * pow(1 - N.Dot(B), 2) ) ) * std::exp(Exponent); // NF + FF

        // Add this contribution
        SumE += ThisEw;

      }

      TVector3DC ThisSumE = SumE * Particle.GetTrajectoryInterpolated().GetDeltaTInclusiveToLevel(iLevel);
      if (PolarizationVector.Mag2() > 0.001) {
        ThisSumE = ThisSumE.Dot(PolarizationVector) * PolarizationVector;
      }
      ThisMag = ThisSumE.Dot( ThisSumE.CC() ).real();

      Result_Precision = fabs(ThisMag - LastMag) / LastMag;
      if (iLevel > 8 && Result_Precision < Precision && MaxDPhase < TOSCARSSR::Pi()) {
        Result_Level = iLevel;
        break;
      }

      LastMag = ThisMag;
    }

    if (Result_Level == -1) {
      Spectrum.SetNotConverged(i);
    }

    // Multiply by constant factor
    SumE *= C0 * Particle.GetTrajectoryInterpolated().GetDeltaTInclusiveToLevel(LastLevel);

    // Correcr for polarization
    if (PolarizationVector.Mag2() > 0.001) {
      SumE = SumE.Dot(PolarizationVector) * PolarizationVector;
    }

    // Set the flux for this frequency / energy point
    // Add to container
    switch (ReturnQuantity) {
      case 1:
        Spectrum.AddToFlux(i, Result_Precision * Weight);
        break;
      case 2:
        Spectrum.AddToFlux(i, ((double) Result_Level) * Weight);
        break;
      default:
        Spectrum.AddToFlux(i, C2 *  SumE.Dot( SumE.CC() ).real() * Weight);
        break;
    }
  }


  // Set done to true
  Done = true;

  return;
}





void OSCARSSR::CalculateSpectrumThreads (TParticleA& Particle,
                                         TVector3D const& Obs,
                                         TSpectrumContainer& Spectrum,
                                         int const NThreads,
                                         std::string const& Polarization,
                                         double const Angle,
                                         TVector3D const& HorizontalDirection,
                                         TVector3D const& PropogationDirection,
                                         double const Precision,
                                         int    const MaxLevel,
                                         int    const MaxLevelExtended,
                                         double const Weight,
                                         int    const ReturnQuantity)
{
  // Calculates spectrum for the given particle and observation point
  // in units of [photons / second / 0.001% BW / mm^2]
  //
  // Surface - Observation Point

  // Calculate trajectory if it doesn't exist
  if (Particle.GetTrajectory().GetNPoints() == 0) {
    this->CalculateTrajectory(Particle);
  }

  // Vector container for threads
  std::vector<std::thread> Threads;

  // Number of points in spectrum
  size_t const NPoints = Spectrum.GetNPoints();

  // How many threads to start in the first for loop
  size_t const NThreadsActual = NPoints > (size_t) NThreads ? NThreads : NPoints;

  // Keep track of which threads are finished and re-joined
  bool *Done = new bool[NThreadsActual];
  bool *Joined = new bool[NThreadsActual];

  // Number per thread plus remainder to be added to first threads
  size_t const NPerThread = NPoints / NThreadsActual;
  size_t const NRemainder = NPoints % NThreadsActual;

  // Start threads and keep in vector
  for (size_t it = 0; it != NThreadsActual; ++it) {

    // First and last points for each thread
    size_t const iFirst = it < NRemainder ? NPerThread * it + it: NPerThread * it + NRemainder;
    size_t const iLast  = it < NRemainder ? iFirst + NPerThread : iFirst + NPerThread - 1;

    // Set Done and joined to false for this thread
    Done[it] = false;
    Joined[it] = false;

    // Start thread for these points
    Threads.push_back(std::thread(&OSCARSSR::CalculateSpectrumPoints,
                                  this,
                                  std::ref(Particle),
                                  std::ref(Obs),
                                  std::ref(Spectrum),
                                  it,
                                  NThreadsActual,
                                  std::ref(Done[it]),
                                  std::ref(Polarization),
                                  Angle,
                                  std::ref(HorizontalDirection),
                                  std::ref(PropogationDirection),
                                  Precision,
                                  MaxLevel,
                                  MaxLevelExtended,
                                  Weight,
                                  ReturnQuantity));

  }

  // Are all of the threads finished or not?  Continue loop until all come back.
  bool AllThreadsFinished = false;
  size_t NThreadsFinished = 0;
  while (!AllThreadsFinished) {


    // So as to not use the current thread at 100%
    std::this_thread::sleep_for(std::chrono::milliseconds(10));

    // Check all threads
    for (size_t it = 0; it != NThreadsActual; ++it) {

      if (Done[it] && !Joined[it]) {
        Threads[it].join();
        Joined[it] = true;
        ++NThreadsFinished;
      }
    }

    // If the number finished is equal to the number of points total then we're done
    if (NThreadsFinished == NThreadsActual) {
      AllThreadsFinished = true;
    }
  }


  // Clear all threads
  Threads.clear();

  // Delete my arrays, I hate new for this purpose
  delete [] Done;
  delete [] Joined;

  return;
}






void OSCARSSR::CalculateSpectrumGPU (TParticleA& Particle,
                                     TVector3D const& ObservationPoint,
                                     TSpectrumContainer& Spectrum,
                                     std::string const& Polarization,
                                     double const Angle,
                                     TVector3D const& HorizontalDirection,
                                     TVector3D const& PropogationDirection,
                                     int const NParticles,
                                     std::vector<int> GPUVector,
                                     double const Precision,
                                     int    const MaxLevel,
                                     int    const MaxLevelExtended,
                                     int    const ReturnQuantity)
{
  // If you compile for Cuda use the GPU in this function, else throw

  // If GPUVector is empty assume you want to use ALL GPUs available
  if (GPUVector.size() == 0) {
    int const NGPUAvailable = this->CheckGPU();
    for (int i = 0; i < NGPUAvailable; ++i) {
      GPUVector.push_back(i);
    }
  }

  #ifdef CUDA
  // Check that the GPU exists
  if (this->CheckGPU() < 1) {
    throw std::invalid_argument("You are requesting the GPU, but none were found");
  }

  // Set the level to stop at if requested, but not above the hard limit
  int const LevelStopMemory = MaxLevel >= -1  && MaxLevel <= TParticleA::kMaxTrajectoryLevel ? MaxLevel : TParticleA::kMaxTrajectoryLevel;
  int const LevelStopWithExtended = MaxLevelExtended > LevelStopMemory ? MaxLevelExtended : LevelStopMemory;


  return OSCARSSR_Cuda_CalculateSpectrumGPU (*this,
                                             ObservationPoint,
                                             Spectrum,
                                             Polarization,
                                             Angle,
                                             HorizontalDirection,
                                             PropogationDirection,
                                             NParticles,
                                             GPUVector,
                                             Precision,
                                             LevelStopWithExtended,
                                             ReturnQuantity);
  #else
  throw std::invalid_argument("GPU functionality not compiled into this binary distribution");
  #endif

  return;
}


















void OSCARSSR::AddToSpectrum (TSpectrumContainer const& S, double const Weight)
{
  // Check if spectrum exists yet or not.  In not, create it
  if (fSpectrum.GetNPoints() == 0) {
    for (size_t i = 0; i != S.GetNPoints(); ++i) {
      fSpectrum.AddPoint(S.GetEnergy(i), S.GetFlux(i) * Weight);
    }
  } else if (fSpectrum.GetNPoints() == S.GetNPoints()) {
    for (size_t i = 0; i != S.GetNPoints(); ++i) {
      fSpectrum.AddToFlux(i, S.GetFlux(i) * Weight);
    }
  } else {
    throw std::out_of_range("spectra dimensions do not match");
  }

  return;
}






void OSCARSSR::AddToFlux (T3DScalarContainer const& F, double const Weight)
{
  // Check if spectrum exists yet or not.  In not, create it
  if (fFlux.GetNPoints() == 0) {
    for (size_t i = 0; i != F.GetNPoints(); ++i) {
      fFlux.AddPoint(F.GetPoint(i).GetX(), F.GetPoint(i).GetV() * Weight);
    }
  } else if (fFlux.GetNPoints() == F.GetNPoints()) {
    for (size_t i = 0; i != F.GetNPoints(); ++i) {
      fFlux.AddToPoint(i, F.GetPoint(i).GetV() * Weight);
    }
  } else {
    throw std::out_of_range("spectra dimensions do not match");
  }

  return;
}






void OSCARSSR::AddToPowerDensity (T3DScalarContainer const& P, double const Weight)
{
  // Check if spectrum exists yet or not.  In not, create it
  if (fPowerDensity.GetNPoints() == 0) {
    for (size_t i = 0; i != P.GetNPoints(); ++i) {
      fPowerDensity.AddPoint(P.GetPoint(i).GetX(), P.GetPoint(i).GetV() * Weight);
    }
  } else if (fPowerDensity.GetNPoints() == P.GetNPoints()) {
    for (size_t i = 0; i != P.GetNPoints(); ++i) {
      fPowerDensity.AddToPoint(i, P.GetPoint(i).GetV() * Weight);
    }
  } else {
    throw std::out_of_range("spectra dimensions do not match");
  }

  return;
}






TSpectrumContainer const& OSCARSSR::GetSpectrum () const
{
  return fSpectrum;
}




void OSCARSSR::ClearSpectrum ()
{
  // Clear the contents of the particle beam container
  fSpectrum.Clear();

  return;
}




T3DScalarContainer const& OSCARSSR::GetFlux () const
{
  return fFlux;
}




void OSCARSSR::ClearFlux ()
{
  // Clear the contents
  fFlux.Clear();

  return;
}




T3DScalarContainer const& OSCARSSR::GetPowerDensity () const
{
  return fPowerDensity;
}




void OSCARSSR::ClearPowerDensity ()
{
  // Clear the contents
  fPowerDensity.Clear();

  return;
}












void OSCARSSR::CalculatePowerDensity (TParticleA& Particle,
                                      TSurfacePoints const& Surface,
                                      T3DScalarContainer& PowerDensityContainer,
                                      bool const Directional,
                                      double const Precision,
                                      int    const MaxLevel,
                                      int    const MaxLevelExtended,
                                      double const Weight,
                                      int    const ReturnQuantity)
{
  // Calculates the single particle spectrum at a given observation point
  // in units of [photons / second / 0.001% BW / mm^2]
  //
  // Particle - Particle, contains trajectory (or if not, calculate it)
  // Surface - Observation Point

  // Calculate trajectory if it doesn't exist
  if (Particle.GetTrajectory().GetNPoints() == 0) {
    this->CalculateTrajectory(Particle);
  }

  // Extra inpts for calculation
  bool Done = false;
  size_t const iFirst = 0;
  size_t const iLast = Surface.GetNPoints() - 1;

  // Calculate the power density
  CalculatePowerDensityPoints(Particle,
                              Surface,
                              PowerDensityContainer,
                              iFirst,
                              iLast,
                              Done,
                              Directional,
                              Precision,
                              MaxLevel,
                              MaxLevelExtended,
                              Weight,
                              ReturnQuantity);

  return;
}




void OSCARSSR::CalculatePowerDensity (TSurfacePoints const& Surface,
                                      T3DScalarContainer& PowerDensityContainer,
                                      int const Dimension,
                                      bool const Directional,
                                      double const Precision,
                                      int    const MaxLevel,
                                      int    const MaxLevelExtended,
                                      int const NParticles,
                                      int const NThreads,
                                      int const GPU,
                                      int const NGPU,
                                      std::vector<int> VGPU,
                                      int const ReturnQuantity)
{
  // Calculates the power density in units of [W / mm^2]
  // THIS is the ENTRY POINT typically
  //

  // Check that particle has been set yet.  If fType is "" it has not been set yet
  if (fParticle.GetType() == "") {
    try {
      this->SetNewParticle();
    } catch (std::exception e) {
      throw std::out_of_range("no beam defined");
    }
  }

  // Number of threads to possibly use
  int const NThreadsToUse = NThreads < 1 ? fNThreadsGlobal : NThreads;
  if (NThreadsToUse <= 0) {
    throw std::out_of_range("NThreads or NThreadsGlobal must be >= 1");
  }

  // Should we use the GPU or not?
  int const NGPUAvailable = this->CheckGPU();
  bool const UseGPU = GPU == 0 ? false : this->GetUseGPUGlobal() && (NGPUAvailable > 0) ? true : (GPU == 1) && (NGPUAvailable > 0);

  // If use GPU, let's check the GPU vector, or construct one if we must
  std::vector<int> GPUVector;
  for (std::vector<int>::const_iterator it = VGPU.begin(); it != VGPU.end(); ++it) {
    GPUVector.push_back(*it);
  }
  if (GPUVector.size() == 0) {
    for (int i = 0; i < NGPU; ++i) {
      GPUVector.push_back(i);
    }
  }
  if (NGPU != -1 && NGPU < (int) GPUVector.size()) {
    GPUVector.resize(NGPU);
  }




  PowerDensityContainer.Clear();
  // Set the output flux container and with correct dimensions
  if (Dimension == 3) {
    for (size_t i = 0; i != Surface.GetNPoints(); ++i) {
      PowerDensityContainer.AddPoint(Surface.GetPoint(i).GetPoint(), 0);
    }
  } else if (Dimension == 2) {
    for (size_t i = 0; i != Surface.GetNPoints(); ++i) {
      PowerDensityContainer.AddPoint( TVector3D(Surface.GetX1(i), Surface.GetX2(i), 0), 0);
    }
  } else if (Dimension == 1) {
    // UPDATE: Needs something for 1D
    if (Surface.GetNPoints() == 1) {
      PowerDensityContainer.AddPoint( Surface.GetPoint(0).GetPoint(), 0 );
    } else {
      TVector3D const x1 = Surface.GetPoint(0).GetPoint();
      TVector3D const x2 = Surface.GetPoint( Surface.GetNPoints() - 1 ).GetPoint();

      TVector3D const C = x1 + (x2 - x1) / 2;
      for (size_t i = 0; i != Surface.GetNPoints(); ++i) {
        double const Sign = i < Surface.GetNPoints() / 2 ? -1. : 1.;

        PowerDensityContainer.AddPoint( TVector3D(Sign * (Surface.GetPoint(i).GetPoint() - C).Mag(), 0, 0), 0 );
      }
    }
  } else {
    throw std::out_of_range("Wrong dimension");
  }




  // Which cpmpute method will we use, gpu, multi-thread, or single-thread
  if (UseGPU) {
    // Send to GPU function
    CalculatePowerDensityGPU(Surface,
                             PowerDensityContainer,
                             NParticles,
                             GPUVector,
                             Directional,
                             Precision,
                             MaxLevel,
                             MaxLevelExtended,
                             ReturnQuantity);
  } else {
    if (NParticles == 0) {
      if (NThreadsToUse == 1) {
        this->CalculatePowerDensity(fParticle,
                                    Surface,
                                    PowerDensityContainer,
                                    Directional,
                                    Precision,
                                    MaxLevel,
                                    MaxLevelExtended,
                                    1,
                                    ReturnQuantity);
      } else {
        this->CalculatePowerDensityThreads(fParticle,
                                           Surface,
                                           PowerDensityContainer,
                                           NThreadsToUse,
                                           Directional,
                                           Precision,
                                           MaxLevel,
                                           MaxLevelExtended,
                                           1,
                                           ReturnQuantity);
      }
    } else {
      // Weight this by the number of particles
      double const Weight = 1.0 / (double) NParticles;

      // Loop over particles
      for (int i = 0; i != NParticles; ++i) {

        // Set a new random particle
        this->SetNewParticle();
        this->CalculateTrajectory();

        if (NThreadsToUse == 1) {
          this->CalculatePowerDensity(fParticle,
                                      Surface,
                                      PowerDensityContainer,
                                      Directional,
                                      Precision,
                                      MaxLevel,
                                      MaxLevelExtended,
                                      Weight,
                                      ReturnQuantity);
        } else {
          this->CalculatePowerDensityThreads(fParticle,
                                             Surface,
                                             PowerDensityContainer,
                                             NThreadsToUse,
                                             Directional,
                                             Precision,
                                             MaxLevel,
                                             MaxLevelExtended,
                                             Weight,
                                             ReturnQuantity);
        }
      }
    }
  }

  return;
}









void OSCARSSR::CalculatePowerDensityPoints (TParticleA& Particle,
                                            TSurfacePoints const& Surface,
                                            T3DScalarContainer& PowerDensityContainer,
                                            size_t const iFirst,
                                            size_t const iLast,
                                            bool& Done,
                                            bool const Directional,
                                            double const Precision,
                                            int    const MaxLevel,
                                            int    const MaxLevelExtended,
                                            double const Weight,
                                            int    const ReturnQuantity)
{
  // Calculates the single particle power density in a range of points
  // in units of [watts / second / mm^2]

  // Check you are not requesting a level above the maximum
  if (MaxLevel > TParticleA::kMaxTrajectoryLevel) {
    std::cerr << "WARNING: MaxLevel > TParticleA::kMaxTrajectoryLevel.  Setting MaxLevel to TParticleA::kMaxTrajectoryLevel" << std::endl;
  }

  // Set the level to stop at if requested, but not above the hard limit
  int const LevelStopMemory = MaxLevel >= -1  && MaxLevel <= TParticleA::kMaxTrajectoryLevel ? MaxLevel : TParticleA::kMaxTrajectoryLevel;
  int const LevelStopWithExtended = MaxLevelExtended > LevelStopMemory ? MaxLevelExtended : LevelStopMemory;

  // Variables for parts of the numerator and denominator of power density equation
  TVector3D Numerator;
  double Denominator;

  bool const HasNormal = Surface.HasNormal();

  // Extended trajectory (not using memory for storage of arrays
  TParticleTrajectoryInterpolatedPoints TE;

  // Alternative outputs
  double Result_Precision = -1;
  int    Result_Level     = -1;

  // Loop over all points in the spectrum container
  for (size_t i = iFirst; i <= iLast; ++i) {

    // Obs point
    TVector3D const Obs    = Surface.GetPoint(i).GetPoint();
    TVector3D const Normal = Surface.GetPoint(i).GetNormal();

    double ThisSum = -1;
    double LastSum = -1;
    int    LastLevel = 0;

    // Summing for this power density
    double Sum = 0;

    for (int iLevel = 0; iLevel <= LevelStopWithExtended; ++iLevel) {
      LastLevel = iLevel;

      // Keep track of Beta for precision
      TVector3D Last_Beta(0, 0, 0);
      double BetaDiffMax = -1;

      // Grab the Trajectory (using memory arrays) if below level threshold, else set NULL
      TParticleTrajectoryPoints const& TM = Particle.GetTrajectoryLevel(iLevel <= LevelStopMemory ? iLevel : 0);
      if (iLevel > LevelStopMemory) {
        TE = Particle.GetTrajectoryExtendedLevel(iLevel);
      }

      // Number of points in the trajectory
      size_t const NTPoints = iLevel <= LevelStopMemory ? TM.GetNPoints() : TE.GetNPoints();

      // Loop over trajectory points
      for (int iT = 0; iT != NTPoints; ++iT) {

        TParticleTrajectoryPoint const& PP = (iLevel <= LevelStopMemory ? TM.GetPoint(iT) : TE.GetTrajectoryPoint(iT));

        // Get position, Beta, and Acceleration (over c)
        TVector3D const& X = PP.GetX();
        TVector3D const& B = PP.GetB();
        TVector3D const& AoverC = PP.GetAoverC();

        double const BetaDiff = (B - Last_Beta).Mag();
        if (iT > 0 && BetaDiff > BetaDiffMax) {
          BetaDiffMax = BetaDiff;
        }
        Last_Beta = B;

        // Define the three normal vectors.  N1 is in the direction of propogation,
        // N2 and N3 are in a plane perpendicular to N1
        TVector3D const N1 = (Obs - X).UnitVector();
        TVector3D const N2 = N1.Orthogonal().UnitVector();
        TVector3D const N3 = N1.Cross(N2).UnitVector();

        // For computing non-normally incidence
        double const N1DotNormal = HasNormal ? N1.Dot(Normal) : 1;
        // UPDATE: URGENT: Check this
        if (Directional && N1DotNormal <= 0) {
          continue;
        }

        // Compute Numerator and denominator
        Numerator = N1.Cross( ( (N1 - B).Cross((AoverC)) ) );
        Denominator = pow(1 - (B).Dot(N1), 5);

        // Add contributions from both N2 and N3
        Sum += pow(Numerator.Dot(N2), 2) / Denominator / (Obs - X).Mag2() * N1DotNormal;
        Sum += pow(Numerator.Dot(N3), 2) / Denominator / (Obs - X).Mag2() * N1DotNormal;

      }

      double const ThisSum = Sum * Particle.GetTrajectoryInterpolated().GetDeltaTInclusiveToLevel(iLevel);

      Result_Precision = fabs(ThisSum - LastSum) / LastSum;
      if (iLevel > 8 && Result_Precision < Precision && BetaDiffMax < 2. / (Particle.GetGamma())) {
        Result_Level = iLevel;
        break;
      } else if (iLevel > 8 && ThisSum == LastSum) {
        // The assumption here is that zero is last and now
        Result_Level = iLevel;
        Result_Precision = 0;
        break;
      }

      LastSum = ThisSum;
    }

    // If a point does not converge mark it
    if (Result_Level == -1) {
      PowerDensityContainer.SetNotConverged(i);
    }

    Sum *= fabs(Particle.GetQ() * Particle.GetCurrent()) / (16 * TOSCARSSR::Pi2() * TOSCARSSR::Epsilon0() * TOSCARSSR::C()) * Particle.GetTrajectoryInterpolated().GetDeltaTInclusiveToLevel(LastLevel);

    // m^2 to mm^2
    Sum /= 1e6;

    if (!Directional) {
      if (Sum < 0) {
        Sum *= -1;
      }
    }

    // Add to container
    switch (ReturnQuantity) {
      case 1:
        PowerDensityContainer.AddToPoint(i, Result_Precision * Weight);
        break;
      case 2:
        PowerDensityContainer.AddToPoint(i, ((double) Result_Level) * Weight);
        break;
      default:
        PowerDensityContainer.AddToPoint(i, Sum * Weight);
        break;
    }

  } // POINTS


  // Set done to true
  Done = true;

  return;
}




void OSCARSSR::CalculatePowerDensityThreads (TParticleA& Particle,
                                             TSurfacePoints const& Surface,
                                             T3DScalarContainer& PowerDensityContainer,
                                             int const NThreads,
                                             bool const Directional,
                                             double const Precision,
                                             int    const MaxLevel,
                                             int    const MaxLevelExtended,
                                             double const Weight,
                                             int    const ReturnQuantity)
{
  // Calculates the single particle power density on surface
  // in units of [watts / second / mm^2]
  //
  // Surface - Observation Points

  // Calculate trajectory if it doesn't exist
  if (Particle.GetTrajectory().GetNPoints() == 0) {
    this->CalculateTrajectory(Particle);
  }



  // Vector for storing threads to rejoin
  std::vector<std::thread> Threads;

  // Number of points in spectrum
  size_t const NPoints = Surface.GetNPoints();

  // How many threads to start in the first for loop
  size_t const NThreadsActual = NPoints > (size_t) NThreads ? NThreads : NPoints;

  // Keep track of which threads are finished and re-joined
  bool *Done = new bool[NThreadsActual];
  bool *Joined = new bool[NThreadsActual];

  // Number per thread plus remainder to be added to first threads
  size_t const NPerThread = NPoints / NThreadsActual;
  size_t const NRemainder = NPoints % NThreadsActual;

  // Start threads and keep in vector
  for (size_t it = 0; it != NThreadsActual; ++it) {

    // First and last points for each thread
    size_t const iFirst = it < NRemainder ? NPerThread * it + it: NPerThread * it + NRemainder;
    size_t const iLast  = it < NRemainder ? iFirst + NPerThread : iFirst + NPerThread - 1;

    // Set Done and joined to false for this thread
    Done[it] = false;
    Joined[it] = false;

    // Start thread for these points
    Threads.push_back(std::thread(&OSCARSSR::CalculatePowerDensityPoints,
                                  this,
                                  std::ref(Particle),
                                  std::ref(Surface),
                                  std::ref(PowerDensityContainer),
                                  iFirst,
                                  iLast,
                                  std::ref(Done[it]),
                                  Directional,
                                  Precision,
                                  MaxLevel,
                                  MaxLevelExtended,
                                  Weight,
                                  ReturnQuantity));
  }

  // Are all of the threads finished or not?  Continue loop until all come back.
  bool AllThreadsFinished = false;
  size_t NThreadsFinished = 0;
  while (!AllThreadsFinished) {

    // So as to not use the current thread at 100%
    std::this_thread::sleep_for(std::chrono::milliseconds(10));

    // Check all threads
    for (size_t it = 0; it != NThreadsActual; ++it) {

      if (Done[it] && !Joined[it]) {
        Threads[it].join();
        Joined[it] = true;
        ++NThreadsFinished;
      }
    }

    // If the number finished is equal to the number of points total then we're done
    if (NThreadsFinished == NThreadsActual) {
      AllThreadsFinished = true;
    }
  }


  // Clear all threads
  Threads.clear();

  // Delete my arrays, I hate new for this purpose
  delete [] Done;
  delete [] Joined;

  return;
}










void OSCARSSR::CalculatePowerDensityGPU (TSurfacePoints const& Surface,
                                         T3DScalarContainer& PowerDensityContainer,
                                         int const NParticles,
                                         std::vector<int> GPUVector,
                                         bool const Directional,
                                         double const Precision,
                                         int    const MaxLevel,
                                         int    const MaxLevelExtended,
                                         int    const ReturnQuantity)
{
  // If you compile for Cuda use the GPU in this function, else throw

  // If GPUVector is empty assume you want to use ALL GPUs available
  if (GPUVector.size() == 0) {
    int const NGPUAvailable = this->CheckGPU();
    for (int i = 0; i < NGPUAvailable; ++i) {
      GPUVector.push_back(i);
    }
  }

  #ifdef CUDA
  // Check that the GPU exists
  if (this->CheckGPU() < 1) {
    throw std::invalid_argument("You are requesting the GPU, but none were found");
  }

  // Set the level to stop at if requested, but not above the hard limit
  int const LevelStopMemory = MaxLevel >= -1  && MaxLevel <= TParticleA::kMaxTrajectoryLevel ? MaxLevel : TParticleA::kMaxTrajectoryLevel;
  int const LevelStopWithExtended = MaxLevelExtended > LevelStopMemory ? MaxLevelExtended : LevelStopMemory;


  return OSCARSSR_Cuda_CalculatePowerDensityGPU (*this,
                                                 Surface,
                                                 PowerDensityContainer,
                                                 NParticles,
                                                 GPUVector,
                                                 Precision,
                                                 LevelStopWithExtended,
                                                 ReturnQuantity);
  #else
  throw std::invalid_argument("GPU functionality not compiled into this binary distribution");
  #endif

  return;
}





double OSCARSSR::CalculateTotalPower ()
{
  // UPDATE: commet

  // Check that particle has been set yet.  If fType is "" it has not been set yet
  if (fParticle.GetType() == "") {
    try {
      this->SetNewParticle();
    } catch (std::exception e) {
      throw std::out_of_range("no beam defined");
    }
  }

  return this->CalculateTotalPower(fParticle);
}





double OSCARSSR::CalculateTotalPower (TParticleA& Particle)
{
  // Calculate total power out

  // Check that particle has been set yet.  If fType is "" it has not been set yet
  if (Particle.GetType() == "") {
    throw std::out_of_range("no particle defined");
  }

  // Grab the Trajectory
  TParticleTrajectoryPoints& T = Particle.GetTrajectory();

  // Calculate trajectory if it doesn't exist
  if (Particle.GetTrajectory().GetNPoints() == 0) {
    this->CalculateTrajectory(Particle);
  }


  // Number of points in Trajectory
  size_t const NTPoints = T.GetNPoints();

  // Timestep from trajectory
  double const DeltaT = T.GetDeltaT();

  // For summing total power
  double TotalPower = 0;


  // Loop over all points in trajectory
  for (size_t i = 0; i != NTPoints; ++i) {
    TVector3D const& B = T.GetB(i);
    TVector3D const& AoverC = T.GetAoverC(i);

    TotalPower += (AoverC.Mag2() - (B.Cross(AoverC)).Mag2()) * DeltaT;

  }

  // Undulators, Wigglers and their applications, p42
  TotalPower *= fabs(Particle.GetQ() * Particle.GetCurrent()) * pow(Particle.GetGamma(), 6) / (6 * TOSCARSSR::Pi() * TOSCARSSR::Epsilon0() * TOSCARSSR::C());


  return TotalPower;
}









void OSCARSSR::CalculateFlux (TParticleA& Particle,
                              TSurfacePoints const& Surface,
                              double const Energy_eV,
                              T3DScalarContainer& FluxContainer,
                              std::string const& Polarization,
                              double const Angle,
                              TVector3D const& HorizontalDirection,
                              TVector3D const& PropogationDirection,
                              double const Precision,
                              int    const MaxLevel,
                              int    const MaxLevelExtended,
                              double const Weight,
                              int    const ReturnQuantity)
{
  // Calculates the single particle flux
  // in units of [photons / second / 0.001% BW / mm^2]

  // Calculate trajectory if it doesn't exist
  if (Particle.GetTrajectory().GetNPoints() == 0) {
    this->CalculateTrajectory(Particle);
  }

  // Extra inpts for calculation
  bool Done = false;
  size_t const iFirst = 0;
  size_t const iLast = Surface.GetNPoints() - 1;

  // Calculate the flux
  CalculateFluxPoints(Particle,
                      Surface,
                      Energy_eV,
                      FluxContainer,
                      iFirst,
                      iLast,
                      Done,
                      Polarization,
                      Angle,
                      HorizontalDirection,
                      PropogationDirection,
                      Precision,
                      MaxLevel,
                      MaxLevelExtended,
                      Weight,
                      ReturnQuantity);

  return;
}





void OSCARSSR::CalculateFlux (TSurfacePoints const& Surface,
                              double const Energy_eV,
                              T3DScalarContainer& FluxContainer,
                              std::string const& Polarization,
                              double const Angle,
                              TVector3D const& HorizontalDirection,
                              TVector3D const& PropogationDirection,
                              int const NParticles,
                              int const NThreads,
                              int const GPU,
                              int const NGPU,
                              std::vector<int> VGPU,
                              double const Precision,
                              int    const MaxLevel,
                              int    const MaxLevelExtended,
                              int    const Dimension,
                              int    const ReturnQuantity)
{
  // Calculate flux on surface
  // THIS is the ENTRY POINT typically

  // Check that particle has been set yet.  If fType is "" it has not been set yet
  if (fParticle.GetType() == "") {
    try {
      this->SetNewParticle();
    } catch (std::exception e) {
      throw std::out_of_range("no beam defined");
    }
  }

  // Number of threads to possibly use
  int const NThreadsToUse = NThreads < 1 ? fNThreadsGlobal : NThreads;
  if (NThreadsToUse <= 0) {
    throw std::out_of_range("NThreads or NThreadsGlobal must be >= 1");
  }

  // Should we use the GPU or not?
  int const NGPUAvailable = this->CheckGPU();
  bool const UseGPU = GPU == 0 ? false : this->GetUseGPUGlobal() && (NGPUAvailable > 0) ? true : (GPU == 1) && (NGPUAvailable > 0);

  // If use GPU, let's check the GPU vector, or construct one if we must
  std::vector<int> GPUVector;
  for (std::vector<int>::const_iterator it = VGPU.begin(); it != VGPU.end(); ++it) {
    GPUVector.push_back(*it);
  }
  if (GPUVector.size() == 0) {
    for (int i = 0; i < NGPU; ++i) {
      GPUVector.push_back(i);
    }
  }
  if (NGPU != -1 && NGPU < (int) GPUVector.size()) {
    GPUVector.resize(NGPU);
  }


  if (Dimension == 3) {
    for (size_t i = 0; i != Surface.GetNPoints(); ++i) {
      FluxContainer.AddPoint(Surface.GetPoint(i).GetPoint(), 0);
    }
  } else if (Dimension == 2) {
    for (size_t i = 0; i != Surface.GetNPoints(); ++i) {
      FluxContainer.AddPoint( TVector3D(Surface.GetX1(i), Surface.GetX2(i), 0), 0);
    }
  } else {
    throw std::out_of_range("wROng dimension");
  }

  // Check polarization
  if (Polarization == "all" ||
      Polarization == "linear-horizontal" ||
      Polarization == "linear-vertical"   ||
      Polarization == "linear"            ||
      Polarization == "circular-left"     ||
      Polarization == "circular-right") {
    // Do nothing
  } else {
    throw std::invalid_argument("Polarization requested not recognized");
  }


  // Which cpmpute method will we use, gpu, multi-thread, or single-thread
  if (UseGPU) {
    // Send to GPU function
    CalculateFluxGPU(Surface,
                     Energy_eV,
                     FluxContainer,
                     Polarization,
                     Angle,
                     HorizontalDirection,
                     PropogationDirection,
                     NParticles,
                     GPUVector,
                     Precision,
                     MaxLevel,
                     MaxLevelExtended,
                     ReturnQuantity);
  } else {
    if (NParticles == 0) {
      if (NThreadsToUse == 1) {
        this->CalculateFlux(fParticle,
                            Surface,
                            Energy_eV,
                            FluxContainer,
                            Polarization,
                            Angle,
                            HorizontalDirection,
                            PropogationDirection,
                            Precision,
                            MaxLevel,
                            MaxLevelExtended,
                            1,
                            ReturnQuantity);
      } else {
        this->CalculateFluxThreads(fParticle,
                                   Surface,
                                   Energy_eV,
                                   FluxContainer,
                                   Polarization,
                                   Angle,
                                   HorizontalDirection,
                                   PropogationDirection,
                                   NThreadsToUse,
                                   Precision,
                                   MaxLevel,
                                   MaxLevelExtended,
                                   1,
                                   ReturnQuantity);
      }
    } else {
      // Weight this by the number of particles
      double const Weight = 1.0 / (double) NParticles;

      // Loop over particles
      for (int i = 0; i != NParticles; ++i) {

        // Set a new random particle
        this->SetNewParticle();
        this->CalculateTrajectory();

        if (NThreadsToUse == 1) {
          this->CalculateFlux(fParticle,
                              Surface,
                              Energy_eV,
                              FluxContainer,
                              Polarization,
                              Angle,
                              HorizontalDirection,
                              PropogationDirection,
                              Precision,
                              MaxLevel,
                              MaxLevelExtended,
                              Weight,
                              ReturnQuantity);
        } else {
          this->CalculateFluxThreads(fParticle,
                                     Surface,
                                     Energy_eV,
                                     FluxContainer,
                                     Polarization,
                                     Angle,
                                     HorizontalDirection,
                                     PropogationDirection,
                                     NThreadsToUse,
                                     Precision,
                                     MaxLevel,
                                     MaxLevelExtended,
                                     Weight,
                                     ReturnQuantity);
        }
      }
    }
  }


  return;
}






void OSCARSSR::CalculateFluxPoints (TParticleA& Particle,
                                    TSurfacePoints const& Surface,
                                    double const Energy_eV,
                                    T3DScalarContainer& FluxContainer,
                                    size_t const iFirst,
                                    size_t const iLast,
                                    bool& Done,
                                    std::string const& Polarization,
                                    double const Angle,
                                    TVector3D const& HorizontalDirection,
                                    TVector3D const& PropogationDirection,
                                    double const Precision,
                                    int    const MaxLevel,
                                    int    const MaxLevelExtended,
                                    double const Weight,
                                    int    const ReturnQuantity)
{
  // Calculates the single particle flux at a given observation point
  // in units of [photons / second / 0.001% BW / mm^2]
  //
  // Particle - the Particle.. with a Trajectory structure hopefully

  // Check number of points
  //if (NTPoints < 1) {
  //  throw std::length_error("no points in trajectory.  Is particle or beam defined?");
  //}

  // Check you are not requesting a level above the maximum
  if (MaxLevel > TParticleA::kMaxTrajectoryLevel) {
    std::cerr << "WARNING: MaxLevel > TParticleA::kMaxTrajectoryLevel.  Setting MaxLevel to TParticleA::kMaxTrajectoryLevel" << std::endl;
  }

  // Set the level to stop at if requested, but not above the hard limit
  int const LevelStopMemory = MaxLevel >= -1  && MaxLevel <= TParticleA::kMaxTrajectoryLevel ? MaxLevel : TParticleA::kMaxTrajectoryLevel;
  int const LevelStopWithExtended = MaxLevelExtended > LevelStopMemory ? MaxLevelExtended : LevelStopMemory;

  // Constant C0 for calculation
  double const C0 = Particle.GetQ() / (TOSCARSSR::FourPi() * TOSCARSSR::C() * TOSCARSSR::Epsilon0() * TOSCARSSR::Sqrt2Pi());

  // Constant for flux calculation at the end
  double const C2 = TOSCARSSR::FourPi() * Particle.GetCurrent() / (TOSCARSSR::H() * fabs(Particle.GetQ()) * TOSCARSSR::Mu0() * TOSCARSSR::C()) * 1e-6 * 0.001;

  // Imaginary "i" and complxe 1+0i
  std::complex<double> const I(0, 1);
  std::complex<double> const One(1, 0);

  // Photon vertical direction and positive and negative helicity
  TVector3D const VerticalDirection = PropogationDirection.Cross(HorizontalDirection).UnitVector();
  TVector3DC const Positive = 1. / sqrt(2) * (TVector3DC(HorizontalDirection) + VerticalDirection * I );
  TVector3DC const Negative = 1. / sqrt(2) * (TVector3DC(HorizontalDirection) - VerticalDirection * I );

  TVector3DC PolarizationVector(0, 0, 0);
  if (Polarization == "all") {
    // Do Nothing
  } else if (Polarization == "linear-horizontal") {
    PolarizationVector = HorizontalDirection;
  } else if (Polarization == "linear-vertical") {
    PolarizationVector = VerticalDirection;
  } else if (Polarization == "linear") {
    TVector3D PolarizationAngle = HorizontalDirection;
    PolarizationAngle.RotateSelf(Angle, PropogationDirection);
    PolarizationVector = PolarizationAngle;
  } else if (Polarization == "circular-left") {
    PolarizationVector = Positive;
  } else if (Polarization == "circular-right") {
    PolarizationVector = Negative;
  } else {
    // Throw invalid argument if polarization is not recognized
    throw std::invalid_argument("Polarization requested not recognized");
  }


  // Angular frequency
  double const Omega = TOSCARSSR::EvToAngularFrequency(Energy_eV);;

  // Extended trajectory (not using memory for storage of arrays
  TParticleTrajectoryInterpolatedPoints TE;

  // Alternative outputs
  double Result_Precision = -1;
  int    Result_Level     = -1;

  // Loop over all points in the spectrum container
  for (size_t i = iFirst; i <= iLast; ++i) {

    // Obs point
    TVector3D ObservationPoint = Surface.GetPoint(i).GetPoint();

    // Electric field summation in frequency space
    TVector3DC SumE(0, 0, 0);

    double ThisMag = -1;
    double LastMag = -1;
    double ThisPhase = -1;
    double LastPhase = -1;
    double MaxDPhase = 0;
    int    LastLevel = 0;

    bool PhaseSatisfied = true;
    for (int iLevel = 0; iLevel <= LevelStopWithExtended; ++iLevel) {
      LastLevel = iLevel;

      // Grab the Trajectory (using memory arrays) if below level threshold, else set NULL
      TParticleTrajectoryPoints const& TM = Particle.GetTrajectoryLevel(iLevel <= LevelStopMemory ? iLevel : 0);
      if (iLevel > LevelStopMemory) {
        TE = Particle.GetTrajectoryExtendedLevel(iLevel);
      }

      // Number of points in the trajectory
      size_t const NTPoints = iLevel <= LevelStopMemory ? TM.GetNPoints() : TE.GetNPoints();


      MaxDPhase = 0;
      // Loop over trajectory points
      for (int iT = 0; iT != NTPoints; ++iT) {

        TParticleTrajectoryPoint const& PP = (iLevel <= LevelStopMemory ? TM.GetPoint(iT) : TE.GetTrajectoryPoint(iT));

        // Get position, Beta, and Acceleration (over c)
        TVector3D const& X = PP.GetX();
        TVector3D const& B = PP.GetB();
        TVector3D const& AoverC = PP.GetAoverC();
        double    const  Time = iLevel <= LevelStopMemory ? TM.GetT(iT) : TE.GetT(iT);


        // Define R and unit vector in direction of R, and D (distance to observer)
        TVector3D const R = ObservationPoint - X;
        TVector3D const N = R.UnitVector();
        double const D = R.Mag();

        // Exponent in transformed field
        ThisPhase = -Omega * (Time + D / TOSCARSSR::C());
        double const PhaseTestValue = fabs(ThisPhase - LastPhase);
        if (iT != 0 && PhaseTestValue > MaxDPhase) {
          MaxDPhase = PhaseTestValue;
        }
        LastPhase = ThisPhase;
        std::complex<double> Exponent(0, ThisPhase);

        TVector3DC const ThisEw = ( ( (1 - (B).Mag2()) * (N - B) ) / ( D * D * (pow(1 - N.Dot(B), 2)) )
            + ( N.Cross( (N - B).Cross(AoverC) ) ) / ( D * pow(1 - N.Dot(B), 2) ) ) * std::exp(Exponent); // NF + FF

        // Add this contribution
        SumE += ThisEw;

      }

      TVector3DC ThisSumE = SumE * Particle.GetTrajectoryInterpolated().GetDeltaTInclusiveToLevel(iLevel);
      if (PolarizationVector.Mag2() > 0.001) {
        ThisSumE = ThisSumE.Dot(PolarizationVector) * PolarizationVector;
      }

      ThisMag = ThisSumE.Dot( ThisSumE.CC() ).real();

      Result_Precision = fabs(ThisMag - LastMag) / LastMag;
      if (iLevel > 8 && Result_Precision < Precision && MaxDPhase < TOSCARSSR::Pi()) {
        Result_Level = iLevel;
        break;
      }

      LastMag = ThisMag;
    }

    if (Result_Level == -1) {
      FluxContainer.SetNotConverged(i);
    }


    // Multiply by constant factor
    SumE *= C0 * Particle.GetTrajectoryInterpolated().GetDeltaTInclusiveToLevel(LastLevel);

    // Correcr for polarization
    if (PolarizationVector.Mag2() > 0.001) {
      SumE = SumE.Dot(PolarizationVector) * PolarizationVector;
    }


    // Add to container
    switch (ReturnQuantity) {
      case 1:
        FluxContainer.AddToPoint(i, Result_Precision * Weight);
        break;
      case 2:
        FluxContainer.AddToPoint(i, ((double) Result_Level) * Weight);
        break;
      default:
        FluxContainer.AddToPoint(i, C2 *  SumE.Dot( SumE.CC() ).real() * Weight);
        break;
    }
  } // POINTS

  // Set done to true
  Done = true;

  return;
}





void OSCARSSR::CalculateFluxThreads (TParticleA& Particle,
                                     TSurfacePoints const& Surface,
                                     double const Energy_eV,
                                     T3DScalarContainer& FluxContainer,
                                     std::string const& Polarization,
                                     double const Angle,
                                     TVector3D const& HorizontalDirection,
                                     TVector3D const& PropogationDirection,
                                     int const NThreads,
                                     double const Precision,
                                     int    const MaxLevel,
                                     int    const MaxLevelExtended,
                                     double const Weight,
                                     int    const ReturnQuantity)
{
  // Calculates the single particle flux on a surface
  // in units of [photons / second / 0.001% BW / mm^2]
  //
  // Surface - Observation Point

  // Calculate the trajectory only if it doesn't exist yet
  if (Particle.GetTrajectory().GetNPoints() == 0) {
    this->CalculateTrajectory(Particle);
  }

  // Vector for storing threads to rejoin
  std::vector<std::thread> Threads;

  // Number of points in spectrum
  size_t const NPoints = Surface.GetNPoints();

  // How many threads to start in the first for loop
  size_t const NThreadsActual = NPoints > (size_t) NThreads ? NThreads : NPoints;

  // Keep track of which threads are finished and re-joined
  bool *Done = new bool[NThreadsActual];
  bool *Joined = new bool[NThreadsActual];

  // Number per thread plus remainder to be added to first threads
  size_t const NPerThread = NPoints / NThreadsActual;
  size_t const NRemainder = NPoints % NThreadsActual;

  // Start threads and keep in vector
  for (size_t it = 0; it != NThreadsActual; ++it) {

    // First and last points for each thread
    size_t const iFirst = it < NRemainder ? NPerThread * it + it: NPerThread * it + NRemainder;
    size_t const iLast  = it < NRemainder ? iFirst + NPerThread : iFirst + NPerThread - 1;

    // Set Done and joined to false for this thread
    Done[it] = false;
    Joined[it] = false;

    // Start thread for these points
    Threads.push_back(std::thread(&OSCARSSR::CalculateFluxPoints,
                                  this,
                                  std::ref(Particle),
                                  std::ref(Surface),
                                  Energy_eV,
                                  std::ref(FluxContainer),
                                  iFirst,
                                  iLast,
                                  std::ref(Done[it]),
                                  Polarization,
                                  Angle,
                                  HorizontalDirection,
                                  PropogationDirection,
                                  Precision,
                                  MaxLevel,
                                  MaxLevelExtended,
                                  Weight,
                                  ReturnQuantity));
  }

  // Are all of the threads finished or not?  Continue loop until all come back.
  bool AllThreadsFinished = false;
  size_t NThreadsFinished = 0;
  while (!AllThreadsFinished) {

    // So as to not use the current thread at 100%
    std::this_thread::sleep_for(std::chrono::milliseconds(10));

    // Check all threads
    for (size_t it = 0; it != NThreadsActual; ++it) {

      if (Done[it] && !Joined[it]) {
        Threads[it].join();
        Joined[it] = true;
        ++NThreadsFinished;
      }
    }

    // If the number finished is equal to the number of points total then we're done
    if (NThreadsFinished == NThreadsActual) {
      AllThreadsFinished = true;
    }
  }

  // Clear all threads
  Threads.clear();

  // Delete my arrays, I hate new for this purpose
  delete [] Done;
  delete [] Joined;

  return;
}





void OSCARSSR::CalculateFluxGPU (TSurfacePoints const& Surface,
                                 double const Energy_eV,
                                 T3DScalarContainer& FluxContainer,
                                 std::string const& Polarization,
                                 double const Angle,
                                 TVector3D const& HorizontalDirection,
                                 TVector3D const& PropogationDirection,
                                 int const NParticles,
                                 std::vector<int> GPUVector,
                                 double const Precision,
                                 int    const MaxLevel,
                                 int    const MaxLevelExtended,
                                 int    const ReturnQuantity)
{
  // If you compile for Cuda use the GPU in this function, else throw

  // If GPUVector is empty assume you want to use ALL GPUs available
  if (GPUVector.size() == 0) {
    int const NGPUAvailable = this->CheckGPU();
    for (int i = 0; i < NGPUAvailable; ++i) {
      GPUVector.push_back(i);
    }
  }

  #ifdef CUDA
  // Check that the GPU exists
  if (this->CheckGPU() < 1) {
    throw std::invalid_argument("You are requesting the GPU, but none were found");
  }

  // Set the level to stop at if requested, but not above the hard limit
  int const LevelStopMemory = MaxLevel >= -1  && MaxLevel <= TParticleA::kMaxTrajectoryLevel ? MaxLevel : TParticleA::kMaxTrajectoryLevel;
  int const LevelStopWithExtended = MaxLevelExtended > LevelStopMemory ? MaxLevelExtended : LevelStopMemory;

  return OSCARSSR_Cuda_CalculateFluxGPU (*this,
                                         Surface,
                                         Energy_eV,
                                         FluxContainer,
                                         Polarization,
                                         Angle,
                                         HorizontalDirection,
                                         PropogationDirection,
                                         NParticles,
                                         GPUVector,
                                         Precision,
                                         LevelStopWithExtended,
                                         ReturnQuantity);
  #else
  throw std::invalid_argument("GPU functionality not compiled into this binary distribution");
  #endif

  return;
}





void OSCARSSR::CalculateElectricFieldTimeDomain (TVector3D const& Observer, T3DScalarContainer& XYZT)
{
  // Check that particle has been set yet.  If fType is "" it has not been set yet
  if (fParticle.GetType() == "") {
    try {
      this->SetNewParticle();
    } catch (std::exception e) {
      throw std::out_of_range("no beam defined");
    }
  }

  this->CalculateElectricFieldTimeDomain(Observer, XYZT, fParticle);
  return;
}





void OSCARSSR::CalculateElectricFieldTimeDomain (TVector3D const& Observer, T3DScalarContainer& XYZT,  TParticleA& Particle)
{
  // Check that particle has been set yet.  If fType is "" it has not been set yet
  if (Particle.GetType() == "") {
    throw std::out_of_range("no beam defined");
  }

  // Calculate trajectory if it doesn't exist
  if (Particle.GetTrajectory().GetNPoints() == 0) {
    this->CalculateTrajectory(Particle);
  }


  // Grab the Trajectory
  TParticleTrajectoryPoints& T = Particle.GetTrajectory();
  size_t const NTPoints = T.GetNPoints();

  double const DeltaT = T.GetDeltaT();

  double const C0 = Particle.GetQ() / (TOSCARSSR::FourPi() * TOSCARSSR::Epsilon0());
  // Loop over trajectory points
  for (size_t iT = 0; iT != NTPoints; ++iT) {

    // Get position, Beta, and Acceleration (over c)
    TVector3D const& X = T.GetX(iT);
    TVector3D const& B = T.GetB(iT);
    TVector3D const& AoverC = T.GetAoverC(iT);

    // Define R and unit vector in direction of R, and D (distance to observer)
    TVector3D const R = Observer - X;
    TVector3D const N = R.UnitVector();
    double const D = R.Mag();

    // Synchrotron Radiation, Philip john Duke, p74

    double    const      Mult = C0 * pow(1.0 / (1.0 - N.Dot(B)), 3);
    TVector3D const NearField = ((1.0 - B.Mag2() ) * (N - B)) / R.Mag2();
    TVector3D const  FarField = (1.0 / TOSCARSSR::C()) * (N.Cross(  (N - B).Cross(AoverC))  ) / R.Mag();

    TVector3D const EField = Mult * (NearField + FarField);
    double    const Time = iT * DeltaT + D / TOSCARSSR::C();

    XYZT.AddPoint(EField, Time);
  }

  //return TVector3D(Mult * (NearField + FarField));
  return;
}





TFieldContainer const& OSCARSSR::GetBFieldContainer () const
{
  // Return the bfield container
  return fBFieldContainer;
}





TFieldContainer const& OSCARSSR::GetEFieldContainer () const
{
  // Return the efield container
  return fEFieldContainer;
}




TDriftVolumeContainer const& OSCARSSR::GetDriftVolumeContainer () const
{
  // Return the drift volume container
  return fDriftVolumeContainer;
}
