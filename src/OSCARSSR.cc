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

#include "TVector3DC.h"
#include "TField3D_Grid.h"
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
                                 std::vector<double> const& Scaling)
{
  // Add a magnetic field from a file to the field container

  // Format string all upper-case (just in case you like to type L.C.).
  std::string FormatUpperCase = Format;
  std::transform(FormatUpperCase.begin(), FormatUpperCase.end(), FormatUpperCase.begin(), ::toupper);

  // Check that the format name is correct
  if ( (FormatUpperCase == "OSCARS" || FormatUpperCase == "SRW" || FormatUpperCase == "SPECTRA") ||
       (FormatUpperCase.size() > 8 && std::string(FormatUpperCase.begin(), FormatUpperCase.begin() + 8) == std::string("OSCARS1D"))) {

    this->fBFieldContainer.AddField( new TField3D_Grid(FileName, Format, Rotations, Translation, Scaling) );

  } else {
    throw std::invalid_argument("Incorrect format in format string");
  }

  // Set the derivs function accordingly
  this->SetDerivativesFunction();

  return;
}



void OSCARSSR::AddMagneticFieldInterpolated (std::vector<std::pair<double, std::string> > const& Mapping,
                                             std::string const Format,
                                             double const Parameter,
                                             TVector3D const& Rotations,
                                             TVector3D const& Translation,
                                             std::vector<double> const& Scaling)
{
  // Add a magnetic field from a file to the field container

  // Format string all upper-case (just in case you like to type L.C.).
  std::string FormatUpperCase = Format;
  std::transform(FormatUpperCase.begin(), FormatUpperCase.end(), FormatUpperCase.begin(), ::toupper);


  // Check that the format name is correct
  if ( (FormatUpperCase == "OSCARS"  || FormatUpperCase == "SRW" || FormatUpperCase == "SPECTRA") ||
     (FormatUpperCase.size() > 8 && std::string(FormatUpperCase.begin(), FormatUpperCase.begin() + 8) == std::string("OSCARS1D")) ) {

    this->fBFieldContainer.AddField( new TField3D_Grid(Mapping, Format, Parameter, Rotations, Translation, Scaling) );

  } else {
    throw std::invalid_argument("Incorrect format in format string");
  }

  // Set the derivs function accordingly
  this->SetDerivativesFunction();

  return;
}




void OSCARSSR::AddMagneticField (TField* Field)
{
  // Add a magnetic field from a file to the field container

  this->fBFieldContainer.AddField(Field);

  // Set the derivs function accordingly
  this->SetDerivativesFunction();

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




double OSCARSSR::GetBx (double const X, double const Y, double const Z) const
{
  // Return summed Bx from container
  return this->fBFieldContainer.GetFx(X, Y, Z);
}





double OSCARSSR::GetBy (double const X, double const Y, double const Z) const
{
  // Return summed By from container
  return this->fBFieldContainer.GetFy(X, Y, Z);
}





double OSCARSSR::GetBz (double const X, double const Y, double const Z) const
{
  // Return summed Bx from container
  return this->fBFieldContainer.GetFz(X, Y, Z);
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










void OSCARSSR::AddElectricField (std::string const FileName, std::string const Format, TVector3D const& Rotations, TVector3D const& Translation, std::vector<double> const& Scaling)
{
  // Add a electric field from a file to the field container
  this->fEFieldContainer.AddField( new TField3D_Grid(FileName, Format, Rotations, Translation) );

  // Set the derivs function accordingly
  this->SetDerivativesFunction();

  return;
}


void OSCARSSR::AddElectricField (TField* F)
{
  // Add a electric field from a file to the field container
  this->fEFieldContainer.AddField(F);

  // Set the derivs function accordingly
  this->SetDerivativesFunction();

  return;
}



void OSCARSSR::ClearElectricFields ()
{
  this->fEFieldContainer.Clear();

  // Set the derivs function accordingly
  this->SetDerivativesFunction();

  return;
}




void OSCARSSR::WriteField (std::string const& BorE, std::string const& OutFileName, std::string const& OutFormat, TVector2D const& XLim, int const NX, TVector2D const& YLim, int const NY, TVector2D const& ZLim, int const NZ, std::string const& Comment)
{
  if (BorE == "B") {
    fBFieldContainer.WriteToFile(OutFileName, OutFormat, XLim, NX, YLim, NY, ZLim, NZ, Comment);
  } else if (BorE == "E") {
    fEFieldContainer.WriteToFile(OutFileName, OutFormat, XLim, NX, YLim, NY, ZLim, NZ, Comment);
  } else {
    std::cerr << "Write failure because not B or E" << std::endl;
    throw;
  }

  return;
}




double OSCARSSR::GetEx (double const X, double const Y, double const Z) const
{
  // Return summed Ex from container
  return this->fEFieldContainer.GetFx(X, Y, Z);
}





double OSCARSSR::GetEy (double const X, double const Y, double const Z) const
{
  // Return summed Ey from container
  return this->fEFieldContainer.GetFy(X, Y, Z);
}





double OSCARSSR::GetEz (double const X, double const Y, double const Z) const
{
  // Return summed Ez from container
  return this->fEFieldContainer.GetFz(X, Y, Z);
}




TVector3D OSCARSSR::GetE (double const X, double const Y, double const Z) const
{
  // Return summed E from container
  return this->fEFieldContainer.GetF(X, Y, Z);
}




TVector3D OSCARSSR::GetE (TVector3D const& X) const
{
  // Return summed E from container
  return this->fEFieldContainer.GetF(X);
}
























void OSCARSSR::AddParticleBeam (std::string const& Type, std::string const& Name, TVector3D const& X0, TVector3D const& V0, double const Energy_GeV, double const T0, double const Current, double const Weight, double const Charge, double const Mass)
{
  // Add a particle beam
  // Type        - The name of the particle type that you want to use
  // Name        - A user specified 'name' for this beam
  // X0          - Initial position in X,Y,Z
  // V0          - A vector pointing in the direction of the velocity of arbitrary magnitude
  // Energy_GeV  - Energy of particle beam in GeV
  // T0          - Time of initial conditions, specified in units of m (for v = c)
  // Current     - Beam current in Amperes
  // Weight      - Relative weight to give this beam when randomly sampling
  // Charge      - Charge of custom particle
  // Mass        - Mass of custom particle

  fParticleBeamContainer.AddNewParticleBeam(Type, Name, X0, V0, Energy_GeV, T0, Current, Weight, Charge, Mass);
  return;
}




void OSCARSSR::AddParticleBeam (std::string const& Beam, std::string const& Name, double const Weight)
{
  // Add a particle beam
  // Beam - The name of the predefined particle beam to add

  fParticleBeamContainer.AddNewParticleBeam(Beam, Name, Weight);
  return;
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

  return;
}




void OSCARSSR::ClearParticleBeams ()
{
  // Clear the contents of the particle beam container
  fParticleBeamContainer.Clear();

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
  if (OSCARSSR_Cuda_GetDeviceCount() > 0) {
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
    return OSCARSSR_Cuda_GetDeviceCount();
  #endif
  return -1;
}




void OSCARSSR::SetNThreadsGlobal (int const N)
{
  fNThreadsGlobal = N;
  return;
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

  // UPDATE: Check this hits correct points in trajectory and correct number of points
  // UPDATE: Do a detailed check of forward and back prop of X, V, and A

  // Check that CTStart is not after T0 of particle
  if (this->GetCTStart() > P.GetT0()) {
    std::cerr << "ERROR: start time is greater than T0" << std::endl;
    throw std::out_of_range("start time is greater than T0");
  }

  // Check that CTStart and CTStop are not the same (probably not defined if this is the case)
  if (this->GetCTStart() >= this->GetCTStop()) {
    std::cerr << "ERROR: start time is >= stop time" << std::endl;
    throw std::out_of_range("start time is greater than stop time.  check that set_ctstartstop is set");
  }

  // Check that particle has been set yet.  If fType is "" it has not been set yet
  if (P.GetType() == "") {
    throw std::out_of_range("particle not initialized.  make sure you have a particle or beam defined");
  }

  // Clear any current trajectory
  P.GetTrajectory().Clear();


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

  // Set delta T for the trajectory
  ParticleTrajectory.SetDeltaT(DeltaT);

  // Loop over points in the forward direction
  for (int i = 0; i != NPointsForward; ++i) {

    // This time
    double t = P.GetT0() + DeltaT * i;


    // Add this point to the trajectory
    ParticleTrajectory.AddPoint(x[0], x[2], x[4], x[1] / TOSCARSSR::C(), x[3] / TOSCARSSR::C(), x[5] / TOSCARSSR::C(), dxdt[1] / TOSCARSSR::C(), dxdt[3] / TOSCARSSR::C(), dxdt[5] / TOSCARSSR::C());

    // Propogate
    (this->*fDerivativesFunction)(t, x, dxdt, P);
    //this->DerivativesB(t, x, dxdt, P);
    RK4(x, dxdt, N, t, DeltaT, x, P);
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
  for (int i = 0; i != NPointsBackward; ++i) {

    // This time
    double t = P.GetT0() + DeltaTReversed * (i + 1);

    // Propogate backward in time!
    (this->*fDerivativesFunction)(t, x, dxdt, P);
    //this->DerivativesB(t, x, dxdt, P);
    RK4(x, dxdt, N, t, DeltaTReversed, x, P);

    // Add the point to the trajectory
    ParticleTrajectory.AddPoint(x[0], x[2], x[4], x[1] / TOSCARSSR::C(), x[3] / TOSCARSSR::C(), x[5] / TOSCARSSR::C(), dxdt[1] / TOSCARSSR::C(), dxdt[3] / TOSCARSSR::C(), dxdt[5] / TOSCARSSR::C());
  }

  // Re-Reverse the trajectory to be in the proper time order
  ParticleTrajectory.ReverseArrays();

  return;
}




TParticleTrajectoryPoints const& OSCARSSR::GetTrajectory ()
{
  // Get the trajectory for *the* current particle in fParticle

  return fParticle.GetTrajectory();
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

  dxdt[0] = x[1];
  dxdt[1] = P.GetQoverMGamma() * E.GetX();
  dxdt[2] = x[3];                                                        
  dxdt[3] = P.GetQoverMGamma() * E.GetY();
  dxdt[4] = x[5];                                                        
  dxdt[5] = P.GetQoverMGamma() * E.GetZ();

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

  dxdt[0] = x[1];
  dxdt[1] = P.GetQoverMGamma() * (-x[5] * B.GetY() + x[3] * B.GetZ());
  dxdt[2] = x[3];                                                                                            
  dxdt[3] = P.GetQoverMGamma() * ( x[5] * B.GetX() - x[1] * B.GetZ());
  dxdt[4] = x[5];                                                                                            
  dxdt[5] = P.GetQoverMGamma() * ( x[1] * B.GetY() - x[3] * B.GetX());

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

  dxdt[0] = x[1];
  dxdt[1] = P.GetQoverMGamma() * (E.GetX() - x[5] * B.GetY() + x[3] * B.GetZ());
  dxdt[2] = x[3];                                                                                            
  dxdt[3] = P.GetQoverMGamma() * (E.GetY() + x[5] * B.GetX() - x[1] * B.GetZ());
  dxdt[4] = x[5];                                                                                            
  dxdt[5] = P.GetQoverMGamma() * (E.GetZ() + x[1] * B.GetY() - x[3] * B.GetX());

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

  dxdt[0] = x[1];
  dxdt[1] = P.GetQoverMGamma() * (-x[5] * B.GetY() + x[3] * B.GetZ());
  dxdt[2] = x[3];                                                                                            
  dxdt[3] = P.GetQoverMGamma() * ( x[5] * B.GetX() - x[1] * B.GetZ());
  dxdt[4] = x[5];                                                                                            
  dxdt[5] = P.GetQoverMGamma() * ( x[1] * B.GetY() - x[3] * B.GetX());

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





void OSCARSSR::CalculateSpectrumGPU (TParticleA& Particle,
                                     TVector3D const& ObservationPoint,
                                     TSpectrumContainer& Spectrum,
                                     std::string const& Polarization,
                                     double const Angle,
                                     TVector3D const& HorizontalDirection,
                                     TVector3D const& PropogationDirection,
                                     double const Weight)
{
  // Check that particle has been set yet.  If fType is "" it has not been set yet
  if (Particle.GetType() == "") {
    try {
      this->SetNewParticle();
    } catch (std::exception e) {
      throw std::out_of_range("no beam defined");
    }
  }

  // Calculate trajectory
  this->CalculateTrajectory(Particle);

  #ifdef CUDA
  // Check that the GPU exists
  if (this->CheckGPU() < 1) {
    throw std::invalid_argument("You are requesting the GPU, but none were found");
  }

  return OSCARSSR_Cuda_CalculateSpectrumGPU (Particle,
                                             ObservationPoint,
                                             Spectrum,
                                             Polarization,
                                             Angle,
                                             HorizontalDirection,
                                             PropogationDirection,
                                             Weight);
  #else
  throw std::invalid_argument("GPU functionality not compiled into this binary distribution");
  #endif

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
                                  int const GPU)
{
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
    std::cerr << "NThreads or NThreadsGlobal must be >= 1" << std::endl;
    throw;
  }

  // Should we use the GPU or not?
  bool const UseGPU = GPU == 0 ? false : this->GetUseGPUGlobal() && (this->CheckGPU() > 0) ? true : false;

  // GPU will outrank NThreads...
  if (NParticles == 0) {
    if (UseGPU == 0) {
      if (NThreadsToUse == 1) {
        this->CalculateSpectrum(fParticle,
                                ObservationPoint,
                                Spectrum,
                                Polarization,
                                Angle,
                                HorizontalDirection,
                                PropogationDirection,
                                1);
      } else {
        this->CalculateSpectrumThreads(fParticle,
                                       ObservationPoint,
                                       Spectrum,
                                       NThreadsToUse,
                                       Polarization,
                                       Angle,
                                       HorizontalDirection,
                                       PropogationDirection,
                                       1);
      }
    } else if (UseGPU == 1) {
      this->CalculateSpectrumGPU(fParticle,
                                 ObservationPoint,
                                 Spectrum,
                                 Polarization,
                                 Angle,
                                 HorizontalDirection,
                                 PropogationDirection,
                                 1);
    }
  } else {

    // Weight this by the number of particles
    double const Weight = 1.0 / (double) NParticles;

    // Loop over particles
    for (int i = 0; i != NParticles; ++i) {

      // Set a new random particle
      this->SetNewParticle();

      // GPU will outrank NThreads...
      if (UseGPU == 0) {
        if (NThreadsToUse == 1) {
          this->CalculateSpectrum(fParticle,
                                  ObservationPoint,
                                  Spectrum,
                                  Polarization,
                                  Angle,
                                  HorizontalDirection,
                                  PropogationDirection,
                                  Weight);
        } else {
          this->CalculateSpectrumThreads(fParticle,
                                         ObservationPoint,
                                         Spectrum,
                                         NThreadsToUse,
                                         Polarization,
                                         Angle,
                                         HorizontalDirection,
                                         PropogationDirection,
                                         Weight);
        }
      } else if (UseGPU == 1) {
        this->CalculateSpectrumGPU(fParticle,
                                   ObservationPoint,
                                   Spectrum,
                                   Polarization,
                                   Angle,
                                   HorizontalDirection,
                                   PropogationDirection,
                                   Weight);
      }
    }
  }

  return;
}




void OSCARSSR::CalculateSpectrumPoint (TParticleA& Particle,
                                       TVector3D const& ObservationPoint,
                                       TSpectrumContainer& Spectrum,
                                       int const i,
                                       bool& Done,
                                       std::string const& Polarization,
                                       double const Angle,
                                       TVector3D const& HorizontalDirection,
                                       TVector3D const& PropogationDirection,
                                       double const Weight)
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

  // Grab the Trajectory
  TParticleTrajectoryPoints& T = Particle.GetTrajectory();


  // Time step.  Expecting it to be constant throughout calculation
  double const DeltaT = T.GetDeltaT();


  // Number of points in the trajectory
  size_t const NTPoints = T.GetNPoints();

  if (NTPoints < 1) {
    throw std::length_error("no points in trajectory.  Is particle or beam defined?");
  }


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


  // Angular frequency
  double const Omega = Spectrum.GetAngularFrequency(i);

  // Constant for field calculation
  std::complex<double> ICoverOmega = I * TOSCARSSR::C() / Omega;

  // Constant for calculation
  std::complex<double> const C1(0, C0 * Omega);

  // Electric field summation in frequency space
  TVector3DC SumE(0, 0, 0);

  // Loop over all points in trajectory
  for (int iT = 0; iT != NTPoints; ++iT) {

    // Particle position
    TVector3D const& X = T.GetX(iT);

    // Particle "Beta" (velocity over speed of light)
    TVector3D const& B = T.GetB(iT);

    // Vector pointing from particle to observer
    TVector3D const R = ObservationPoint - X;

    // Unit vector pointing from particl to observer
    TVector3D const N = R.UnitVector();

    // Distance from particle to observer
    double const D = R.Mag();

    // Exponent for fourier transformed field
    std::complex<double> Exponent(0, Omega * (DeltaT * iT + D / TOSCARSSR::C()));

    // Sum in fourier transformed field (integral)
    SumE += (TVector3DC(B) - (N * ( One + (ICoverOmega / (D))))) / D * std::exp(Exponent);
  }

  // Multiply field by Constant C1 and time step
  SumE *= C1 * DeltaT;

  // If a polarization is specified, calculate it
  if (Polarization == "all") {
    // Do nothing, it is already ALL
  } else if (Polarization == "linear-horizontal") {
    SumE = SumE.Dot(HorizontalDirection) * HorizontalDirection;
  } else if (Polarization == "linear-vertical") {
    SumE = SumE.Dot(VerticalDirection) * VerticalDirection;
  } else if (Polarization == "linear") {
    TVector3D PolarizationAngle = HorizontalDirection;
    PolarizationAngle.RotateSelf(Angle, PropogationDirection);
    SumE = SumE.Dot(PolarizationAngle) * PolarizationAngle;
  } else if (Polarization == "circular-left") {
    SumE = SumE.Dot(Positive.CC()) * Positive;
  } else if (Polarization == "circular-right") {
    SumE = SumE.Dot(Negative.CC()) * Negative;
  } else {
    // Throw invalid argument if polarization is not recognized
    throw std::invalid_argument("Polarization requested not recognized");
  }

  // Set the flux for this frequency / energy point
  Spectrum.AddToFlux(i, C2 *  SumE.Dot( SumE.CC() ).real() * Weight);

  // Noify that this calculation is finished
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
                                         double const Weight)
{
  // Calculates spectrum for the given particle and observation point
  // in units of [photons / second / 0.001% BW / mm^2]
  //
  // Surface - Observation Point

  // Check that particle has been set yet.  If fType is "" it has not been set yet
  if (Particle.GetType() == "") {
    try {
      this->SetNewParticle();
    } catch (std::exception e) {
      throw std::out_of_range("no beam defined");
    }
  }

  // Calculate trajectory before we fanout into threads
  this->CalculateTrajectory(Particle);

  // Calculate the trajectory from scratch
  this->CalculateTrajectory(Particle);

  // Vector container for threads
  std::vector<std::thread> Threads;

  // Number of points in spectrum
  size_t const NPoints = Spectrum.GetNPoints();

  // How many threads to start in the first for loop
  size_t const NThreadsActual = NPoints > NThreads ? NThreads : NPoints;

  // Keep track of which threads are finished and started
  bool Done[NThreadsActual];
  bool Submitted[NThreadsActual];

  // Number of points for which a thread was started so far
  size_t NPointsStarted = 0;

  // Start threads and keep in vector
  for (size_t io = 0; io != NThreadsActual; ++io) {

    // Set Done to false for this thread
    Done[io] = false;

    // Start thread for this point
    Threads.push_back(std::thread(&OSCARSSR::CalculateSpectrumPoint,
                                  this,
                                  std::ref(Particle),
                                  std::ref(Obs),
                                  std::ref(Spectrum),
                                  (int) io,
                                  std::ref(Done[io]),
                                  std::ref(Polarization),
                                  Angle,
                                  std::ref(HorizontalDirection),
                                  std::ref(PropogationDirection),
                                  Weight));

    // Incement the number of points started so far
    ++NPointsStarted;

    // Set the submitted flag for this thread to true
    Submitted[io] = true;
  }

  // Keep track of how many threads have completed
  size_t NThreadsFinished = 0;

  // Are all of the threads finished or not?  Continue loop until all come back.
  // UPDATE: Could think about applying a timeout feature here
  bool AllThreadsFinished = false;
  while (!AllThreadsFinished) {

    // So as to not use the current thread at 100%
    std::this_thread::sleep_for(std::chrono::milliseconds(10));

    // Check all threads
    for (size_t it = 0; it != NThreadsActual; ++it) {

      // If it hasn't reported done leave it alone, otherwise see about creating another thread
      if (!Done[it]) {
        continue;
      } else if (Submitted[it] && Done[it]) {

        // Join the thread
        Threads[it].join();

        // Set the submitted bit to false (so it's open)
        Submitted[it] = false;

        // Increment the number of finished threads
        ++NThreadsFinished;

        // If the number of points started is less than the total number start another point
        if (NPointsStarted < NPoints) {

          // Set done bit to false for this thread
          Done[it] = false;

          // Index of the point of interest
          int const ThisPoint = (int) NPointsStarted;

          // Create a new thread for this point
          Threads[it] = std::thread(&OSCARSSR::CalculateSpectrumPoint,
                                    this,
                                    std::ref(Particle),
                                    std::ref(Obs),
                                    std::ref(Spectrum),
                                    ThisPoint,
                                    std::ref(Done[it]),
                                    std::ref(Polarization),
                                    Angle,
                                    std::ref(HorizontalDirection),
                                    std::ref(PropogationDirection),
                                    Weight);

          // Increment the number of points that have been started
          ++NPointsStarted;

          // Set the submitted bit to true
          Submitted[it] = true;
        }
      }
    }

    // If the number finished is equal to the number of points total then we're done
    if (NThreadsFinished == NPoints) {
      AllThreadsFinished = true;
    }
  }

  // Clear all threads
  Threads.clear();

  return;
}

















void OSCARSSR::CalculateSpectrum (TParticleA& Particle,
                                  TVector3D const& ObservationPoint,
                                  TSpectrumContainer& Spectrum,
                                  std::string const& Polarization,
                                  double const Angle,
                                  TVector3D const& HorizontalDirection,
                                  TVector3D const& PropogationDirection,
                                  double const Weight)
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

  // Calculate trajectory
  this->CalculateTrajectory(Particle);

  // Grab the Trajectory
  TParticleTrajectoryPoints& T = Particle.GetTrajectory();


  // Time step.  Expecting it to be constant throughout calculation
  double const DeltaT = T.GetDeltaT();


  // Number of points in the trajectory
  size_t const NTPoints = T.GetNPoints();

  if (NTPoints < 1) {
    throw std::length_error("no points in trajectory.  Is particle or beam defined?");
  }

  // Number of points in the spectrum container
  size_t const NEPoints = Spectrum.GetNPoints();

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


  // Loop over all points in the spectrum container
  for (size_t i = 0; i != NEPoints; ++i) {

    // Angular frequency
    double const Omega = Spectrum.GetAngularFrequency(i);

    // Constant for field calculation
    std::complex<double> ICoverOmega = I * TOSCARSSR::C() / Omega;

    // Constant for calculation
    std::complex<double> const C1(0, C0 * Omega);

    // Electric field summation in frequency space
    TVector3DC SumE(0, 0, 0);

    // Loop over all points in trajectory
    for (int iT = 0; iT != NTPoints; ++iT) {

      // Particle position
      TVector3D const& X = T.GetX(iT);

      // Particle "Beta" (velocity over speed of light)
      TVector3D const& B = T.GetB(iT);

      // Vector pointing from particle to observer
      TVector3D const R = ObservationPoint - X;

      // Unit vector pointing from particl to observer
      TVector3D const N = R.UnitVector();

      // Distance from particle to observer
      double const D = R.Mag();

      // Exponent for fourier transformed field
      std::complex<double> Exponent(0, Omega * (DeltaT * iT + D / TOSCARSSR::C()));

      // Sum in fourier transformed field (integral)
      SumE += (TVector3DC(B) - (N * ( One + (ICoverOmega / (D))))) / D * std::exp(Exponent);
    }

    // Multiply field by Constant C1 and time step
    SumE *= C1 * DeltaT;

    // If a polarization is specified, calculate it
    if (Polarization == "all") {
      // Do nothing, it is already ALL
    } else if (Polarization == "linear-horizontal") {
      SumE = SumE.Dot(HorizontalDirection) * HorizontalDirection;
    } else if (Polarization == "linear-vertical") {
      SumE = SumE.Dot(VerticalDirection) * VerticalDirection;
    } else if (Polarization == "linear") {
      TVector3D PolarizationAngle = HorizontalDirection;
      PolarizationAngle.RotateSelf(Angle, PropogationDirection);
      SumE = SumE.Dot(PolarizationAngle) * PolarizationAngle;
    } else if (Polarization == "circular-left") {
      SumE = SumE.Dot(Positive.CC()) * Positive;
    } else if (Polarization == "circular-right") {
      SumE = SumE.Dot(Negative.CC()) * Negative;
    } else {
      // Throw invalid argument if polarization is not recognized
      throw std::invalid_argument("Polarization requested not recognized");
    }


    // Set the flux for this frequency / energy point
    Spectrum.AddToFlux(i, C2 *  SumE.Dot( SumE.CC() ).real() * Weight);
  }



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
    std::cerr << "ERROR: incorrect dimension in spectrum" << std::endl;
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
    std::cerr << "ERROR: incorrect dimension in spectrum" << std::endl;
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
    std::cerr << "ERROR: incorrect dimension in spectrum" << std::endl;
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
                                      double const Weight)
{
  // Calculates the single particle spectrum at a given observation point
  // in units of [photons / second / 0.001% BW / mm^2]
  //
  // Particle - Particle, contains trajectory (or if not, calculate it)
  // Surface - Observation Point


  // Check that particle has been set yet.  If fType is "" it has not been set yet
  if (Particle.GetType() == "") {
    throw std::out_of_range("no particle defined");
  }

  // Calculate the trajectory from scratch
  this->CalculateTrajectory(Particle);

  // Grab the Trajectory
  TParticleTrajectoryPoints& T = Particle.GetTrajectory();

  // Number of points in Trajectory
  size_t const NTPoints = T.GetNPoints();

  // Timestep from trajectory
  double const DeltaT = T.GetDeltaT();

  // Variables for parts of the numerator and denominator of power density equation
  TVector3D Numerator;
  double Denominator;

  //std::cout << "Directional: " << Directional << std::endl;

  // Loop over all points in the given surface
  for (size_t io = 0; io < Surface.GetNPoints(); ++io) {

    // Get the observation point (on the surface, and its "normal"
    TVector3D const Obs = Surface.GetPoint(io).GetPoint();
    TVector3D const Normal = Surface.GetPoint(io).GetNormal();


    // For summing power contributions
    double Sum = 0;

    // Loop over all points in the trajectory
    for (int iT = 0; iT != NTPoints ; ++iT) {

      // Get current position, Beta, and Acceleration(over c)
      TVector3D const& X = T.GetX(iT);
      TVector3D const& B = T.GetB(iT);
      TVector3D const& AoverC = T.GetAoverC(iT);

      // Define the three normal vectors.  N1 is in the direction of propogation,
      // N2 and N3 are in a plane perpendicular to N1
      TVector3D const N1 = (Obs - X).UnitVector();
      TVector3D const N2 = N1.Orthogonal().UnitVector();
      TVector3D const N3 = N1.Cross(N2).UnitVector();

      // For computing non-normally incidence
      double const N1DotNormal = N1.Dot(Normal);
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

    // Undulators, Wigglers and their applications, p42
    Sum *= fabs(Particle.GetQ() * Particle.GetCurrent()) / (16 * TOSCARSSR::Pi2() * TOSCARSSR::Epsilon0() * TOSCARSSR::C()) * DeltaT;

    // m^2 to mm^2
    Sum /= 1e6;

    // Weighting
    Sum *= Weight;


    // If you don't care about the direction of the normal vector
    // UPDATE: Check
    if (!Directional) {
      if (Sum < 0) {
        Sum *= -1;
      }
    }

    // Add to container
    PowerDensityContainer.AddToPoint(io, Sum);
  }

  return;
}




void OSCARSSR::CalculatePowerDensity (TSurfacePoints const& Surface,
                                      T3DScalarContainer& PowerDensityContainer,
                                      int const Dimension,
                                      bool const Directional,
                                      int const NParticles,
                                      int const NThreads,
                                      int const GPU)
{
  // Calculates the power density
  // in units of [W / mm^2]
  //
  // UPDATE: inputs


  // How many threads to use.
  int const NThreadsToUse = NThreads < 1 ? fNThreadsGlobal : NThreads;
  if (NThreadsToUse <= 0) {
    std::cerr << "NThreads or NThreadsGlobal must be >= 1" << std::endl;
    throw;
  }

  // Should we use the GPU or not?
  bool const UseGPU = GPU == 0 ? false : this->GetUseGPUGlobal() && (this->CheckGPU() > 0) ? true : false;


  // Check that particle has been set yet.  If fType is "" it has not been set yet
  if (fParticle.GetType() == "") {
    try {
      this->SetNewParticle();
    } catch (std::exception e) {
      throw std::out_of_range("no beam defined");
    }
  }

  if (Dimension == 3) {
    for (int i = 0; i != Surface.GetNPoints(); ++i) {
      PowerDensityContainer.AddPoint(Surface.GetPoint(i).GetPoint(), 0);
    }
  } else if (Dimension == 2) {
    for (int i = 0; i != Surface.GetNPoints(); ++i) {
      PowerDensityContainer.AddPoint( TVector3D(Surface.GetX1(i), Surface.GetX2(i), 0), 0);
    }
  } else {
    std::cerr << "Wrong dimension" << std::endl;
    throw;
  }

  // GPU will outrank NThreads...
  if (NParticles == 0) {
    // nparticles = 0 has special meaning.  Set new ideal particle.  Of course if you have multiple
    // beams defined this is a bit strange as it will pick randomly
    // UPDATE: Consider adding "beam" as an input
    this->SetNewParticle("", "ideal");

    if (UseGPU == 0) {
      if (NThreadsToUse == 1) {
        this->CalculatePowerDensity(fParticle, Surface, PowerDensityContainer, Directional, 1);
      } else {
        this->CalculatePowerDensityThreads(fParticle, Surface, PowerDensityContainer, NThreadsToUse, Directional, 1);
      }
    } else if (UseGPU == 1) {
      this->CalculatePowerDensityGPU(fParticle, Surface, PowerDensityContainer, Directional, 1);
    }
  } else {
    double const Weight = 1.0 / (double) NParticles;
    for (int i = 0; i != NParticles; ++i) {
      this->SetNewParticle();
      if (UseGPU == 0) {
        if (NThreadsToUse == 1) {
          this->CalculatePowerDensity(fParticle, Surface, PowerDensityContainer, Directional, Weight);
        } else {
          this->CalculatePowerDensityThreads(fParticle, Surface, PowerDensityContainer,  NThreadsToUse, Directional, Weight);
        }
      } else if (UseGPU == 1) {
        this->CalculatePowerDensityGPU(fParticle, Surface, PowerDensityContainer, Directional, Weight);
      }
    }
  }

  return;
}






void OSCARSSR::CalculatePowerDensityPoint (TParticleA& Particle,
                                           TSurfacePoints const& Surface,
                                           T3DScalarContainer& PowerDensityContainer,
                                           size_t const io,
                                           bool& Done,
                                           bool const Directional,
                                           double const Weight)
{
  // Calculates the single particle spectrum at a given observation point
  // in units of [photons / second / 0.001% BW / mm^2]
  //
  // Particle - Particle, contains trajectory (or if not, calculate it)
  // Surface - Observation Point


  // Check that particle has been set yet.  If fType is "" it has not been set yet
  if (Particle.GetType() == "") {
    throw std::out_of_range("no particle defined");
  }


  // Grab the Trajectory
  TParticleTrajectoryPoints& T = Particle.GetTrajectory();

  // Number of points in Trajectory
  size_t const NTPoints = T.GetNPoints();

  // Timestep from trajectory
  double const DeltaT = T.GetDeltaT();

  // Variables for parts of the numerator and denominator of power density equation
  TVector3D Numerator;
  double Denominator;

  //std::cout << "Directional: " << Directional << std::endl;


  TVector3D const Obs = Surface.GetPoint(io).GetPoint();
  TVector3D const Normal = Surface.GetPoint(io).GetNormal();


  // For summing power contributions
  double Sum = 0;

  // Loop over all points in the trajectory
  for (int iT = 0; iT != NTPoints ; ++iT) {

    // Get current position, Beta, and Acceleration(over c)
    TVector3D const& X = T.GetX(iT);
    TVector3D const& B = T.GetB(iT);
    TVector3D const& AoverC = T.GetAoverC(iT);

    // Define the three normal vectors.  N1 is in the direction of propogation,
    // N2 and N3 are in a plane perpendicular to N1
    TVector3D const N1 = (Obs - X).UnitVector();
    TVector3D const N2 = N1.Orthogonal().UnitVector();
    TVector3D const N3 = N1.Cross(N2).UnitVector();

    // Speed up here if you want
    //if (B.Angle(N1) > 5. / (Particle.GetGamma())) {
    //  continue;
    //}

    // For computing non-normally incidence
    double const N1DotNormal = N1.Dot(Normal);
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
  // Undulators, Wigglers and their applications, p42
  Sum *= fabs(Particle.GetQ() * Particle.GetCurrent()) / (16 * TOSCARSSR::Pi2() * TOSCARSSR::Epsilon0() * TOSCARSSR::C()) * DeltaT;

  // m^2 to mm^2
  Sum /= 1e6;

  // Weight this event
  Sum *= Weight;

  // If you don't care about the direction of the normal vector
  // UPDATE: Check
  if (!Directional) {
    if (Sum < 0) {
      Sum *= -1;
    }
  }

  // Add to container
  PowerDensityContainer.AddToPoint(io, Sum);

  // Set the done bit
  Done = true;

  return;
}








void OSCARSSR::CalculatePowerDensityThreads (TParticleA& Particle,
                                             TSurfacePoints const& Surface,
                                             T3DScalarContainer& PowerDensityContainer,
                                             int const NThreads,
                                             bool const Directional,
                                             double const Weight)
{
  // Calculates the single particle spectrum at a given observation point
  // in units of [photons / second / 0.001% BW / mm^2]
  //
  // Surface - Observation Point

  // Check that particle has been set yet.  If fType is "" it has not been set yet
  if (Particle.GetType() == "") {
    try {
      this->SetNewParticle();
    } catch (std::exception e) {
      throw std::out_of_range("no beam defined");
    }
  }

  // Calculate the trajectory from scratch
  this->CalculateTrajectory(Particle);

  std::vector<std::thread> Threads;

  // Number of points in spectrum
  size_t const NPoints = Surface.GetNPoints();

  // How many threads to start in the first for loop
  size_t const NThreadsActual = (size_t) NPoints > NThreads ? NThreads : NPoints;

  // Keep track of which threads are finished and started
  bool Done[NThreadsActual];
  bool Submitted[NThreadsActual];

  // Number of points for which a thread was started so far
  size_t NPointsStarted = 0;

  // Start threads and keep in vector
  for (size_t it = 0; it != NThreadsActual; ++it) {

    // Set Done to false for this thread
    Done[it] = false;

    // Start thread for this point
    Threads.push_back(std::thread(&OSCARSSR::CalculatePowerDensityPoint,
                                  this,
                                  std::ref(Particle),
                                  std::ref(Surface),
                                  std::ref(PowerDensityContainer),
                                  (int) it,
                                  std::ref(Done[it]),
                                  Directional,
                                  Weight));

    // Incement the number of points started so far
    ++NPointsStarted;

    // Set the submitted flag for this thread to true
    Submitted[it] = true;
  }


  // Keep track of how many threads have completed
  size_t NThreadsFinished = 0;

  // Are all of the threads finished or not?  Continue loop until all come back.
  // UPDATE: Could think about applying a timeout feature here
  bool AllThreadsFinished = false;
  while (!AllThreadsFinished) {

    // So as to not use the current thread at 100%
    std::this_thread::sleep_for(std::chrono::milliseconds(10));

    // Check all threads
    for (size_t it = 0; it != NThreadsActual; ++it) {

      // If it hasn't reported done leave it alone, otherwise see about creating another thread
      if (!Done[it]) {
        continue;
      } else if (Submitted[it] && Done[it]) {

        // Join the thread
        Threads[it].join();

        // Set the submitted bit to false (so it's open)
        Submitted[it] = false;

        // Increment the number of finished threads
        ++NThreadsFinished;

        // If the number of points started is less than the total number start another point
        if (NPointsStarted < NPoints) {

          // Set done bit to false for this thread
          Done[it] = false;

          // Index of the point of interest
          int const ThisPoint = (int) NPointsStarted;

          // Create a new thread for this point
          Threads[it] = std::thread(&OSCARSSR::CalculatePowerDensityPoint,
                                    this,
                                    std::ref(Particle),
                                    std::ref(Surface),
                                    std::ref(PowerDensityContainer),
                                    ThisPoint,
                                    std::ref(Done[it]),
                                    Directional,
                                    Weight);

          // Increment the number of points that have been started
          ++NPointsStarted;

          // Set the submitted bit to true
          Submitted[it] = true;
        }
      }
    }

    // If the number finished is equal to the number of points total then we're done
    if (NThreadsFinished == NPoints) {
      AllThreadsFinished = true;
    }
  }

  // Clear all threads
  Threads.clear();

  return;
}










void OSCARSSR::CalculatePowerDensityGPU (TParticleA& Particle,
                                         TSurfacePoints const& Surface,
                                         T3DScalarContainer& PowerDensityContainer,
                                         bool const Directional,
                                         double const Weight)
{
  // If you compile for Cuda use the GPU in this function, else throw

  // Check that particle has been set yet.  If fType is "" it has not been set yet
  if (Particle.GetType() == "") {
    throw std::out_of_range("no particle defined");
  }

  // Calculate the trajectory from scratch
  this->CalculateTrajectory(Particle);

  #ifdef CUDA
  // Check that the GPU exists
  if (this->CheckGPU() < 1) {
    throw std::invalid_argument("You are requesting the GPU, but none were found");
  }

  return OSCARSSR_Cuda_CalculatePowerDensityGPU (Particle, Surface, PowerDensityContainer, Directional, Weight);
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

  // Calculate Trajectory from scratch
  this->CalculateTrajectory(Particle);

  // Number of points in Trajectory
  size_t const NTPoints = T.GetNPoints();

  // Timestep from trajectory
  double const DeltaT = T.GetDeltaT();

  // For summing total power
  double TotalPower = 0;


  // Loop over all points in trajectory
  for (int i = 0; i != NTPoints; ++i) {
    TVector3D const& B = T.GetB(i);
    TVector3D const& AoverC = T.GetAoverC(i);

    TotalPower += (AoverC.Mag2() - (B.Cross(AoverC)).Mag2()) * DeltaT;

  }

  // Undulators, Wigglers and their applications, p42
  TotalPower *= fabs(Particle.GetQ() * Particle.GetCurrent()) * pow(Particle.GetGamma(), 6) / (6 * TOSCARSSR::Pi() * TOSCARSSR::Epsilon0() * TOSCARSSR::C());


  return TotalPower;
}










void OSCARSSR::CalculateFluxPoint (TParticleA& Particle,
                                   TSurfacePoints const& Surface,
                                   double const Energy_eV,
                                   T3DScalarContainer& FluxContainer,
                                   size_t const i,
                                   bool& Done,
                                   std::string const& Polarization,
                                   double const Angle,
                                   TVector3D const& HorizontalDirection,
                                   TVector3D const& PropogationDirection,
                                   int const Dimension,
                                   double const Weight)
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

  // Grab the Trajectory
  TParticleTrajectoryPoints& T = Particle.GetTrajectory();

  // Time step.  Expecting it to be constant throughout calculation
  double const DeltaT = T.GetDeltaT();


  // Number of points in the trajectory
  size_t const NTPoints = T.GetNPoints();

  if (NTPoints < 1) {
    throw std::length_error("no points in trajectory.  Is particle or beam defined?");
  }

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


  // Angular frequency
  double const Omega = TOSCARSSR::EvToAngularFrequency(Energy_eV);;

  // Constant for field calculation
  std::complex<double> ICoverOmega = I * TOSCARSSR::C() / Omega;

  // Constant for calculation
  std::complex<double> const C1(0, C0 * Omega);

  // Obs point
  TVector3D ObservationPoint = Surface.GetPoint(i).GetPoint();

  // Electric field summation in frequency space
  TVector3DC SumE(0, 0, 0);

  // Loop over all points in trajectory
  for (int iT = 0; iT != NTPoints; ++iT) {

    // Particle position
    TVector3D const& X = T.GetX(iT);

    // Particle "Beta" (velocity over speed of light)
    TVector3D const& B = T.GetB(iT);

    // Vector pointing from particle to observer
    TVector3D const R = ObservationPoint - X;

    // Unit vector pointing from particl to observer
    TVector3D const N = R.UnitVector();

    // Distance from particle to observer
    double const D = R.Mag();

    // Exponent for fourier transformed field
    std::complex<double> Exponent(0, Omega * (DeltaT * iT + D / TOSCARSSR::C()));

    // Sum in fourier transformed field (integral)
    SumE += (TVector3DC(B) - (N * ( One + (ICoverOmega / (D))))) / D * std::exp(Exponent);
  }


  // Multiply field by Constant C1 and time step
  SumE *= C1 * DeltaT;

  // If a polarization is specified, calculate it
  if (Polarization == "all") {
    // Do nothing, it is already ALL
  } else if (Polarization == "linear-horizontal") {
    SumE = SumE.Dot(HorizontalDirection) * HorizontalDirection;
  } else if (Polarization == "linear-vertical") {
    SumE = SumE.Dot(VerticalDirection) * VerticalDirection;
  } else if (Polarization == "linear") {
    TVector3D PolarizationAngle = HorizontalDirection;
    PolarizationAngle.RotateSelf(Angle, PropogationDirection);
    SumE = SumE.Dot(PolarizationAngle) * PolarizationAngle;
  } else if (Polarization == "circular-left") {
    SumE = SumE.Dot(Positive.CC()) * Positive;
  } else if (Polarization == "circular-right") {
    SumE = SumE.Dot(Negative.CC()) * Negative;
  } else {
    // Throw invalid argument if polarization is not recognized
    throw std::invalid_argument("Polarization requested not recognized");
  }

  // Set the flux for this frequency / energy point
  double const ThisFlux = C2 *  SumE.Dot( SumE.CC() ).real() * Weight;

  // Add point to flux container
  FluxContainer.AddToPoint(i, ThisFlux);

  // Set done to true
  Done = true;

  return;
}




void OSCARSSR::CalculateFlux1 (TParticleA& Particle, TSurfacePoints const& Surface, double const Energy_eV, T3DScalarContainer& FluxContainer)
{
  // Here for historical reasons.
  //
  // Calculates the single particle spectrum at a given observation point
  // in units of [photons / second / 0.001% BW / mm^2]
  //
  // T - Trajectory of particle
  // Surface - Surface of Observation Points
  // Current - beam current
  // Energy - beam energy in eV

  // Check that particle has been set yet.  If fType is "" it has not been set yet
  if (Particle.GetType() == "") {
    throw std::out_of_range("no particle defined");
  }

  // UPDATE: Dimension is hard coded
  int const Dimension = 2;


  // Calculate trajectory
  this->CalculateTrajectory(Particle);

  // Grab the Trajectory
  TParticleTrajectoryPoints& T = Particle.GetTrajectory();

  // Number of points in trajectory
  size_t const NTPoints = T.GetNPoints();

  // Time step size
  double const DeltaT = T.GetDeltaT();

  // Constand C0
  double const C0 = Particle.GetQ() / (4 * TOSCARSSR::Pi() * TOSCARSSR::C() * TOSCARSSR::Epsilon0() * TOSCARSSR::Sqrt2Pi());


  // Constant for flux calculation at the end
  double const C2 = TOSCARSSR::FourPi() * Particle.GetCurrent() / (TOSCARSSR::H() * fabs(Particle.GetQ()) * TOSCARSSR::Mu0() * TOSCARSSR::C()) * 1e-6 * 0.001;

  // Complex number i
  std::complex<double> const I(0, 1);

  // UPDATE: use common function
  double const Omega = Energy_eV * TOSCARSSR::TwoPi() / 4.1357e-15;

  // Constant C1 (complex)
  //std::complex<double> const C1(0, C0 * Omega);

  // Loop over all surface points
  for (size_t ip = 0; ip != Surface.GetNPoints(); ++ip) {

    // Observation point
    TVector3D const Obs = Surface.GetPoint(ip).GetPoint();

    // Sum E-field
    TVector3DC SumE(0, 0, 0);

    // Loop over trajectory points
    for (int iT = 0; iT != NTPoints; ++iT) {

      // Get position, Beta, and Acceleration (over c)
      TVector3D const& X = T.GetX(iT);
      TVector3D const& B = T.GetB(iT);
      TVector3D const& AoverC = T.GetAoverC(iT);

      // Define R and unit vector in direction of R, and D (distance to observer)
      TVector3D const R = Obs - X;
      TVector3D const N = R.UnitVector();
      double const D = R.Mag();

      // Exponent in transformed field
      std::complex<double> Exponent(0, -Omega * (DeltaT * iT + D / TOSCARSSR::C()));

      // UPDATE: with better field avoiding A
      // TVector3DC const ThisEw = ( N.Cross( (N - B).Cross(AoverC) ) ) / ( D * pow(1 - N.Dot(B), 2) ) * std::exp(Exponent) * DeltaT; // FF only

      TVector3DC const ThisEw = ( ( (1 - (B).Mag2()) * (N - B) ) / ( D * D * (pow(1 - N.Dot(B), 2)) )
          + ( N.Cross( (N - B).Cross(AoverC) ) ) / ( D * pow(1 - N.Dot(B), 2) ) ) * std::exp(Exponent) * DeltaT; // NF + FF

      // Add this contribution
      SumE += ThisEw;

    }

    // Multiply by constant factor
    SumE *= C0;


    double const ThisFlux = C2 * SumE.Dot( SumE.CC() ).real();
  }

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
                              int const Dimension)
{
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
    std::cerr << "NThreads or NThreadsGlobal must be >= 1" << std::endl;
    throw;
  }

  // Should we use the GPU or not?
  bool const UseGPU = GPU == 0 ? false : this->GetUseGPUGlobal() && (this->CheckGPU() > 0) ? true : false;

  if (Dimension == 3) {
    for (int i = 0; i != Surface.GetNPoints(); ++i) {
      FluxContainer.AddPoint(Surface.GetPoint(i).GetPoint(), 0);
    }
  } else if (Dimension == 2) {
    for (int i = 0; i != Surface.GetNPoints(); ++i) {
      FluxContainer.AddPoint( TVector3D(Surface.GetX1(i), Surface.GetX2(i), 0), 0);
    }
  } else {
    std::cerr << "wROng dimension" << std::endl;
    throw;
  }

  // GPU will outrank NThreads...
  if (NParticles == 0) {
    if (UseGPU == 0) {
      this->CalculateFluxThreads(fParticle,
                                 Surface,
                                 Energy_eV,
                                 FluxContainer,
                                 Polarization,
                                 Angle,
                                 HorizontalDirection,
                                 PropogationDirection,
                                 NThreadsToUse,
                                 Dimension,
                                 1);
    } else if (UseGPU == 1) {
      this->CalculateFluxGPU(fParticle,
                             Surface,
                             Energy_eV,
                             FluxContainer,
                             Polarization,
                             Angle,
                             HorizontalDirection,
                             PropogationDirection,
                             Dimension,
                             1);
    }
  } else {

    // Weight this by the number of particles
    double const Weight = 1.0 / (double) NParticles;

    // Loop over particles
    for (int i = 0; i != NParticles; ++i) {

      // Set a new random particle
      this->SetNewParticle();

      // GPU will outrank NThreads...
      if (UseGPU == 0) {
        this->CalculateFluxThreads(fParticle,
                                   Surface,
                                   Energy_eV,
                                   FluxContainer,
                                   Polarization,
                                   Angle,
                                   HorizontalDirection,
                                   PropogationDirection,
                                   NThreadsToUse,
                                   Dimension,
                                   Weight);
      } else if (UseGPU == 1) {
        this->CalculateFluxGPU(fParticle,
                               Surface,
                               Energy_eV,
                               FluxContainer,
                               Polarization,
                               Angle,
                               HorizontalDirection,
                               PropogationDirection,
                               Dimension,
                               Weight);
      }
    }
  }

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
                                     int const Dimension,
                                     double const Weight)
{
  // Calculates the single particle spectrum at a given observation point
  // in units of [photons / second / 0.001% BW / mm^2]
  //
  // Surface - Observation Point

  // Check that particle has been set yet.  If fType is "" it has not been set yet
  if (Particle.GetType() == "") {
    try {
      this->SetNewParticle();
    } catch (std::exception e) {
      throw std::out_of_range("no beam defined");
    }
  }

  // Calculate the trajectory from scratch
  this->CalculateTrajectory(Particle);

  std::vector<std::thread> Threads;

  // Number of points in spectrum
  size_t const NPoints = Surface.GetNPoints();

  // How many threads to start in the first for loop
  size_t const NThreadsActual = (size_t) NPoints > NThreads ? NThreads : NPoints;

  // Keep track of which threads are finished and started
  bool Done[NThreadsActual];
  bool Submitted[NThreadsActual];

  // Number of points for which a thread was started so far
  size_t NPointsStarted = 0;

  // Start threads and keep in vector
  for (size_t it = 0; it != NThreadsActual; ++it) {

    // Set Done to false for this thread
    Done[it] = false;

    // Start thread for this point
    Threads.push_back(std::thread(&OSCARSSR::CalculateFluxPoint,
                                  this,
                                  std::ref(Particle),
                                  std::ref(Surface),
                                  Energy_eV,
                                  std::ref(FluxContainer),
                                  (int) it,
                                  std::ref(Done[it]),
                                  Polarization,
                                  Angle,
                                  HorizontalDirection,
                                  PropogationDirection,
                                  Dimension,
                                  Weight));

    // Incement the number of points started so far
    ++NPointsStarted;

    // Set the submitted flag for this thread to true
    Submitted[it] = true;
  }


  // Keep track of how many threads have completed
  size_t NThreadsFinished = 0;

  // Are all of the threads finished or not?  Continue loop until all come back.
  // UPDATE: Could think about applying a timeout feature here
  bool AllThreadsFinished = false;
  while (!AllThreadsFinished) {

    // So as to not use the current thread at 100%
    std::this_thread::sleep_for(std::chrono::milliseconds(10));

    // Check all threads
    for (size_t it = 0; it != NThreadsActual; ++it) {

      // If it hasn't reported done leave it alone, otherwise see about creating another thread
      if (!Done[it]) {
        continue;
      } else if (Submitted[it] && Done[it]) {

        // Join the thread
        Threads[it].join();

        // Set the submitted bit to false (so it's open)
        Submitted[it] = false;

        // Increment the number of finished threads
        ++NThreadsFinished;

        // If the number of points started is less than the total number start another point
        if (NPointsStarted < NPoints) {

          // Set done bit to false for this thread
          Done[it] = false;

          // Index of the point of interest
          int const ThisPoint = (int) NPointsStarted;

          // Create a new thread for this point
          Threads[it] = std::thread(&OSCARSSR::CalculateFluxPoint,
                                    this,
                                    std::ref(Particle),
                                    std::ref(Surface),
                                    Energy_eV,
                                    std::ref(FluxContainer),
                                    ThisPoint,
                                    std::ref(Done[it]),
                                    Polarization,
                                    Angle,
                                    HorizontalDirection,
                                    PropogationDirection,
                                    Dimension,
                                    Weight);

          // Increment the number of points that have been started
          ++NPointsStarted;

          // Set the submitted bit to true
          Submitted[it] = true;
        }
      }
    }

    // If the number finished is equal to the number of points total then we're done
    if (NThreadsFinished == NPoints) {
      AllThreadsFinished = true;
    }

  }

  // Clear all threads
  Threads.clear();

  return;
}





void OSCARSSR::CalculateFluxGPU (TSurfacePoints const& Surface,
                                 double const Energy_eV,
                                 T3DScalarContainer& FluxContainer,
                                 std::string const& Polarization,
                                 double const Angle,
                                 TVector3D const& HorizontalDirection,
                                 TVector3D const& PropogationDirection,
                                 int const Dimension,
                                 double const Weight)
{
  // If you compile for Cuda use the GPU in this function, else throw

  // Check that particle has been set yet.  If fType is "" it has not been set yet
  if (fParticle.GetType() == "") {
    try {
      this->SetNewParticle();
    } catch (std::exception e) {
      throw std::out_of_range("no beam defined");
    }
  }

  this->CalculateFluxGPU(fParticle, Surface, Energy_eV, FluxContainer, Polarization, Angle, HorizontalDirection, PropogationDirection, Dimension, Weight);

  return;
}





void OSCARSSR::CalculateFluxGPU (TParticleA& Particle,
                                 TSurfacePoints const& Surface,
                                 double const Energy_eV,
                                 T3DScalarContainer& FluxContainer,
                                 std::string const& Polarization,
                                 double const Angle,
                                 TVector3D const& HorizontalDirection,
                                 TVector3D const& PropogationDirection,
                                 int const Dimension,
                                 double const Weight)
{
  // If you compile for Cuda use the GPU in this function, else throw

  // Add points to flux container
  for (size_t i = 0; i != Surface.GetNPoints(); ++i) {
    FluxContainer.AddPoint(Surface.GetPoint(i).GetPoint(), 0);
  }

  // Check that particle has been set yet.  If fType is "" it has not been set yet
  if (Particle.GetType() == "") {
    throw std::out_of_range("no particle defined");
  }

  // Calculate the trajectory from scratch
  this->CalculateTrajectory(Particle);

  #ifdef CUDA
  // Check that the GPU exists
  if (this->CheckGPU() < 1) {
    throw std::invalid_argument("You are requesting the GPU, but none were found");
  }

  return OSCARSSR_Cuda_CalculateFluxGPU(Particle, Surface, Energy_eV, FluxContainer, Polarization, Angle, HorizontalDirection, PropogationDirection, Dimension, Weight);
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

  // Calculate Trajectory
  this->CalculateTrajectory(Particle);

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
