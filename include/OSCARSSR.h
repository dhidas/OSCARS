#ifndef GUARD_OSCARSSR_h
#define GUARD_OSCARSSR_h
////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Wed May 18 20:34:55 EDT 2016
//
// Class to contain all elements and functions for radiation
// simulation.  THIS is the c++ API
//
////////////////////////////////////////////////////////////////////

#define OSCARSSR_VMAJOR 1
#define OSCARSSR_VMINOR 34
#define OSCARSSR_REVISION 0

#include "TOSCARSSR.h"

#include <string>

#include "OSCARSSR_Cuda.h"
#include "TFieldContainer.h"
#include "TParticleBeamContainer.h"
#include "TSurfacePoints.h"
#include "TSpectrumContainer.h"
#include "T3DScalarContainer.h"
#include "TRandomA.h"


class OSCARSSR
{
  // This class is meant to be the main interface to the simulation,
  // also, from all extensions

  public:
    OSCARSSR ();
    ~OSCARSSR ();

    static std::string GetVersionString ();

    // Functions related to the magnetic field
    void AddMagneticField (std::string const FileName,
                           std::string const Format,
                           TVector3D const& R = TVector3D(0, 0, 0),
                           TVector3D const& D = TVector3D(0, 0, 0),
                           std::vector<double> const& S = std::vector<double>());

    void AddMagneticFieldInterpolated (std::vector<std::pair<double, std::string> > const& Mapping,
                                       std::string const Format,
                                       double const Parameter,
                                       TVector3D const& Rotations = TVector3D(0, 0, 0),
                                       TVector3D const& Translation = TVector3D(0, 0, 0),
                                       std::vector<double> const& Scaling = std::vector<double>());
    void AddMagneticField (TField*);
    void ClearMagneticFields ();

    void AddElectricField (std::string const, std::string const, TVector3D const& R = TVector3D(0, 0, 0), TVector3D const& D = TVector3D(0, 0, 0), std::vector<double> const& S = std::vector<double>());
    void AddElectricField (TField*);
    void ClearElectricFields ();

    void WriteField (std::string const&, std::string const&, std::string const&, TVector2D const&, int const, TVector2D const&, int const, TVector2D const&, int const, std::string const& Comment = "");


    double    GetBx (double const, double const, double const) const;
    double    GetBy (double const, double const, double const) const;
    double    GetBz (double const, double const, double const) const;
    TVector3D GetB  (double const, double const, double const) const;
    TVector3D GetB  (TVector3D const&) const;

    double    GetEx (double const, double const, double const) const;
    double    GetEy (double const, double const, double const) const;
    double    GetEz (double const, double const, double const) const;
    TVector3D GetE  (double const, double const, double const) const;
    TVector3D GetE  (TVector3D const&) const;


    // Functions related to the particle beam(s)
    void AddParticleBeam (std::string const&, std::string const&, TVector3D const&, TVector3D const&, double const, double const, double const, double const, double const Charge = 0, double const Mass = 0);
    TParticleBeam& GetParticleBeam (std::string const&);
    size_t GetNParticleBeams () const;
    TParticleA GetNewParticle ();
    TParticleA const&  GetCurrentParticle () const;
    void SetNewParticle ();
    void SetNewParticle (std::string const&, std::string const&);
    void ClearParticleBeams ();


    // Functions related to Trajectory
    void CalculateTrajectory ();
    void CalculateTrajectory (TParticleA&);
    TParticleTrajectoryPoints const& GetTrajectory ();

    void SetNPointsTrajectory (size_t const);
    void SetNPointsPerMeter (size_t const);
    void SetCTStartStop (double const, double const);

    size_t GetNPointsTrajectory () const;
    double GetCTStart () const;
    double GetCTStop  () const;

    // Global threads and GPU settings
    void SetUseGPUGlobal (int const);
    int  CheckGPU () const;
    void SetNThreadsGlobal (int const);

    // Random seed setting and random numbers
    void SetSeed (int const) const;
    double GetRandomNormal () const;
    double GetRandomUniform () const;

    void CalculateSpectrumGPU (TParticleA&,
                               TVector3D const& ObservationPoint,
                               TSpectrumContainer& Spectrum,
                               std::string const& Polarization = "all",
                               double const Angle = 0,
                               TVector3D const& HorizontalDirection = TVector3D(0, 0, 0),
                               TVector3D const& PropogationDirection = TVector3D(0, 0, 0),
                               double const Weight = 1,
                               std::string const OutFileName = "");

    void CalculateSpectrum (TVector3D const& ObservationPoint,
                            TSpectrumContainer& Spectrum,
                            std::string const& Polarization = "all",
                            double const Angle = 0,
                            TVector3D const& HorizontalDirection = TVector3D(0, 0, 0),
                            TVector3D const& PropogationDirection = TVector3D(0, 0, 0),
                            int const NParticles = 0,
                            int const NThreads = 0,
                            int const GPU = 0);

    void CalculateSpectrumThreads (TParticleA& Particle,
                                   TVector3D const& Obs,
                                   TSpectrumContainer& Spectrum,
                                   int const NThreads,
                                   std::string const& Polarization = "all",
                                   double const Angle = 0,
                                   TVector3D const& HorizontalDirection = TVector3D(0, 0, 0),
                                   TVector3D const& PropogationDirection = TVector3D(0, 0, 0),
                                   double const Weight = 1,
                                   std::string const& OutFileName = "");

    void CalculateSpectrumPoint (TParticleA& ,
                                 TVector3D const&,
                                 TSpectrumContainer&,
                                 int const i,
                                 std::string const& Polarization = "all",
                                 double const Angle = 0,
                                 TVector3D const& HorizontalDirection = TVector3D(0, 0, 0),
                                 TVector3D const& PropogationDirection = TVector3D(0, 0, 0),
                                 double const Weight = 1);

    void CalculateSpectrum (TParticleA&,
                            TVector3D const&,
                            TSpectrumContainer&,
                            std::string const& Polarization = "all",
                            double const Angle = 0,
                            TVector3D const& HorizontalDirection = TVector3D(0, 0, 0),
                            TVector3D const& PropogationDirection = TVector3D(0, 0, 0),
                            double const Weight = 1);


    void AddToSpectrum (TSpectrumContainer const&, double const Weight = 1);
    void AddToFlux (T3DScalarContainer const&, double const Weight = 1);
    void AddToPowerDensity (T3DScalarContainer const&, double const Weight = 1);

    TSpectrumContainer const& GetSpectrum () const;
    void ClearSpectrum ();
    T3DScalarContainer const& GetFlux () const;
    void ClearFlux ();
    T3DScalarContainer const& GetPowerDensity () const;
    void ClearPowerDensity ();

    // Power Density calculation
    void CalculatePowerDensity (TParticleA&, TSurfacePoints const&, int const Dimension = 3, bool const Directional = true, double const Weight = 1, std::string const& OutFileName = "");
    void CalculatePowerDensity (TParticleA&, TSurfacePoints const&, T3DScalarContainer&, int const Dimension = 3, bool const Directional = true, double const Weight = 1, std::string const& OutFileName = "");
    void CalculatePowerDensity (TSurfacePoints const&, T3DScalarContainer&, int const Dimension = 3, bool const Directional = true, int const NParticles = 0, std::string const& OutFileName = "", int const NThreads = 0, int const GPU = 0);
    void CalculatePowerDensityGPU (TParticleA&, TSurfacePoints const&, T3DScalarContainer&, int const Dimension = 3, bool const Directional = true, double const Weight = 1, std::string const& OutFileName = "");
    void CalculatePowerDensityGPU (TSurfacePoints const&, T3DScalarContainer&, int const Dimension = 3, bool const Directional = true, double const Weight = 1, std::string const& OutFileName = "");
    void CalculatePowerDensityThreads (TParticleA&, TSurfacePoints const&, T3DScalarContainer&, int const, int const Dimension = 3, bool const Directional = true, double const Weight = 1, std::string const& OutFileName = "");
    void CalculatePowerDensityPoint (TParticleA&, TSurfacePoints const&, T3DScalarContainer&, size_t const, int const Dimension = 3, bool const Directional = true, double const Weight = 1);
    double CalculateTotalPower ();
    double CalculateTotalPower (TParticleA&);


    // Flux Calculations
    void CalculateFlux2   (TParticleA& Particle,
                           TSurfacePoints const& Surface,
                           double const Energy_eV,
                           T3DScalarContainer& FluxContainer,
                           std::string const& Polarization = "all",
                           double const Angle = 0,
                           TVector3D const& HorizontalDirection = TVector3D(0, 0, 0),
                           TVector3D const& PropogationDirection = TVector3D(0, 0, 0),
                           int const Dimension = 3,
                           double const Weight = 1);

    void CalculateFluxPoint (TParticleA& Particle,
                             TSurfacePoints const& Surface,
                             double const Energy_eV,
                             T3DScalarContainer& FluxContainer,
                             size_t const i,
                             std::string const& Polarization = "all",
                             double const Angle = 0,
                             TVector3D const& HorizontalDirection = TVector3D(0, 0, 0),
                             TVector3D const& PropogationDirection = TVector3D(0, 0, 0),
                             int const Dimension = 3,
                             double const Weight = 1);

    void CalculateFlux    (TParticleA&,
                           TSurfacePoints const&,
                           double const,
                           T3DScalarContainer&,
                           std::string const& Polarization = "all",
                           double const Angle = 0,
                           TVector3D const& HorizontalDirection = TVector3D(0, 0, 0),
                           TVector3D const& PropogationDirection = TVector3D(0, 0, 0),
                           int const Dimension = 3,
                           double const Weight = 1,
                           std::string const& OutFileName = "");

    void CalculateFlux1   (TParticleA&, TSurfacePoints const&, double const, T3DScalarContainer&, std::string const& OutFileName = "");
    void CalculateFlux    (TSurfacePoints const& Surface,
                           double const Energy_eV,
                           T3DScalarContainer& FluxContainer,
                           std::string const& Polarization = "all",
                           double const Angle = 0,
                           TVector3D const& HorizontalDirection = TVector3D(0, 0, 0),
                           TVector3D const& PropogationDirection = TVector3D(0, 0, 0),
                           int const NParticles = 0,
                           int const NThreads = 0,
                           int const GPU = 0,
                           int const Dimension = 3,
                           std::string const& OutFileName = "");

    void CalculateFluxThreads (TParticleA& Particle,
                               TSurfacePoints const& Surface,
                               double const Energy_eV,
                               T3DScalarContainer& FluxContainer,
                               std::string const& Polarization = "all",
                               double const Angle = 0,
                               TVector3D const& HorizontalDirection = TVector3D(0, 0, 0),
                               TVector3D const& PropogationDirection = TVector3D(0, 0, 0),
                               int const NThreads = 0,
                               int const Dimension = 3,
                               double const Weight = 1,
                               std::string const& OutFileName = "");

    void CalculateFluxGPU (TParticleA& Particle,
                           TSurfacePoints const& Surface,
                           double const Energy_eV,
                           T3DScalarContainer& FluxContainer,
                           std::string const& Polarization = "all",
                           double const Angle = 0,
                           TVector3D const& HorizontalDirection = TVector3D(0, 0, 0),
                           TVector3D const& PropogationDirection = TVector3D(0, 0, 0),
                           int const Dimension = 3,
                           double const Weight = 1,
                           std::string const& OutFileName = "");

    void CalculateFluxGPU (TSurfacePoints const&,
                           double const,
                           T3DScalarContainer&,
                           std::string const& Polarization = "all",
                           double const Angle = 0,
                           TVector3D const& HorizontalDirection = TVector3D(0, 0, 0),
                           TVector3D const& PropogationDirection = TVector3D(0, 0, 0),
                           int const Dimension = 3,
                           double const Weight = 1,
                           std::string const& OutFileName = "");

    // Electric Field Calculations
    void CalculateElectricFieldTimeDomain (TVector3D const& Observer, T3DScalarContainer&);
    void CalculateElectricFieldTimeDomain (TVector3D const& Observer, T3DScalarContainer&, TParticleA& Particle);

  private:
    TFieldContainer  fBFieldContainer;
    TFieldContainer  fEFieldContainer;

    TParticleBeamContainer fParticleBeamContainer;

    void SetDerivativesFunction ();

    void DerivativesE (double t, double x[], double dxdt[], TParticleA const&);
    void DerivativesB (double t, double x[], double dxdt[], TParticleA const&);
    void DerivativesEB (double t, double x[], double dxdt[], TParticleA const&);
    void Derivatives (double t, double x[], double dxdt[], TParticleA const&);
    void RK4 (double y[], double dydx[], int n, double x, double h, double yout[], TParticleA const&);


    double fCTStart;
    double fCTStop;
    size_t fNPointsTrajectory;
    size_t fNPointsPerMeter;


    // Current particle for calculations and rel parameters
    TParticleA fParticle;
    double fCurrent;

    // Spectrum container
    TSpectrumContainer fSpectrum;
    T3DScalarContainer fFlux;
    T3DScalarContainer fPowerDensity;

    // Global thread and GPU settings
    int fNThreadsGlobal;
    bool fUseGPUGlobal;

    // Function pointer for which function to use in the RK4 propogation
    void (OSCARSSR::*fDerivativesFunction)(double, double*, double*, TParticleA const&);

};








#endif
