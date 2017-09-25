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


#include "TOSCARSSR.h"

#include <string>

#include "OSCARSSR_Cuda.h"
#include "TFieldContainer.h"
#include "TParticleBeamContainer.h"
#include "TDriftVolumeContainer.h"
#include "TSurfacePoints.h"
#include "TSpectrumContainer.h"
#include "T3DScalarContainer.h"
#include "TParticleTrajectoryInterpolated.h"
#include "TRandomA.h"


class OSCARSSR
{
  // This class is meant to be the main interface to the simulation,
  // also, from all extensions

  public:
    OSCARSSR ();
    ~OSCARSSR ();

    // Functions related to the magnetic field
    void AddMagneticField (std::string const FileName,
                           std::string const Format,
                           TVector3D const& R = TVector3D(0, 0, 0),
                           TVector3D const& D = TVector3D(0, 0, 0),
                           std::vector<double> const& S = std::vector<double>(),
                           std::string const& Name = "");

    void AddMagneticFieldInterpolated (std::vector<std::pair<double, std::string> > const& Mapping,
                                       std::string const Format,
                                       double const Parameter,
                                       TVector3D const& Rotations = TVector3D(0, 0, 0),
                                       TVector3D const& Translation = TVector3D(0, 0, 0),
                                       std::vector<double> const& Scaling = std::vector<double>(),
                                       std::string const& Name = "");

    void AddMagneticField (TField*);

    void RemoveMagneticField (std::string const& Name);

    void ClearMagneticFields ();

    void AddElectricField (std::string const FileName,
                           std::string const Format,
                           TVector3D const& Rotations = TVector3D(0, 0, 0),
                           TVector3D const& Translation = TVector3D(0, 0, 0),
                           std::vector<double> const& Scaling = std::vector<double>(),
                           std::string const& Name = "");

    void AddElectricFieldInterpolated (std::vector<std::pair<double, std::string> > const& Mapping,
                                       std::string const Format,
                                       double const Parameter,
                                       TVector3D const& Rotations = TVector3D(0, 0, 0),
                                       TVector3D const& Translation = TVector3D(0, 0, 0),
                                       std::vector<double> const& Scaling = std::vector<double>(),
                                       std::string const& Name = "");

    void AddElectricField (TField*);

    void RemoveElectricField (std::string const& Name);

    void ClearElectricFields ();

    void WriteField (std::string const& BorE,
                     std::string const& OutFileName,
                     std::string const& OutFormat,
                     TVector2D const& XLim,
                     int const NX,
                     TVector2D const& YLim,
                     int const NY,
                     TVector2D const& ZLim,
                     int const NZ,
                     std::string const& Comment);

    void WriteFieldBinary (std::string const& BorE,
                           std::string const& OutFileName,
                           std::string const& OutFormat,
                           TVector2D const& XLim,
                           int const NX,
                           TVector2D const& YLim,
                           int const NY,
                           TVector2D const& ZLim,
                           int const NZ,
                           std::string const& Comment,
                           int const Version);


    TVector3D GetB  (double const, double const, double const) const;
    TVector3D GetB  (TVector3D const&) const;

    TVector3D GetE  (double const, double const, double const) const;
    TVector3D GetE  (TVector3D const&) const;


    // Functions related to the particle beam(s)
    TParticleBeam& AddParticleBeam (std::string const& Type,
                                    std::string const& Name,
                                    TVector3D const& X0,
                                    TVector3D const& V0,
                                    double const Energy_GeV,
                                    double const T0,
                                    double const Current,
                                    double const Weight,
                                    double const Charge = 0,
                                    double const Mass = 0);

    TParticleBeam& AddParticleBeam (std::string const& Beam,
                                    std::string const& Name,
                                    double const Weight = 1);

    TParticleBeamContainer& GetParticleBeamContainer () 
    {
      return fParticleBeamContainer;
    }

    TParticleBeam& GetParticleBeam (std::string const&);
    size_t GetNParticleBeams () const;
    TParticleA GetNewParticle ();
    TParticleA const&  GetCurrentParticle () const;
    void SetNewParticle ();
    void SetNewParticle (std::string const&, std::string const&);
    void ClearParticleBeams ();

    void SetEmittance (std::string const& Beam,
                       TVector2D const& Emittance);

    void SetTwissParameters (std::string const& Beam,
                             TVector2D const& Beta,
                             TVector2D const& Alpha,
                             TVector2D const& Gamma,
                             TVector3D const& Lattice_Reference,
                             bool const HasReferencePoint);


    // Functions related to drift spaces
    void AddDriftVolume (TDriftVolume*);
    void RemoveDriftVolume (std::string const& Name);
    void ClearDriftVolumes ();

    // Functions related to Trajectory
    void CorrectTrajectory ();
    void CalculateTrajectory ();
    void CalculateTrajectory (TParticleA&);
    TParticleTrajectoryPoints const& GetTrajectory ();
    void ClearTrajectory ();

    void SetNPointsTrajectory (size_t const);
    void SetNPointsPerMeterTrajectory (size_t const);
    void SetCTStartStop (double const, double const);

    size_t GetNPointsTrajectory () const;
    double GetCTStart () const;
    double GetCTStop  () const;

    // Global threads and GPU settings
    bool SetUseGPUGlobal (int const);
    int  GetUseGPUGlobal () const;
    int  CheckGPU () const;
    std::string GetGPUInfo (int const) const;
    void SetNThreadsGlobal (int const);
    int  GetNThreadsGlobal () const;

    // Random seed setting and random numbers
    void SetSeed (int const) const;
    double GetRandomNormal () const;
    double GetRandomUniform () const;

    void CalculateSpectrum (TParticleA& Particle,
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
                            int    const ReturnQuantity);

    void CalculateSpectrum (TVector3D const& ObservationPoint,
                            TSpectrumContainer& Spectrum,
                            std::string const& Polarization = "all",
                            double const Angle = 0,
                            TVector3D const& HorizontalDirection = TVector3D(0, 0, 0),
                            TVector3D const& PropogationDirection = TVector3D(0, 0, 0),
                            int const NParticles = 0,
                            int const NThreads = 0,
                            int const GPU = 0,
                            int const NGPU = -1,
                            std::vector<int> VGPU = std::vector<int>(),
                            double const Precision = 0.01,
                            int    const MaxLevel = -2,
                            int    const MaxLevelExtended = 0,
                            int    const ReturnQuantity = 0);

    void CalculateSpectrumPoints (TParticleA& Particle,
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
                                  int    const ReturnQuantity);

    void CalculateSpectrumThreads (TParticleA& Particle,
                                   TVector3D const& Obs,
                                   TSpectrumContainer& Spectrum,
                                   int const NThreads,
                                   std::string const& Polarization = "all",
                                   double const Angle = 0,
                                   TVector3D const& HorizontalDirection = TVector3D(0, 0, 0),
                                   TVector3D const& PropogationDirection = TVector3D(0, 0, 0),
                                   double const Precision = 0.01,
                                   int    const MaxLevel = -2,
                                   int    const MaxLevelExtended = 0,
                                   double const Weight = 1,
                                   int    const ReturnQuantity = 0);

    void CalculateSpectrumGPU (TParticleA& Particle,
                               TVector3D const& ObservationPoint,
                               TSpectrumContainer& Spectrum,
                               std::string const& Polarization = "all",
                               double const Angle = 0,
                               TVector3D const& HorizontalDirection = TVector3D(0, 0, 0),
                               TVector3D const& PropogationDirection = TVector3D(0, 0, 0),
                               int const NParticles = 0,
                               std::vector<int> GPUVector = std::vector<int>(),
                               double const Precision = 0.01,
                               int    const MaxLevel = -2,
                               int    const MaxLevelExtended = 0,
                               int    const ReturnQuantity = 0);

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
    void CalculatePowerDensity (TParticleA& Particle,
                                TSurfacePoints const& Surface,
                                T3DScalarContainer& PowerDensityContainer,
                                bool const Directional,
                                double const Precision,
                                int    const MaxLevel,
                                int    const MaxLevelExtended,
                                double const Weight,
                                int    const ReturnQuantity);

    void CalculatePowerDensity (TSurfacePoints const& Surface,
                                T3DScalarContainer& PowerDensityContainer,
                                int const Dimension,
                                bool const Directional,
                                double const Precision,
                                int    const MaxLevel,
                                int    const MaxLevelExtended,
                                int const NParticles,
                                int const NThreads,
                                int const GPU,
                                int const NGPU = -1,
                                std::vector<int> VGPU = std::vector<int>(),
                                int const ReturnQuantity = 0);

    void CalculatePowerDensityPoints (TParticleA& Particle,
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
                                      int    const ReturnQuantity);

    void CalculatePowerDensityThreads (TParticleA& Particle,
                                       TSurfacePoints const& Surface,
                                       T3DScalarContainer& PowerDensityContainer,
                                       int const NThreads,
                                       bool const Directional,
                                       double const Precision,
                                       int    const MaxLevel,
                                       int    const MaxLevelExtended,
                                       double const Weight,
                                       int    const ReturnQuantity);

    void CalculatePowerDensityGPU (TSurfacePoints const& Surface,
                                   T3DScalarContainer& PowerDensityContainer,
                                   int const NParticles,
                                   std::vector<int> GPUVector,
                                   bool const Directional,
                                   double const Precision,
                                   int    const MaxLevel,
                                   int    const MaxLevelExtended,
                                   int    const ReturnQuantity);

    double CalculateTotalPower ();
    double CalculateTotalPower (TParticleA&);


    // Flux Calculations
    void CalculateFlux (TParticleA& Particle,
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
                        int    const ReturnQuantity);

    void CalculateFlux (TSurfacePoints const& Surface,
                        double const Energy_eV,
                        T3DScalarContainer& FluxContainer,
                        std::string const& Polarization = "all",
                        double const Angle = 0,
                        TVector3D const& HorizontalDirection = TVector3D(0, 0, 0),
                        TVector3D const& PropogationDirection = TVector3D(0, 0, 0),
                        int const NParticles = 0,
                        int const NThreads = 0,
                        int const GPU = 0,
                        int const NGPU = -1,
                        std::vector<int> VGPU = std::vector<int>(),
                        double const Precision = 0.01,
                        int    const MaxLevel = -2,
                        int    const MaxLevelExtended = 0,
                        int    const Dimension = 3,
                        int    const ReturnQuantity = 0);

    void CalculateFluxPoints (TParticleA& Particle,
                              TSurfacePoints const& Surface,
                              double const Energy_eV,
                              T3DScalarContainer& FluxContainer,
                              size_t const iFirst,
                              size_t const iLast,
                              bool& Done,
                              std::string const& Polarization = "all",
                              double const Angle = 0,
                              TVector3D const& HorizontalDirection = TVector3D(0, 0, 0),
                              TVector3D const& PropogationDirection = TVector3D(0, 0, 0),
                              double const Precision = 0.01,
                              int    const MaxLevel = -2,
                              int    const MaxLevelExtended = 0,
                              double const Weight = 1,
                              int    const ReturnQuantity = 0);

    void CalculateFluxThreads (TParticleA& Particle,
                               TSurfacePoints const& Surface,
                               double const Energy_eV,
                               T3DScalarContainer& FluxContainer,
                               std::string const& Polarization = "all",
                               double const Angle = 0,
                               TVector3D const& HorizontalDirection = TVector3D(0, 0, 0),
                               TVector3D const& PropogationDirection = TVector3D(0, 0, 0),
                               int const NThreads = 0,
                               double const Precision = 0.01,
                               int    const MaxLevel = -2,
                               int    const MaxLevelExtended = 0,
                               double const Weight = 1,
                               int    const ReturnQuantity = 0);

    void CalculateFluxGPU (TSurfacePoints const& Surface,
                           double const Energy_eV,
                           T3DScalarContainer& FluxContainer,
                           std::string const& Polarization = "all",
                           double const Angle = 0,
                           TVector3D const& HorizontalDirection = TVector3D(0, 0, 0),
                           TVector3D const& PropogationDirection = TVector3D(0, 0, 0),
                           int const NParticles = 0,
                           std::vector<int> GPUVector = std::vector<int>(),
                           double const Precision = 0.01,
                           int    const MaxLevel = -2,
                           int    const MaxLevelExtended = 0,
                           int    const ReturnQuantity = 0);


    void CalculateFluxGPUNew (TSurfacePoints const& Surface,
                           double const Energy_eV,
                           T3DScalarContainer& FluxContainer,
                           std::string const& Polarization = "all",
                           double const Angle = 0,
                           TVector3D const& HorizontalDirection = TVector3D(0, 0, 0),
                           TVector3D const& PropogationDirection = TVector3D(0, 0, 0),
                           int const NParticles = 0,
                           std::vector<int> GPUVector = std::vector<int>(),
                           double const Precision = 0.01,
                           int    const MaxLevel = -2,
                           int    const MaxLevelExtended = 0,
                           int    const ReturnQuantity = 0);

    // Electric Field Calculations
    void CalculateElectricFieldTimeDomain (TVector3D const& Observer, T3DScalarContainer&);
    void CalculateElectricFieldTimeDomain (TVector3D const& Observer, T3DScalarContainer&, TParticleA& Particle);

    TFieldContainer const& GetBFieldContainer () const;
    TFieldContainer const& GetEFieldContainer () const;
    TDriftVolumeContainer const& GetDriftVolumeContainer () const;

  private:
    TFieldContainer  fBFieldContainer;
    TFieldContainer  fEFieldContainer;

    TParticleBeamContainer fParticleBeamContainer;

    TDriftVolumeContainer fDriftVolumeContainer;

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
