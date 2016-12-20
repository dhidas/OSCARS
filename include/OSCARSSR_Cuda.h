#ifndef GUARD_OSCARSSR_Cuda_h
#define GUARD_OSCARSSR_Cuda_h
////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Wed Aug 17 08:29:54 EDT 2016
//
////////////////////////////////////////////////////////////////////

#include "OSCARSSR.h"

#include <cmath>
#include <complex>
#include <fstream>
#include <sstream>
#include <string>

#include "TVector3DC.h"
#include "TSpectrumContainer.h"
#include "TParticleA.h"
#include "TSurfacePoints.h"
#include "T3DScalarContainer.h"

extern "C" int  OSCARSSR_Cuda_GetDeviceCount ();

extern "C" void OSCARSSR_Cuda_CalculateFluxGPU (TParticleA& Particle,
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

extern "C" void OSCARSSR_Cuda_CalculateSpectrumGPU (TParticleA& Particle,
                                                    TVector3D const& ObservationPoint,
                                                    TSpectrumContainer& Spectrum,
                                                    std::string const& Polarization = "all",
                                                    double const Angle = 0,
                                                    TVector3D const& HorizontalDirection = TVector3D(0, 0, 0),
                                                    TVector3D const& PropogationDirection = TVector3D(0, 0, 0),
                                                    double const Weight = 1);

extern "C" void OSCARSSR_Cuda_CalculatePowerDensityGPU (TParticleA& Particle, TSurfacePoints const& Surface, T3DScalarContainer& PowerDensityContainer, int const Dimension, bool const Directional, double const Weight, std::string const& OutFileName = "");




















#endif
