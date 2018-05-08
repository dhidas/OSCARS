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

class OSCARSSR;


extern "C" int  OSCARSSR_Cuda_GetDeviceCount ();

std::string OSCARSSR_Cuda_GetDeviceProperties (int const);

extern "C" void OSCARSSR_Cuda_CalculateFluxGPU (OSCARSSR& OSR,
                                                TSurfacePoints const& Surface,
                                                double const Energy_eV,
                                                T3DScalarContainer& FluxContainer,
                                                std::string const& PolarizationIn = "all",
                                                double const Angle = 0,
                                                TVector3D const& HorizontalDirection = TVector3D(0, 0, 0),
                                                TVector3D const& PropogationDirection = TVector3D(0, 0, 0),
                                                int const NParticles = 0,
                                                std::vector<int> const& GPUVector = std::vector<int>(),
                                                double const Precision = 0.01,
                                                int const MaxLevel = 25,
                                                int    const ReturnQuantity = 0);

extern "C" void OSCARSSR_Cuda_CalculateSpectrumGPU (OSCARSSR& OSR,
                                                    TVector3D const& ObservationPoint,
                                                    TSpectrumContainer& Spectrum,
                                                    std::string const& PolarizationIn,
                                                    double const Angle,
                                                    TVector3D const& HorizontalDirection,
                                                    TVector3D const& PropogationDirection,
                                                    int const NParticles,
                                                    std::vector<int> const& GPUVector,
                                                    double const Precision,
                                                    int const MaxLevel,
                                                    int const ReturnQuantity);

extern "C" void OSCARSSR_Cuda_CalculatePowerDensityGPU (OSCARSSR& OSR,
                                                        TSurfacePoints const& Surface,
                                                        T3DScalarContainer& PowerDensityContainer,
                                                        int const NParticles,
                                                        std::vector<int> const& GPUVector,
                                                        double const Precision,
                                                        int const MaxLevel,
                                                        int const ReturnQuantity = 0);





#endif
