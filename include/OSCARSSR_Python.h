#ifndef GUARD_OSCARSSR_Python_h
#define GUARD_OSCARSSR_Python_h
////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Tue Jul 12 08:36:18 EDT 2016
//
// This is a header file for the OSCARSSR_Python 'OSCARSSR' module
//
////////////////////////////////////////////////////////////////////

// Include Python.h first!
#include <Python.h>

#include "OSCARSSR.h"

// The python OSCARSSR object
typedef struct {
  // Define the OSCARSSRObject struct which contains the class I want
  PyObject_HEAD
  OSCARSSR* obj;
} OSCARSSRObject;




static void OSCARSSR_dealloc(OSCARSSRObject* self);
static PyObject* OSCARSSR_new (PyTypeObject* type, PyObject* args, PyObject* kwds);
static PyObject* OSCARSSR_Pi (OSCARSSRObject* self, PyObject* arg);
static PyObject* OSCARSSR_Qe (OSCARSSRObject* self, PyObject* arg);
static PyObject* OSCARSSR_Me (OSCARSSRObject* self, PyObject* arg);
static PyObject* OSCARSSR_C (OSCARSSRObject* self, PyObject* arg);
static PyObject* OSCARSSR_Random (OSCARSSRObject* self, PyObject* arg);
static PyObject* OSCARSSR_RandomNormal (OSCARSSRObject* self, PyObject* arg);
static PyObject* OSCARSSR_SetSeed (OSCARSSRObject* self, PyObject* arg);
static PyObject* OSCARSSR_SetGPUGlobal (OSCARSSRObject* self, PyObject* arg);
static PyObject* OSCARSSR_CheckGPU (OSCARSSRObject* self, PyObject* arg);
static PyObject* OSCARSSR_SetNThreadsGlobal (OSCARSSRObject* self, PyObject* arg);
static PyObject* OSCARSSR_MandelbrotSet (OSCARSSRObject* self, PyObject* arg);
static PyObject* OSCARSSR_GetCTStart (OSCARSSRObject* self);
static PyObject* OSCARSSR_GetCTStop (OSCARSSRObject* self);
static PyObject* OSCARSSR_SetCTStartStop (OSCARSSRObject* self, PyObject* args);
static PyObject* OSCARSSR_GetNPointsTrajectory (OSCARSSRObject* self);
static PyObject* OSCARSSR_SetNPointsTrajectory (OSCARSSRObject* self, PyObject* arg);
static PyObject* OSCARSSR_SetNPointsPerMeterTrajectory (OSCARSSRObject* self, PyObject* arg);
static PyObject* OSCARSSR_AddMagneticField (OSCARSSRObject* self, PyObject* args, PyObject* keywds);
static PyObject* OSCARSSR_AddMagneticFieldInterpolated (OSCARSSRObject* self, PyObject* args, PyObject* keywds);
static PyObject* OSCARSSR_AddMagneticFieldFunction (OSCARSSRObject* self, PyObject* args, PyObject* keywds);
static PyObject* OSCARSSR_AddMagneticFieldGaussian (OSCARSSRObject* self, PyObject* args, PyObject* keywds);
static PyObject* OSCARSSR_AddMagneticFieldUniform (OSCARSSRObject* self, PyObject* args, PyObject* keywds);
static PyObject* OSCARSSR_AddMagneticFieldIdealUndulator (OSCARSSRObject* self, PyObject* args, PyObject* keywds);
static PyObject* OSCARSSR_AddMagneticFieldQuadrupole (OSCARSSRObject* self, PyObject* args, PyObject* keywds);
static PyObject* OSCARSSR_RemoveMagneticField (OSCARSSRObject* self, PyObject* args, PyObject* keywds);
static PyObject* OSCARSSR_GetBField (OSCARSSRObject* self, PyObject* args, PyObject* keywds);
static PyObject* OSCARSSR_ClearMagneticFields (OSCARSSRObject* self);
static PyObject* OSCARSSR_PrintMagneticFields (OSCARSSRObject* self);
static PyObject* OSCARSSR_AddElectricField (OSCARSSRObject* self, PyObject* args, PyObject* keywds);
static PyObject* OSCARSSR_AddElectricFieldFunction (OSCARSSRObject* self, PyObject* args, PyObject* keywds);
static PyObject* OSCARSSR_AddElectricFieldGaussian (OSCARSSRObject* self, PyObject* args, PyObject* keywds);
static PyObject* OSCARSSR_AddElectricFieldUniform (OSCARSSRObject* self, PyObject* args, PyObject* keywds);
static PyObject* OSCARSSR_AddElectricFieldIdealUndulator (OSCARSSRObject* self, PyObject* args, PyObject* keywds);
static PyObject* OSCARSSR_GetEField (OSCARSSRObject* self, PyObject* args, PyObject* keywds);
static PyObject* OSCARSSR_ClearElectricFields (OSCARSSRObject* self);
static PyObject* OSCARSSR_PrintElectricFields (OSCARSSRObject* self);
static PyObject* OSCARSSR_WriteMagneticField (OSCARSSRObject* self, PyObject* args, PyObject* keywds);
static PyObject* OSCARSSR_WriteElectricField (OSCARSSRObject* self, PyObject* args, PyObject* keywds);
static PyObject* OSCARSSR_SetParticleBeam (OSCARSSRObject* self, PyObject* args, PyObject* keywds);
static PyObject* OSCARSSR_AddParticleBeam (OSCARSSRObject* self, PyObject* args, PyObject* keywds);
static PyObject* OSCARSSR_SetParticleBeamSize (OSCARSSRObject* self, PyObject* args, PyObject* keywds);
static PyObject* OSCARSSR_ClearParticleBeams (OSCARSSRObject* self);
static PyObject* OSCARSSR_PrintParticleBeams (OSCARSSRObject* self);
static PyObject* OSCARSSR_SetTwissParameters (OSCARSSRObject* self, PyObject* args, PyObject* keywds);
static PyObject* OSCARSSR_GetEmittance (OSCARSSRObject* self, PyObject* args, PyObject* keywds);
static PyObject* OSCARSSR_GetTwissBetaX0 (OSCARSSRObject* self, PyObject* args, PyObject* keywds);
static PyObject* OSCARSSR_GetTwissAlphaX0 (OSCARSSRObject* self, PyObject* args, PyObject* keywds);
static PyObject* OSCARSSR_GetTwissGammaX0 (OSCARSSRObject* self, PyObject* args, PyObject* keywds);
static PyObject* OSCARSSR_SetNewParticle (OSCARSSRObject* self, PyObject* args, PyObject* keywds);
static PyObject* OSCARSSR_GetBeamX0 (OSCARSSRObject* self, PyObject* args, PyObject* keywds);
static PyObject* OSCARSSR_GetBeamU0 (OSCARSSRObject* self, PyObject* args, PyObject* keywds);
static PyObject* OSCARSSR_GetBeamHorizontalDirection (OSCARSSRObject* self, PyObject* args, PyObject* keywds);
static PyObject* OSCARSSR_GetBeamVerticalDirection (OSCARSSRObject* self, PyObject* args, PyObject* keywds);
static PyObject* OSCARSSR_GetParticleX0 (OSCARSSRObject* self);
static PyObject* OSCARSSR_GetParticleBeta0 (OSCARSSRObject* self);
static PyObject* OSCARSSR_GetParticleE0 (OSCARSSRObject* self);
static PyObject* OSCARSSR_AddDriftVolume_Box (OSCARSSRObject* self, PyObject* args, PyObject* keywds);
static PyObject* OSCARSSR_RemoveDriftVolume (OSCARSSRObject* self, PyObject* args, PyObject* keywds);
static PyObject* OSCARSSR_ClearDriftVolumes (OSCARSSRObject* self);
static PyObject* OSCARSSR_PrintDriftVolumes (OSCARSSRObject* self);
static PyObject* OSCARSSR_CalculateTrajectory (OSCARSSRObject* self, PyObject* args, PyObject* keywds);
static PyObject* OSCARSSR_GetTrajectory (OSCARSSRObject* self);
static PyObject* OSCARSSR_SetTrajectory (OSCARSSRObject* self, PyObject* args, PyObject* keywds);
static PyObject* OSCARSSR_CalculateSpectrum (OSCARSSRObject* self, PyObject* args, PyObject* keywds);
static PyObject* OSCARSSR_CalculateTotalPower (OSCARSSRObject* self, PyObject* args, PyObject *keywds);
static PyObject* OSCARSSR_CalculatePowerDensity (OSCARSSRObject* self, PyObject* args, PyObject *keywds);
static PyObject* OSCARSSR_CalculatePowerDensityRectangle (OSCARSSRObject* self, PyObject* args, PyObject *keywds);
static PyObject* OSCARSSR_CalculateFlux (OSCARSSRObject* self, PyObject* args, PyObject *keywds);
static PyObject* OSCARSSR_CalculateFluxRectangle (OSCARSSRObject* self, PyObject* args, PyObject *keywds);
static PyObject* OSCARSSR_WriteSpectrum (OSCARSSRObject* self, PyObject* args, PyObject *keywds);
static PyObject* OSCARSSR_AverageSpectra (OSCARSSRObject* self, PyObject* args, PyObject *keywds);
static PyObject* OSCARSSR_AddToSpectrum (OSCARSSRObject* self, PyObject* args, PyObject *keywds);
static PyObject* OSCARSSR_GetSpectrum (OSCARSSRObject* self);
static PyObject* OSCARSSR_AverageT3DScalars (OSCARSSRObject* self, PyObject* args, PyObject *keywds);
static PyObject* OSCARSSR_AddToFlux (OSCARSSRObject* self, PyObject* args, PyObject *keywds);
static PyObject* OSCARSSR_GetFlux (OSCARSSRObject* self);
static PyObject* OSCARSSR_AddToPowerDensity (OSCARSSRObject* self, PyObject* args, PyObject *keywds);
static PyObject* OSCARSSR_GetPowerDensity (OSCARSSRObject* self);
static PyObject* OSCARSSR_CalculateElectricFieldTimeDomain (OSCARSSRObject* self, PyObject* args, PyObject *keywds);
static PyObject* OSCARSSR_PrintGPU (OSCARSSRObject* self);
static PyObject* OSCARSSR_PrintNThreads (OSCARSSRObject* self);
static PyObject* OSCARSSR_PrintTrajectory (OSCARSSRObject* self);
static PyObject* OSCARSSR_PrintAll (OSCARSSRObject* self);
static PyObject* OSCARSSR_COUT (OSCARSSRObject* self, PyObject* args, PyObject *keywds);
static PyObject* OSCARSSR_CERR (OSCARSSRObject* self, PyObject* args, PyObject *keywds);


static PyObject* OSCARSSR_Fake (OSCARSSRObject* self, PyObject* args, PyObject *keywds);
static PyObject* OSCARSSR_GetT3DScalarAsList (T3DScalarContainer const& C);

TSpectrumContainer OSCARSSR_GetSpectrumFromList (PyObject* List);
T3DScalarContainer OSCARSSR_GetT3DScalarContainerFromList (PyObject* List);























#endif
