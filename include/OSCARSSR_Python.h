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
static TVector3D OSCARSSR_ListAsTVector3D (PyObject* List);
static PyObject* OSCARSSR_TVector3DAsList (TVector3D const& V);
static PyObject* OSCARSSR_Pi (OSCARSSRObject* self, PyObject* arg);
static PyObject* OSCARSSR_GetCTStart (OSCARSSRObject* self);
static PyObject* OSCARSSR_GetCTStop (OSCARSSRObject* self);
static PyObject* OSCARSSR_SetCTStartStop (OSCARSSRObject* self, PyObject* args);
static PyObject* OSCARSSR_GetNPointsTrajectory (OSCARSSRObject* self);
static PyObject* OSCARSSR_SetNPointsTrajectory (OSCARSSRObject* self, PyObject* arg);
static PyObject* OSCARSSR_AddMagneticField (OSCARSSRObject* self, PyObject* args, PyObject* keywds);
static PyObject* OSCARSSR_AddMagneticFieldFunction (OSCARSSRObject* self, PyObject* args);
static PyObject* OSCARSSR_AddMagneticFieldGaussian (OSCARSSRObject* self, PyObject* args, PyObject* keywds);
static PyObject* OSCARSSR_ClearMagneticFields (OSCARSSRObject* self);
static PyObject* OSCARSSR_GetBField (OSCARSSRObject* self, PyObject* args);
static PyObject* OSCARSSR_AddParticleBeam (OSCARSSRObject* self, PyObject* args, PyObject* keywds);
static PyObject* OSCARSSR_SetNewParticle (OSCARSSRObject* self, PyObject* args, PyObject* keywds);
static PyObject* OSCARSSR_ClearParticleBeams (OSCARSSRObject* self);
static PyObject* OSCARSSR_CalculateTrajectory (OSCARSSRObject* self);
static PyObject* OSCARSSR_GetTrajectory (OSCARSSRObject* self);
static PyObject* OSCARSSR_GetSpectrum (OSCARSSRObject* self);
static PyObject* OSCARSSR_CalculateSpectrum (OSCARSSRObject* self, PyObject* args, PyObject* keywds);
static PyObject* OSCARSSR_CalculateTotalPower (OSCARSSRObject* self);
static PyObject* OSCARSSR_CalculatePowerDensityRectangle (OSCARSSRObject* self, PyObject* args, PyObject *keywds);
static PyObject* OSCARSSR_CalculateFluxRectangle (OSCARSSRObject* self, PyObject* args, PyObject *keywds);























#endif
