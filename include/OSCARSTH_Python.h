#ifndef GUARD_OSCARSTH_Python_h
#define GUARD_OSCARSTH_Python_h
////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Tue Jul 12 08:36:18 EDT 2016
//
// This is a header file for the OSCARSTH_Python 'OSCARSTH' module
//
////////////////////////////////////////////////////////////////////

// Include Python.h first!
#include <Python.h>

#include "OSCARSTH.h"

// The python OSCARSTH object
typedef struct {
  // Define the OSCARSTHObject struct which contains the class I want
  PyObject_HEAD
  OSCARSTH* obj;
} OSCARSTHObject;




static void OSCARSTH_dealloc(OSCARSTHObject* self);
static PyObject* OSCARSTH_new (PyTypeObject* type, PyObject* args, PyObject* kwds);
static PyObject* OSCARSTH_Version (OSCARSTHObject* self, PyObject* arg);
static PyObject* OSCARSTH_UndulatorK (OSCARSTHObject* self, PyObject* args, PyObject* keywds);
static PyObject* OSCARSTH_UndulatorBField (OSCARSTHObject* self, PyObject* args, PyObject* keywds);
static PyObject* OSCARSTH_UndulatorPeriod (OSCARSTHObject* self, PyObject* args, PyObject* keywds);
static PyObject* OSCARSTH_DipoleSpectrum (OSCARSTHObject* self, PyObject* args, PyObject* keywds);
static PyObject* OSCARSTH_DipoleSpectrumPoint (OSCARSTHObject* self, PyObject* args, PyObject* keywds);
static PyObject* OSCARSTH_DipoleCriticalEnergy (OSCARSTHObject* self, PyObject* args, PyObject* keywds);
static PyObject* OSCARSTH_DipoleCriticalWavelength (OSCARSTHObject* self, PyObject* args, PyObject* keywds);
static PyObject* OSCARSTH_DipoleBrightness (OSCARSTHObject* self, PyObject* args, PyObject* keywds);
static PyObject* OSCARSTH_UndulatorFluxOnAxis (OSCARSTHObject* self, PyObject* args, PyObject* keywds);
static PyObject* OSCARSTH_UndulatorBrightness (OSCARSTHObject* self, PyObject* args, PyObject* keywds);
static PyObject* OSCARSTH_UndulatorEnergyHarmonic (OSCARSTHObject* self, PyObject* args, PyObject* keywds);
static PyObject* OSCARSTH_WigglerSpectrum (OSCARSTHObject* self, PyObject* args, PyObject* keywds);
static PyObject* OSCARSTH_WigglerFluxRectangle (OSCARSTHObject* self, PyObject* args, PyObject* keywds);
static PyObject* OSCARSTH_BesselJ (OSCARSTHObject* self, PyObject* args, PyObject* keywds);
static PyObject* OSCARSTH_BesselK (OSCARSTHObject* self, PyObject* args, PyObject* keywds);
static PyObject* OSCARSTH_SetParticleBeam (OSCARSTHObject* self, PyObject* args, PyObject* keywds);
static PyObject* OSCARSTH_AddParticleBeam (OSCARSTHObject* self, PyObject* args, PyObject* keywds);
static PyObject* OSCARSTH_PrintAll (OSCARSTHObject* self);
static PyObject* OSCARSTH_Fake (OSCARSTHObject* self, PyObject* args, PyObject *keywds);





















#endif
