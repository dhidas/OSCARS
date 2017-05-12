////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Tue Jan 31 17:38:57 EST 2017
//
// This is the python-C extension which allows access to the c++
// class OSCARSTH.
//
////////////////////////////////////////////////////////////////////

// Include Python.h first!
#include <Python.h>

#include "OSCARSTH_Python.h"

#include "OSCARSTH.h"
#include "Version.h"

#include "TOMATH.h"

#include "TOSCARSSR.h"
#include "TVector2D.h"

#include <iostream>
#include <vector>
#include <algorithm>
#include <cstring>
#include <stdexcept>
#include <sstream>






static TVector2D OSCARSTH_ListAsTVector2D (PyObject* List)
{
  TVector2D V;
  if (PyList_Size(List) == 2) {
    Py_INCREF(List);
    V.SetXY(PyFloat_AsDouble(PyList_GetItem(List, 0)),
             PyFloat_AsDouble(PyList_GetItem(List, 1)));
    Py_DECREF(List);
  } else {
    throw std::length_error("number of elements not 2");
  }

  // Return the python list
  return V;
}



static PyObject* OSCARSTH_TVector2DAsList (TVector2D const& V)
{
  // Turn a TVector3D into a list (like a vector)

  // Create a python list
  PyObject *PList = PyList_New(0);

  PyList_Append(PList, Py_BuildValue("f", V.GetX()));
  PyList_Append(PList, Py_BuildValue("f", V.GetY()));

  // Return the python list
  return PList;
}









static void OSCARSTH_dealloc(OSCARSTHObject* self)
{
  // Python needs to know how to deallocate things in the struct

  delete self->obj;
  //self->ob_type->tp_free((PyObject*) self);
  Py_TYPE(self)->tp_free((PyObject*) self);
}




static PyObject* OSCARSTH_new (PyTypeObject* type, PyObject* args, PyObject* kwds)
{
  // Python needs to know how to create things in this struct

  OSCARSTHObject* self = (OSCARSTHObject*) type->tp_alloc(type, 0);
  if (self != NULL) {

    // Create the new object for self
    self->obj = new OSCARSTH();
  }

  // Return myself
  return (PyObject*) self;
}





static int th_init(OSCARSTHObject* self, PyObject* args, PyObject* kwds)
{
  return 0;
}
















const char* DOC_OSCARSTH_UndulatorK = "Get the undulator K parameter vaule";
static PyObject* OSCARSTH_UndulatorK (OSCARSTHObject* self, PyObject* args, PyObject* keywds)
{
  // Return the undulator K parameter given peak bfield and period

  // Require 2 arguments
  double BFieldMax = 0;
  double Period = 0;

  // Input variables and parsing
  static char *kwlist[] = {"bfield", "period", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, keywds, "dd", kwlist, &BFieldMax, &Period)) {
    return NULL;
  }


  // Return the deflection parameter K
  return Py_BuildValue("d", self->obj->UndulatorK(BFieldMax, Period));
}





const char* DOC_OSCARSTH_UndulatorBField = "Get the undulator BFieldMax parameter value";
static PyObject* OSCARSTH_UndulatorBField (OSCARSTHObject* self, PyObject* args, PyObject* keywds)
{
  // Return the undulator BFieldMax parameter give the K and period

  // Require 2 arguments
  double K = 0;
  double Period = 0;

  // Input variables and parsing
  static char *kwlist[] = {"K", "period", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, keywds, "dd", kwlist, &K, &Period)) {
    return NULL;
  }


  // Return the BField Max value
  return Py_BuildValue("d", self->obj->UndulatorBField(K, Period));
}





const char* DOC_OSCARSTH_UndulatorPeriod = "Get the undulator Period parameter vaule";
static PyObject* OSCARSTH_UndulatorPeriod (OSCARSTHObject* self, PyObject* args, PyObject* keywds)
{
  // Return the undulator Period parameter given peak bfield and K

  // Require 2 arguments
  double BFieldMax = 0;
  double K = 0;

  // Input variables and parsing
  static char *kwlist[] = {"bfield", "K", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, keywds, "dd", kwlist, &BFieldMax, &K)) {
    return NULL;
  }


  // Return the undulator Period
  return Py_BuildValue("d", self->obj->UndulatorPeriod(BFieldMax, K));
}





const char* DOC_OSCARSTH_DipoleSpectrum = "Get the spectrum from ideal dipole field";
static PyObject* OSCARSTH_DipoleSpectrum (OSCARSTHObject* self, PyObject* args, PyObject* keywds)
{
  // Return a list of points corresponding to the flux in a given energy range for a given vertical angle.
  // This approximation assumes that the particle beam is perpendicular to the magnetic field

  // Require 2 arguments
  double BField = 0;
  double BeamEnergy = 0;
  double Angle = 0;
  double Energy_eV = 0;
  //PyObject* List_EnergyRange = PyList_New(0);

  // Input variables and parsing
  static char *kwlist[] = {"bfield", "beam_energy_GeV", "angle", "energy_eV", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, keywds, "dddd", kwlist, &BField, &BeamEnergy, &Angle, &Energy_eV)) {
    return NULL;
  }

  // Check that beam energy makes sense
  if (BeamEnergy <= 0) {
    PyErr_SetString(PyExc_ValueError, "'beam_energy_GeV' must be > 0");
    return NULL;
  }

  //TVector2D const EnergyRange = OSCARSTH_ListAsTVector2D(List_EnergyRange);
  //if (EnergyRange[0] >= EnergyRange[1] || EnergyRange[0] <= 1 || EnergyRange[1] <= 0) {
  //  PyErr_SetString(PyExc_ValueError, "'energy_range_eV' is incorrect");
  //  return NULL;
  //}

  // Calculate the spectrum
  double const Result = self->obj->DipoleSpectrum(BField, BeamEnergy, Angle, Energy_eV);

  return Py_BuildValue("d", Result);
  // Must return python object None in a special way
  //Py_INCREF(Py_None);
  //return Py_None;
}






const char* DOC_OSCARSTH_DipoleCriticalEnergy = "Get the critical energy for bending magnet in [eV]";
static PyObject* OSCARSTH_DipoleCriticalEnergy (OSCARSTHObject* self, PyObject* args, PyObject* keywds)
{
  // Return a list of points corresponding to the flux in a given energy range for a given vertical angle.
  // This approximation assumes that the particle beam is perpendicular to the magnetic field

  // Require 2 arguments
  double BField = 0;
  double BeamEnergy = 0;

  // Input variables and parsing
  static char *kwlist[] = {"bfield", "beam_energy_GeV", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, keywds, "dd", kwlist, &BField, &BeamEnergy)) {
    return NULL;
  }

  // Check that beam energy makes sense
  if (BeamEnergy <= 0) {
    PyErr_SetString(PyExc_ValueError, "'beam_energy_GeV' must be > 0");
    return NULL;
  }

  // Calculate the spectrum
  double const Result = self->obj->DipoleCriticalEnergy(BField, BeamEnergy);

  return Py_BuildValue("d", Result);
}






const char* DOC_OSCARSTH_DipoleBrightness = "Get the spectrum from ideal dipole field";
static PyObject* OSCARSTH_DipoleBrightness (OSCARSTHObject* self, PyObject* args, PyObject* keywds)
{
  // Return a list of points corresponding to the flux in a given energy range for a given vertical angle.
  // This approximation assumes that the particle beam is perpendicular to the magnetic field

  // Require 2 arguments
  double BField = 0;
  double BeamEnergy = 0;
  double Angle = 0;
  double Energy_eV = 0;
  //PyObject* List_EnergyRange = PyList_New(0);

  // Input variables and parsing
  static char *kwlist[] = {"bfield", "beam_energy_GeV", "angle", "energy_eV", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, keywds, "|dddd", kwlist, &BField, &BeamEnergy, &Angle, &Energy_eV)) {
    return NULL;
  }

  // Check that beam energy makes sense
  //if (BeamEnergy <= 0) {
  //  PyErr_SetString(PyExc_ValueError, "'beam_energy_GeV' must be > 0");
  //  return NULL;
  //}

  //TVector2D const EnergyRange = OSCARSTH_ListAsTVector2D(List_EnergyRange);
  //if (EnergyRange[0] >= EnergyRange[1] || EnergyRange[0] <= 1 || EnergyRange[1] <= 0) {
  //  PyErr_SetString(PyExc_ValueError, "'energy_range_eV' is incorrect");
  //  return NULL;
  //}

  // Calculate the spectrum
  double const Result = self->obj->DipoleBrightness();

  return Py_BuildValue("d", Result);
  // Must return python object None in a special way
  //Py_INCREF(Py_None);
  //return Py_None;
}







const char* DOC_OSCARSTH_UndulatorFlux = "Get the flux for an ideal undulator";
static PyObject* OSCARSTH_UndulatorFlux (OSCARSTHObject* self, PyObject* args, PyObject* keywds)
{
  // Return the flux [gamma/s/mrad^2/0.1%bw] at a given angle
  // This approximation assumes that the particle beam is perpendicular to the magnetic field

  // Require 2 arguments
  double BField = 0;
  double NPeriods = 0;
  double Period = 0;
  double BeamEnergy = 0;
  double AngleV = 0;
  double AngleH = 0;
  double Energy_eV = 0;

  // Input variables and parsing
  static char *kwlist[] = {"bfield", "period", "nperiods", "beam_energy_GeV", "angle_vertical", "angle_horizontal", "energy_eV", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, keywds, "ddddddd", kwlist, &BField, &Period, &NPeriods, &BeamEnergy, &AngleV, &AngleH, &Energy_eV)) {
    return NULL;
  }

  // Check that beam energy makes sense
  if (BeamEnergy <= 0) {
    PyErr_SetString(PyExc_ValueError, "'beam_energy_GeV' must be > 0");
    return NULL;
  }

  // Check period and nperiods
  if (Period <= 0 || NPeriods <= 0) {
    PyErr_SetString(PyExc_ValueError, "'period' and 'nperiods' must be > 0");
    return NULL;
  }

  // Check that photon energy is > 0
  if (Energy_eV <= 0) {
    PyErr_SetString(PyExc_ValueError, "'energy_eV' must be > 0");
    return NULL;
  }



  // Calculate the spectrum
  double const Result = self->obj->UndulatorFlux(BField, Period, NPeriods, BeamEnergy, AngleV, AngleH, Energy_eV);

  return Py_BuildValue("d", Result);
}

















const char* DOC_OSCARSTH_UndulatorFluxOnAxis = "Get the on-axis flux for an ideal undulator given K for a specific harmonic";
static PyObject* OSCARSTH_UndulatorFluxOnAxis (OSCARSTHObject* self, PyObject* args, PyObject* keywds)
{
  // Return the flux [gamma/s/mrad^2/0.1%bw] at a given K for a given harmonic

  // Require some arguments
  int    NPeriods = 0;
  double Period = 0;
  int    Harmonic = 0;
  double K = 0;
  double BField = 0;

  // Input variables and parsing
  static char *kwlist[] = {"period", "nperiods", "harmonic", "K", "bfield", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, keywds, "dii|dd", kwlist, 
                                                           &Period,
                                                           &NPeriods,
                                                           &Harmonic,
                                                           &K,
                                                           &BField
                                                           )) {
    return NULL;
  }

  // Check that bfield is > 0 or K > 0
  if (BField <= 0 && K <= 0) {
    PyErr_SetString(PyExc_ValueError, "'bfield' or 'K' must be > 0");
    return NULL;
  }

  // Check that bfield and K not both defined
  if (BField > 0 && K > 0) {
    PyErr_SetString(PyExc_ValueError, "use 'bfield' or 'K' but not both");
    return NULL;
  }


  // Check period and nperiods
  if (Period <= 0 || NPeriods <= 0) {
    PyErr_SetString(PyExc_ValueError, "'period' and 'nperiods' must be > 0");
    return NULL;
  }

  // Check that photon energy is > 0
  if (Harmonic <= 0) {
    PyErr_SetString(PyExc_ValueError, "'harmonic' must be > 0");
    return NULL;
  }


  // Calculate the spectrum
  TVector2D Result;

  if (K > 0) {
    Result = self->obj->UndulatorFluxOnAxisK(K, Period, NPeriods, Harmonic);
  } else if (BField > 0) {
    Result = self->obj->UndulatorFluxOnAxisB(BField, Period, NPeriods, Harmonic);
  } else {
    PyErr_SetString(PyExc_ValueError, "not B nor K.  please report this bug");
    return NULL;
  }

  return OSCARSTH_TVector2DAsList(Result);
}









const char* DOC_OSCARSTH_UndulatorFluxWeak = "Undulator flux weak approx";
static PyObject* OSCARSTH_UndulatorFluxWeak (OSCARSTHObject* self, PyObject* args, PyObject* keywds)
{
  // Return the flux [gamma/s/mrad^2/0.1%bw] at a given K for a given harmonic

  // Require some arguments
  double K = 0;
  int    NPeriods = 0;
  double Period = 0;
  double BeamEnergy = 0;
  int    Harmonic = 0;

  // Input variables and parsing
  static char *kwlist[] = {"k", "period", "nperiods", "beam_energy_GeV", "harmonic", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, keywds, "ddidi", kwlist, 
                                                          &K,
                                                          &Period,
                                                          &NPeriods,
                                                          &BeamEnergy,
                                                          &Harmonic)) {
    std::cout << "K             " << K << std::endl;
    std::cout << "Period:       " << Period << std::endl;
    std::cout << "NPeriods:     " << NPeriods << std::endl;
    std::cout << "BeamEnergy:   " << BeamEnergy << std::endl;
    std::cout << "Harmonic      " << Harmonic << std::endl;
    std::cout << "Cannot parse for some reason" << std::endl;
    return NULL;
  }

  // Check that K is > 0
  if (K <= 0) {
    PyErr_SetString(PyExc_ValueError, "'k' must be > 0");
    return NULL;
  }


  // Check that beam energy makes sense
  if (BeamEnergy <= 0) {
    PyErr_SetString(PyExc_ValueError, "'beam_energy_GeV' must be > 0");
    return NULL;
  }

  // Check period and nperiods
  if (Period <= 0 || NPeriods <= 0) {
    PyErr_SetString(PyExc_ValueError, "'period' and 'nperiods' must be > 0");
    return NULL;
  }

  // Check that photon energy is > 0
  if (Harmonic <= 0) {
    PyErr_SetString(PyExc_ValueError, "'harmonic' must be > 0");
    return NULL;
  }

  // Calculate the spectrum
  double const Result = self->obj->UndulatorFluxWeak(K, Period, NPeriods, BeamEnergy, Harmonic);

  return Py_BuildValue("d", Result);
}









const char* DOC_OSCARSTH_UndulatorBrightness = "Undulator brightness calculation";
static PyObject* OSCARSTH_UndulatorBrightness (OSCARSTHObject* self, PyObject* args, PyObject* keywds)
{
  // Return the brightness [gamma/s/mrad^2/mm^2/0.1%bw] at a given K for a given harmonic

  // Require some arguments
  double Period   = 0;
  int    NPeriods = 0;
  int    Harmonic = 0;
  double BField   = 0;
  double K        = 0;

  // Input variables and parsing
  static char *kwlist[] = {"period", "nperiods", "harmonic", "bfield", "K", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, keywds, "dii|dd", kwlist, 
                                                         &Period,
                                                         &NPeriods,
                                                         &Harmonic,
                                                         &BField,
                                                         &K
                                                         )) {
    return NULL;
  }


  // Check that bfield is > 0 or K > 0
  if (BField <= 0 && K <= 0) {
    PyErr_SetString(PyExc_ValueError, "'bfield' or 'K' must be > 0");
    return NULL;
  }

  // Check that bfield and K not both defined
  if (BField > 0 && K > 0) {
    PyErr_SetString(PyExc_ValueError, "use 'bfield' or 'K' but not both");
    return NULL;
  }


  // Check period and nperiods
  if (Period <= 0 || NPeriods <= 0) {
    PyErr_SetString(PyExc_ValueError, "'period' and 'nperiods' must be > 0");
    return NULL;
  }

  // Check that photon energy is > 0
  if (Harmonic <= 0) {
    PyErr_SetString(PyExc_ValueError, "'harmonic' must be > 0");
    return NULL;
  }


  TVector2D Result;

  if (K > 0) {
    Result = self->obj->UndulatorBrightnessK(K, Period, NPeriods, Harmonic);
  } else if (BField > 0) {
    Result = self->obj->UndulatorBrightnessB(BField, Period, NPeriods, Harmonic);
  } else {
    PyErr_SetString(PyExc_ValueError, "not B nor K.  please report this bug");
    return NULL;
  }


  return OSCARSTH_TVector2DAsList(Result);
}









const char* DOC_OSCARSTH_UndulatorEnergyHarmonic = "Undulator photon energy at particular harmonic";
static PyObject* OSCARSTH_UndulatorEnergyHarmonic (OSCARSTHObject* self, PyObject* args, PyObject* keywds)
{
  // Return the brightness [gamma/s/mrad^2/mm^2/0.1%bw] at a given K for a given harmonic

  // Require some arguments
  double Period = 0;
  int    Harmonic = 0;
  double BeamEnergy_GeV = 0;
  double BField = 0;
  double K = 0;

  // Input variables and parsing
  static char *kwlist[] = {"period", "harmonic", "bfield", "K", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, keywds, "di|dd", kwlist, 
                                                           &Period,
                                                           &Harmonic,
                                                           &BField,
                                                           &K
                                                           )) {
    return NULL;
  }

  // Check that K is > 0
  if (BField <= 0 && K <= 0) {
    PyErr_SetString(PyExc_ValueError, "'bfield' or 'K' must be > 0");
    return NULL;
  }


  // Check period and nperiods
  if (Period <= 0) {
    PyErr_SetString(PyExc_ValueError, "'period' must be > 0");
    return NULL;
  }

  // Check that photon energy is > 0
  if (Harmonic <= 0) {
    PyErr_SetString(PyExc_ValueError, "'harmonic' must be > 0");
    return NULL;
  }


  double Result = 0;
  if (BField > 0) {
    Result = self->obj->UndulatorEnergyAtHarmonicB(BField, Period, Harmonic);
  } else if (K > 0) {
    Result = self->obj->UndulatorEnergyAtHarmonicK(K, Period, Harmonic);
  } else {
    PyErr_SetString(PyExc_ValueError, "Something is wrong with bfield and K.  Repot this bug.");
    return NULL;
  }

  return Py_BuildValue("d", Result);
}












































const char* DOC_OSCARSTH_BesselJ = "Get the value of the modified bessel function J_nu for integer nu";
static PyObject* OSCARSTH_BesselJ (OSCARSTHObject* self, PyObject* args, PyObject* keywds)
{
  // Return Bessel J_nu(x)

  // Require 2 arguments
  int    Nu = 0;
  double X  = 0;

  // Input variables and parsing
  static char *kwlist[] = {"nu", "x", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, keywds, "id", kwlist, &Nu, &X)) {
    return NULL;
  }


  // Return the internal OSCARSTH number constant pi
  return Py_BuildValue("d", TOMATH::BesselJ(Nu, X));
}










const char* DOC_OSCARSTH_BesselK = "Get the value of the modified bessel function K_nu";
static PyObject* OSCARSTH_BesselK (OSCARSTHObject* self, PyObject* args, PyObject* keywds)
{
  // Return Bessel K_nu(x)

  // Require 2 arguments
  double Nu = 0;
  double X = 0;

  // Input variables and parsing
  static char *kwlist[] = {"nu", "x", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, keywds, "dd", kwlist, &Nu, &X)) {
    return NULL;
  }


  // Return the internal OSCARSTH number constant pi
  return Py_BuildValue("d", TOMATH::BesselK(Nu, X));
}











const char* DOC_OSCARSTH_SetParticleBeam = "Set the particle beam parameters";
static PyObject* OSCARSTH_SetParticleBeam (OSCARSTHObject* self, PyObject* args, PyObject* keywds)
{
  // Add a particle beam to the experiment

  // Lists and variables some with initial values
  char const* Type                       = "";
  char const* Name                       = "";
  double      Energy_GeV                 = 0;
  double      Sigma_Energy_GeV           = 0;
  double      T0                         = 0;
  double      Current                    = 0;
  double      Weight                     = 1;
  double      Mass                       = 0;
  double      Charge                     = 0;
  char const* Beam                       = "";
  PyObject*   List_Position              = PyList_New(0);
  PyObject*   List_Direction             = PyList_New(0);
  PyObject*   List_Rotations             = PyList_New(0);
  PyObject*   List_Translation           = PyList_New(0);
  PyObject*   List_Horizontal_Direction  = PyList_New(0);
  PyObject*   List_Beta                  = PyList_New(0);
  PyObject*   List_Emittance             = PyList_New(0);
  PyObject*   List_Lattice_Reference     = PyList_New(0);

  TVector3D Position(0, 0, 0);
  TVector3D Direction;
  TVector3D Rotations(0, 0, 0);
  TVector3D Translation(0, 0, 0);
  TVector3D Horizontal_Direction;
  TVector2D Beta(0, 0);
  TVector2D Emittance(0, 0);
  TVector3D Lattice_Reference(0, 0, 0);


  // Input variables and parsing
  static char *kwlist[] = {"type", "name", "energy_GeV", "d0", "x0", "beam", "sigma_energy_GeV", "t0", "current", "weight", "rotations", "translation", "horizontal_direction", "beta", "emittance", "lattice_reference", "mass", "charge", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, keywds, "|ssdOOsddddOOOOOOdd", kwlist,
                                                                        &Type,
                                                                        &Name,
                                                                        &Energy_GeV,
                                                                        &List_Direction,
                                                                        &List_Position,
                                                                        &Beam,
                                                                        &Sigma_Energy_GeV,
                                                                        &T0,
                                                                        &Current,
                                                                        &Weight,
                                                                        &List_Rotations,
                                                                        &List_Translation,
                                                                        &List_Horizontal_Direction,
                                                                        &List_Beta,
                                                                        &List_Emittance,
                                                                        &List_Lattice_Reference,
                                                                        &Mass,
                                                                        &Charge)) {
    return NULL;
  }


  // Are you asking for one of the predefined beams?
  bool const HasPredefinedBeam = std::strlen(Beam) != 0 ? true : false;

  // Check if beam is defined (for predefined beams)
  if (HasPredefinedBeam) {
    try {
      self->obj->SetParticleBeam(Beam);
    } catch (...) {
      PyErr_SetString(PyExc_ValueError, "Error in predefined beam name");
      return NULL;
    }
  }


  if (Sigma_Energy_GeV == 0) {
    // Do nothing.  zero energy diff is alright
  } else if (Sigma_Energy_GeV < 0) {
    PyErr_SetString(PyExc_ValueError, "'sigma_energy_GeV' cannot be less than zero");
    return NULL;
  } else {
    // Change predefined beam accordingly
    if (HasPredefinedBeam) {
      self->obj->GetParticleBeam().SetSigmaEnergyGeV(Sigma_Energy_GeV);
    }
  }


  // Check for Beta in the input
  if (PyList_Size(List_Beta) != 0) {
    try {
      Beta = OSCARSTH_ListAsTVector2D(List_Beta);
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'beta'");
      return NULL;
    }
  }


  // Check for Emittance in the input
  if (PyList_Size(List_Emittance) != 0) {
    try {
      Emittance = OSCARSTH_ListAsTVector2D(List_Emittance);
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'emittance'");
      return NULL;
    }
  }

  // Check type
  if (std::string(Type) == "") {
    Type = "electron";
  }

  // Add the particle beam
  if (std::strlen(Beam) == 0) {
    try {
      if (std::string(Type) != "electron") {
        PyErr_SetString(PyExc_ValueError, "type must be electron");
        return NULL;
      } else {
        self->obj->SetParticleBeam(Energy_GeV, Current, Beta, Emittance, Sigma_Energy_GeV);
      }
    } catch (std::invalid_argument e) {
      PyErr_SetString(PyExc_ValueError, "invalid argument in adding particle beam.  possibly 'name' already exists");
      return NULL;
    }
  }

  if (T0 != 0) {
    self->obj->GetParticleBeam().SetT0(T0);
  }

  // Must return python object None in a special way
  Py_INCREF(Py_None);
  return Py_None;
}






const char* DOC_OSCARSTH_Print = "Print beam information";
static PyObject* OSCARSTH_Print (OSCARSTHObject* self)
{
  // Print all particle beams stored in OSCARSSR

  // Out string stream for printing beam information
  std::ostringstream ostream;
  ostream << "oscars.th object:\n" << self->obj->GetParticleBeam() << std::endl;

  // Python printing
  PyObject* sys = PyImport_ImportModule("sys");
  PyObject* s_out = PyObject_GetAttrString(sys, "stdout");
  std::string Message = ostream.str();
  PyObject_CallMethod(s_out, "write", "s", Message.c_str());

  // Must return python object None in a special way
  Py_INCREF(Py_None);
  return Py_None;
}



































static PyMethodDef OSCARSTH_methods[] = {
  // We must tell python about the function we allow access as well as give them nice
  // python names, and tell python the method of input parameters.

  {"undulator_K",                                (PyCFunction) OSCARSTH_UndulatorK,                              METH_VARARGS | METH_KEYWORDS,                  DOC_OSCARSTH_UndulatorK},
  {"undulator_bfield",                           (PyCFunction) OSCARSTH_UndulatorBField,                         METH_VARARGS | METH_KEYWORDS,                  DOC_OSCARSTH_UndulatorBField},
  {"undulator_period",                           (PyCFunction) OSCARSTH_UndulatorPeriod,                         METH_VARARGS | METH_KEYWORDS,                  DOC_OSCARSTH_UndulatorPeriod},

  {"dipole_spectrum",                            (PyCFunction) OSCARSTH_DipoleSpectrum,                          METH_VARARGS | METH_KEYWORDS,                  DOC_OSCARSTH_DipoleSpectrum},
  {"dipole_critical_energy",                     (PyCFunction) OSCARSTH_DipoleCriticalEnergy,                    METH_VARARGS | METH_KEYWORDS,                  DOC_OSCARSTH_DipoleCriticalEnergy},
  {"dipole_brightness",                          (PyCFunction) OSCARSTH_DipoleBrightness,                        METH_VARARGS | METH_KEYWORDS,                  DOC_OSCARSTH_DipoleBrightness},
  {"undulator_flux",                             (PyCFunction) OSCARSTH_UndulatorFlux,                           METH_VARARGS | METH_KEYWORDS,                  DOC_OSCARSTH_UndulatorFlux},
  {"undulator_flux_onaxis",                      (PyCFunction) OSCARSTH_UndulatorFluxOnAxis,                     METH_VARARGS | METH_KEYWORDS,                  DOC_OSCARSTH_UndulatorFluxOnAxis},
  {"undulator_flux_weak",                        (PyCFunction) OSCARSTH_UndulatorFluxWeak,                       METH_VARARGS | METH_KEYWORDS,                  DOC_OSCARSTH_UndulatorFluxWeak},
  {"undulator_brightness",                       (PyCFunction) OSCARSTH_UndulatorBrightness,                     METH_VARARGS | METH_KEYWORDS,                  DOC_OSCARSTH_UndulatorBrightness},
  {"undulator_energy_harmonic",                  (PyCFunction) OSCARSTH_UndulatorEnergyHarmonic,                 METH_VARARGS | METH_KEYWORDS,                  DOC_OSCARSTH_UndulatorEnergyHarmonic},
  {"bessel_j",                                   (PyCFunction) OSCARSTH_BesselJ,                                 METH_VARARGS | METH_KEYWORDS,                  DOC_OSCARSTH_BesselK},
  {"bessel_k",                                   (PyCFunction) OSCARSTH_BesselK,                                 METH_VARARGS | METH_KEYWORDS,                  DOC_OSCARSTH_BesselK},

  {"set_particle_beam",                          (PyCFunction) OSCARSTH_SetParticleBeam,                         METH_VARARGS | METH_KEYWORDS,                  DOC_OSCARSTH_SetParticleBeam},
  {"add_particle_beam",                          (PyCFunction) OSCARSTH_SetParticleBeam,                         METH_VARARGS | METH_KEYWORDS,                  DOC_OSCARSTH_SetParticleBeam},

  {"print",                                      (PyCFunction) OSCARSTH_Print,                                   METH_NOARGS,                                   DOC_OSCARSTH_Print},
  {"print_particle_beams",                       (PyCFunction) OSCARSTH_Print,                                   METH_NOARGS,                                   DOC_OSCARSTH_Print},

  {NULL}  /* Sentinel */
};


#if PY_MAJOR_VERSION >= 3
PyMODINIT_FUNC PyInit_th(void);
#else
PyMODINIT_FUNC initth(OSCARSTHObject* self, PyObject* args, PyObject* kwds);
#endif





#if PY_MAJOR_VERSION >= 3
static PyTypeObject OSCARSTHType = {
  PyVarObject_HEAD_INIT(NULL, 0)
  "th",            /* tp_name */
  sizeof(OSCARSTHObject),       /* tp_basicsize */
  0,                          /* tp_itemsize */
  (destructor)OSCARSTH_dealloc, /* tp_dealloc */
  0,                          /* tp_print */
  0,                          /* tp_getattr */
  0,                          /* tp_setattr */
  0,                          /* tp_reserved */
  0,                          /* tp_repr */
  0,                          /* tp_as_number */
  0,                          /* tp_as_sequence */
  0,                          /* tp_as_mapping */
  0,                          /* tp_hash  */
  0,                          /* tp_call */
  0,                          /* tp_str */
  0,                          /* tp_getattro */
  0,                          /* tp_setattro */
  0,                          /* tp_as_buffer */
  Py_TPFLAGS_DEFAULT |
  Py_TPFLAGS_BASETYPE,        /* tp_flags */
  "oscars th class",           /* tp_doc */
  0,                          /* tp_traverse */
  0,                          /* tp_clear */
  0,                          /* tp_richcompare */
  0,                          /* tp_weaklistoffset */
  0,                          /* tp_iter */
  0,                          /* tp_iternext */
  OSCARSTH_methods,             /* tp_methods */
  0,                          /* tp_members */
  0,                          /* tp_getset */
  0,                          /* tp_base */
  0,                          /* tp_dict */
  0,                          /* tp_descr_get */
  0,                          /* tp_descr_set */
  0,                          /* tp_dictoffset */
  0,      /* tp_init */
  0,                          /* tp_alloc */
  OSCARSTH_new,                 /* tp_new */
};
#else
static PyTypeObject OSCARSTHType = {
  // The python object.  Fully defined elsewhere.  only put here what you need,
  // otherwise default values

  PyObject_HEAD_INIT(NULL)
  0,                                        /* ob_size */
  "th",                                 /* tp_name */
  sizeof(OSCARSTHObject),                         /* tp_basicsize */
  0,                                        /* tp_itemsize */
  (destructor) OSCARSTH_dealloc,                 /* tp_dealloc */
  0,                                        /* tp_print */
  0,                                        /* tp_getattr */
  0,                                        /* tp_setattr */
  0,                                        /* tp_compare */
  0,                                        /* tp_repr */
  0,                                        /* tp_as_number */
  0,                                        /* tp_as_sequence */
  0,                                        /* tp_as_mapping */
  0,                                        /* tp_hash */
  0,                                        /* tp_call */
  0,                                        /* tp_str */
  0,                                        /* tp_getattro */
  0,                                        /* tp_setattro */
  0,                                        /* tp_as_buffer */
  Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, /* tp_flags */
  "oscars th class",                              /* tp_doc */
  0,                                        /* tp_traverse */
  0,                                        /* tp_clear */
  0,                                        /* tp_richcompare */
  0,                                        /* tp_weaklistoffset */
  0,                                        /* tp_iter */
  0,                                        /* tp_iternext */
  OSCARSTH_methods,                             /* tp_methods */
  0,                                        /* tp_members */
  0,                                        /* tp_getset */
  0,                                        /* tp_base */
  0,                                        /* tp_dict */
  0,                                        /* tp_descr_get */
  0,                                        /* tp_descr_set */
  0,                                        /* tp_dictoffset */
  0,                                        /* tp_init */
  0,                                        /* tp_alloc */
  OSCARSTH_new,                                  /* tp_new */
};
#endif




static PyMethodDef module_methods[] = {
  // I do not need
  {NULL}  /* Sentinel */
};


#if PY_MAJOR_VERSION >= 3
static PyModuleDef OSCARSTHmodule = {
  PyModuleDef_HEAD_INIT,
  "th",
  "OSCARSTH module extension.",
  -1,
  NULL, NULL, NULL, NULL, NULL
};
#endif


#if PY_MAJOR_VERSION >= 3
PyMODINIT_FUNC PyInit_th(void)
#else
PyMODINIT_FUNC initth(OSCARSTHObject* self, PyObject* args, PyObject* kwds)
#endif
{
  if (PyType_Ready(&OSCARSTHType) < 0) {
#if PY_MAJOR_VERSION >= 3
    return NULL;
#else
    return;
#endif
  }

#if PY_MAJOR_VERSION >= 3
  PyObject* m = PyModule_Create(&OSCARSTHmodule);
#else
  PyObject *m = Py_InitModule("oscars.th", OSCARSTH_methods);
#endif
  if (m == NULL) {
#if PY_MAJOR_VERSION >= 3
    return NULL;
#else
    return;
#endif
  }

  Py_INCREF(&OSCARSTHType);
  PyModule_AddObject(m, "th", (PyObject *)&OSCARSTHType);

  // Print copyright notice
  PyObject* sys = PyImport_ImportModule( "sys");
  PyObject* s_out = PyObject_GetAttrString(sys, "stdout");
  std::string Message = "OSCARS v" + OSCARS::GetVersionString() + " - Open Source Code for Advanced Radiation Simulation\nBrookhaven National Laboratory, Upton NY, USA\nhttp://oscars.bnl.gov\noscars@bnl.gov\n";
  PyObject_CallMethod(s_out, "write", "s", Message.c_str());

#if PY_MAJOR_VERSION >= 3
  return m;
#else
  return;
#endif
}



