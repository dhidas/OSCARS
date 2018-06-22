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

#include "OSCARSPY.h"
#include "OSCARSTH.h"
#include "Version.h"

#include "TOMATH.h"

#include "TOSCARSSR.h"
#include "TVector2D.h"
#include "TSurfacePoints_Rectangle.h"

#include <iostream>
#include <vector>
#include <algorithm>
#include <cstring>
#include <stdexcept>
#include <sstream>






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





//static int th_init(OSCARSTHObject* self, PyObject* args, PyObject* kwds)
//{
//  return 0;
//}








const char* DOC_OSCARSTH_Version = R"docstring(
version()

Version ID

Returns
-------
version : str
)docstring";
static PyObject* OSCARSTH_Version (OSCARSTHObject* self, PyObject* arg)
{
  return Py_BuildValue("s", OSCARSPY::GetVersionString().c_str());
}













const char* DOC_OSCARSTH_UndulatorK = R"docstring(
undulator_K(bfield, period)

Get the undulator K parameter vaule

Parameters
----------
bfield : float
    Peak magnetic field in [T]

period : float
    Magnetic period in [m]

Returns
-------
K : float
    Undulator deflectrion parameter
)docstring";
static PyObject* OSCARSTH_UndulatorK (OSCARSTHObject* self, PyObject* args, PyObject* keywds)
{
  // Return the undulator K parameter given peak bfield and period

  // Require 2 arguments
  double BFieldMax = 0;
  double Period = 0;

  // Input variables and parsing
  static const char *kwlist[] = {"bfield",
                                 "period",
                                 NULL};

  if (!PyArg_ParseTupleAndKeywords(args, keywds, "dd",
                                   const_cast<char **>(kwlist),
                                   &BFieldMax,
                                   &Period)) {
    return NULL;
  }


  // Return the deflection parameter K
  return Py_BuildValue("d", self->obj->UndulatorK(BFieldMax, Period));
}





const char* DOC_OSCARSTH_UndulatorBField = R"docstring(
undulator_bfield(K, period)

Get the undulator BFieldMax parameter value given K and the period

Parameters
----------
K : float
    Undulator deflection parameters

period : float
    Magnetic period in [m]

Returns
-------
bfield : float
    Peak bfield
)docstring";
static PyObject* OSCARSTH_UndulatorBField (OSCARSTHObject* self, PyObject* args, PyObject* keywds)
{
  // Return the undulator BFieldMax parameter give the K and period

  // Require 2 arguments
  double K = 0;
  double Period = 0;

  // Input variables and parsing
  static const char *kwlist[] = {"K",
                                 "period",
                                 NULL};

  if (!PyArg_ParseTupleAndKeywords(args, keywds, "dd",
                                   const_cast<char **>(kwlist),
                                   &K,
                                   &Period)) {
    return NULL;
  }

  // Return the BField Max value
  return Py_BuildValue("d", self->obj->UndulatorBField(K, Period));
}





const char* DOC_OSCARSTH_UndulatorPeriod = R"docstring(
undulator_period(bfield, K)

Get the undulator Period parameter vaule given this bfield and K

Parameters
----------

K : float
    Undulator deflection parameters

bfield : float
    Peak magnetic field in [T]

Returns
-------
period : float
    Magnetic period
)docstring";
static PyObject* OSCARSTH_UndulatorPeriod (OSCARSTHObject* self, PyObject* args, PyObject* keywds)
{
  // Return the undulator Period parameter given peak bfield and K

  // Require 2 arguments
  double BFieldMax = 0;
  double K = 0;

  // Input variables and parsing
  static const char *kwlist[] = {"bfield",
                                 "K",
                                 NULL};

  if (!PyArg_ParseTupleAndKeywords(args, keywds, "dd",
                                   const_cast<char **>(kwlist),
                                   &BFieldMax,
                                   &K)) {
    return NULL;
  }


  // Return the undulator Period
  return Py_BuildValue("d", self->obj->UndulatorPeriod(BFieldMax, K));
}





const char* DOC_OSCARSTH_DipoleSpectrum = R"docstring(
dipole_spectrum(bfield [, energy_range_eV, energy_points_eV, energy_eV, angle_integrated, angle_range, angle_points, angle, npoints, minimum, ofile, bofile])

Get the spectrum from ideal dipole field.  One can calculate the spectrum as a function of photon energy at a given angle or the angular dependence for a particular energy.

Parameters
----------
bfield : float
    Magnetic field of the dipole in [T]

energy_range_eV : list
    [min, max] photon energy of interest in [eV]

energy_points_eV : list
    List of energy points in [eV]

energy_eV : float
    Photon energy of interest

angle_integrated : boolean
    Set to true if you want the vertical angle integrated

angle_range : list
    [min, max] angle of interest in [rad]

angle_points : list
    List of angle points in [rad]

angle : float
    Angle of interest in [rad]

npoints : int
    Number of points for _range requests

ofile : str
    Output file name

bofile : str
    Binary output file name

Returns
-------
flux : list
    A list of flux values for given points [[p0, f0], [p1, f1], ...].  The points can be in energy or angle depending on the input parameters.  Output is in [photons / s / mrad^2 / 0.1\%bw]

Examples
--------
Calculate the spectrum in a given energy range

    >>> oth.dipole_spectrum(bfield=0.4, energy_range_eV=[10, 30000], npoints=1000)

Calculate the spectrum in a given energy range at a vertical angle of 0.0005 [rad]

    >>> oth.dipole_spectrum(bfield=0.4, energy_range_eV=[10, 30000], npoints=1000, angle=0.0005)

Calculate the flux in a given angle range at an energy of 2394 [eV]

    >>> oth.dipole_spectrum(bfield=0.4, angle_range=[-0.0005, 0.0005], npoints=1000, energy_eV=2394)
)docstring";
static PyObject* OSCARSTH_DipoleSpectrum (OSCARSTHObject* self, PyObject* args, PyObject* keywds)
{
  // Return a list of points corresponding to the flux in a given energy range for a given vertical angle.
  // This approximation assumes that the particle beam is perpendicular to the magnetic field

  // Require 2 arguments
  double      BField               = 0;
  PyObject*   List_EnergyRange_eV  = PyList_New(0);
  PyObject*   List_EnergyPoints_eV = PyList_New(0);
  double      Energy_eV            = 0;
  bool        AngleIntegrated      = false;
  PyObject*   List_AngleRange      = PyList_New(0);
  PyObject*   List_AnglePoints     = PyList_New(0);
  double      Angle                = 0;
  int         NPoints              = 0;
  const char* OutFileNameText      = "";
  const char* OutFileNameBinary    = "";

  // Input variable list
  static const char *kwlist[] = {"bfield",
                                 "energy_range_eV",
                                 "energy_points_eV",
                                 "energy_eV",
                                 "angle_integrated",
                                 "angle_range",
                                 "angle_points",
                                 "angle",
                                 "npoints",
                                 "ofile",
                                 "bofile",
                                 NULL};

  // Parse inputs
  if (!PyArg_ParseTupleAndKeywords(args, keywds, "d|OOdpOOdiss",
                                   const_cast<char **>(kwlist),
                                   &BField,
                                   &List_EnergyRange_eV,
                                   &List_EnergyPoints_eV,
                                   &Energy_eV,
                                   &AngleIntegrated,
                                   &List_AngleRange,
                                   &List_AnglePoints,
                                   &Angle,
                                   &NPoints,
                                   &OutFileNameText,
                                   &OutFileNameBinary)) {
    return NULL;
  }

  // CHeck if beam is ok
  if (!self->obj->CheckBeam()) {
    PyErr_SetString(PyExc_ValueError, "particle beam not correctly defined");
    return NULL;
  }

  // Must have a bfield
  if (fabs(BField) == 0) {
    PyErr_SetString(PyExc_ValueError, "'bfield' must not be zero");
    return NULL;
  }

  // Must have some combination that makes sense
  if (PyList_Size(List_EnergyRange_eV) != 0 && PyList_Size(List_AngleRange) != 0) {
    PyErr_SetString(PyExc_ValueError, "can only specify 'energy_range_eV' or 'angle_range', but not both");
    return NULL;
  } else if (PyList_Size(List_EnergyPoints_eV) != 0 && PyList_Size(List_AnglePoints) != 0) {
    PyErr_SetString(PyExc_ValueError, "cannot specify both energy and angle lists");
    return NULL;
  } else if (PyList_Size(List_EnergyRange_eV) != 0 && NPoints == 0) {
    PyErr_SetString(PyExc_ValueError, "must specify 'npoints' > 0");
    return NULL;
  } else if (PyList_Size(List_AngleRange) != 0 && NPoints == 0) {
    PyErr_SetString(PyExc_ValueError, "must specify 'npoints' > 0");
    return NULL;
  }

  // Container for spectrum
  TSpectrumContainer SpectrumContainer;

  // Grap points if they are there
  std::vector<double> VEnergyPoints_eV;
  for (int i = 0; i < PyList_Size(List_EnergyPoints_eV); ++i) {
    VEnergyPoints_eV.push_back(PyFloat_AsDouble(PyList_GetItem(List_EnergyPoints_eV, i)));
  }
  std::vector<double> VAnglePoints;
  for (int i = 0; i < PyList_Size(List_AnglePoints); ++i) {
    VAnglePoints.push_back(PyFloat_AsDouble(PyList_GetItem(List_AnglePoints, i)));
  }

  // Select on inputs and calculate
  if (PyList_Size(List_EnergyRange_eV) != 0 && NPoints > 0) {
    TVector2D const EnergyRange_eV = OSCARSPY::ListAsTVector2D(List_EnergyRange_eV);
    SpectrumContainer.Init(NPoints, EnergyRange_eV[0], EnergyRange_eV[1]);
    if (AngleIntegrated) {
      self->obj->DipoleSpectrumEnergyAngleIntegrated(BField, SpectrumContainer);
    } else {
      self->obj->DipoleSpectrumEnergy(BField, SpectrumContainer, Angle);
    }
  } else if (PyList_Size(List_EnergyPoints_eV) != 0) {
    SpectrumContainer.Init(VEnergyPoints_eV);
  } else if (PyList_Size(List_AngleRange) != 0 && NPoints > 0) {
    TVector2D const AngleRange = OSCARSPY::ListAsTVector2D(List_AngleRange);
    SpectrumContainer.Init(NPoints, AngleRange[0], AngleRange[1]);
    self->obj->DipoleSpectrumAngle(BField, SpectrumContainer, Energy_eV);
  } else if (PyList_Size(List_AnglePoints) != 0 && Energy_eV > 0) {
    SpectrumContainer.Init(VAnglePoints);
  } else if (Energy_eV > 0 && (NPoints == 0 || NPoints == 1)) {
    SpectrumContainer.Init(1, Energy_eV, Energy_eV);
    if (AngleIntegrated) {
      self->obj->DipoleSpectrumEnergyAngleIntegrated(BField, SpectrumContainer);
    } else {
      self->obj->DipoleSpectrumEnergy(BField, SpectrumContainer, Angle);
    }
  } else {
    PyErr_SetString(PyExc_ValueError, "Incorrect combination of or missing input parameters.  Please see documentation for this function");
    return NULL;
  }



  TVector2D EnergyRange_eV;
  TVector2D AngleRange;

  // Check for energy range
  if (PyList_Size(List_EnergyRange_eV) != 0) {
    try {
      EnergyRange_eV = OSCARSPY::ListAsTVector2D(List_EnergyRange_eV);
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'energy_range_eV'");
      return NULL;
    }
  }

  // Check for angle range
  if (PyList_Size(List_AngleRange) != 0) {
    try {
      AngleRange = OSCARSPY::ListAsTVector2D(List_AngleRange);
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'angle_range'");
      return NULL;
    }
  }

  // Write to file if output is requested
  if (std::string(OutFileNameText) != "") {
    SpectrumContainer.WriteToFileText(OutFileNameText);
  }
  if (std::string(OutFileNameBinary) != "") {
    SpectrumContainer.WriteToFileBinary(OutFileNameBinary);
  }

  // Return the spectrum
  return OSCARSPY::GetSpectrumAsList(SpectrumContainer);
}





/*
const char* DOC_OSCARSTH_DipoleSpectrumPoint = "Get the spectrum from ideal dipole field.";
static PyObject* OSCARSTH_DipoleSpectrumPoint (OSCARSTHObject* self, PyObject* args, PyObject* keywds)
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
  static const char *kwlist[] = {"bfield", "beam_energy_GeV", "angle", "energy_eV", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, keywds, "dddd", const_cast<char **>(kwlist), &BField, &BeamEnergy, &Angle, &Energy_eV)) {
    return NULL;
  }

  // CHeck if beam is ok
  if (!self->obj->CheckBeam()) {
    PyErr_SetString(PyExc_ValueError, "particle beam not correctly defined");
    return NULL;
  }

  // Check that beam energy makes sense
  if (BeamEnergy <= 0) {
    PyErr_SetString(PyExc_ValueError, "'beam_energy_GeV' must be > 0");
    return NULL;
  }

  //TVector2D const EnergyRange = OSCARSPY::ListAsTVector2D(List_EnergyRange);
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
*/





const char* DOC_OSCARSTH_DipoleCriticalEnergy = R"docstring(
dipole_critical_energy(bfield)

Get the critical energy for bending magnet in [eV]

Parameters
----------
bfield : float
    Magnetic field of dipole in [T]

Returns
-------
energy : float
    Dipole critical energy in [eV]
)docstring";
static PyObject* OSCARSTH_DipoleCriticalEnergy (OSCARSTHObject* self, PyObject* args, PyObject* keywds)
{
  // Require 1 arguments
  double BField = 0;

  // Input variables and parsing
  static const char *kwlist[] = {"bfield",
                                 NULL};

  if (!PyArg_ParseTupleAndKeywords(args, keywds, "d",
                                   const_cast<char **>(kwlist),
                                   &BField)) {
    return NULL;
  }

  // CHeck if beam is ok
  if (!self->obj->CheckBeam()) {
    PyErr_SetString(PyExc_ValueError, "particle beam not correctly defined");
    return NULL;
  }

  // Calculate the spectrum
  double const Result = self->obj->DipoleCriticalEnergy(BField);

  return Py_BuildValue("d", Result);
}






const char* DOC_OSCARSTH_DipoleCriticalWavelength = R"docstring(
dipole_critical_wavelength(bfield)

Get the critical wavelength for bending magnet in [m]

Parameters
----------
bfield : float
    Magnetic field of dipole in [T]

Returns
-------
energy : float
    Dipole critical wavelength in [m]
)docstring";
static PyObject* OSCARSTH_DipoleCriticalWavelength (OSCARSTHObject* self, PyObject* args, PyObject* keywds)
{

  // Require 1 arguments
  double BField = 0;

  // Input variables and parsing
  static const char *kwlist[] = {"bfield",
                                 NULL};

  if (!PyArg_ParseTupleAndKeywords(args, keywds, "d",
                                   const_cast<char **>(kwlist),
                                   &BField)) {
    return NULL;
  }

  // CHeck if beam is ok
  if (!self->obj->CheckBeam()) {
    PyErr_SetString(PyExc_ValueError, "particle beam not correctly defined");
    return NULL;
  }

  // Calculate the spectrum
  double const Result = TOSCARSSR::EvToWavelength(self->obj->DipoleCriticalEnergy(BField));

  return Py_BuildValue("d", Result);
}





const char* DOC_OSCARSTH_DipoleBrightness = "Not Implemented";
static PyObject* OSCARSTH_DipoleBrightness (OSCARSTHObject* self, PyObject* args, PyObject* keywds)
{
  // Return a list of points corresponding to the flux in a given energy range for a given vertical angle.
  // This approximation assumes that the particle beam is perpendicular to the magnetic field

  // Require 2 arguments
  double BField = 0;
  double Energy_eV = 0;
  PyObject* List_EnergyRange_eV = 0x0;
  int    NPoints = 0;

  // Input variables and parsing
  static const char *kwlist[] = {"bfield",
                                 "energy_eV",
                                 "energy_range_eV",
                                 "npoints",
                                 NULL};
  if (!PyArg_ParseTupleAndKeywords(args, keywds, "d|dOi",
                                                 const_cast<char **>(kwlist),
                                                 &BField,
                                                 &Energy_eV,
                                                 &List_EnergyRange_eV,
                                                 &NPoints)) {
    return NULL;
  }

  // Container for spectrum
  TSpectrumContainer SpectrumContainer;

  TVector2D const EnergyRange_eV = List_EnergyRange_eV != 0x0 ? OSCARSPY::ListAsTVector2D(List_EnergyRange_eV) : TVector2D(0, 0);
  if (List_EnergyRange_eV != 0x0 && PyList_Size(List_EnergyRange_eV) != 0 && NPoints > 0) {
    if (EnergyRange_eV[0] >= EnergyRange_eV[1] || EnergyRange_eV[0] <= 1 || EnergyRange_eV[1] <= 0) {
      PyErr_SetString(PyExc_ValueError, "'energy_range_eV' is incorrect");
      return NULL;
    } else {
      SpectrumContainer.Init(NPoints, EnergyRange_eV[0], EnergyRange_eV[1]);
    }
  } else if (Energy_eV > 0) {
    SpectrumContainer.Init(1, Energy_eV, Energy_eV);
  } else {
    PyErr_SetString(PyExc_ValueError, "Incorrect combination of or missing input parameters.  Please see documentation for this function");
    return NULL;
  }


  self->obj->DipoleBrightness(BField, SpectrumContainer);

  // Return the spectrum
  return OSCARSPY::GetSpectrumAsList(SpectrumContainer);
}
















const char* DOC_OSCARSTH_UndulatorFlux = R"docstring(
undulator_flux(period, nperiods, harmonic, [, bfield_range, K_range, npoints, bfield_points, K_points, minimum, ofile, bofile])

Get the on-axis flux for an ideal undulator given K for a specific harmonic.  Should specify either K or bfield, but not both.  You *must* have previously defined a beam including the emittance and beta function.

Parameters
----------
period : float
    Undulator period length [m]

nperiods : int
    Number of periods

harmonic : int
    Harmonic number of interest

bfield_range : list
    [min, max] of bfield range of interest

K_range : list
    [min, max] of K range of interest

npoints : int
    number of points to use when bfield_range or K_range is specified

bfield_points : list
    List of bfield points to calculate flux at

K_points : list
    List of K points to calculate flux at

minimum : float
    Any flux below the minimum will not be included in the return list

ofile : str
    Output file name

bofile : str
    Binary output file name

Returns
-------
[energy_eV, flux]s : list[[float, float], ...]
    Photon energy [eV] and flux [photons/s/0.1%bw/mrad^2] for the given parameters

Examples
--------
Get the harmonic peak spectrum for a undulator with a period of 0.050 [m] having 41 periods and a bfield ranging from 0.05 to 0.8 [T] from the 1st harmonic

    >>> oth.undulator_flux(period=0.050, nperiods=41, harmonic=1, bfield_range=[0.05, 0.8], npoints=1000)

Get the photon energy and flux for a undulator with a period of 0.050 [m] having 41 periods and K value of 1.8674577 from the 1st harmonic

    >>> oth.undulator_flux(period=0.050, nperiods=41, harmonic=1, K=1.8674577)
)docstring";
static PyObject* OSCARSTH_UndulatorFlux (OSCARSTHObject* self, PyObject* args, PyObject* keywds)
{
  // Return a list of points corresponding to the flux in a given energy range for a given vertical angle.
  // This approximation assumes that the particle beam is perpendicular to the magnetic field

  // Require 2 arguments
  double       Period            = 0;
  int          NPeriods          = 0;
  int          Harmonic          = 0;
  PyObject*    List_BFieldRange  = PyList_New(0);
  PyObject*    List_KRange       = PyList_New(0);
  int          NPoints           = 0;
  PyObject*    List_BFieldPoints = PyList_New(0);
  PyObject*    List_KPoints      = PyList_New(0);
  double       Minimum           = 0;
  const char*  OutFileNameText   = "";
  const char*  OutFileNameBinary = "";

  // Input variable list
  static const char *kwlist[] = {"period",
                                 "nperiods",
                                 "harmonic",
                                 "bfield_range",
                                 "K_range",
                                 "npoints",
                                 "bfield_points",
                                 "K_points",
                                 "minimum",
                                 "ofile",
                                 "bofile",
                                 NULL};

  // Parse inputs
  if (!PyArg_ParseTupleAndKeywords(args, keywds, "dii|OOiOOdss",
                                   const_cast<char **>(kwlist),
                                   &Period,
                                   &NPeriods,
                                   &Harmonic,
                                   &List_BFieldRange,
                                   &List_KRange,
                                   &NPoints,
                                   &List_BFieldPoints,
                                   &List_KPoints,
                                   &Minimum,
                                   &OutFileNameText,
                                   &OutFileNameBinary)) {
    return NULL;
  }

  // CHeck if beam is ok
  if (!self->obj->CheckBeam()) {
    PyErr_SetString(PyExc_ValueError, "particle beam not correctly defined");
    return NULL;
  }

  // Check harmonic number is > 0
  if (Harmonic <= 0) {
    PyErr_SetString(PyExc_ValueError, "'harmonic' must be > 0");
    return NULL;
  }

  // Check period
  if (Period <= 0) {
    PyErr_SetString(PyExc_ValueError, "'period' must be > 0");
    return NULL;
  }

  // Check nperiods
  if (NPeriods <= 0) {
    PyErr_SetString(PyExc_ValueError, "'nperiod' must be > 0");
    return NULL;
  }

  // Check not overlapping definitions
  int const SizeSumLists =  PyList_Size(List_BFieldRange)
                             + PyList_Size(List_KRange)
                             + PyList_Size(List_BFieldPoints)
                             + PyList_Size(List_KPoints);
  if (!(PyList_Size(List_BFieldRange)  == SizeSumLists ||
        PyList_Size(List_KRange)       == SizeSumLists ||
        PyList_Size(List_BFieldPoints) == SizeSumLists ||
        PyList_Size(List_KPoints)      == SizeSumLists)) {
    PyErr_SetString(PyExc_ValueError, "May only specify one of: 'bfield_range', 'K_range', 'bfield_points', 'K_points'");
    return NULL;
  }

  // Container for spectrum
  TSpectrumContainer SpectrumContainer;

  if (Harmonic % 2 != 1) {
    return OSCARSPY::GetSpectrumAsList(SpectrumContainer);
  }

  TVector2D Range;
  std::vector<double> Points;

  // Init container based on inputs
  if (PyList_Size(List_BFieldRange) > 0 && NPoints > 1) {
    try {
      Range = OSCARSPY::ListAsTVector2D(List_BFieldRange);
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'bfield_range'");
      return NULL;
    }

    // Add each point
    for (int i = 0; i < NPoints; ++i) {
      double BField = Range[0] + (Range[1] - Range[0]) / (double) (NPoints - 1) * (double) i;
      TVector2D const Result = self->obj->UndulatorFluxB(BField, Period, NPeriods, Harmonic);
      if (Result[1] >= Minimum) {
        SpectrumContainer.AddPoint(Result[0], Result[1]);
      }
    }

  } else if (PyList_Size(List_KRange) > 0 && NPoints > 1) {
    try {
      Range = OSCARSPY::ListAsTVector2D(List_KRange);
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'K_range'");
      return NULL;
    }

    // Add each point
    for (int i = 0; i < NPoints; ++i) {
      double K = Range[0] + (Range[1] - Range[0]) / (double) (NPoints - 1) * (double) i;
      TVector2D const Result = self->obj->UndulatorFluxK(K, Period, NPeriods, Harmonic);
      if (Result[1] >= Minimum) {
        SpectrumContainer.AddPoint(Result[0], Result[1]);
      }
    }

  } else if (PyList_Size(List_BFieldPoints) > 0) {
    for (int i = 0; i < PyList_Size(List_BFieldPoints); ++i) {
      Points.push_back(PyFloat_AsDouble(PyList_GetItem(List_BFieldPoints, i)));
    }

    // Add each point
    for (size_t i = 0; i < Points.size(); ++i) {
      TVector2D const Result = self->obj->UndulatorFluxB(Points[i], Period, NPeriods, Harmonic);
      if (Result[1] >= Minimum) {
        SpectrumContainer.AddPoint(Result[0], Result[1]);
      }
    }

  } else if (PyList_Size(List_KPoints) > 0) {
    for (int i = 0; i < PyList_Size(List_KPoints); ++i) {
      Points.push_back(PyFloat_AsDouble(PyList_GetItem(List_KPoints, i)));
    }

    // Add each point
    for (size_t i = 0; i < Points.size(); ++i) {
      TVector2D const Result = self->obj->UndulatorFluxK(Points[i], Period, NPeriods, Harmonic);
      if (Result[1] >= Minimum) {
        SpectrumContainer.AddPoint(Result[0], Result[1]);
      }
    }

  } else {
    PyErr_SetString(PyExc_ValueError, "Incorrect input format in input, possibly check that npoints > 1?");
    return NULL;
  }

  // Don't need to hold on to this any longer
  Points.clear();

  // Write to file if output is requested
  if (std::string(OutFileNameText) != "") {
    SpectrumContainer.WriteToFileText(OutFileNameText);
  }
  if (std::string(OutFileNameBinary) != "") {
    SpectrumContainer.WriteToFileBinary(OutFileNameBinary);
  }

  // Return the spectrum
  return OSCARSPY::GetSpectrumAsList(SpectrumContainer);
}












const char* DOC_OSCARSTH_UndulatorFluxOnAxis = R"docstring(
undulator_flux_onaxis(period, nperiods, harmonic, [, bfield_range, K_range, npoints, bfield_points, K_points, minimum, ofile, bofile])

Get the on-axis flux for an ideal undulator given K for a specific harmonic.  Should specify either K or bfield, but not both.  You *must* have previously defined a beam.

Parameters
----------
period : float
    Undulator period length [m]

nperiods : int
    Number of periods

harmonic : int
    Harmonic number of interest

bfield_range : list
    [min, max] of bfield range of interest

K_range : list
    [min, max] of K range of interest

npoints : int
    number of points to use when bfield_range or K_range is specified

bfield_points : list
    List of bfield points to calculate flux at

K_points : list
    List of K points to calculate flux at

minimum : float
    Any flux below the minimum will not be included in the return list

ofile : str
    Output file name

bofile : str
    Binary output file name

Returns
-------
[energy_eV, flux]s : list[[float, float], ...]
    Photon energy [eV] and flux [photons/s/0.1%bw/mrad^2] for the given parameters

Examples
--------
Get the harmonic peak spectrum for a undulator with a period of 0.050 [m] having 41 periods and a bfield ranging from 0.05 to 0.8 [T] from the 1st harmonic

    >>> oth.undulator_flux_onaxis(period=0.050, nperiods=41, harmonic=1, bfield_range=[0.05, 0.8], npoints=1000)

Get the photon energy and flux for a undulator with a period of 0.050 [m] having 41 periods and K value of 1.8674577 from the 1st harmonic

    >>> oth.undulator_flux_onaxis(period=0.050, nperiods=41, harmonic=1, K=1.8674577)
)docstring";
static PyObject* OSCARSTH_UndulatorFluxOnAxis (OSCARSTHObject* self, PyObject* args, PyObject* keywds)
{
  // Return a list of points corresponding to the flux in a given energy range for a given vertical angle.
  // This approximation assumes that the particle beam is perpendicular to the magnetic field

  // Require 2 arguments
  double       Period            = 0;
  int          NPeriods          = 0;
  int          Harmonic          = 0;
  PyObject*    List_BFieldRange  = PyList_New(0);
  PyObject*    List_KRange       = PyList_New(0);
  int          NPoints           = 0;
  PyObject*    List_BFieldPoints = PyList_New(0);
  PyObject*    List_KPoints      = PyList_New(0);
  double       Minimum           = 0;
  const char*  OutFileNameText   = "";
  const char*  OutFileNameBinary = "";

  // Input variable list
  static const char *kwlist[] = {"period",
                                 "nperiods",
                                 "harmonic",
                                 "bfield_range",
                                 "K_range",
                                 "npoints",
                                 "bfield_points",
                                 "K_points",
                                 "minimum",
                                 "ofile",
                                 "bofile",
                                 NULL};

  // Parse inputs
  if (!PyArg_ParseTupleAndKeywords(args, keywds, "dii|OOiOOdss",
                                   const_cast<char **>(kwlist),
                                   &Period,
                                   &NPeriods,
                                   &Harmonic,
                                   &List_BFieldRange,
                                   &List_KRange,
                                   &NPoints,
                                   &List_BFieldPoints,
                                   &List_KPoints,
                                   &Minimum,
                                   &OutFileNameText,
                                   &OutFileNameBinary)) {
    return NULL;
  }

  // CHeck if beam is ok
  if (!self->obj->CheckBeam()) {
    PyErr_SetString(PyExc_ValueError, "particle beam not correctly defined");
    return NULL;
  }

  // Check harmonic number is > 0
  if (Harmonic <= 0) {
    PyErr_SetString(PyExc_ValueError, "'harmonic' must be > 0");
    return NULL;
  }

  // Check period
  if (Period <= 0) {
    PyErr_SetString(PyExc_ValueError, "'period' must be > 0");
    return NULL;
  }

  // Check nperiods
  if (NPeriods <= 0) {
    PyErr_SetString(PyExc_ValueError, "'nperiod' must be > 0");
    return NULL;
  }

  // Check not overlapping definitions
  int const SizeSumLists =  PyList_Size(List_BFieldRange)
                             + PyList_Size(List_KRange)
                             + PyList_Size(List_BFieldPoints)
                             + PyList_Size(List_KPoints);
  if (!(PyList_Size(List_BFieldRange)  == SizeSumLists ||
        PyList_Size(List_KRange)       == SizeSumLists ||
        PyList_Size(List_BFieldPoints) == SizeSumLists ||
        PyList_Size(List_KPoints)      == SizeSumLists)) {
    PyErr_SetString(PyExc_ValueError, "May only specify one of: 'bfield_range', 'K_range', 'bfield_points', 'K_points'");
    return NULL;
  }

  // Container for spectrum
  TSpectrumContainer SpectrumContainer;

  if (Harmonic % 2 != 1) {
    return OSCARSPY::GetSpectrumAsList(SpectrumContainer);
  }

  TVector2D Range;
  std::vector<double> Points;

  // Init container based on inputs
  if (PyList_Size(List_BFieldRange) > 0 && NPoints > 1) {
    try {
      Range = OSCARSPY::ListAsTVector2D(List_BFieldRange);
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'bfield_range'");
      return NULL;
    }

    // Add each point
    for (int i = 0; i < NPoints; ++i) {
      double BField = Range[0] + (Range[1] - Range[0]) / (double) (NPoints - 1) * (double) i;
      TVector2D const Result = self->obj->UndulatorFluxOnAxisB(BField, Period, NPeriods, Harmonic);
      if (Result[1] >= Minimum) {
        SpectrumContainer.AddPoint(Result[0], Result[1]);
      }
    }

  } else if (PyList_Size(List_KRange) > 0 && NPoints > 1) {
    try {
      Range = OSCARSPY::ListAsTVector2D(List_KRange);
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'K_range'");
      return NULL;
    }

    // Add each point
    for (int i = 0; i < NPoints; ++i) {
      double K = Range[0] + (Range[1] - Range[0]) / (double) (NPoints - 1) * (double) i;
      TVector2D const Result = self->obj->UndulatorFluxOnAxisK(K, Period, NPeriods, Harmonic);
      if (Result[1] >= Minimum) {
        SpectrumContainer.AddPoint(Result[0], Result[1]);
      }
    }

  } else if (PyList_Size(List_BFieldPoints) > 0) {
    for (int i = 0; i < PyList_Size(List_BFieldPoints); ++i) {
      Points.push_back(PyFloat_AsDouble(PyList_GetItem(List_BFieldPoints, i)));
    }

    // Add each point
    for (size_t i = 0; i < Points.size(); ++i) {
      TVector2D const Result = self->obj->UndulatorFluxOnAxisB(Points[i], Period, NPeriods, Harmonic);
      if (Result[1] >= Minimum) {
        SpectrumContainer.AddPoint(Result[0], Result[1]);
      }
    }

  } else if (PyList_Size(List_KPoints) > 0) {
    for (int i = 0; i < PyList_Size(List_KPoints); ++i) {
      Points.push_back(PyFloat_AsDouble(PyList_GetItem(List_KPoints, i)));
    }

    // Add each point
    for (size_t i = 0; i < Points.size(); ++i) {
      TVector2D const Result = self->obj->UndulatorFluxOnAxisK(Points[i], Period, NPeriods, Harmonic);
      if (Result[1] >= Minimum) {
        SpectrumContainer.AddPoint(Result[0], Result[1]);
      }
    }

  } else {
    PyErr_SetString(PyExc_ValueError, "Incorrect input format in input, possibly check that npoints > 1?");
    return NULL;
  }

  // Don't need to hold on to this any longer
  Points.clear();

  // Write to file if output is requested
  if (std::string(OutFileNameText) != "") {
    SpectrumContainer.WriteToFileText(OutFileNameText);
  }
  if (std::string(OutFileNameBinary) != "") {
    SpectrumContainer.WriteToFileBinary(OutFileNameBinary);
  }

  // Return the spectrum
  return OSCARSPY::GetSpectrumAsList(SpectrumContainer);
}





const char* DOC_OSCARSTH_UndulatorBrightness = R"docstring(
undulator_brightness(period, nperiods, harmonic, [, bfield_range, K_range, npoints, bfield_points, K_points, minimum, ofile, bofile])

Get the brightness for an ideal undulator given K for a specific harmonic.  Should specify either K or bfield, but not both.  You must have previously defined a beam, including the beta and emittance values.

Parameters
----------
period : float
    Undulator period length [m]

nperiods : int
    Number of periods

harmonic : int
    Harmonic number of interest

bfield_range : list
    [min, max] of bfield range of interest

K_range : list
    [min, max] of K range of interest

npoints : int
    number of points to use when bfield_range or K_range is specified

bfield_points : list
    List of bfield points to calculate brightness at

K_points : list
    List of K points to calculate brightness at

minimum : float
    Any brightness below the minimum will not be included in the return list

ofile : str
    Output file name

bofile : str
    Binary output file name

Returns
-------
[energy_eV, brightness]s : list[[float, float], ...]
    Photon energy [eV] and brightness [photons/s/0.1%bw/mrad^2/mm] for the given parameters

Examples
--------
Get the harmonic peak spectrum for a undulator with a period of 0.050 [m] having 41 periods and a bfield ranging from 0.05 to 0.8 [T] from the 1st harmonic

    >>> oth.undulator_flux_onaxis(period=0.050, nperiods=41, harmonic=1, bfield_range=[0.05, 0.8], npoints=1000)

Get the photon energy and flux for a undulator with a period of 0.050 [m] having 41 periods and K value of 1.8674577 from the 1st harmonic

    >>> oth.undulator_flux_onaxis(period=0.050, nperiods=41, harmonic=1, K=1.8674577)
)docstring";
static PyObject* OSCARSTH_UndulatorBrightness (OSCARSTHObject* self, PyObject* args, PyObject* keywds)
{
  // Return a list of points corresponding to the brightness.

  // Require 2 arguments
  double      Period            = 0;
  int         NPeriods          = 0;
  int         Harmonic          = 0;
  PyObject*   List_BFieldRange  = PyList_New(0);
  PyObject*   List_KRange       = PyList_New(0);
  int         NPoints           = 0;
  PyObject*   List_BFieldPoints = PyList_New(0);
  PyObject*   List_KPoints      = PyList_New(0);
  double      Minimum           = 0;
  const char* OutFileNameText   = "";
  const char* OutFileNameBinary = "";

  // Input variable list
  static const char *kwlist[] = {"period",
                                 "nperiods",
                                 "harmonic",
                                 "bfield_range",
                                 "K_range",
                                 "npoints",
                                 "bfield_points",
                                 "K_points",
                                 "minimum",
                                 "ofile",
                                 "bofile",
                                 NULL};

  // Parse inputs
  if (!PyArg_ParseTupleAndKeywords(args, keywds, "dii|OOiOOdss",
                                   const_cast<char **>(kwlist),
                                   &Period,
                                   &NPeriods,
                                   &Harmonic,
                                   &List_BFieldRange,
                                   &List_KRange,
                                   &NPoints,
                                   &List_BFieldPoints,
                                   &List_KPoints,
                                   &Minimum,
                                   &OutFileNameText,
                                   &OutFileNameBinary)) {
    return NULL;
  }

  // CHeck if beam is ok
  if (!self->obj->CheckBeam()) {
    PyErr_SetString(PyExc_ValueError, "particle beam not correctly defined");
    return NULL;
  }

  // Check harmonic number is > 0
  if (Harmonic <= 0) {
    PyErr_SetString(PyExc_ValueError, "'harmonic' must be > 0");
    return NULL;
  }

  // Check period
  if (Period <= 0) {
    PyErr_SetString(PyExc_ValueError, "'period' must be > 0");
    return NULL;
  }

  // Check nperiods
  if (NPeriods <= 0) {
    PyErr_SetString(PyExc_ValueError, "'nperiod' must be > 0");
    return NULL;
  }

  // Check not overlapping definitions
  int const SizeSumLists =  PyList_Size(List_BFieldRange)
                            + PyList_Size(List_KRange)
                            + PyList_Size(List_BFieldPoints)
                            + PyList_Size(List_KPoints);
  if (!(PyList_Size(List_BFieldRange)  == SizeSumLists ||
        PyList_Size(List_KRange)       == SizeSumLists ||
        PyList_Size(List_BFieldPoints) == SizeSumLists ||
        PyList_Size(List_KPoints)      == SizeSumLists)) {
    PyErr_SetString(PyExc_ValueError, "May only specify one of: 'bfield_range', 'K_range', 'bfield_points', 'K_points'");
    return NULL;
  }

  // Container for spectrum
  TSpectrumContainer SpectrumContainer;

  TVector2D Range;
  std::vector<double> Points;

  // Init container based on inputs
  try {
    if (PyList_Size(List_BFieldRange) > 0 && NPoints > 1) {
      try {
        Range = OSCARSPY::ListAsTVector2D(List_BFieldRange);
      } catch (std::length_error e) {
        PyErr_SetString(PyExc_ValueError, "Incorrect format in 'bfield_range'");
        return NULL;
      }

      // Add each point
      for (int i = 0; i < NPoints; ++i) {
        double BField = Range[0] + (Range[1] - Range[0]) / (double) (NPoints - 1) * (double) i;
        TVector2D const Result = self->obj->UndulatorBrightnessB(BField, Period, NPeriods, Harmonic);
        if (Result[1] >= Minimum) {
          SpectrumContainer.AddPoint(Result[0], Result[1]);
        }
      }

    } else if (PyList_Size(List_KRange) > 0 && NPoints > 1) {
      try {
        Range = OSCARSPY::ListAsTVector2D(List_KRange);
      } catch (std::length_error e) {
        PyErr_SetString(PyExc_ValueError, "Incorrect format in 'K_range'");
        return NULL;
      }

      // Add each point
      for (int i = 0; i < NPoints; ++i) {
        double K = Range[0] + (Range[1] - Range[0]) / (double) (NPoints - 1) * (double) i;
        TVector2D const Result = self->obj->UndulatorBrightnessK(K, Period, NPeriods, Harmonic);
        if (Result[1] >= Minimum) {
          SpectrumContainer.AddPoint(Result[0], Result[1]);
        }
      }

    } else if (PyList_Size(List_BFieldPoints) > 0) {
      for (int i = 0; i < PyList_Size(List_BFieldPoints); ++i) {
        Points.push_back(PyFloat_AsDouble(PyList_GetItem(List_BFieldPoints, i)));
      }

      // Add each point
      for (size_t i = 0; i < Points.size(); ++i) {
        TVector2D const Result = self->obj->UndulatorBrightnessB(Points[i], Period, NPeriods, Harmonic);
        if (Result[1] >= Minimum) {
          SpectrumContainer.AddPoint(Result[0], Result[1]);
        }
      }

    } else if (PyList_Size(List_KPoints) > 0) {
      for (int i = 0; i < PyList_Size(List_KPoints); ++i) {
        Points.push_back(PyFloat_AsDouble(PyList_GetItem(List_KPoints, i)));
      }

      // Add each point
      for (size_t i = 0; i < Points.size(); ++i) {
        TVector2D const Result = self->obj->UndulatorBrightnessK(Points[i], Period, NPeriods, Harmonic);
        if (Result[1] >= Minimum) {
          SpectrumContainer.AddPoint(Result[0], Result[1]);
        }
      }

    } else {
      PyErr_SetString(PyExc_ValueError, "Incorrect input format in input, possibly check that npoints > 1?");
      return NULL;
    }
  } catch (std::invalid_argument e) {
    PyErr_SetString(PyExc_ValueError, e.what());
    return NULL;
  }

  // Don't need to hold on to this any longer
  Points.clear();

  // Write to file if output is requested
  if (std::string(OutFileNameText) != "") {
    SpectrumContainer.WriteToFileText(OutFileNameText);
  }
  if (std::string(OutFileNameBinary) != "") {
    SpectrumContainer.WriteToFileBinary(OutFileNameBinary);
  }

  // Return the spectrum
  return OSCARSPY::GetSpectrumAsList(SpectrumContainer);
}





const char* DOC_OSCARSTH_UndulatorEnergyHarmonic = R"docstring(
undulator_energy_harmonic(period, harmonic [, K, bfield])

Get the undulator photon energy at particular harmonic.  Specify either K or bfield, but not both.

period : float
    Undulator period length [m]

harmonic : int
    Harmonic number of interest

K : float
    Undulator deflection parameter

bfield : float
    Magnetic field of interest [T]

Returns
-------
energy : float
    Photon energy [eV]

Examples
--------
Find the energy of the 5th harmonic for an undulator of period 0.050 [m] and magnetic field of 0.4 [T]

    >>> oth.undulator_energy_harmonic(period=0.050, harmonic=5, bfield=0.4)

Find the energy of the 5th harmonic for an undulator of period 0.050 [m] and K value of 1.8674577

    >>> oth.undulator_energy_harmonic(period=0.050, harmonic=5, K=1.8674577)
)docstring";
static PyObject* OSCARSTH_UndulatorEnergyHarmonic (OSCARSTHObject* self, PyObject* args, PyObject* keywds)
{
  // Return the brightness [gamma/s/mrad^2/mm^2/0.1%bw] at a given K for a given harmonic

  // Require some arguments
  double Period = 0;
  int    Harmonic = 0;
  double K = 0;
  double BField = 0;

  // Input variables and parsing
  static const char *kwlist[] = {"period",
                                 "harmonic",
                                 "K",
                                 "bfield",
                                 NULL};

  if (!PyArg_ParseTupleAndKeywords(args, keywds, "di|dd",
                                   const_cast<char **>(kwlist),
                                   &Period,
                                   &Harmonic,
                                   &K,
                                   &BField
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





const char* DOC_OSCARSTH_WigglerSpectrum = R"docstring(
wiggler_spectrum(bfield, period, length [, energy_range_eV, energy_points_eV, energy_eV, angle_integrated, angle_range, angle_points, angle, npoints, minimum, ofile, bofile])

Get the spectrum from ideal dipole field.  One can calculate the spectrum as a function of photon energy at a given angle or the angular dependence for a particular energy.

Parameters
----------
bfield : float
    Magnetic field of the dipole in [T]

period : float
    Length of magnetic period in [m]

length : float
    Length of device in [m]

energy_range_eV : list
    [min, max] photon energy of interest in [eV]

energy_points_eV : list
    List of energy points in [eV]

energy_eV : float
    Photon energy of interest

angle_integrated : boolean
    Set to true if you want the vertical angle integrated

angle_range : list
    [min, max] angle of interest in [rad]

angle_points : list
    List of angle points in [rad]

angle : float
    Angle of interest in [rad]

npoints : int
    Number of points for _range requests

ofile : str
    Output file name

bofile : str
    Binary output file name

Returns
-------
flux : list
    A list of flux values for given points [[p0, f0], [p1, f1], ...].  The points can be in energy or angle depending on the input parameters.  Output is in [photons / s / mrad^2 / 0.1\%bw]

Examples
--------
Calculate the spectrum in a given energy range

    >>> oth.wiggler_spectrum(bfield=0.4, period=0.100, length=7.0, energy_range_eV=[10, 30000], npoints=1000)
)docstring";
static PyObject* OSCARSTH_WigglerSpectrum (OSCARSTHObject* self, PyObject* args, PyObject* keywds)
{
  // Return a list of points corresponding to the flux in a given energy range for a given vertical angle.
  // This approximation assumes that the particle beam is perpendicular to the magnetic field

  // Require 2 arguments
  double      BField               = 0;
  double      Period               = 0;
  double      Length               = 0;
  PyObject*   List_EnergyRange_eV  = PyList_New(0);
  PyObject*   List_EnergyPoints_eV = PyList_New(0);
  double      Energy_eV            = 0;
  bool        AngleIntegrated      = false;
  PyObject*   List_AngleRange      = PyList_New(0);
  PyObject*   List_AnglePoints     = PyList_New(0);
  double      Angle                = 0;
  int         NPoints              = 500;
  const char* OutFileNameText      = "";
  const char* OutFileNameBinary    = "";

  // Input variable list
  static const char *kwlist[] = {"bfield",
                                 "period",
                                 "length",
                                 "energy_range_eV",
                                 "energy_points_eV",
                                 "energy_eV",
                                 "angle_integrated",
                                 "angle_range",
                                 "angle_points",
                                 "angle",
                                 "npoints",
                                 "ofile",
                                 "bofile",
                                 NULL};

  // Parse inputs
  if (!PyArg_ParseTupleAndKeywords(args, keywds, "ddd|OOdpOOdiss",
                                   const_cast<char **>(kwlist),
                                   &BField,
                                   &Period,
                                   &Length,
                                   &List_EnergyRange_eV,
                                   &List_EnergyPoints_eV,
                                   &Energy_eV,
                                   &AngleIntegrated,
                                   &List_AngleRange,
                                   &List_AnglePoints,
                                   &Angle,
                                   &NPoints,
                                   &OutFileNameText,
                                   &OutFileNameBinary)) {
    return NULL;
  }

  // CHeck if beam is ok
  if (!self->obj->CheckBeam()) {
    PyErr_SetString(PyExc_ValueError, "particle beam not correctly defined");
    return NULL;
  }

  // Must have a bfield
  if (fabs(BField) == 0) {
    PyErr_SetString(PyExc_ValueError, "'bfield' must not be zero");
    return NULL;
  }

  // Must have a period
  if (Period <= 0) {
    PyErr_SetString(PyExc_ValueError, "'period' must not be > 0");
    return NULL;
  }

  // Must have a length greater than or equal to 1/2 period
  if (Length < Period / 2.) {
    PyErr_SetString(PyExc_ValueError, "'Length' must be at least one half-period");
    return NULL;
  }

  // Number of periods as a double
  double const NPeriods = (double) (int) (Length / Period);


  // Must have some combination that makes sense
  if (PyList_Size(List_EnergyRange_eV) != 0 && PyList_Size(List_AngleRange) != 0) {
    PyErr_SetString(PyExc_ValueError, "can only specify 'energy_range_eV' or 'angle_range', but not both");
    return NULL;
  } else if (PyList_Size(List_EnergyPoints_eV) != 0 && PyList_Size(List_AnglePoints) != 0) {
    PyErr_SetString(PyExc_ValueError, "cannot specify both energy and angle lists");
    return NULL;
  } else if (PyList_Size(List_EnergyRange_eV) != 0 && NPoints == 0) {
    PyErr_SetString(PyExc_ValueError, "must specify 'npoints' > 0");
    return NULL;
  } else if (PyList_Size(List_AngleRange) != 0 && NPoints == 0) {
    PyErr_SetString(PyExc_ValueError, "must specify 'npoints' > 0");
    return NULL;
  }

  // Container for spectrum
  TSpectrumContainer SpectrumContainer;

  // Grap points if they are there
  std::vector<double> VEnergyPoints_eV;
  for (int i = 0; i < PyList_Size(List_EnergyPoints_eV); ++i) {
    VEnergyPoints_eV.push_back(PyFloat_AsDouble(PyList_GetItem(List_EnergyPoints_eV, i)));
  }
  std::vector<double> VAnglePoints;
  for (int i = 0; i < PyList_Size(List_AnglePoints); ++i) {
    VAnglePoints.push_back(PyFloat_AsDouble(PyList_GetItem(List_AnglePoints, i)));
  }

  // Select on inputs and calculate
  if (PyList_Size(List_EnergyRange_eV) != 0 && NPoints > 0) {
    TVector2D const EnergyRange_eV = OSCARSPY::ListAsTVector2D(List_EnergyRange_eV);
    SpectrumContainer.Init(NPoints, EnergyRange_eV[0], EnergyRange_eV[1]);
    if (AngleIntegrated) {
      self->obj->DipoleSpectrumEnergyAngleIntegrated(BField, SpectrumContainer);
    } else {
      self->obj->DipoleSpectrumEnergy(BField, SpectrumContainer, Angle);
    }
  } else if (PyList_Size(List_EnergyPoints_eV) != 0) {
    SpectrumContainer.Init(VEnergyPoints_eV);
  } else if (PyList_Size(List_AngleRange) != 0 && NPoints > 0) {
    TVector2D const AngleRange = OSCARSPY::ListAsTVector2D(List_AngleRange);
    SpectrumContainer.Init(NPoints, AngleRange[0], AngleRange[1]);
    self->obj->DipoleSpectrumAngle(BField, SpectrumContainer, Energy_eV);
  } else if (PyList_Size(List_AnglePoints) != 0 && Energy_eV > 0) {
    SpectrumContainer.Init(VAnglePoints);
  } else {
    PyErr_SetString(PyExc_ValueError, "Incorrect combination of or missing input parameters.  Please see documentation for this function");
    return NULL;
  }

  // Insert factor for Wiggler
  SpectrumContainer.Scale(2. * NPeriods);


  TVector2D EnergyRange_eV;
  TVector2D AngleRange;

  // Check for energy range
  if (PyList_Size(List_EnergyRange_eV) != 0) {
    try {
      EnergyRange_eV = OSCARSPY::ListAsTVector2D(List_EnergyRange_eV);
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'energy_range_eV'");
      return NULL;
    }
  }

  // Check for angle range
  if (PyList_Size(List_AngleRange) != 0) {
    try {
      AngleRange = OSCARSPY::ListAsTVector2D(List_AngleRange);
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'angle_range'");
      return NULL;
    }
  }

  // Write to file if output is requested
  if (std::string(OutFileNameText) != "") {
    SpectrumContainer.WriteToFileText(OutFileNameText);
  }
  if (std::string(OutFileNameBinary) != "") {
    SpectrumContainer.WriteToFileBinary(OutFileNameBinary);
  }

  // Return the spectrum
  return OSCARSPY::GetSpectrumAsList(SpectrumContainer);
}





const char* DOC_OSCARSTH_WigglerFluxRectangle = R"docstring(
wiggler_flux_rectangle(plane, energy_eV, period, nperiods, width, npoints [, bfield, K, ofile, bofile])

Get the flux for an ideal wiggler according to R. P. Walker XXX XXX.  Should specify either K or bfield, but not both.  You *must* have previously defined a beam.

Parameters
----------
period : float
    Magnetic period length [m]

nperiods : int
    Number of periods

harmonic : int
    Harmonic number of interest

bfield_range : list
    [min, max] of bfield range of interest

K_range : list
    [min, max] of K range of interest

npoints : int
    number of points to use when bfield_range or K_range is specified

bfield_points : list
    List of bfield points to calculate flux at

K_points : list
    List of K points to calculate flux at

minimum : float
    Any flux below the minimum will not be included in the return list

ofile : str
    Output file name

bofile : str
    Binary output file name

Returns
-------
[energy_eV, flux]s : list[[float, float], ...]
    Photon energy [eV] and flux [photons/s/0.1%bw/mrad^2] for the given parameters

)docstring";
static PyObject* OSCARSTH_WigglerFluxRectangle (OSCARSTHObject* self, PyObject* args, PyObject* keywds)
{
  // Return a list of points corresponding to the flux in a given energy range for a given vertical angle.
  // This approximation assumes that the particle beam is perpendicular to the magnetic field

  // Require 2 arguments
  char const*  SurfacePlane      = "";
  double       Energy_eV         = 0;
  double       Period            = 0;
  int          NPeriods          = 0;
  PyObject*    List_NPoints      = PyList_New(0);
  PyObject*    List_Width        = PyList_New(0);
  PyObject*    List_X0X1X2       = PyList_New(0);
  double       BField            = 0;
  double       K                 = 0;
  PyObject*    List_Rotations    = PyList_New(0);
  PyObject*    List_Translation  = PyList_New(0);
  int          NormalDirection   = 0;
  int          Dim               = 2;
  int          NThreads          = 0;
  int          GPU               = -1;
  const char*  OutFileNameText   = "";
  const char*  OutFileNameBinary = "";

  size_t NX1 = 0;
  size_t NX2 = 0;

  // Input variable list
  static const char *kwlist[] = {
                                 "plane",
                                 "energy_eV",
                                 "period",
                                 "nperiods",
                                 "npoints",
                                 "width",
                                 "x0x1x2",
                                 "bfield",
                                 "K",
                                 "rotations",
                                 "translation",
                                 "normal",
                                 "dim",
                                 "nthreads",
                                 "gpu",
                                 "ofile",
                                 "bofile",
                                 NULL};

  // Parse inputs
  if (!PyArg_ParseTupleAndKeywords(args, keywds, "sddiO|OOddOOiiiiss",
                                   const_cast<char **>(kwlist),
                                   &SurfacePlane,
                                   &Energy_eV,
                                   &Period,
                                   &NPeriods,
                                   &List_NPoints,
                                   &List_Width,
                                   &List_X0X1X2,
                                   &BField,
                                   &K,
                                   &List_Rotations,
                                   &List_Translation,
                                   &NormalDirection,
                                   &Dim,
                                   &NThreads,
                                   &GPU,
                                   &OutFileNameText,
                                   &OutFileNameBinary)) {
    return NULL;
  }

  // CHeck if beam is ok
  if (!self->obj->CheckBeam()) {
    PyErr_SetString(PyExc_ValueError, "particle beam not correctly defined");
    return NULL;
  }

  // Check period
  if (Period <= 0) {
    PyErr_SetString(PyExc_ValueError, "'period' must be > 0");
    return NULL;
  }

  // Check nperiods
  if (NPeriods <= 0) {
    PyErr_SetString(PyExc_ValueError, "'nperiod' must be > 0");
    return NULL;
  }

  TVector2D Width;
  try {
    Width = OSCARSPY::ListAsTVector2D(List_Width);
  } catch (...) {
    PyErr_SetString(PyExc_ValueError, "'width' has incorrect format");
    return NULL;
  }

  if (PyList_Size(List_NPoints) == 2) {
    // NPoints in [m]
    NX1 = PyLong_AsSsize_t(PyList_GetItem(List_NPoints, 0));
    NX2 = PyLong_AsSsize_t(PyList_GetItem(List_NPoints, 1));
  } else {
    PyErr_SetString(PyExc_ValueError, "'npoints' must be [int, int]");
    return NULL;
  }

  // Check BField and K
  if (!((BField != 0) ^ (K != 0))) {
    std::cout << BField << std::endl;
    std::cout << K << std::endl;
    PyErr_SetString(PyExc_ValueError, "Must specify one and only one of: 'bfield' or 'K'");
    return NULL;
  }

  // Vectors for rotations and translations.  Default to 0
  TVector3D Rotations(0, 0, 0);
  TVector3D Translation(0, 0, 0);

  // Check for Rotations in the input
  if (PyList_Size(List_Rotations) != 0) {
    try {
      Rotations = OSCARSPY::ListAsTVector3D(List_Rotations);
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'rotations'");
      return NULL;
    }
  }


  // Check for Translation in the input
  if (PyList_Size(List_Translation) != 0) {
    try {
      Translation = OSCARSPY::ListAsTVector3D(List_Translation);
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'translation'");
      return NULL;
    }
  }

  // Check normal
  if (abs(NormalDirection) > 1) {
    PyErr_SetString(PyExc_ValueError, "'normal' must be -1, 0, or 1");
    return NULL;
  }

  // Check dim
  if (Dim != 2 && Dim != 3) {
    PyErr_SetString(PyExc_ValueError, "'dim' must be 2 or 3");
    return NULL;
  }

  // Check GPU parameter
  if (GPU != 0 && GPU != 1 && GPU != -1) {
    PyErr_SetString(PyExc_ValueError, "'gpu' must be 0 or 1");
    return NULL;
  }

  // Check NThreads parameter
  if (NThreads < 0) {
    PyErr_SetString(PyExc_ValueError, "'nthreads' must be > 0");
    return NULL;
  }

  // Check you are not trying to use threads and GPU
  if (NThreads > 0 && GPU == 1) {
    PyErr_SetString(PyExc_ValueError, "gpu is 1 and nthreads > 0.  Both are not currently allowed.");
    return NULL;
  }


  // The rectangular surface object we'll use
  TSurfacePoints_Rectangle Surface;

  // If you are requesting a simple surface plane, check that you have widths
  if (std::strlen(SurfacePlane) != 0 && Width[0] > 0 && Width[1] > 0) {
    try {
      Surface.Init(SurfacePlane, (int) NX1, (int) NX2, Width[0], Width[1], Rotations, Translation, NormalDirection);
    } catch (std::invalid_argument e) {
      PyErr_SetString(PyExc_ValueError, e.what());
      return NULL;
    }
  } else if (PyList_Size(List_X0X1X2) != 0) {
    std::vector<TVector3D> X0X1X2;
    if (PyList_Size(List_X0X1X2) == 3) {
      for (int i = 0; i != 3; ++i) {
        PyObject* List_X = PyList_GetItem(List_X0X1X2, i);

        try {
          X0X1X2.push_back(OSCARSPY::ListAsTVector3D(List_X));
        } catch (std::length_error e) {
          PyErr_SetString(PyExc_ValueError, "Incorrect format in 'x0x1x2'");
          return NULL;
        }
      }
    } else {
      PyErr_SetString(PyExc_ValueError, "'x0x1x2' must have 3 XYZ points defined correctly");
      return NULL;
    }

    for (std::vector<TVector3D>::iterator it = X0X1X2.begin(); it != X0X1X2.end(); ++it) {
      it->RotateSelfXYZ(Rotations);
      *it += Translation;
    }

    // UPDATE: Check for orthogonality
    Surface.Init((int) NX1, (int) NX2, X0X1X2[0], X0X1X2[1], X0X1X2[2], NormalDirection);
  }

  // Container for Flux
  T3DScalarContainer FluxContainer;

  //self->obj->CalculateFluxWiggler(Surface, Energy_eV, FluxContainer, Polarization, Angle, HorizontalDirection, PropogationDirection, NParticles, NThreads, GPU, Dim);

  if (BField != 0) {
    self->obj->WigglerFluxB(BField, Period, NPeriods, Surface, Energy_eV, FluxContainer, NThreads, GPU);
  } else if (K != 0) {
    self->obj->WigglerFluxK(K,      Period, NPeriods, Surface, Energy_eV, FluxContainer, NThreads, GPU);
  }

  // Write to file if output is requested
  if (std::string(OutFileNameText) != "") {
    FluxContainer.WriteToFileText(OutFileNameText, Dim);
  }
  if (std::string(OutFileNameBinary) != "") {
    FluxContainer.WriteToFileBinary(OutFileNameBinary, Dim);
  }

  // Build the output list of: [[[x, y, z], Flux], [...]]
  // Create a python list
  PyObject *PList = PyList_New(0);

  size_t const NPoints = FluxContainer.GetNPoints();

  for (size_t i = 0; i != NPoints; ++i) {
    T3DScalar P = FluxContainer.GetPoint(i);

    // Inner list for each point
    PyObject *PList2 = PyList_New(0);


    // Add position and value to list
    PyList_Append(PList2, OSCARSPY::TVector3DAsList(P.GetX()));
    PyList_Append(PList2, Py_BuildValue("f", P.GetV()));
    PyList_Append(PList, PList2);

  }

  return PList;
}












































const char* DOC_OSCARSTH_BesselJ = R"docstring(
bessel_j(nu, x)

Get the value of the bessel function of the first kind J_nu for integer nu

Parameters
----------
nu : int

x : float

Returns
-------
besselj : float
    Value of the bessel function for nu at x
)docstring";
static PyObject* OSCARSTH_BesselJ (OSCARSTHObject* self, PyObject* args, PyObject* keywds)
{
  // Return Bessel J_nu(x)

  // Require 2 arguments
  int    Nu = 0;
  double X  = 0;

  // Input variables and parsing
  static const char *kwlist[] = {"nu",
                                 "x",
                                 NULL};

  if (!PyArg_ParseTupleAndKeywords(args, keywds, "id",
                                   const_cast<char **>(kwlist),
                                   &Nu,
                                   &X)) {
    return NULL;
  }


  // Return the internal OSCARSTH number constant pi
  return Py_BuildValue("d", TOMATH::BesselJ(Nu, X));
}










const char* DOC_OSCARSTH_BesselK = R"docstring(
bessel_k(nu, x)

Get the value of the modified bessel function of the second kind K_nu

Parameters
----------
nu : float

x : float

Returns
-------
besselk : float
    Value of the bessel function for nu at x
)docstring";
static PyObject* OSCARSTH_BesselK (OSCARSTHObject* self, PyObject* args, PyObject* keywds)
{
  // Return Bessel K_nu(x)

  // Require 2 arguments
  double Nu = 0;
  double X = 0;

  // Input variables and parsing
  static const char *kwlist[] = {"nu",
                                 "x",
                                 NULL};

  if (!PyArg_ParseTupleAndKeywords(args, keywds, "dd",
                                   const_cast<char **>(kwlist),
                                   &Nu,
                                   &X)) {
    return NULL;
  }


  // Return the internal OSCARSTH number constant pi
  return Py_BuildValue("d", TOMATH::BesselK(Nu, X));
}







const char* DOC_OSCARSTH_SetParticleBeam = R"docstring(
set_particle_beam([, type, name, energy_GeV, d0, x0, beam, sigma_energy_GeV, t0, current, weight, rotations, translation, horizontal_direction, beta, alpha, gamma, emittance, eta, lattice_reference, mass, charge])

This function is the same as add_particle_beam(), but it clears all particle beams before the 'add'.
)docstring";
static PyObject* OSCARSTH_SetParticleBeam (OSCARSTHObject* self, PyObject* args, PyObject* keywds)
{
  // Clear all particle beams, add this beam, and set a new particle

  self->obj->ClearParticleBeams();

  PyObject* ret = OSCARSTH_AddParticleBeam(self, args, keywds);
  if (ret == NULL) {
    return ret;
  }

  //self->obj->SetNewParticle("", "ideal");

  return ret;
}




const char* DOC_OSCARSTH_AddParticleBeam = R"docstring(
add_particle_beam([, type, name, energy_GeV, d0, x0, beam, sigma_energy_GeV, t0, current, weight, rotations, translation, horizontal_direction, beta, alpha, gamma, emittance, eta, lattice_reference, mass, charge])

Add a particle beam to the OSCARS object with a name given by *name*.  There is no limit to the number of different particle beams one can add.  They are added with a *weight* which is by default 1.  The weight is used in random sampling when asking for a new particle, for example in set_new_particle().  If the *beam* parameter is given you only need to specify *name* and *x0*.

Supported particle types for *type* are:
    * electron

Parameters
----------

type : str
    One of the built-in types of particle beams, or 'custom'.  If you use custom you must also specify *mass* and *charge*.

name : str
    User identified of this beam

energy_GeV : float
    Beam energy in [GeV]

d0 : list
    Vector [float, float, float] representing the default direction of the beam at the initial position.  The normalization of this vector does not matter.

x0 : list
    Coordinates of the initial position in [m] [x, y, z]

beam : str
    Name of predefined beam

sigma_energy_GeV : float
    Beam energy spread in [GeV]

t0 : float
    Initial time in [m] at the initial_position.  Time here is in terms of ct.

current : float
    Beam current in [A].  If this parameter is 0, the current defaults to the single particle charge

weight : float
    Weight to give this beam for random sampling when there are multiple beams

rotations : list, optional
    3-element list representing rotations around x, y, and z axes: [:math:`\theta_x, \theta_y, \theta_z`]

translation : list, optional
    3-element list representing a translation in space [x, y, z]

horizontal_direction : list
    A list representing the *horizontal* beam direction.  This must be perpendicular to the beam *direction*.  The vertical direction is defined internally as :math:`[direction] \times [horizontal\_direction]`.

beta : list
    Values of the horizontal and vertical beta funtion at the point *lattice_center* [beta_x, beta_y]

alpha : list
    Values of the horizontal and vertical alpha funtion at the point *lattice_center* [alpha_x, alpha_y]

gamma : list
    Values of the horizontal and vertical gamma funtion at the point *lattice_center* [gamma_x, gamma_y]

emittance : list
    values of the horizontal and vertical emittance [emittance_x, emittance_y]

eta : list
    Values of the horizontal and vertical dispersion at the point *lattice_center* [eta_x, eta_y].  Currently eta_y is ignored and dispersion is only used in oscars.th calculations

lattice_reference : list
    Coordinates of the lattice center [x, y, z] (must be on-axis with respect to the beam)

mass : float
    mass of a *custom* particle.  This is only used if *type* = 'custom'.  Must be non-zero.

charge : float
    Charge of a *custom* particle.  This is only used if *type* = 'custom'.  Must be non-zero.

Returns
-------
None

Currently the predefined beam types are:
* NSLSII
* NSLSII-ShortStraight
* NSLSII-LongStraight

Examples
--------
Add an electron beam with 0.500 [A] current at an initial position of [0, 0, 0] in the Y direction with energy of 3 [GeV]

    >>> oth.add_particle_beam(type='electron', energy_GeV=3, current=0.500)

Add a predefined beam for the NSLSII short straight section

    >>> oth.add_particle_beam(beam='NSLSII-ShortStraight')
)docstring";
static PyObject* OSCARSTH_AddParticleBeam (OSCARSTHObject* self, PyObject* args, PyObject* keywds)
{
  // Add a particle beam to the experiment

  // Lists and variables some with initial values
  char const* Type                       = "electron";
  char const* Name                       = "";
  double      Energy_GeV                 = -1;
  double      Sigma_Energy_GeV           = 0;
  double      T0                         = 0;
  double      Current                    = 0;
  double      Weight                     = 1;
  double      Mass                       = 0;
  double      Charge                     = 0;
  char const* Beam                       = "";
  PyObject*   List_Position              = 0x0;
  PyObject*   List_Direction             = 0x0;
  PyObject*   List_Rotations             = 0x0;
  PyObject*   List_Translation           = 0x0;
  PyObject*   List_Horizontal_Direction  = 0x0;
  PyObject*   List_Beta                  = 0x0;
  PyObject*   List_Alpha                 = 0x0;
  PyObject*   List_Gamma                 = 0x0;
  PyObject*   List_Emittance             = 0x0;
  PyObject*   List_Eta                   = 0x0;
  PyObject*   List_Lattice_Reference     = 0x0;

  TVector3D Position(0, 0, 0);
  TVector3D Direction(0, 0, 1);
  TVector3D Rotations(0, 0, 0);
  TVector3D Translation(0, 0, 0);
  TVector3D Horizontal_Direction;
  TVector2D Beta(0, 0);
  TVector2D Alpha(0, 0);
  TVector2D Gamma(0, 0);
  TVector2D Emittance(0, 0);
  TVector3D Lattice_Reference(0, 0, 0);
  TVector2D Eta(0, 0);


  // Input variables and parsing
  static const char *kwlist[] = {"type",
                                 "name",
                                 "energy_GeV",
                                 "d0",
                                 "x0",
                                 "beam",
                                 "sigma_energy_GeV",
                                 "t0",
                                 "current",
                                 "weight",
                                 "rotations",
                                 "translation",
                                 "horizontal_direction",
                                 "beta",
                                 "alpha",
                                 "gamma",
                                 "emittance",
                                 "eta",
                                 "lattice_reference",
                                 "mass",
                                 "charge",
                                 NULL};

  if (!PyArg_ParseTupleAndKeywords(args, keywds, "|ssdOOsddddOOOOOOOOOdd",
                                   const_cast<char **>(kwlist),
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
                                   &List_Alpha,
                                   &List_Gamma,
                                   &List_Emittance,
                                   &List_Eta,
                                   &List_Lattice_Reference,
                                   &Mass,
                                   &Charge)) {
    return NULL;
  }

  // Clear particle beam (only one beam in TH)
  self->obj->ClearParticleBeams();

  // Are you asking for one of the predefined beams?
  bool const HasPredefinedBeam = std::strlen(Beam) != 0 ? true : false;

  TParticleBeam* ThisBeam = 0x0;

  // Check if beam is defined (for predefined beams)
  if (HasPredefinedBeam) {
    try {
      ThisBeam = &(self->obj->AddParticleBeam(Beam, Name, Weight));
    } catch (...) {
      PyErr_SetString(PyExc_ValueError, "Error in predefined beam name / definition");
      return NULL;
    }

    if (Energy_GeV >= 0) {
      ThisBeam->SetE0(Energy_GeV);
    }

  }


  // default is 0
  if (Energy_GeV == -1) {
    Energy_GeV = 0;
  }

  // Initial position
  if (List_Position != 0x0) {
    try {
      Position = OSCARSPY::ListAsTVector3D(List_Position);
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'x0'");
      return NULL;
    }

  }

  // Energy check
  if (Energy_GeV < 0) {
    PyErr_SetString(PyExc_ValueError, "We do not currently support negative energy beams.  Please change 'energy_GeV' to >= 0");
    return NULL;
  }

  // Initial direction
  if (List_Direction != 0x0) {
    try {
      Direction = OSCARSPY::ListAsTVector3D(List_Direction);
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'd0'");
      return NULL;
    }
  }

  // Check for Rotations in the input
  if (List_Rotations != 0x0) {
    try {
      Rotations = OSCARSPY::ListAsTVector3D(List_Rotations);
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'rotations'");
      return NULL;
    }
  }


  // Check for Translation in the input
  if (List_Translation != 0x0) {
    try {
      Translation = OSCARSPY::ListAsTVector3D(List_Translation);
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'translation'");
      return NULL;
    }
  }


  // Check for Horizontal_Direction in the input
  if (List_Horizontal_Direction != 0x0) {
    try {
      Horizontal_Direction = OSCARSPY::ListAsTVector3D(List_Horizontal_Direction);
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'horizontal_direction'");
      return NULL;
    }
  } else {
    Horizontal_Direction = -Direction.Orthogonal();
  }
  Horizontal_Direction = Horizontal_Direction.UnitVector();


  // Rotate beam parameters
  Position.RotateSelfXYZ(Rotations);
  Direction.RotateSelfXYZ(Rotations);
  Position += Translation;

  // Add the particle beam
  if (std::strlen(Beam) == 0) {
    try {
      if (std::string(Type) == "custom") {
        if (Mass == 0 || Charge == 0) {
          PyErr_SetString(PyExc_ValueError, "'mass' or 'charge' is zero");
          return NULL;
        }
        // UPDATE: for custom beams
        ThisBeam = &(self->obj->AddParticleBeam(Type, Name, Position, Direction, Energy_GeV, T0, Current, Weight, Charge, Mass));
      } else {
        ThisBeam = &(self->obj->AddParticleBeam(Type, Name, Position, Direction, Energy_GeV, T0, Current, Weight));
      }
    } catch (std::invalid_argument e) {
      PyErr_SetString(PyExc_ValueError, "invalid argument in adding particle beam.  possibly 'name' already exists");
      return NULL;
    }
  }

  // Change predefined beam accordingly
  if (HasPredefinedBeam) {
    if (List_Direction != 0x0) {
      ThisBeam->SetU0(Direction);
    }
    ThisBeam->SetX0(Position);
  }

  // Set horizontal direction
  ThisBeam->SetHorizontalDirection(Horizontal_Direction);

  // UPDATE
  // Check for Emittance in the input
  if (List_Emittance != 0x0) {
    try {
      Emittance = OSCARSPY::ListAsTVector2D(List_Emittance);
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'emittance'");
      return NULL;
    }

    // Set emittance and beam distribution
    ThisBeam->SetEmittance(Emittance);
    ThisBeam->SetBeamDistribution(TParticleBeam::kBeamDistribution_Gaussian);
  } else {
    // Beam distribution to filament
    if (std::strlen(Beam) == 0) {
      ThisBeam->SetBeamDistribution(TParticleBeam::kBeamDistribution_Filament);
    }
  }

  // Check for no beam
  if (ThisBeam == 0x0) {
    throw std::logic_error("logic error. I do not see a beam at this point.  please report");
  }

  if (Sigma_Energy_GeV == 0) {
    // Do nothing.  zero energy diff is alright
  } else if (Sigma_Energy_GeV < 0) {
    PyErr_SetString(PyExc_ValueError, "'sigma_energy_GeV' cannot be less than zero");
    return NULL;
  } else {
    ThisBeam->SetSigmaEnergyGeV(Sigma_Energy_GeV);
  }

  // Check for beta, alpha, gammain the input
  int HasBAG = 0x0;
  if (List_Beta != 0x0) {
    try {
      Beta = OSCARSPY::ListAsTVector2D(List_Beta);
      HasBAG |= 0x4;
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'beta'");
      return NULL;
    }
  }
  if (List_Alpha != 0x0) {
    try {
      Alpha = OSCARSPY::ListAsTVector2D(List_Alpha);
      HasBAG |= 0x2;
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'alpha'");
      return NULL;
    }
  }
  if (List_Gamma != 0x0) {
    try {
      Gamma = OSCARSPY::ListAsTVector2D(List_Gamma);
      HasBAG |= 0x1;
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'gamma'");
      return NULL;
    }
  }

  // Check for Lattice reference in the input
  bool HasReferencePoint = false;
  if (List_Lattice_Reference != 0x0) {
    try {
      Lattice_Reference = OSCARSPY::ListAsTVector3D(List_Lattice_Reference);
      HasReferencePoint = true;
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'lattice_reference'");
      return NULL;
    }
  }
  ThisBeam->SetTwissLatticeReference(Lattice_Reference);

  // Set correct twiss parameters
  switch (HasBAG) {
    case 0x7:
      ThisBeam->SetTwissParameters(Beta, Alpha, Gamma, Lattice_Reference, HasReferencePoint);
      break;
    case 0x6:
      ThisBeam->SetTwissBetaAlpha(Beta, Alpha, Lattice_Reference, HasReferencePoint);
      break;
    case 0x5:
      ThisBeam->SetTwissBetaGamma(Beta, Gamma, Lattice_Reference, HasReferencePoint);
      break;
    case 0x3:
      ThisBeam->SetTwissAlphaGamma(Alpha, Gamma, Lattice_Reference, HasReferencePoint);
      break;
    case 0x4:
      ThisBeam->SetTwissBetaAlpha(Beta, TVector2D(0, 0), Lattice_Reference, HasReferencePoint);
      break;
    default:
      break;
  }

  if (List_Eta != 0x0) {
    try {
      Eta = OSCARSPY::ListAsTVector2D(List_Eta);
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'gamma'");
      return NULL;
    }
  }
  ThisBeam->SetEta(Eta);


  if (T0 != 0) {
    ThisBeam->SetT0(T0);
  }

  // Should set weight on input
  //self->obj->GetParticleBeam(Name).SetWeight(Weight);

  // Must return python object None in a special way
  Py_INCREF(Py_None);
  return Py_None;
}



















const char* DOC_OSCARSTH_PrintParticleBeams = R"docstring(
print_particle_beams()

Print information about what is stored in the th object

Returns
-------
None
)docstring";
const char* DOC_OSCARSTH_PrintAll = R"docstring(
print_all()

Print information about what is stored in the th object

Returns
-------
None
)docstring";
static PyObject* OSCARSTH_PrintAll (OSCARSTHObject* self)
{
  // Print all particle beams stored

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

































static PyObject* OSCARSTH_Fake (OSCARSTHObject* self, PyObject* args, PyObject *keywds)
{
    PyErr_SetString(PyExc_RuntimeError, "You must create an object to use this function: oth = oscars.th.th()");
    return NULL;
}


static PyMethodDef OSCARSTH_methods[] = {
  // We must tell python about the function we allow access as well as give them nice
  // python names, and tell python the method of input parameters.

  {"version",                                    (PyCFunction) OSCARSTH_Version,                                 METH_NOARGS,                                   DOC_OSCARSTH_Version},
  {"undulator_K",                                (PyCFunction) OSCARSTH_UndulatorK,                              METH_VARARGS | METH_KEYWORDS,                  DOC_OSCARSTH_UndulatorK},
  {"undulator_bfield",                           (PyCFunction) OSCARSTH_UndulatorBField,                         METH_VARARGS | METH_KEYWORDS,                  DOC_OSCARSTH_UndulatorBField},
  {"undulator_period",                           (PyCFunction) OSCARSTH_UndulatorPeriod,                         METH_VARARGS | METH_KEYWORDS,                  DOC_OSCARSTH_UndulatorPeriod},

  {"dipole_spectrum",                            (PyCFunction) OSCARSTH_DipoleSpectrum,                          METH_VARARGS | METH_KEYWORDS,                  DOC_OSCARSTH_DipoleSpectrum},
  //{"dipole_spectrum_point",                      (PyCFunction) OSCARSTH_DipoleSpectrumPoint,                       METH_VARARGS | METH_KEYWORDS,                  DOC_OSCARSTH_DipoleSpectrumPoint},
  {"dipole_critical_energy",                     (PyCFunction) OSCARSTH_DipoleCriticalEnergy,                    METH_VARARGS | METH_KEYWORDS,                  DOC_OSCARSTH_DipoleCriticalEnergy},
  {"dipole_critical_wavelength",                 (PyCFunction) OSCARSTH_DipoleCriticalWavelength,                METH_VARARGS | METH_KEYWORDS,                  DOC_OSCARSTH_DipoleCriticalWavelength},
  {"dipole_brightness",                          (PyCFunction) OSCARSTH_DipoleBrightness,                        METH_VARARGS | METH_KEYWORDS,                  DOC_OSCARSTH_DipoleBrightness},

  {"undulator_flux",                             (PyCFunction) OSCARSTH_UndulatorFlux,                           METH_VARARGS | METH_KEYWORDS,                  DOC_OSCARSTH_UndulatorFlux},
  {"undulator_flux_onaxis",                      (PyCFunction) OSCARSTH_UndulatorFluxOnAxis,                     METH_VARARGS | METH_KEYWORDS,                  DOC_OSCARSTH_UndulatorFluxOnAxis},
  {"undulator_brightness",                       (PyCFunction) OSCARSTH_UndulatorBrightness,                     METH_VARARGS | METH_KEYWORDS,                  DOC_OSCARSTH_UndulatorBrightness},
  {"undulator_energy_harmonic",                  (PyCFunction) OSCARSTH_UndulatorEnergyHarmonic,                 METH_VARARGS | METH_KEYWORDS,                  DOC_OSCARSTH_UndulatorEnergyHarmonic},

  {"wiggler_spectrum",                           (PyCFunction) OSCARSTH_WigglerSpectrum,                         METH_VARARGS | METH_KEYWORDS,                  DOC_OSCARSTH_WigglerSpectrum},
  {"wiggler_flux",                               (PyCFunction) OSCARSTH_WigglerFluxRectangle,                    METH_VARARGS | METH_KEYWORDS,                  DOC_OSCARSTH_WigglerFluxRectangle},

  {"bessel_j",                                   (PyCFunction) OSCARSTH_BesselJ,                                 METH_VARARGS | METH_KEYWORDS,                  DOC_OSCARSTH_BesselJ},
  {"bessel_k",                                   (PyCFunction) OSCARSTH_BesselK,                                 METH_VARARGS | METH_KEYWORDS,                  DOC_OSCARSTH_BesselK},

  {"set_particle_beam",                          (PyCFunction) OSCARSTH_SetParticleBeam,                         METH_VARARGS | METH_KEYWORDS,                  DOC_OSCARSTH_SetParticleBeam},
  {"add_particle_beam",                          (PyCFunction) OSCARSTH_SetParticleBeam,                         METH_VARARGS | METH_KEYWORDS,                  DOC_OSCARSTH_SetParticleBeam},

  {"print_all",                                  (PyCFunction) OSCARSTH_PrintAll,                                METH_NOARGS,                                   DOC_OSCARSTH_PrintAll},
  {"print_particle_beams",                       (PyCFunction) OSCARSTH_PrintAll,                                METH_NOARGS,                                   DOC_OSCARSTH_PrintParticleBeams},

  {NULL}  /* Sentinel */
};


static PyMethodDef OSCARSTH_methods_fake[] = {
  // We must tell python about the function we allow access as well as give them nice
  // python names, and tell python the method of input parameters.

  {"version",                                    (PyCFunction) OSCARSTH_Version, METH_NOARGS,                                DOC_OSCARSTH_Version},
  {"undulator_K",                                (PyCFunction) OSCARSTH_Fake, METH_VARARGS | METH_KEYWORDS,                  DOC_OSCARSTH_UndulatorK},
  {"undulator_bfield",                           (PyCFunction) OSCARSTH_Fake, METH_VARARGS | METH_KEYWORDS,                  DOC_OSCARSTH_UndulatorBField},
  {"undulator_period",                           (PyCFunction) OSCARSTH_Fake, METH_VARARGS | METH_KEYWORDS,                  DOC_OSCARSTH_UndulatorPeriod},

  {"dipole_spectrum",                            (PyCFunction) OSCARSTH_Fake, METH_VARARGS | METH_KEYWORDS,                  DOC_OSCARSTH_DipoleSpectrum},
  //{"dipole_spectrum_point",                      (PyCFunction) OSCARSTH_Fake, METH_VARARGS | METH_KEYWORDS,                  DOC_OSCARSTH_DipoleSpectrumPoint},
  {"dipole_critical_energy",                     (PyCFunction) OSCARSTH_Fake, METH_VARARGS | METH_KEYWORDS,                  DOC_OSCARSTH_DipoleCriticalEnergy},
  {"dipole_critical_wavelength",                 (PyCFunction) OSCARSTH_Fake, METH_VARARGS | METH_KEYWORDS,                  DOC_OSCARSTH_DipoleCriticalWavelength},
  {"dipole_brightness",                          (PyCFunction) OSCARSTH_Fake, METH_VARARGS | METH_KEYWORDS,                  DOC_OSCARSTH_DipoleBrightness},

  {"undulator_flux",                             (PyCFunction) OSCARSTH_Fake, METH_VARARGS | METH_KEYWORDS,                  DOC_OSCARSTH_UndulatorFlux},
  {"undulator_flux_onaxis",                      (PyCFunction) OSCARSTH_Fake, METH_VARARGS | METH_KEYWORDS,                  DOC_OSCARSTH_UndulatorFluxOnAxis},
  {"undulator_brightness",                       (PyCFunction) OSCARSTH_Fake, METH_VARARGS | METH_KEYWORDS,                  DOC_OSCARSTH_UndulatorBrightness},
  {"undulator_energy_harmonic",                  (PyCFunction) OSCARSTH_Fake, METH_VARARGS | METH_KEYWORDS,                  DOC_OSCARSTH_UndulatorEnergyHarmonic},

  {"wiggler_spectrum",                           (PyCFunction) OSCARSTH_Fake, METH_VARARGS | METH_KEYWORDS,                  DOC_OSCARSTH_WigglerSpectrum},
  {"wiggler_flux",                               (PyCFunction) OSCARSTH_Fake, METH_VARARGS | METH_KEYWORDS,                  DOC_OSCARSTH_WigglerFluxRectangle},

  {"bessel_j",                                   (PyCFunction) OSCARSTH_Fake, METH_VARARGS | METH_KEYWORDS,                  DOC_OSCARSTH_BesselK},
  {"bessel_k",                                   (PyCFunction) OSCARSTH_Fake, METH_VARARGS | METH_KEYWORDS,                  DOC_OSCARSTH_BesselK},

  {"set_particle_beam",                          (PyCFunction) OSCARSTH_Fake, METH_VARARGS | METH_KEYWORDS,                  DOC_OSCARSTH_SetParticleBeam},
  {"add_particle_beam",                          (PyCFunction) OSCARSTH_Fake, METH_VARARGS | METH_KEYWORDS,                  DOC_OSCARSTH_SetParticleBeam},

  {"print_all",                                  (PyCFunction) OSCARSTH_Fake, METH_NOARGS,                                   DOC_OSCARSTH_PrintAll},
  {"print_particle_beams",                       (PyCFunction) OSCARSTH_Fake, METH_NOARGS,                                   DOC_OSCARSTH_PrintParticleBeams},

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




//static PyMethodDef module_methods[] = {
//  // I do not need
//  {NULL}  /* Sentinel */
//};


#if PY_MAJOR_VERSION >= 3
static PyModuleDef OSCARSTHmodule = {
  PyModuleDef_HEAD_INIT,
  "th",
  "OSCARSTH module extension.",
  -1,
  OSCARSTH_methods_fake,
  NULL, NULL, NULL
};
#endif




#if PY_MAJOR_VERSION >= 3
PyMODINIT_FUNC PyInit_th(void)
{
  if (PyType_Ready(&OSCARSTHType) < 0) {
    return NULL;
  }
  PyObject* m = PyModule_Create(&OSCARSTHmodule);
  if (m == NULL) {
    return NULL;
  }

  Py_INCREF(&OSCARSTHType);
  PyModule_AddObject(m, "th", (PyObject *)&OSCARSTHType);

  // Print copyright notice
  PyObject* sys = PyImport_ImportModule( "sys");
  PyObject* s_out = PyObject_GetAttrString(sys, "stdout");
  std::string Message = "OSCARS v" + OSCARSPY::GetVersionString() + " - Open Source Code for Advanced Radiation Simulation\nBrookhaven National Laboratory, Upton NY, USA\nhttp://oscars.bnl.gov\noscars@bnl.gov\n";
  PyObject_CallMethod(s_out, "write", "s", Message.c_str());

  return m;
}
#else
PyMODINIT_FUNC initth(OSCARSTHObject* self, PyObject* args, PyObject* kwds)
{
  if (PyType_Ready(&OSCARSTHType) < 0) {
    return;
  }
  PyObject *m = Py_InitModule("oscars.th", OSCARSTH_methods);
  if (m == NULL) {
    return;
  }

  Py_INCREF(&OSCARSTHType);
  PyModule_AddObject(m, "th", (PyObject *)&OSCARSTHType);

  // Print copyright notice
  PyObject* sys = PyImport_ImportModule( "sys");
  PyObject* s_out = PyObject_GetAttrString(sys, "stdout");
  std::string Message = "OSCARS v" + OSCARSPY::GetVersionString() + " - Open Source Code for Advanced Radiation Simulation\nBrookhaven National Laboratory, Upton NY, USA\nhttp://oscars.bnl.gov\noscars@bnl.gov\n";
  PyObject_CallMethod(s_out, "write", "s", Message.c_str());

  return;
}
#endif







