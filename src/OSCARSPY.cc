////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Wed May 31 08:40:08 EDT 2017
//
// This is a namespace for common c-python functions used throughout
// the oscars c code.
//
////////////////////////////////////////////////////////////////////

#include "OSCARSPY.h"


#include <stdexcept>

namespace OSCARSPY {


char* GetAsString (PyObject* S)
{
  // Return the correct string version depending on the python version

  #if PY_MAJOR_VERSION >= 3
  return PyUnicode_AsUTF8(S);
  #else
  return PyString_AsString(S);
  #endif
}




PyObject* GetSpectrumAsList (TSpectrumContainer const& Spectrum)
{
  // Get the spectrum as a list format for python output

  // Create a python list
  PyObject *PList = PyList_New(0);

  // Number of points in trajectory calculation
  size_t NSPoints = Spectrum.GetNPoints();

  // Loop over all points in trajectory
  for (int iS = 0; iS != NSPoints; ++iS) {
    // Create a python list for X and Beta
    PyObject *PList2 = PyList_New(0);

    // Add position and Beta to list
    PyList_Append(PList2, Py_BuildValue("f", Spectrum.GetEnergy(iS)));
    PyList_Append(PList2, Py_BuildValue("f", Spectrum.GetFlux(iS)));
    PyList_Append(PList, PList2);
  }

  // Return the python list
  return PList;
}




TVector2D ListAsTVector2D (PyObject* List)
{
  // Get a list as a TVector2D

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




TVector3D ListAsTVector3D (PyObject* List)
{
  // Get a list as a TVector2D

  TVector3D V;
  if (PyList_Size(List) == 3) {
    Py_INCREF(List);
    V.SetXYZ(PyFloat_AsDouble(PyList_GetItem(List, 0)),
             PyFloat_AsDouble(PyList_GetItem(List, 1)),
             PyFloat_AsDouble(PyList_GetItem(List, 2)));
    Py_DECREF(List);
  } else {
    throw std::length_error("number of elements not 3");
  }

  // Return the python list
  return V;
}




PyObject* TVector2DAsList (TVector2D const& V)
{
  // Turn a TVector3D into a list (like a vector)

  // Create a python list
  PyObject *PList = PyList_New(0);

  PyList_Append(PList, Py_BuildValue("f", V.GetX()));
  PyList_Append(PList, Py_BuildValue("f", V.GetY()));

  // Return the python list
  return PList;
}




PyObject* TVector3DAsList (TVector3D const& V)
{
  // Turn a TVector3D into a list (like a vector)

  // Create a python list
  PyObject *PList = PyList_New(0);

  PyList_Append(PList, Py_BuildValue("f", V.GetX()));
  PyList_Append(PList, Py_BuildValue("f", V.GetY()));
  PyList_Append(PList, Py_BuildValue("f", V.GetZ()));

  // Return the python list
  return PList;
}










}