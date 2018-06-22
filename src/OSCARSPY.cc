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
#include "Python.h"

#include "OSCARSPY.h"

#include <stdexcept>

namespace OSCARSPY {

std::string GetVersionString ()
{
  // Get the version string based off of the compiler defines
  char ver[200];
  if (OSCARS_RELEASE == NULL) {
    sprintf(ver, "%i.%i.%i", OSCARS_VMAJOR, OSCARS_VMINOR, OSCARS_REVISION);
  } else {
    sprintf(ver, "%i.%i.%i.%s", OSCARS_VMAJOR, OSCARS_VMINOR, OSCARS_REVISION, OSCARS_RELEASE);
  }
  return std::string(ver);
}





void PyPrint_stderr (std::string const& s)
{
  PyObject* sys = PyImport_ImportModule( "sys");
  PyObject* s_out = PyObject_GetAttrString(sys, "stderr");
  PyObject_CallMethod(s_out, "write", "s", s.c_str());

  return;
}




void PyPrint_stdout (std::string const& s)
{
  PyObject* sys = PyImport_ImportModule( "sys");
  PyObject* s_out = PyObject_GetAttrString(sys, "stdout");
  PyObject_CallMethod(s_out, "write", "s", s.c_str());

  return;
}




#if PY_MAJOR_VERSION < 3
char* GetAsString (PyObject* S)
{
  return PyString_AsString(S);
}
#elif PY_MAJOR_VERSION == 3
#if PY_MINOR_VERSION <= 6
char* GetAsString (PyObject* S)
{
  return PyUnicode_AsUTF8(S);
}
#else
const char* GetAsString (PyObject* S)
{
  return PyUnicode_AsUTF8(S);
}
#endif
#endif






#if PY_MAJOR_VERSION == 3
#if PY_MINOR_VERSION >= 7
const char* GetVersionOfModule (std::string const& ModuleName)
#else
char* GetVersionOfModule (std::string const& ModuleName)
#endif
#elif PY_MAJOR_VERSION < 3
char* GetVersionOfModule (std::string const& ModuleName)
#endif
{
  PyObject* pkg_resources = PyImport_ImportModule("pkg_resources");
  if (pkg_resources == NULL) {
    throw std::invalid_argument("cannot import pkg_resources");
  }
  PyObject* dist = PyObject_CallMethod(pkg_resources, "get_distribution", "s", ModuleName.c_str());
  if (dist == NULL) {
    throw std::invalid_argument("cannot call get_distribution with this argument");
  }
  PyObject* ver = PyObject_GetAttrString(dist, "version");
  if (ver == NULL) {
    throw std::invalid_argument("cannot find version");
  }

  return GetAsString(ver);
}




PyObject* GetSpectrumAsList (TSpectrumContainer const& Spectrum)
{
  // Get the spectrum as a list format for python output

  // Create a python list
  PyObject *List = PyList_New(0);

  PyObject* Value;

  // Number of points in trajectory calculation
  size_t NSPoints = Spectrum.GetNPoints();

  // Loop over all points in trajectory
  for (size_t iS = 0; iS != NSPoints; ++iS) {
    // Create a python list for X and Beta
    PyObject *List2 = PyList_New(0);

    Value = Py_BuildValue("f", Spectrum.GetEnergy(iS));
    PyList_Append(List2, Value);
    Py_DECREF(Value);

    Value = Py_BuildValue("f", Spectrum.GetFlux(iS));
    PyList_Append(List2, Value);
    Py_DECREF(Value);

    PyList_Append(List, List2);
    Py_DECREF(List2);
  }

  // Return the python list
  return List;
}




TSpectrumContainer GetSpectrumFromList (PyObject* List)
{
  // Take an input list in spectrum format and convert it to TSpectrumContainer object

  // Increment reference for list
  Py_INCREF(List);

  // Get size of input list
  size_t const NPoints = PyList_Size(List);
  if (NPoints <= 0) {
    throw std::length_error("GetSpectrumFromList reporting no points");
  }

  TSpectrumContainer S;

  for (size_t ip = 0; ip != NPoints; ++ip) {
    PyObject* List_Point = PyList_GetItem(List, ip);
    if (PyList_Size(List_Point) == 2) {
      S.AddPoint(PyFloat_AsDouble(PyList_GetItem(List_Point, 0)), PyFloat_AsDouble(PyList_GetItem(List_Point, 1)));
    } else {
      throw std::length_error("GetSpectrumFromList reporting not 2 points");
    }
  }


  // Increment reference for list
  Py_DECREF(List);

  // Return the object
  return S;
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
  // Get a list as a TVector3D

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




void ListToVectorInt (PyObject* List, std::vector<int>& V)
{
  // Get a list as std::vector

  V.clear();
  V.resize(PyList_Size(List));
  for (int i = 0; i < PyList_Size(List); ++i) {
    V[i] = (int) PyLong_AsLong(PyList_GetItem(List, i));
  }

  return;
}




PyObject* VectorIntToList (std::vector<int>& V)
{
  // Get a vector<int> as a PyObject list

  // Create a python list
  PyObject *List = PyList_New(0);

  PyObject* Value;

  for (std::vector<int>::const_iterator it = V.begin(); it != V.end(); ++it) {
    Value = Py_BuildValue("i", *it);
    PyList_Append(List, Value);
    Py_DECREF(Value);
  }

  return List;
}




PyObject* TVector2DAsList (TVector2D const& V)
{
  // Turn a TVector2D into a list (like a vector)

  // Create a python list
  PyObject *List = PyList_New(0);

  // Value for input which can be Py_DECREF when done with
  PyObject* Value;

  Value = Py_BuildValue("f", V.GetX());
  PyList_Append(List, Value);
  Py_DECREF(Value);

  Value = Py_BuildValue("f", V.GetY());
  PyList_Append(List, Value);
  Py_DECREF(Value);

  // Return the python list
  return List;
}




PyObject* TVector3DAsList (TVector3D const& V)
{
  // Turn a TVector3D into a list (like a vector)

  // Create a python list
  PyObject *List = PyList_New(0);

  // Value for input which can be Py_DECREF when done with
  PyObject* Value;

  Value = Py_BuildValue("f", V.GetX());
  PyList_Append(List, Value);
  Py_DECREF(Value);

  Value = Py_BuildValue("f", V.GetY());
  PyList_Append(List, Value);
  Py_DECREF(Value);

  Value = Py_BuildValue("f", V.GetZ());
  PyList_Append(List, Value);
  Py_DECREF(Value);

  // Return the python list
  return List;
}







T3DScalarContainer GetT3DScalarContainerFromList (PyObject* List)
{
  // Take an input list and convert it to T3DScalarContainer object

  // Increment reference for list
  Py_INCREF(List);

  // Get size of input list
  size_t const NPoints = PyList_Size(List);
  if (NPoints <= 0) {
    throw std::length_error("GetT3DScalarContainerFromList reporting no points");
  }

  T3DScalarContainer F;

  for (size_t ip = 0; ip != NPoints; ++ip) {
    PyObject* List_Point = PyList_GetItem(List, ip);
    if (PyList_Size(List_Point) == 2) {
      F.AddPoint(OSCARSPY::ListAsTVector3D(PyList_GetItem(List_Point, 0)), PyFloat_AsDouble(PyList_GetItem(List_Point, 1)));
    } else {
      throw std::length_error("GetT3DScalarContainerFromList reporting not 2 points");
    }
  }


  // Increment reference for list
  Py_DECREF(List);

  // Return the object
  return F;
}










}
