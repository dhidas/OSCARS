#ifndef GUARD_OSCARSPY_h
#define GUARD_OSCARSPY_h
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

#include <Python.h>

#include "Version.h"
#include "TVector2D.h"
#include "TVector3D.h"
#include "TSpectrumContainer.h"
#include "T3DScalarContainer.h"



namespace OSCARSPY {

  std::string GetVersionString ();

  void PyPrint_stderr (std::string const&);
  void PyPrint_stdout (std::string const&);

#if PY_MAJOR_VERSION == 3
#if PY_MINOR_VERSION >= 7
  const char* GetAsString (PyObject* S);
#else
  char* GetAsString (PyObject* S);
#endif
#elif PY_MAJOR_VERSION < 3
  char* GetAsString (PyObject* S);
#endif
  char* GetVersionOfModule (std::string const&);

  PyObject* GetSpectrumAsList (TSpectrumContainer const& Spectrum);
  TSpectrumContainer GetSpectrumFromList (PyObject* List);

  TVector2D ListAsTVector2D (PyObject* List);
  TVector3D ListAsTVector3D (PyObject* List);

  void      ListToVectorInt (PyObject* List, std::vector<int>& V);
  PyObject* VectorIntToList (std::vector<int>& V);

  PyObject* TVector2DAsList (TVector2D const& V);
  PyObject* TVector3DAsList (TVector3D const& V);

  T3DScalarContainer GetT3DScalarContainerFromList (PyObject* List);

}


#endif
