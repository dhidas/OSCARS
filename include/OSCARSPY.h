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



namespace OSCARSPY {

  std::string GetVersionString ();

  char* GetAsString (PyObject* S);
  char* GetVersionOfModule (std::string const&);

  PyObject* GetSpectrumAsList (TSpectrumContainer const& Spectrum);

  TVector2D ListAsTVector2D (PyObject* List);
  TVector3D ListAsTVector3D (PyObject* List);

  PyObject* TVector2DAsList (TVector2D const& V);
  PyObject* TVector3DAsList (TVector3D const& V);

}


#endif
