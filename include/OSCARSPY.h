#ifndef GUARD_OSCARSPY_h
#define GUARD_OSCARSPY_h
////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Wed May 31 08:40:08 EDT 2017
//
// This is a namespace for common c-python functions used throughout
// the oscars c code.  It also provides some common constants to the
// python layer.
//
////////////////////////////////////////////////////////////////////

#include <Python.h>

#include "Version.h"
#include "TVector2D.h"
#include "TVector3D.h"
#include "TSpectrumContainer.h"
#include "T3DScalarContainer.h"
#include "TOSCARS.h"
#include "TRandomA.h"

// External global random generator
extern TRandomA* gRandomA;

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

#if PY_MAJOR_VERSION == 3
#if PY_MINOR_VERSION >= 7
  const char* GetVersionOfModule (std::string const& ModuleName);
#else
    char* GetVersionOfModule (std::string const& ModuleName);
#endif
#elif PY_MAJOR_VERSION < 3
    char* GetVersionOfModule (std::string const& ModuleName);
#endif

  PyObject* GetSpectrumAsList (TSpectrumContainer const& Spectrum);
  TSpectrumContainer GetSpectrumFromList (PyObject* List);

  TVector2D ListAsTVector2D (PyObject* List);
  TVector3D ListAsTVector3D (PyObject* List);

  void      ListToVectorInt (PyObject* List, std::vector<int>& V);
  PyObject* VectorIntToList (std::vector<int>& V);

  PyObject* TVector2DAsList (TVector2D const& V);
  PyObject* TVector3DAsList (TVector3D const& V);

  T3DScalarContainer GetT3DScalarContainerFromList (PyObject* List);


  // Common fucntions/constants made available in python
  static PyObject* Py_Version ()  { return Py_BuildValue("s", GetVersionString().c_str());}
  static PyObject* Py_Pi()        { return Py_BuildValue("d", TOSCARS::Pi()); }
  static PyObject* Py_Alpha()     { return Py_BuildValue("d", TOSCARS::Alpha()); }
  static PyObject* Py_E()         { return Py_BuildValue("d", TOSCARS::E()); }
  static PyObject* Py_C()         { return Py_BuildValue("d", TOSCARS::C()); }
  static PyObject* Py_H()         { return Py_BuildValue("d", TOSCARS::H()); }
  static PyObject* Py_Hbar()      { return Py_BuildValue("d", TOSCARS::Hbar()); }
  static PyObject* Py_Qe()        { return Py_BuildValue("d", TOSCARS::Qe()); }
  static PyObject* Py_Me()        { return Py_BuildValue("d", TOSCARS::Me()); }
  static PyObject* Py_Epsilon0()  { return Py_BuildValue("d", TOSCARS::Epsilon0()); }
  static PyObject* Py_Mu0()       { return Py_BuildValue("d", TOSCARS::Mu0()); }

  static PyObject* Py_Random ()       { return Py_BuildValue("d", gRandomA->Uniform()); }
  static PyObject* Py_RandomNormal () { return Py_BuildValue("d", gRandomA->Normal()); }
  static PyObject* Py_SetSeed (PyObject* self, PyObject* arg) {
    gRandomA->SetSeed(PyFloat_AsDouble(arg));
    Py_INCREF(Py_None);
    return Py_None;
  }

  static PyObject* Py_COUT (PyObject* self, PyObject* args, PyObject *keywds)
  {
    char const* Out = "";
    static const char *kwlist[] = {"out", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, keywds, "s", const_cast<char **>(kwlist), &Out)) {
      return NULL;
    }
    std::cout << Out << std::endl;

    // Must return python object None in a special way
    Py_INCREF(Py_None);
    return Py_None;
  }
  static PyObject* Py_CERR (PyObject* self, PyObject* args, PyObject *keywds)
  {
    char const* Out = "";
    static const char *kwlist[] = {"out", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, keywds, "s", const_cast<char **>(kwlist), &Out)) {
      return NULL;
    }
    std::cerr << Out << std::endl;

    // Must return python object None in a special way
    Py_INCREF(Py_None);
    return Py_None;
  }

}


#endif
