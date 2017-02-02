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

























#endif
