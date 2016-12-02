#include "TFieldPythonFunction.h"

#include <stdexcept>

// UPDATE: exceptions

TFieldPythonFunction::TFieldPythonFunction (PyObject* Function)
{
  // Constructor takes a python object, which should be a function
  // Increment reference because we're going to keep it..

  Py_INCREF(Function);
  fPythonFunction = Function;

  // Check to see the function is callable
  if (!PyCallable_Check(fPythonFunction)) {
    throw std::invalid_argument("python function not callable");
  }
}



TFieldPythonFunction::~TFieldPythonFunction ()
{
  // When exit, decrement reference since we're done with it
  Py_DECREF(fPythonFunction);
}




TVector3D TFieldPythonFunction::GetF (TVector3D const& X) const
{
  // Get the magnetic field from a python function.

  // For the future
  double T = 0;

  // Check to see the function is callable
  if (!PyCallable_Check(fPythonFunction)) {
    throw;
  }

  // Build the input object for the python function
  PyObject* InputTuple;
  InputTuple = Py_BuildValue("(dddd)", X.GetX(), X.GetY(), X.GetZ(), T);

  // Call python function
  PyObject* OutputTuple = PyEval_CallObject(fPythonFunction, InputTuple);

  // We're done with the input object
  Py_DECREF(InputTuple);

  // If the output is null we didn't get anything
  if (OutputTuple == NULL) {
    throw;
  }

  // Get a python list from output tuple
  PyObject* OutputList;
  if (!PyArg_Parse(OutputTuple, "O!", &PyList_Type, &OutputList)) {
    throw;
  }



  // Observation point from python list
  TVector3D ReturnVector(PyFloat_AsDouble(PyList_GetItem(OutputList, 0)),
                         PyFloat_AsDouble(PyList_GetItem(OutputList, 1)),
                         PyFloat_AsDouble(PyList_GetItem(OutputList, 2)));

  // Decrement object references no longer needed
  Py_DECREF(OutputTuple);
  Py_DECREF(OutputList);

  // Return the magnetic field vector
  return ReturnVector;
}




TVector3D TFieldPythonFunction::GetF (double const X, double const Y, double const Z) const
{
  return this->GetF(TVector3D(X, Y, Z));
}




double TFieldPythonFunction::GetFx (double const X, double const Y, double const Z) const
{
  return this->GetF(X, Y, Z).GetX();
}




double TFieldPythonFunction::GetFy (double const X, double const Y, double const Z) const
{
  return this->GetF(X, Y, Z).GetY();
}




double TFieldPythonFunction::GetFz (double const X, double const Y, double const Z) const
{
  return this->GetF(X, Y, Z).GetZ();
}
