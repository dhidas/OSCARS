#include "TFieldPythonFunction.h"

#include "OSCARSPY.h"

#include <stdexcept>

// UPDATE: exceptions

TFieldPythonFunction::TFieldPythonFunction (PyObject* Function,
                                            TVector3D const& Rotations,
                                            TVector3D const& Translation,
                                            double const TimeOffset,
                                            std::string const& Name)
{
  // Constructor takes a python object, which should be a function
  // Increment reference because we're going to keep it..

  Py_INCREF(Function);
  fPythonFunction = Function;
  this->SetName(Name);
  this->SetScaleFactorMinimumMaximum();

  fRotations = Rotations;
  fTranslation = Translation;
  fTimeOffset = TimeOffset;

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




TVector3D TFieldPythonFunction::GetF (TVector3D const& X, double const T) const
{
  // Get the magnetic field from a python function.

  // Check to see the function is callable
  if (!PyCallable_Check(fPythonFunction)) {
    throw std::runtime_error("For some reason GetF is saying python function not callable.  Please report this.");
  }

  // Translate back into box frame
  TVector3D XNew = X;
  XNew.RotateSelfXYZ(fRotations);
  XNew -= fTranslation;

  // Build the input object for the python function
  PyObject* InputTuple;
  InputTuple = Py_BuildValue("(dddd)", XNew.GetX(), XNew.GetY(), XNew.GetZ(), T + fTimeOffset);

  // Call python function
  PyObject* OutputTuple = PyEval_CallObject(fPythonFunction, InputTuple);

  // We're done with the input object
  Py_DECREF(InputTuple);

  // If the output is null we didn't get anything
  if (OutputTuple == NULL) {
    throw std::logic_error("TFieldPythonFunction::GetF output tuple is NULL.  Please report this.");
  }

  // Get a python list from output tuple
  PyObject* OutputList;
  if (!PyArg_Parse(OutputTuple, "O!", &PyList_Type, &OutputList)) {
    throw std::logic_error("TFieldPythonFunction::GetF cannot get from output tuple.  Please report this.");
  }

  // Rotate field
  TVector3D ReturnVector = OSCARSPY::ListAsTVector3D(OutputList);
  ReturnVector.RotateSelfXYZ(fRotations);

  // Decrement object references no longer needed
  Py_DECREF(OutputTuple);
  Py_DECREF(OutputList);

  // Return the magnetic field vector
  return ReturnVector;
}




TVector3D TFieldPythonFunction::GetF (double const X, double const Y, double const Z, double const T) const
{
  return this->GetF(TVector3D(X, Y, Z), T);
}




TVector3D TFieldPythonFunction::GetRotations () const
{
  return fRotations;
}




TVector3D TFieldPythonFunction::GetTranslation () const
{
  return fTranslation;
}




double TFieldPythonFunction::GetTimeOffset () const
{
  // Return the time offset
  return fTimeOffset;
}








void TFieldPythonFunction::Print (std::ostream& os) const
{
  os << *this << std::endl;
  return;
}

