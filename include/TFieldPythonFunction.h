#ifndef GUARD_TFieldPythonFunction_h
#define GUARD_TFieldPythonFunction_h
////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Mon Jun 27 07:54:27 EDT 2016
//
////////////////////////////////////////////////////////////////////
#include "Python.h"

#include "TField.h"

class TFieldPythonFunction : public TField
{
  public:
    TFieldPythonFunction (PyObject*);
    ~TFieldPythonFunction ();

    double    GetFx (double const, double const, double const) const;
    double    GetFy (double const, double const, double const) const;
    double    GetFz (double const, double const, double const) const;
    TVector3D GetF  (double const, double const, double const) const;
    TVector3D GetF  (TVector3D const&) const;

    void Print (std::ostream& os) const;

  private:
    PyObject* fPythonFunction;

};



inline std::ostream& operator << (std::ostream& os, TFieldPythonFunction const& o)
{
  // For easy printing
  os << "TFieldPythonFunction at address " << &o << "\n";

  return os;
}




#endif
