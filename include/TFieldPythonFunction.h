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
    TFieldPythonFunction (PyObject* Function,
                          TVector3D const& Rotations = TVector3D(0, 0, 0),
                          TVector3D const& Translation = TVector3D(0, 0, 0),
                          double const TimeOffset = 0,
                          std::string const& Name = "");
    ~TFieldPythonFunction ();

    TVector3D GetF  (double const, double const, double const, double const T = 0) const;
    TVector3D GetF  (TVector3D const&, double const T = 0) const;

    TVector3D GetRotations () const;
    TVector3D GetTranslation () const;
    double GetTimeOffset () const;

    void Print (std::ostream& os) const;

  private:
    PyObject* fPythonFunction;

    TVector3D fRotations;
    TVector3D fTranslation;
    double fTimeOffset;

};



inline std::ostream& operator << (std::ostream& os, TFieldPythonFunction const& o)
{
  // For easy printing
  os << "TFieldPythonFunction\n"
     << "  Rotations          " << o.GetRotations() << "\n"
     << "  Translation        " << o.GetTranslation() << "\n"
     << "  TimeOffset         " << o.GetTimeOffset() << "\n"
     << "  Name               " << o.GetName() << "\n"
     << "  at address " << &o << "\n";

  return os;
}




#endif
