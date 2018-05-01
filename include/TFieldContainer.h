#ifndef GUARD_TFieldContainer
#define GUARD_TFieldContainer
////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Thu Sep 22 08:19:53 EDT 2016
//
// This class is meant to be a container class for all fields in
// a given simulation.  It will sum all B contributions and return
// the sum Fx, Fy, Fz, or zero where there is no defined field
//
////////////////////////////////////////////////////////////////////

#include <vector>

#include "TField.h"
#include "TVector2D.h"
#include "TVector3D.h"

class TFieldContainer
{
  public:
    TFieldContainer ();
    TFieldContainer (TField*);
    ~TFieldContainer ();

    void AddField (TField*);
    void RemoveField (std::string const& Name);

    TVector3D GetF  (double const, double const, double const, double const T = 0, std::string const& Name = "") const;
    TVector3D GetF  (TVector3D const&, double const T = 0, std::string const& Name = "") const;

    TField const& GetField (size_t const) const;

    size_t GetNFields () const;

    void      Clear ();

    void WriteToFile (std::string const& OutFileName,
                      std::string const& OutFormat,
                      TVector2D const& XLim,
                      int const NX,
                      TVector2D const& YLim,
                      int const NY,
                      TVector2D const& ZLim,
                      int const NZ,
                      std::string const Comment = "");

    void WriteToFileBinary (std::string const& OutFileName,
                            std::string const& OutFormat,
                            TVector2D const& XLim,
                            int const NX,
                            TVector2D const& YLim,
                            int const NY,
                            TVector2D const& ZLim,
                            int const NZ,
                            std::string const Comment = "",
                            int const Version = 1);

    void WriteToFileBinary_v1 (std::string const& OutFileName,
                               std::string const& OutFormat,
                               TVector2D const& XLim,
                               int const NX,
                               TVector2D const& YLim,
                               int const NY,
                               TVector2D const& ZLim,
                               int const NZ,
                               std::string const Comment = "");

  private:
    std::vector<TField*> fFields;
};


inline std::ostream& operator << (std::ostream& os, TFieldContainer const& o)
{
  // For easy printing
  os << "TFieldContainer has " << o.GetNFields() << " fields" << std::endl;

  size_t const N = o.GetNFields();

  for (size_t i = 0; i != N; ++i) {
    o.GetField(i).Print(os);
  }

  return os;
}




























#endif
