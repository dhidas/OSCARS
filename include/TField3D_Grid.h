#ifndef GUARD_TField3D_Grid_h
#define GUARD_TField3D_Grid_h
////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Tue Sep 20 07:55:00 EDT 2016
//
////////////////////////////////////////////////////////////////////

#include "TField.h"

#include <string>
#include <vector>
#include <array>

#include "TVector3D.h"
#include "TVector2D.h"

class TField3D_Grid : public TField
{

  public:
    TField3D_Grid ();
    TField3D_Grid (std::string         const& InFileName,
                   std::string         const& FileFormat = "OSCARS",
                   TVector3D           const& Rotations = TVector3D(0, 0, 0),
                   TVector3D           const& Translation = TVector3D(0, 0, 0),
                   std::vector<double> const& Scaling = std::vector<double>(),
                   char                const  CommentChar = '#');

    TField3D_Grid (std::vector<std::pair<double, std::string> > Mapping,
                   std::string                           const& FileFormat,
                   double                                const  Parameter,
                   TVector3D                             const& Rotations = TVector3D(0, 0, 0),
                   TVector3D                             const& Translation = TVector3D(0, 0, 0),
                   std::vector<double>                   const& Scaling = std::vector<double>(),
                   char                                  const  CommentChar = '#');

    ~TField3D_Grid ();

    double    GetFx (double const X, double const Y, double const Z) const;
    double    GetFy (double const X, double const Y, double const Z) const;
    double    GetFz (double const X, double const Y, double const Z) const;
    TVector3D GetF  (double const X, double const Y, double const Z) const;
    TVector3D GetF  (TVector3D const& X) const;

    size_t GetIndex (size_t const ix, size_t const iy, size_t const iz) const;

    double GetHeaderValue    (std::string const&) const;
    double GetHeaderValueSRW (std::string const&, const char CommentChar = '#') const;

    void ReadFile (std::string         const& InFileName,
                   TVector3D           const& Rotations = TVector3D(0, 0, 0),
                   TVector3D           const& Translation = TVector3D(0, 0, 0),
                   std::vector<double> const& Scaling = std::vector<double>(),
                   char                const  CommentChar = '#');

    void ReadFile_OSCARS1D  (std::string         const& InFileName,
                             std::string         const& InFormat,
                             TVector3D           const& Rotations = TVector3D(0, 0, 0),
                             TVector3D           const& Translation = TVector3D(0, 0, 0),
                             std::vector<double> const& Scaling = std::vector<double>(),
                             char                const  CommentChar = '#');

    void ReadFile_SRW       (std::string const& InFileName,
                             TVector3D   const& Rotations = TVector3D(0, 0, 0),
                             TVector3D   const& Translation = TVector3D(0, 0, 0),
                             char        const  CommentChar = '#');

    void ReadFile_SPECTRA   (std::string const& InFileName,
                             TVector3D   const& Rotations = TVector3D(0, 0, 0),
                             TVector3D   const& Translation = TVector3D(0, 0, 0),
                             char        const  CommentChar = '#');

    void InterpolateFromFiles (std::vector<std::pair<double, std::string> > const& Mapping,
                               double                                       const  Parameter,
                               TVector3D                                    const& Rotations = TVector3D(0, 0, 0),
                               TVector3D                                    const& Translation = TVector3D(0, 0, 0),
                               std::vector<double>                          const& Scaling = std::vector<double>(),
                               char                                         const  CommentChar = '#');



    void InterpolateFromFiles_OSCARS1D (std::vector<std::pair<double, std::string> > const& Mapping,
                                        double                                       const  Parameter,
                                        std::string                                  const& InFormat,
                                        TVector3D                                    const& Rotations,
                                        TVector3D                                    const& Translation,
                                        std::vector<double>                          const& Scaling,
                                        char                                         const  CommentChar = '#');

    void InterpolateFromFiles_SRW (std::vector<std::pair<double, std::string> > const& Mapping,
                                   double                                       const  Parameter,
                                   TVector3D                                    const& Rotations = TVector3D(0, 0, 0),
                                   TVector3D                                    const& Translation = TVector3D(0, 0, 0),
                                   std::vector<double>                          const& Scaling = std::vector<double>(),
                                   char                                         const  CommentChar = '#');



    TVector3D InterpolateFields (std::vector<double>    const& Parameters,
                                 std::vector<TVector3D> const& Fields,
                                 double                 const  Parameter);

    static bool CompareField1D (std::array<double, 4> const& A,
                                std::array<double, 4> const& B);

    static bool CompareMappingElements (std::pair<double, std::string> const& A,
                                        std::pair<double, std::string> const& B);

    TVector2D GetXRange () const;
    TVector2D GetYRange () const;
    TVector2D GetZRange () const;

    double GetXStep () const;
    double GetYStep () const;
    double GetZStep () const;

    void Print (std::ostream& os) const;


    enum TField3D_Grid_DIMX {
      kDIMX_X,
      kDIMX_Y,
      kDIMX_Z,
      kDIMX_XY,
      kDIMX_XZ,
      kDIMX_YZ,
      kDIMX_XYZ
    };



  private:
    // Dimension and position data
    int    fNX;
    int    fNY;
    int    fNZ;
    double fXStart;
    double fYStart;
    double fZStart;
    double fXStep;
    double fYStep;
    double fZStep;
    double fXStop;
    double fYStop;
    double fZStop;

    bool fHasX;
    bool fHasY;
    bool fHasZ;
    int  fXDIM;

    TField3D_Grid_DIMX fDIMX;

    // Rotations and Translations.  Field is stored rotated, point must be rotated into grid
    // coordinate system before asking for field, but remember field is already rotated.
    TVector3D fRotated;
    TVector3D fTranslation;

    // Field data
    std::vector<TVector3D> fData;

};





inline std::ostream& operator << (std::ostream& os, TField3D_Grid const& o)
{
  // For easy printing
  os << "TField3D_Grid " << "\n"
     << "XRange        " << o.GetXRange() << "\n"
     << "YRange        " << o.GetYRange() << "\n"
     << "ZRange        " << o.GetZRange() << "\n"
     << "XStep         " << o.GetXStep() << "\n"
     << "YStep         " << o.GetYStep() << "\n"
     << "ZStep         " << o.GetZStep() << "\n";

  return os;
}
















#endif

