#ifndef GUARD_TSpectrumContainer_h
#define GUARD_TSpectrumContainer_h
////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Thu May 19 16:05:08 EDT 2016
//
////////////////////////////////////////////////////////////////////

#include <vector>
#include <string>

class TSpectrumContainer
{
  public:
    TSpectrumContainer ();
    TSpectrumContainer (std::vector<double> const&);
    TSpectrumContainer (size_t const, double const, double const);
    ~TSpectrumContainer ();

    void Init (size_t const, double const, double const);
    void Init (std::vector<double> const&);

    void   SetFlux   (size_t const, double const);
    void   SetPoint  (size_t const, double const, double const);
    void   AddPoint  (double const, double const Flux = 0);
    void   AddToFlux (size_t const, double const);
    double GetFlux   (size_t const) const;
    double GetEnergy (size_t const) const;
    double GetAngularFrequency (size_t const) const;
    size_t GetNPoints () const;

    void WriteToFileText (std::string const, std::string const Header = "") const;
    void WriteToFileBinary (std::string const, std::string const Header = "") const;

    void Clear ();
    void AverageFromFilesText (std::vector<std::string> const&);
    void AverageFromFilesBinary (std::vector<std::string> const&);

  private:

    std::vector< std::pair<double, double> > fSpectrumPoints;
    std::vector<double> fCompensation;


};








#endif
