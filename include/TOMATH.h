#ifndef GUARD_TOMATH_h
#define GUARD_TOMATH_h
////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Thu Mar  9 11:26:17 EST 2017
//
////////////////////////////////////////////////////////////////////

#include <vector>
#include <string>


namespace TOMATH
{

// Modified bessel function K
double BesselK  (double const nu, double const x);
double BesselJ0 (double const x);
double BesselJ1 (double const x);
double BesselJ  (int    const nu, double const x);




template <class T> class TSpline1D3
{
  public:
    TSpline1D3 ();
    TSpline1D3 (std::vector<double> const&, std::vector<T> const&);
    ~TSpline1D3 ();


    void Set (std::vector<double> const&, std::vector<T> const&);
    void Set (double const*, T const*, int const);
    T GetValue (double const x) const;
    void Derivatives ();

  private:
    std::vector<double> fX;   // x-values
    std::vector<T> fY;        // known y-values
    std::vector<T> fYPP;      // calculated derivs
};




class TSpline1D3_1d : public TSpline1D3<double>
{
  public:
    TSpline1D3_1d ();
    ~TSpline1D3_1d ();

    void ReadFile (std::string const&);
};


}






#endif
