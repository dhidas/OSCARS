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
#include <iostream>
#include <fstream>
#include <stdexcept>

#include "TVector3D.h"

namespace TOMATH
{

// Modified bessel function K
double BesselK  (double const nu, double const x);
double BesselK_IntegralToInfty (double const nu, double const x);
double BesselJ0 (double const x);
double BesselJ1 (double const x);
double BesselJ  (int    const nu, double const x);




template <class T> class TSpline1D3
{
  // This is a template class for cubic spline interpolation in 1D
  // of any class T where T must have the correct operators defined
  public:
    TSpline1D3 ()
    {
    }
    ~TSpline1D3 ()
    {
    }

    TSpline1D3 (std::vector<double> const& X, std::vector<T> const& Y)
    {
      this->Set(X, Y);
    }

    void Set (std::vector<double> const& X, std::vector<T> const& Y)
    {
      fX.clear();
      fY.clear();
      fYPP.clear();

      if (X.size() != Y.size()) {
        throw std::length_error("TSpline1D3 detected the length of each input is different");
      }

      for (size_t i = 0; i != X.size(); ++i) {
        fX.push_back(X[i]);
        fY.push_back(Y[i]);
      }

      Derivatives();

      return;
    }

    void Set (double const* X, T const* Y, int const N)
    {
      fX.clear();
      fY.clear();
      fYPP.clear();
      if (N <= 0) {
        throw std::length_error("TSpline1D3 believes that N <= 0");
      }

      std::vector<double> XI(N);
      std::vector<T> YI(N);

      for (int i = 0; i < N; ++i) {
        XI[i] = X[i];
        YI[i] = Y[i];
      }

      this->Set(XI, YI);

      return;
    }


    T GetValue (double const x) const
    {
      // Return the Y-value according to spline

      // Find the index in fX before and after x
      //std::vector<double>::const_iterator ilo = std::lower_bound(fX.begin(), fX.end(), x);
      //size_t const klo = std::distance(fX.begin(), ilo) - 0;
      //size_t const khi = klo + 1;

      int klo=0;
      int khi = (int) fX.size() - 1;
      int k;
      while (khi - klo > 1) {
        k = (khi + klo) >> 1;
        if (fX[k] > x) {
          khi = k;
        } else {
          klo = k;
        }
      }

      // Distance between points, check that it isn't zero!
      double const h = fX[khi] - fX[klo];
      if (h == 0) {
        throw std::out_of_range("TSpline1D3 sees that the stepsize h is zero");
      }

      // Fractional distance to the points on either side
      double const a = (fX[khi] - x) / h;
      double const b = (x - fX[klo]) / h;

      // Return the value of Y

      return a * fY[klo] + b * fY[khi] + ((a * a * a - a) * fYPP[klo] + (b * b * b - b) * fYPP[khi]) * (h * h) / 6.;
    }




    T GetDerivative (int const i) const
    {
      // Return the dY/dx-value according to spline

      if (i >= fYPP.size()-1) {
        return this->GetDerivative(i-1);
      }
      return (fY[i+1] - fY[i]) / (fX[i+1] - fX[i]) - fYPP[i] * (fX[i+1] - fX[i]) / 3 - fYPP[i+1] * (fX[i+1] - fX[i]) / 6;
    }

    void Derivatives ()
    {
      // Number in vectors
      int const N = (int) fX.size();
      if (N != (int) fY.size()) {
        throw std::length_error("TSpline1D3 sees that N is not equal to the length of Y");
      }


      if (N < 3) {
        throw std::length_error("TSpline1D3 does not see enough points for interpolation.  Please use >= 3 points");
      }

      // Resize output
      fYPP.resize(N);

      T p;
      T sig;

      std::vector<T> u(N);

      T const dydx0 = (fY[1] - fY[0]) / (fX[1] - fX[0]);

      // Use natural spline
      //if (dydx0 > 0.99e30) {
      //  fYPP[0] = T(0);
      //  u[0] = T(0);
      //} else {
        fYPP[0] = T(-0.5);
        u[0] = (3. / (fX[1] - fX[0])) * ((fY[1] - fY[0]) / (fX[1] - fX[0]) - dydx0);
      //}


      for (int i = 1; i < N-1; ++i) {
        sig = (fX[i] - fX[i-1]) / (fX[i+1] - fX[i-1]);
        p = sig * fYPP[i-1] + T(2);
        fYPP[i] = (sig - 1.) / p;
        u[i] = (fY[i+1] - fY[i]) / (fX[i+1] - fX[i]) - (fY[i] - fY[i-1]) / (fX[i] - fX[i-1]);
        u[i] = (6. * u[i] / (fX[i+1] - fX[i-1]) - sig * u[i-1]) / p;
      }

      // Left as a reminder for now
      double qn;
      T un;
      T const dydxN = (fY[N-1] - fY[N-2]) / (fX[N-1] - fX[N-2]);

      // Use natural spline
      //if (dydxN > 0.99e30) {
      //  qn = 0;
      //  un = T(0);
      //} else {
        qn = 0.5;
        un = (3. / (fX[N-1] - fX[N-2])) * (dydxN - (fY[N-1] - fY[N-2]) / (fX[N-1] - fX[N-2]));
      //}

      fYPP[N-1] = (un - qn * u[N-2]) / (qn * fYPP[N-2] + T(1.));

      for (int k = N-2; k >= 0; --k) {
        fYPP[k] = fYPP[k] * fYPP[k+1] + u[k];
      }

      return;
    }

    void Clear ()
    {
      // Clear all vectors
      fX.clear();
      fY.clear();
      fYPP.clear();
    }

    size_t GetNPoints () const
    {
      return fX.size();
    }

    double GetXStart () const
    {
      return fX.front();
    }

    double GetXStop () const
    {
      return fX.back();
    }

    double GetX (int const i) const
    {
      return fX[i];
    }

    T const& GetY (int const i) const
    {
      return fY[i];
    }

    T const& GetYPP (int const i) const
    {
      return fYPP[i];
    }

  private:
    std::vector<double> fX;   // x-values
    std::vector<T> fY;        // known y-values
    std::vector<T> fYPP;      // calculated derivs
};




class TSpline1D3_1d : public TSpline1D3<double>
{
  public:
    TSpline1D3_1d ()
    {
    }
    ~TSpline1D3_1d ()
    {
    }

    void ReadFile (std::string const& InFileName)
    {
      // Unfinished

      std::vector<double> X;
      std::vector<double> Y;

      double x;
      double y;

      std::ifstream f(InFileName);
      if (!f.is_open()) {
        throw std::ifstream::failure("TSpline1D3_1d cannot open file for writing");
      }

      while (!f.eof()) {
        f >> x >> y;
        if (!f.bad()) {
          this->Set(X, Y);
        }
      }

      return;
    }

};









}






#endif
