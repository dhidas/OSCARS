////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Thu Mar  9 11:26:17 EST 2017
//
////////////////////////////////////////////////////////////////////

#include "TOMATH.h"

#include <cmath>

namespace TOMATH
{

  double BesselK (double const nu, double const x)
  {
    // Compute the modified bessel function according to eqn 15 of:
    //   VACLAV O. KOSTROUN 
    //   SIMPLE NUMERICAL EVALUATION OF MODIFIED BESSEL FUNCTIONS Kv(x) OF FRACTIONAL ORDER AND THE INTEGRAL fxKv(rT)dr7
    //   Nuclear Instruments and Methods 172 (1980) 371-374

    // The interval (0.5 by authros suitable)
    double const h = 0.5;

    // Required precision of the r-th term.  Authors use 1e-5, we'll use higher
    double const epsilon = 1e-15;

    // This is the 0-th term
    double RthTerm = exp(-x) / 2. * h;

    // This is the return value, at the moment only containing the 0-th term
    double K = RthTerm;

    // r is the summation index (ie the r-th term)
    int r = 0;

    // Continue until the r-th term satisfies the precision requirement
    // (Prefer to sum from small to large, but ok)
    while (RthTerm > epsilon) {

      // Increment r
      ++r;

      // Calculate the r-th term
      RthTerm = exp(-x * cosh(r * h)) * cosh(nu * r * h);

      // Add r-th term to return vvalue
      K += RthTerm * h;
    }

    // Return the value of the modified bessel function
    return K;
  }
}


