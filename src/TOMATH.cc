////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Thu Mar  9 11:26:17 EST 2017
//
////////////////////////////////////////////////////////////////////

#include "TOMATH.h"

#include <cmath>
#include <fstream>


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






double BesselK_IntegralToInfty (double const nu, double const x)
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
    RthTerm = exp(-x * cosh(r * h)) * cosh(nu * r * h) / cosh(r * h);

    // Add r-th term to return vvalue
    K += RthTerm * h;
  }

  // Return the value of the modified bessel function
  return K;
}















double BesselJ0 (double const x)
{
  double ax,z;
  double xx,y,result,result1,result2;
  double const p1  = 57568490574.0, p2  = -13362590354.0, p3 = 651619640.7;
  double const p4  = -11214424.18,  p5  = 77392.33017,    p6 = -184.9052456;
  double const p7  = 57568490411.0, p8  = 1029532985.0,   p9 = 9494680.718;
  double const p10 = 59272.64853,   p11 = 267.8532712;

  double const q1  = 0.785398164;
  double const q2  = -0.1098628627e-2,  q3  = 0.2734510407e-4;
  double const q4  = -0.2073370639e-5,  q5  = 0.2093887211e-6;
  double const q6  = -0.1562499995e-1,  q7  = 0.1430488765e-3;
  double const q8  = -0.6911147651e-5,  q9  = 0.7621095161e-6;
  double const q10 =  0.934935152e-7,   q11 = 0.636619772;

  if ((ax=fabs(x)) < 8) {
    y=x*x;
    result1 = p1 + y*(p2 + y*(p3 + y*(p4  + y*(p5  + y*p6))));
    result2 = p7 + y*(p8 + y*(p9 + y*(p10 + y*(p11 + y))));
    result  = result1/result2;
  } else {
    z  = 8/ax;
    y  = z*z;
    xx = ax-q1;
    result1 = 1  + y*(q2 + y*(q3 + y*(q4 + y*q5)));
    result2 = q6 + y*(q7 + y*(q8 + y*(q9 - y*q10)));
    result  = sqrt(q11/ax)*(cos(xx)*result1-z*sin(xx)*result2);
  }
  return result;
}


double BesselJ1 (double const x)
{
  double ax,z;
  double xx,y,result,result1,result2;
  double const p1  = 72362614232.0,  p2  = -7895059235.0, p3 = 242396853.1;
  double const p4  = -2972611.439,   p5  = 15704.48260,   p6 = -30.16036606;
  double const p7  = 144725228442.0, p8  = 2300535178.0,  p9 = 18583304.74;
  double const p10 = 99447.43394,    p11 = 376.9991397;

  double const q1  = 2.356194491;
  double const q2  = 0.183105e-2,     q3  = -0.3516396496e-4;
  double const q4  = 0.2457520174e-5, q5  = -0.240337019e-6;
  double const q6  = 0.04687499995,   q7  = -0.2002690873e-3;
  double const q8  = 0.8449199096e-5, q9  = -0.88228987e-6;
  double const q10 = 0.105787412e-6,  q11 = 0.636619772;

  if ((ax=fabs(x)) < 8) {
    y=x*x;
    result1 = x*(p1 + y*(p2 + y*(p3 + y*(p4  + y*(p5  + y*p6)))));
    result2 = p7    + y*(p8 + y*(p9 + y*(p10 + y*(p11 + y))));
    result  = result1/result2;
  } else {
    z  = 8/ax;
    y  = z*z;
    xx = ax-q1;
    result1 = 1  + y*(q2 + y*(q3 + y*(q4 + y*q5)));
    result2 = q6 + y*(q7 + y*(q8 + y*(q9 + y*q10)));
    result  = sqrt(q11/ax)*(cos(xx)*result1-z*sin(xx)*result2);
    if (x < 0) result = -result;
  }
  return result;
}


#define ACC 40.0
#define BIGNO 1.0e10
#define BIGNI 1.0e-10


double BesselJ (int const nu, double const x)
{
  // Function to calculate the bessel J for integer nu.


  int const n = abs(nu);

  int j,jsum,m;
  float ax,bj,bjm,bjp,sum,tox,ans;

  if (nu == 0) {
    return BesselJ0(x);
  } else if (nu == 1) {
    return BesselJ1(x);
  } else if (nu == -1) {
    return -BesselJ1(x);
  }

  ax = fabs(x);
  if (ax == 0) {
    return 0;
  } else if (ax > (double) n) {
    tox = 2.0 / ax;
    bjm = BesselJ0(ax);
    bj = BesselJ1(ax);
    for (j = 1; j < n; j++) {
      bjp = j * tox * bj - bjm;
      bjm = bj;
      bj = bjp;
    }
    ans = bj;
  } else {
    tox = 2. / ax;
    m = 2 * ((n + (int) sqrt(ACC * n)) / 2);
    jsum = 0;
    bjp = ans = sum = 0.;
    bj = 1.;
    for (j = m; j > 0; j--) {
      bjm = j * tox * bj - bjp;
      bjp = bj;
      bj = bjm;
      if (fabs(bj) > BIGNO) {
        bj *= BIGNI;
        bjp *= BIGNI;
        ans *= BIGNI;
        sum *= BIGNI;
      }
      if (jsum) {
        sum += bj;
      }
      jsum = !jsum;
      if (j == n) {
        ans = bjp;
      }
    }
    sum = 2. * sum - bj;
    ans /= sum;
  }

  if (nu < 0 && (nu & 1)) {
    return x < 0. && (n & 1) ? ans : -ans;
  }

  return x < 0. && (n & 1) ? -ans : ans;
}





























}


