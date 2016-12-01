////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Tue Jul 19 08:13:02 EDT 2016
//
////////////////////////////////////////////////////////////////////

#include "TBField3D_IdealUndulator.h"

#include "TSRS.h"
#include <cmath>


TBField3D_IdealUndulator::TBField3D_IdealUndulator ()
{
  // Constructor
}




TBField3D_IdealUndulator::TBField3D_IdealUndulator (TVector3D const& BField, TVector3D const& Period, int const NPeriods, TVector3D const& Center, double const Phase)
{
  // Typical constructor.  This will generate a magnetic field in the direction specified by BField
  // varying as the sine depending on the period, number of periods and phase, centered at a given location.
  // The length of the period is input as a 3-D vector.  The period is the magnitude of this vector,
  // and the direction is taken to be the axis of variation

  // BField - Peak magnetic field in [T]
  // Period - Magnitude is the period while direction is the axis for that variation
  // NPeriods - Number of periods
  // Center - Where in space the center of this field will be
  // Phase - A phase offset for the sine function given in [rad]

  this->Init(BField, Period, NPeriods, Center, Phase);
}




TBField3D_IdealUndulator::~TBField3D_IdealUndulator ()
{
  // Destruction!
}




void TBField3D_IdealUndulator::Init (TVector3D const& BField, TVector3D const& Period, int const NPeriods, TVector3D const& Center, double const Phase)
{
  // Initialization function.  This will generate a magnetic field in the direction specified by BField
  // varying as the sine depending on the period, number of periods and phase, centered at a given location.
  // The length of the period is input as a 3-D vector.  The period is the magnitude of this vector,
  // and the direction is taken to be the axis of variation

  // BField - Peak magnetic field in [T]
  // Period - Magnitude is the period while direction is the axis for that variation
  // NPeriods - Number of periods
  // Center - Where in space the center of this field will be
  // Phase - A phase offset for the sine function given in [rad]

  fBField   = BField;
  fPeriod   = Period;
  fNPeriods = NPeriods;
  fCenter   = Center;
  fPhase    = Phase;

  fPeriodLength = fPeriod.Mag();
  fPeriodUnitVector = fPeriod.UnitVector();

  // Length is 2 periods longer than NPeriods to account for terminating fields
  fUndulatorLength = fPeriod.Mag() * (fNPeriods + 2);

  return;
}




double TBField3D_IdealUndulator::GetBx (double const X, double const Y, double const Z) const
{
  return this->GetB(X, Y, Z).GetX();
}




double TBField3D_IdealUndulator::GetBy (double const X, double const Y, double const Z) const
{
  return this->GetB(X, Y, Z).GetY();
}




double TBField3D_IdealUndulator::GetBz (double const X, double const Y, double const Z) const
{
  return this->GetB(X, Y, Z).GetZ();
}




TVector3D TBField3D_IdealUndulator::GetB (TVector3D const& X) const
{

  // UPDATE: Check min.max limits here?  maybe not

  // How far are you from the "center" in the correct direction
  double const D = (X - fCenter).Dot( fPeriodUnitVector );

  // Phase shift in length
  double const PhaseShift = fPhase * fPeriod.Mag() / TSRS::TwoPi ();

  // Vector we will return
  TVector3D B(0, 0, 0);

  // Check if we are outside of the NPeriod + termination range
  if (D > fUndulatorLength / 2. + PhaseShift || D < -fUndulatorLength / 2. + PhaseShift) {
    return B;
  }

  // Check if we are within the termination period
  if (D < -fUndulatorLength / 2. + PhaseShift + fPeriodLength || D > fUndulatorLength / 2. + PhaseShift - fPeriodLength) {
    if (D < -fUndulatorLength / 2. + PhaseShift + fPeriodLength / 2. || D > fUndulatorLength / 2. + PhaseShift - fPeriodLength / 2.) {

      return 0.25 * fBField * sin(TSRS::TwoPi() * (D - PhaseShift) / fPeriodLength);
    } 

    return 0.75 * fBField * sin(TSRS::TwoPi() * (D - PhaseShift) / fPeriodLength);
  } 

  return fBField * sin(TSRS::TwoPi() * (D - PhaseShift) / fPeriodLength);
}




TVector3D TBField3D_IdealUndulator::GetB (double const X, double const Y, double const Z) const
{
  return this->GetB(TVector3D(X, Y, Z));
}


