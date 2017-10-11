////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Tue Jul 19 08:13:02 EDT 2016
//
////////////////////////////////////////////////////////////////////

#include "TField3D_IdealUndulator.h"

#include "TOSCARSSR.h"
#include <cmath>


TField3D_IdealUndulator::TField3D_IdealUndulator (std::string const& Name)
{
  // Constructor

  // Set the name and default scale factors
  this->SetName(Name);
  this->SetScaleFactorMinimumMaximum();
}




TField3D_IdealUndulator::TField3D_IdealUndulator (TVector3D const& Field,
                                                  TVector3D const& Period,
                                                  int const NPeriods,
                                                  TVector3D const& Center,
                                                  double const Phase,
                                                  double const Taper,
                                                  std::string const& Name)
{
  // Typical constructor.  This will generate a magnetic field in the direction specified by Field
  // varying as the sine depending on the period, number of periods and phase, centered at a given location.
  // The length of the period is input as a 3-D vector.  The period is the magnitude of this vector,
  // and the direction is taken to be the axis of variation

  // Field - Peak magnetic field in [T]
  // Period - Magnitude is the period while direction is the axis for that variation
  // NPeriods - Number of periods
  // Center - Where in space the center of this field will be
  // Phase - A phase offset for the sine function given in [rad]

  this->Init(Field, Period, NPeriods, Center, Phase, Taper, Name);
}




TField3D_IdealUndulator::~TField3D_IdealUndulator ()
{
  // Destruction!
}




void TField3D_IdealUndulator::Init (TVector3D const& Field,
                                    TVector3D const& Period,
                                    int const NPeriods,
                                    TVector3D const& Center,
                                    double const Phase,
                                    double const Taper,
                                    std::string const& Name)
{
  // Initialization function.  This will generate a magnetic field in the direction specified by Field
  // varying as the sine depending on the period, number of periods and phase, centered at a given location.
  // The length of the period is input as a 3-D vector.  The period is the magnitude of this vector,
  // and the direction is taken to be the axis of variation

  // Field - Peak magnetic field in [T]
  // Period - Magnitude is the period while direction is the axis for that variation
  // NPeriods - Number of periods
  // Center - Where in space the center of this field will be
  // Phase - A phase offset for the sine function given in [rad]
  // Name - Name of this field

  // Set the name and default scale factors
  this->SetName(Name);
  this->SetScaleFactorMinimumMaximum();

  fField   = Field;
  fPeriod   = Period;
  fNPeriods = NPeriods;
  fCenter   = Center;
  fPhase    = Phase;
  fTaper    = Taper;

  fPeriodLength = fPeriod.Mag();
  fPeriodUnitVector = fPeriod.UnitVector();

  // Length is 2 periods longer than NPeriods to account for terminating fields
  fUndulatorLength = fPeriod.Mag() * (fNPeriods + 2);

  return;
}




double TField3D_IdealUndulator::GetFx (double const X, double const Y, double const Z) const
{
  return this->GetF(X, Y, Z).GetX();
}




double TField3D_IdealUndulator::GetFy (double const X, double const Y, double const Z) const
{
  return this->GetF(X, Y, Z).GetY();
}




double TField3D_IdealUndulator::GetFz (double const X, double const Y, double const Z) const
{
  return this->GetF(X, Y, Z).GetZ();
}




TVector3D TField3D_IdealUndulator::GetF (TVector3D const& X) const
{

  // UPDATE: Check min.max limits here?  maybe not

  // How far are you from the "center" in the correct direction
  double const D = (X - fCenter).Dot( fPeriodUnitVector );

  // How much taper correction do we make?
  double const TaperCorrection = 1 + D * fTaper;

  // Phase shift in length
  double const PhaseShift = fPhase * fPeriod.Mag() / TOSCARSSR::TwoPi ();

  // Vector we will return
  TVector3D F(0, 0, 0);

  // Check if we are outside of the NPeriod + termination range
  if (D > fUndulatorLength / 2. + PhaseShift || D < -fUndulatorLength / 2. + PhaseShift) {
    return F;
  }

  // Check if we are within the termination period
  if (D < -fUndulatorLength / 2. + PhaseShift + fPeriodLength || D > fUndulatorLength / 2. + PhaseShift - fPeriodLength) {
    if (D < -fUndulatorLength / 2. + PhaseShift + fPeriodLength / 2. || D > fUndulatorLength / 2. + PhaseShift - fPeriodLength / 2.) {

      return 0.25 * fField * sin(TOSCARSSR::TwoPi() * (D - PhaseShift) / fPeriodLength) * TaperCorrection;
    } 

    return 0.75 * fField * sin(TOSCARSSR::TwoPi() * (D - PhaseShift) / fPeriodLength) * TaperCorrection;
  } 

  return fField * sin(TOSCARSSR::TwoPi() * (D - PhaseShift) / fPeriodLength) * TaperCorrection;
}




TVector3D TField3D_IdealUndulator::GetF (double const X, double const Y, double const Z) const
{
  return this->GetF(TVector3D(X, Y, Z));
}





TVector3D TField3D_IdealUndulator::GetField () const
{
  // Return the peak field
  return fField;
}




TVector3D TField3D_IdealUndulator::GetPeriod () const
{
  // Return the period
  return fPeriod;
}




int TField3D_IdealUndulator::GetNPeriods () const
{
  // Return the peak field
  return fNPeriods;
}




TVector3D TField3D_IdealUndulator::GetCenter () const
{
  // Return the center position
  return fCenter;
}




double TField3D_IdealUndulator::GetPhase () const
{
  // Return the phase shift
  return fPhase;
}




double TField3D_IdealUndulator::GetTaper () const
{
  // Return the taper
  return fTaper;
}




void TField3D_IdealUndulator::Print (std::ostream& os) const
{
  os << *this << std::endl;
  return;
}
