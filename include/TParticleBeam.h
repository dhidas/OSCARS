#ifndef GUARD_TParticleBeam_h
#define GUARD_TParticleBeam_h
////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Tue Jun  7 18:01:19 EDT 2016
//
////////////////////////////////////////////////////////////////////

#include "TParticleA.h"

#include <map>
#include <string>
#include <vector>

#include "TVector3D.h"
#include "TVector2D.h"

class TParticleBeam : public TParticleA
{
  public:
    TParticleBeam ();
    TParticleBeam (std::string const& PredefinedBeamType,
                   std::string const& Name,
                   double const Weight = 1);

    TParticleBeam (std::string const& ParticleType,
                   std::string const& Name,
                   double const Energy,
                   double const Current,
                   double const Weight = 1);

    TParticleBeam (std::string const& ParticleType,
                   std::string const& Name,
                   TVector3D const& X0,
                   TVector3D const& D0,
                   double const Energy,
                   double const Current,
                   double const Weight = 1);

    TParticleBeam (std::string const& ParticleType,
                   std::string const& Name,
                   TVector3D const& X0,
                   TVector3D const& D0,
                   double const Energy,
                   double const T0,
                   double const Current,
                   double const Charge = 0,
                   double const Mass = 0,
                   double const Weight = 1);

    ~TParticleBeam ();

    void SetPredefinedBeam (std::string const& Beam);

    void SetInitialConditions (double const X,
                               double const Y,
                               double const Z,
                               double const Dx,
                               double const Dy,
                               double const Dz,
                               double const E0,
                               double const T);

    void SetInitialConditions (TVector3D const& X,
                               TVector3D const& D,
                               double    const E0,
                               double    const T0);

    void SetTwissBetaAlpha (TVector2D const& Beta,
                            TVector2D const& Alpha,
                            TVector3D const& Reference = TVector3D(0, 0, 0),
                            bool const HasReference = false);

    void SetTwissBetaGamma (TVector2D const& Beta,
                            TVector2D const& Gamma,
                            TVector3D const& Reference = TVector3D(0, 0, 0),
                            bool const HasReference = false);

    void SetTwissAlphaGamma (TVector2D const& Alpha,
                             TVector2D const& Gamma,
                             TVector3D const& Reference = TVector3D(0, 0, 0),
                             bool const HasReference = false);

    void SetEmittance (TVector2D const& Emittance);

    void SetTwissParametersAtX0 ();

    void SetTwissParameters (TVector2D const& Beta,
                             TVector2D const& Alpha,
                             TVector2D const& Gamma,
                             TVector3D const& Reference = TVector3D(0, 0, 0),
                             bool const HasReference = false);

    void SetHorizontalDirection (TVector3D const& HorizontalDirection);

    void SetTwissLatticeReference (TVector3D const& L);

    TVector3D GetTwissLatticeReference () const;


    TVector3D const& GetX0 () const;
    TVector3D const& GetU0 () const;
    double           GetE0 () const;
    double           GetT0 () const;

    std::string const& GetName () const;
    double             GetWeight () const;

    TVector2D GetTwissBeta () const;
    TVector2D GetTwissAlpha () const;
    TVector2D GetTwissGamma () const;
    TVector2D GetTwissBetaX0 () const;
    TVector2D GetTwissAlphaX0 () const;
    TVector2D GetTwissGammaX0 () const;

    TVector2D GetEmittance () const;

    void SetEta (TVector2D const&);
    TVector2D GetEta () const;

    void SetVerticalDirection (TVector3D const& VerticalDirection);

    TVector3D GetVerticalDirection () const;
    TVector3D GetHorizontalDirection () const;
    TVector3D GetSigmaAt () const;

    void SetX0 (TVector3D const&);
    void SetU0 (TVector3D const&);
    void SetE0 (double const);
    void SetT0 (double const);
    void SetSigmaEnergyGeV (double const);
    double GetSigmaEnergyGeV () const;

    void SetName (std::string const&);
    void SetWeight (double const);

    TParticleA GetNewParticle ();
    TParticleA GetNewParticle (std::string const&);

    // Types of beam distributions supported
    enum TParticleBeam_BeamDistribution {
      kBeamDistribution_None,
      kBeamDistribution_Filament,
      kBeamDistribution_Gaussian,
      kBeamDistribution_KV
    };

    void SetBeamDistribution (TParticleBeam_BeamDistribution const D);
    void SetBeamDistribution (std::string const& Name);
    TParticleBeam_BeamDistribution const GetBeamDistribution () const;

    TParticleBeam_BeamDistribution GetBeamDistribution (std::string const& Name) const;
    std::string GetBeamDistributionName () const;

  private:
    std::string fName;
    double      fWeight;

    TVector3D fX0;  // Coordinates of initial conditions [m]
    TVector3D fU0;  // Direction at fX0 (stored as a unit vector) [unitless]
    double    fE0;  // Energy [GeV]
    double    fT0;  // Time at initial conditions [m]

    // Twiss parameters, emittance, and dispersion (at reference given)
    TVector2D fTwissBeta;
    TVector2D fTwissAlpha;
    TVector2D fTwissGamma;
    TVector2D fEmittance;
    TVector2D fEta;
    TVector3D fTwissLatticeReference;

    // Twiss parameters at beam x0
    TVector2D fTwissBetaX0;
    TVector2D fTwissAlphaX0;
    TVector2D fTwissGammaX0;

    // Beam distribution type
    TParticleBeam_BeamDistribution fBeamDistribution;

    // Source size, calculated from Beta and Emittance
    TVector2D fSigmaU;
    TVector2D fSigmaUP;

    // Energy spread (in units of [GeV])
    double    fSigmaEnergyGeV;

    // Reference point for Beta and Emittance (and SigmaU and SigmaUP)
    TVector3D fSigmaAt;

    TVector3D fHorizontalDirection;
    TVector3D fVerticalDirection;
};




inline std::ostream& operator << (std::ostream& os, TParticleBeam const& o)
{
  os << "Name:             " << o.GetName() << "\n"
     << "Weight:           " << o.GetWeight() << "\n"
     << "X0:               " << o.GetX0() << "\n"
     << "U0:               " << o.GetU0() << "\n"
     << "T0:               " << o.GetT0() << " [m]  " << o.GetT0() / TOSCARSSR::C() << " [s]\n"
     << "E0:               " << o.GetE0() << "\n"
     << "SigmaE:           " << o.GetSigmaEnergyGeV() << "\n"
     << "Current           " << o.GetCurrent() << "\n"
     << "Emittance         " << o.GetEmittance() << "\n"
     << "V-direction       " << o.GetVerticalDirection() << "\n"
     << "H-direction       " << o.GetHorizontalDirection() << "\n"
     << "BeamDistribution  " << o.GetBeamDistributionName() << "\n"
     << "TwissBeta         " << o.GetTwissBeta() << "\n"
     << "TwissAlpha        " << o.GetTwissAlpha() << "\n"
     << "TwissGamma        " << o.GetTwissGamma() << "\n"
     << "Twiss Lattice Ref " << o.GetTwissLatticeReference() << "\n"
     << "TwissBetaX0       " << o.GetTwissBetaX0() << "\n"
     << "TwissAlphaX0      " << o.GetTwissAlphaX0() << "\n"
     << "TwissGammaX0      " << o.GetTwissGammaX0() << "\n"
     << "Eta               " << o.GetEta() << "\n";

  return os;
}

















#endif
