#ifndef GUARD_OSCARSTH_h
#define GUARD_OSCARSTH_h
////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Thu Feb  2 11:55:24 EST 2017 (minus a day or so)
//
////////////////////////////////////////////////////////////////////

#include "TVector2D.h"
#include "TOMATH.h"
#include "TParticleBeamContainer.h"
#include "TSpectrumContainer.h"
#include "TSurfacePoints.h"
#include "T3DScalarContainer.h"


class OSCARSTH
{
  public:
    OSCARSTH ();
    ~OSCARSTH ();

    double UndulatorK (double const BFieldMax,
                       double const Period) const;

    double UndulatorBField (double const K,
                            double const Period) const;

    double UndulatorPeriod (double const BField,
                            double const K) const;

    double DipoleCriticalEnergy (double const BField) const;

    void DipoleSpectrumEnergy (double const BField, 
                               TSpectrumContainer& Spectrum,
                               double const Angle) const;

    void DipoleSpectrumAngle (double const BField, 
                              TSpectrumContainer& Spectrum,
                              double const Energy_eV) const;

    void DipoleSpectrumEnergyAngleIntegrated (double const BField, 
                                              TSpectrumContainer& Spectrum) const;

    double DipoleSpectrum (double const BField,
                           double const BeamEnergy_GeV,
                           double const Angle,
                           double const Energy_eV) const;

    double DipoleSpectrumAngleIntegrated (double const BField,
                                          double const BeamEnergy_GeV,
                                          double const Energy_eV) const;

    double UndulatorFlux (double const BField,
                          double const Period,
                          double const NPeriods,
                          double const BeamEnergy,
                          double const AngleV,
                          double const AngleH,
                          double const Energy_eV) const;


    TVector2D UndulatorFluxOnAxisK (double const K,
                                    double const Period,
                                    double const NPeriods,
                                    int    const Harmonic) const;

    TVector2D UndulatorFluxOnAxisB (double const BField,
                                    double const Period,
                                    double const NPeriods,
                                    int    const Harmonic) const;

    TVector2D UndulatorFluxB (double const BField,
                              double const Period,
                              int    const NPeriods,
                              int    const Harmonic
                              ) const;

    TVector2D UndulatorFluxK (double const K,
                              double const Period,
                              int    const NPeriods,
                              int    const N
                              ) const;

    double UndulatorFluxWeak (double const K,
                              double const Period,
                              double const NPeriods,
                              double const BeamEnergy,
                              int const Harmonic) const;

    void DipoleBrightness (double const BField,
                           TSpectrumContainer& SpectrumContainer) const;

    double UndulatorEnergyAtHarmonicK (double const K,
                                       double const Period,
                                       int    const Harmonic) const;

    double UndulatorEnergyAtHarmonicB (double const BField,
                                       double const Period,
                                       int    const Harmonic) const;

    TVector2D UndulatorBrightnessK (double const BField,
                                    double const Period,
                                    int    const NPeriods,
                                    int    const N) const;

    TVector2D UndulatorBrightnessB (double const BField,
                                    double const Period,
                                    int    const NPeriods,
                                    int    const N) const;


    void WigglerFluxK (double         const  K,
                       double         const  Period,
                       int            const  NPeriods,
                       TSurfacePoints const& Surface,
                       double         const  Energy_eV,
                       T3DScalarContainer&   FluxContainer) const;

    void WigglerFluxK (double         const  K,
                       double         const  Period,
                       int            const  NPeriods,
                       TSurfacePoints const& Surface,
                       double         const  Energy_eV,
                       T3DScalarContainer&   FluxContainer,
                       int            const  NThreads,
                       int            const  GPU) const;

    void WigglerFluxB (double         const  K,
                       double         const  Period,
                       int            const  NPeriods,
                       TSurfacePoints const& Surface,
                       double         const  Energy_eV,
                       T3DScalarContainer&   FluxContainer,
                       int            const  NThreads,
                       int            const  GPU) const;

    void WigglerFluxKPoints (double         const  K,
                             double         const  Period,
                             int            const  NPeriods,
                             TSurfacePoints const& Surface,
                             double         const  Energy_eV,
                             T3DScalarContainer&   FluxContainer,
                             size_t const iFirst,
                             size_t const iLast,
                             bool& Done
                            ) const;





    // Functions related to the particle beam
    /*
    TParticleBeam& SetParticleBeam (std::string const& Beam,
                                    std::string const& Name = "default_name");

    TParticleBeam& SetParticleBeam (double const Energy_GeV,
                                    double const Current,
                                    TVector2D const& Beta = TVector2D(0, 0),
                                    TVector2D const& Emittance = TVector2D(0, 0),
                                    double const SigmaEnergyGeV = 0,
                                    TVector2D const& Eta = TVector2D(0, 0),
                                    std::string const& Name = "default_name");
    */


    // Functions related to the particle beam(s)
    TParticleBeam& AddParticleBeam (std::string const& Type,
                                    std::string const& Name,
                                    TVector3D const& X0,
                                    TVector3D const& V0,
                                    double const Energy_GeV,
                                    double const T0,
                                    double const Current,
                                    double const Weight,
                                    double const Charge = 0,
                                    double const Mass = 0);

    TParticleBeam& AddParticleBeam (std::string const& Beam,
                                    std::string const& Name,
                                    double const Weight = 1);

    void ClearParticleBeams ();
    TParticleBeam& GetParticleBeam ();

    bool CheckBeam () const;


    // Global threads and GPU settings
    bool SetUseGPUGlobal (int const);
    int  GetUseGPUGlobal () const;
    int  CheckGPU () const;
    void SetNThreadsGlobal (int const);

  private:
    TParticleBeam fParticleBeam;
    TParticleBeamContainer fParticleBeamContainer;

    // Global thread and GPU settings
    int fNThreadsGlobal;
    bool fUseGPUGlobal;

};














#endif
