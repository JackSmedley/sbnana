// BH - 2023
// HEAVILY based on NuMI flux syst, and thanks to Tony Wood for discussing the right histograms to use

#pragma once

#include "sbnana/CAFAna/Core/Var.h"

#include <string>

class TH1;

namespace ana
{

  //===============================================================
  // 05/19/25 JK) A base ppfx correction from "2023-07-31_out_450.37_7991.98_79512.66_QEL11.root" which does not include G3Chase or G4Update
  class NuMIPpfxFluxWeight
  {
  public:
    NuMIPpfxFluxWeight();
    ~NuMIPpfxFluxWeight();
    mutable TH1* fWeight[2][2][2]; // [fhc/rhc][nue/numu][nu/nubar]

    static NuMIPpfxFluxWeight& Instance();

  protected:
    std::string fFluxFilePath;
  };

  // set up to use the flux weight
  extern const Var kGetNuMIFluxWeight;
  extern const TruthVar kGetTruthNuMIFluxWeight;

  //===============================================================
  // 05/19/25 JK) This one still uses "2023-07-31_out_450.37_7991.98_79512.66_QEL11.root" for the ppfx correction,
  //              but provides a CV correction to G3Chase and G4Update from "g3Chase_weights_rewritten.root"
  class NuMIPpfxFluxWeightG3Chase
  {
  public:
    // FOR NOW WE ONLY HAVE FHC...
    // parent idx --> pi+/-, K+/-, mu, K0L
    NuMIPpfxFluxWeightG3Chase();
    ~NuMIPpfxFluxWeightG3Chase();

    double GetWeightFromSRTrueInt(const caf::SRTrueInteractionProxy* nu, const bool applyKaonRW) const;
    unsigned int ParentPDGToIdx(int pdg) const;

    mutable TH1* fWeight[2][2][2]; // [fhc/rhc][nue/numu][nu/nubar]
    mutable TH1* fWeightG3Chase[2][2][4]; // [nue/numu][nu/nubar][parent idx]
    mutable TH1* fWeightG4Kaon[2]; // [nu/nubar]

    static NuMIPpfxFluxWeightG3Chase& Instance();

  protected:
    std::string fFluxFilePath;
    std::string fFluxFilePathG3Chase;
    std::string fFluxFilePathG4Kaon;
  };

  // set up to use the flux weight with G3Chase param controlling the concrete
  extern const Var kGetNuMIFluxWeightG3Chase;
  extern const TruthVar kGetTruthNuMIFluxWeightG3Chase;

  extern const Var kGetNuMIFluxWeightUpdated;
  extern const TruthVar kGetTruthNuMIFluxWeightUpdated;


  //===============================================================
  // 05/19/25 JK) This one now uses "2024-10-03_out_450.37_7991.98_79512.66.root",
  //              a newer file than "2023-07-31_out_450.37_7991.98_79512.66_QEL11.root"
  //              Here the ppfx correction and uncertainties are all re-calculated with G3chase and G4update
  //              The CV correction here includes both G3Chase and G4Update
  class NuMIPpfxFluxWeightG4Update
  {
  public:
    // FOR NOW WE ONLY HAVE FHC...
    // parent idx --> pi+/-, K+/-, mu, K0L
    NuMIPpfxFluxWeightG4Update();
    ~NuMIPpfxFluxWeightG4Update();

    double GetWeightFromSRTrueInt(const caf::SRTrueInteractionProxy* nu) const;
    unsigned int ParentPDGToIdx(int pdg) const;

    mutable TH1* fWeight[2][2][2]; // [fhc/rhc][nue/numu][nu/nubar]
    mutable TH1* fWeightG4Update[2][2][2][4]; // [fhc/rhc][nue/numu][nu/nubar][parent pid (pipm/kpm/k0l/mu)]

    static NuMIPpfxFluxWeightG4Update& Instance();

  protected:
    std::string fFluxFilePath;
  };

  //static const NuMIPpfxFluxWeightG4Update FluxWeightNuMIG4Update;
  extern const Var kGetNuMIFluxWeightG4Update;
  extern const TruthVar kGetTruthNuMIFluxWeightG4Update;

  //===============================================================
  // 05/19/25 JK) This one uses "2025-04-08_out_450.37_7991.98_79512.66.root",
  //              to provide a CV correction of beam width setting 1.5mm/1.4mm;
  //              NuMI2023 reprocessing simulation was done with CV beam width of 1.4mm,
  //              but the data is more like 1.5mm.
  class NuMIBeamWidthCorrection
  {
  public:
    NuMIBeamWidthCorrection();
    ~NuMIBeamWidthCorrection();

    double GetWeightFromSRTrueInt(const caf::SRTrueInteractionProxy* nu) const;
    unsigned int ParentPDGToIdx(int pdg) const;

    // PPFX correction
    mutable TH1* fWeight[2][2][2]; // [fhc/rhc][nue/numu][nu/nubar]
    // Additional CV correction from NuMI reproc-to-PPFXCalculationNominal of this file
    mutable TH1* fWeightCVCorr[2][2][2][4]; // [fhc/rhc][nue/numu][nu/nubar][parent pid (pipm/kpm/k0l/mu)]

    static NuMIBeamWidthCorrection& Instance();

  protected:
    std::string fFluxFilePath;
  };

  extern const Var kGetNuMIBeamWidthCorrection;
  extern const TruthVar kGetTruthNuMIBeamWidthCorrection;





}
