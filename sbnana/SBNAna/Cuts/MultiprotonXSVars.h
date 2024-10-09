#ifndef MULTIPROTONXSVARS_H
#define MULTIPROTONXSVARS_H

//###########  Vars and Cuts for the ICARUS NuMI 1mu>1p0pi XS analysis 
#include <iostream>

#include "sbnana/CAFAna/Core/Cut.h"
#include "sbnana/CAFAna/Core/Var.h"
#include "sbnana/CAFAna/Core/Var.h"
#include "sbnana/CAFAna/Core/MultiVar.h"
#include "sbnana/SBNAna/Cuts/TruthCuts.h"
#include "sbnana/SBNAna/Cuts/NuMIXSecSysts.h"
#include "sbnana/SBNAna/Vars/NuMIFlux.h"
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

#include "TVector3.h"

namespace ana{

/////////////////////////////////////////////////
// SLICE LEVEL

//! Product of PPFX weight, single pion production reweight from MINERvA data, and track data-driven track split reweight
extern const Var kTotalCVWeight;

extern const Var kOne;

extern const Var kZero;

extern const Var kIsFHC;

extern const Var kIsRHC;

extern const Cut isInFV_Vtx;

extern const Cut isInAV_Vtx;

extern const Cut kFV;

/////////////////////////////////////////////////
// Slice Vars and Cuts (that aren't already in another header or in the fiducial volume stuff above)
/////////////////////////////////////////////////

extern const Cut kIsNC;

extern const Cut kIsCC;

extern const Var kRecoMuonIdx;

extern const Cut kHasMuon;

extern const Var kRecoMuonPNew;

extern const Cut kRecoMuonContained;

extern const Var kRecoProtonIdx;

extern const Var kRecoProtonP;

extern const Cut kRecoProtonIsTrueProton;

extern const Var kProtonMult_Reco_wStubs;

extern const Cut kHasScdyPionTrack;

extern const Cut kThreePrimaryTracks;

extern const SpillCut kTriggeredEvent;

extern const Cut k1mu2p0pi;

extern const SpillCut kSpill1mu2p0pi; 

extern const Cut k1mu2p0piThreshold;

extern const Cut k1mu2pNpi;

extern const Cut k1mu0p0pi;

extern const Cut k1mu0pNpi;

extern const Cut k1mu1p0pi;

extern const Cut k1mu1pNpi;

extern const Cut k1mu3p0pi;

extern const Cut k1mu3pNpi;

extern const Cut k1mu2p;

extern const Cut k1mu0p;

extern const Cut k1mu1p;

extern const Cut k1mu3p;

extern const MultiVar kSelectedCPiP;

extern const Cut kIsSignal;

extern const MultiVar kSelectedPi0P;

extern const Var NonSignalPrinting;

extern const Var kNPrimary;

extern const Var kNPrimaryProtons;

extern const Var kNPrimaryPions;

extern const Var kNPrimaryMuons;

extern const Var kNPrimaryNeutrons;

extern const Var kNPrimaryPi0s;

extern const Var kNPrimaryPiMins;

extern const Var kNPrimaryOther;

extern const Var kNPrimaryPFPs;

extern const Cut kFivePrimaryPFPs;

extern const MultiVar kPrimaryShwEnergies;

extern const Var kMinShwEnergy;

extern const MultiVar kPrimaryTrkScores;

extern const MultiVar kPrimaryCollectionSums;

extern const Var kMinCollectionSum;

extern const MultiVar kPrimaryPDGs;

extern const Cut kLeadingProtonThreshold;

extern const Cut kPreselection;

extern const Cut k1mu1p1X;

extern const Var kScndProtonIdx;

extern const Var kThirdProton;

extern const MultiVar kRecoProtonIndices;

extern const Var kSidebandPion;

extern const Cut kThirdRENAMEIFNEEDEDPrimaryContainedTruth;

extern const Cut kMIPLikeTrack;

extern const Cut kProtonLikeTrack;

extern const Cut kScndProtonThreshold;

extern const Cut kTruePion;

extern const Cut kTrueMuon;

extern const Cut kTrueProton;

extern const Cut kMuPiMixing;

extern const Var kRecoMuonThetaNuMI;

extern const Var kRecoMuonPNuMI;

extern const Var kRecoMuonTrackLength;

extern const Var kRecoMuonTrueLength;

extern const Cut kForwardMuon;

extern const Var kNRecoProtons;

extern const Var kNTrueProtons;

extern const Var kNProtonsOverCap;

extern const Var kNSoftProtons;

extern const Var kNPions;

extern const Var kNTruePiPlus;

extern const Var kTruePiPlusP;

extern const Var kNTruePiMinus;

extern const Var kTruePiMinusP;

extern const Var kNPi0s;

extern const Var kNTrueNeutrons;

extern const Var kTrueNeutronP;

extern const Cut kCountPrimaryPFPs;

extern const Cut kHasSidebandPion;

extern const Var kHasSidebandPionVar;

extern const Cut kContainedSidebandHadron;

extern const Cut kThreeRecoProtons;

extern const Cut kPionSidebandPFPCut;

extern const Var kNoExtraMIP;

extern const Cut kNoExtraMIPCut;

extern const Var kExtraMIPIdx;

extern const Var kNoExtraMIPPrimariesOnly;

extern const Var kNoExtraShower;

extern const Cut kNoExtraShowerCut;

extern const Var kExtraShowerIdx;

extern const Var kNoExtraShowerPrimariesOnly;

extern const Cut kHadronicContainment;

extern const Var kHadronicContainmentPrimariesOnly;

extern const Var kFurthestFromVtx;

extern const Var kExtraPrimaryIndex;

extern const Var kExtraPrimaryLinFitLength;

extern const Cut kExtraPrimaryLinFitLengthCut;

extern const Cut kDynamicLengthCut;

extern const Cut kScndProtonCandidate;

extern const Cut kPionSidebandBase;

extern const Cut kPionSideband1p;

extern const Cut kPionSidebandNp;

extern const Cut kSignalOrSideband;

extern const Var kCutType;

extern const Var kNPFPs;

extern const Var kPrintNonprimaryPFPs;

extern const Var kExtraEnergy;

extern const Var kVertexX;

extern const Var kVertexY;

extern const Var kVertexZ;

extern const Var kDistToEdgeX;

extern const Var kDistToEdgeY;

extern const Var kDistToEdgeZ;

extern const Var kRecoMuonTruePDG;

extern const Var kRecoProtonTruePDG;

extern const Var kScndProtonTruePDG;

extern const Var kRecoMuonChi2Muon;

extern const Var kRecoMuonChi2Proton;

extern const Var kRecoMuonChi2Pion;

extern const Var kRecoProtonChi2Muon;

extern const Var kRecoProtonChi2Proton;

extern const Var kRecoProtonChi2Pion;

extern const Var kScndProtonChi2Muon;

extern const Var kScndProtonChi2Proton;

extern const Var kScndProtonChi2Pion;

extern const Var kScndProtonChi2ProtonNEGATIVE;

extern const Var kScndProtonNDaughters;

extern const Var kRecoProtonNDaughters;

extern const Var kTotalNDaughters;

extern const Var kScndProtonLength;

extern const Var kScndProtonTrueLength;

extern const Var kScndProtonLengthResid;

extern const Var kScndProtonPionP;

extern const Var kScndProtonTrueP;

extern const Var kScndProtonPionPResid;

extern const Cut kTruePiPlusUncontained;

extern const Var kRecoProtonEnergyPerLength;

extern const Var kScndProtonEnergyPerLength;

extern const Var kRecoProtonChi2Comp;

extern const Var kScndProtonChi2Comp;

extern const Var kRecoMuonTrackScore;

extern const Var kRecoProtonTrackScore;

extern const Var kScndProtonTrackScore;

extern const Var kRecoProtonTrackLength;

extern const Var kScndProtonTrackLength;

extern const MultiVar kScndProtonDaughterPDG;

extern const MultiVar kScndProtonDaughterTrkScore;

extern const Var kScndProtonDirChange;

extern const Var kScndProtonProtonP;

extern const Var kScndProtonProtonPResid;

extern const Var kRecoMuonEndProcess;

extern const Var kRecoProtonEndProcess;

extern const Var kScndProtonEndProcess;

extern const Var kRecoMuonStartProcess;

extern const Var kRecoProtonStartProcess;

extern const Var kScndProtonStartProcess;

extern const Var kRecoMuonParentPDG;

extern const Var kRecoProtonParentPDG;

extern const Var kScndProtonParentPDG;

extern const Var kRecoProtonChgEndFrac;

extern const Var kScndProtonChgEndFrac;

extern const Var kRecoProtonChgFracSpread;

extern const Var kScndProtonChgFracSpread;

extern const Var kRecoProtonLinFitDiff;

extern const Var kScndProtonLinFitDiff;

extern const Var kRecoProtonLinFitLen;

extern const Var kScndProtonLinFitLen;

extern const Var kRecoProtonLinFitGapLen;

extern const Var kScndProtonLinFitGapLen;

extern const Var kRecoProtonLinFitRMS;

extern const Var kScndProtonLinFitRMS;

extern const Var kRecoProtonOpenAngleDiff;

extern const Var kScndProtonOpenAngleDiff;

extern const Var kRecoProtonPCA2Ratio;

extern const Var kScndProtonPCA2Ratio;

extern const Var kRecoProtonPCA3Ratio;

extern const Var kScndProtonPCA3Ratio;

extern const Var kRecoProtonVtxDist;

extern const Var kScndProtonVtxDist;

extern const Var kRecoMuonShowerLength;

extern const Var kRecoProtonShowerLength;

extern const Var kScndProtonShowerLength;

extern const Var kRecoProtonShowerEnergy;

extern const Var kScndProtonShowerEnergy;

extern const Var kExtraPrimaryShowerLength;

extern const Var kExtraPrimaryShowerEnergy;

extern const Var kExtraPrimaryShowerdEdx;

extern const Var kExtraPrimaryShowerConvGap;

extern const Var kExtraPrimaryShowerDensity;

extern const Var kExtraPrimaryShowerOpenAngle;

extern const Var kExtraPrimaryChi2Muon;

extern const Var kExtraPrimaryChi2Proton;

extern const Var kExtraPrimaryTruePDG;

extern const Var kExtraPrimaryStartProcess;

extern const Var kSplitMuonDistance;

extern const Var kTrueSplitMuon;

extern const Var kTrueNeutrinoEnergy;

extern const Cut kSignalQE;

extern const Cut kSignalMEC;

extern const Cut kSignalRes;

extern const Cut kSignalDIS;

extern const Cut kSignalCoh;

extern const Cut kCCBG;

extern const Var kCCBGPrintout;

extern const Var kPrintPrimaryPFPs;
    
extern const Var kCCBGReason;

extern const Cut kCat0;

extern const Cut kCat1;

extern const Cut kCat2;

extern const Cut kCat3;

extern const Cut kCat4;

extern const Cut kCat5;

extern const Cut kCat6;

extern const Cut kCat7;

extern const Cut kCat8;

extern const Cut kCat9;

extern const Cut kCat10;

extern const Cut kCat11;

extern const Cut kOOFV;

extern const Cut kCCOther;

extern const Var kPrintCCBG;

extern const Var kRecoMuonIsPrimary;

extern const Var kRecoProtonIsPrimary;

extern const Var kScndProtonIsPrimary;

extern const Var kCountTeVParticles;

extern const Var kPrintTeVParticles;

extern const Var kNeutronEnergy;

extern const Var kRecoMuonDistToVtx;

extern const Var kRecoProtonDistToVtx;

extern const Var kScndProtonDistToVtx;

extern const MultiVar kScndProtondEdx;

extern const MultiVar kScndProtonRR;

extern const MultiVar kRecoProtondEdx;

extern const MultiVar kRecoProtonRR;

extern const MultiVar kBothProtonsdEdx;

extern const MultiVar kBothProtonsRR;

extern const Var kTrueNeutrinoDirection;

extern const Var kThirdProtonP;

extern const Var kThirdProtonTruthP;

extern const Var kThirdProtonPResid;

extern const Var kHadronicOpeningAngle;

extern const Var kHadronicOpeningAngleTruth;

extern const Var kHadronicOpeningAngleOldTruth;

extern const Var kHadronicOpeningAngleResid;

extern const Var kMuonHadronAngle;

extern const Var kMuonHadronAngleTruth;

extern const Var kMuonHadronAngleFullTruth;

extern const Var kMuonHadronAngleOldTruth;

extern const Var kMuonHadronAngleResid;

extern const Var kMuonPionAngle;

extern const Var kMuonPionAngleTruth;

extern const Var kSidebandPionThetaNuMI;

extern const Var kSidebandPionTruthThetaNuMI;

extern const Var kRecoMuonTruthThetaNuMI;

extern const Var kRecoMuonThetaNuMIResid;

extern const Var kRecoMuonStartDotEnd;

extern const Var kRecoProtonThetaNuMI;

extern const Var kRecoProtonTruthThetaNuMI;

extern const Var kRecoProtonThetaNuMIResid;

extern const Var kRecoProtonTruthDirDotReco;

extern const Var kRecoProtonTruthDirDotRecoWide;

extern const Var kRecoProtonStartDotEnd;

extern const Var kScndProtonThetaNuMI;

extern const Var kScndProtonTruthThetaNuMI;
    
extern const Var kScndProtonThetaNuMIResid;

extern const Var kScndProtonTruthDirDotReco;

extern const Var kScndProtonTruthDirDotRecoWide;

extern const Var kScndProtonStartDotEnd;

extern const Var kRecoProtonTruthP;

extern const Var kRecoProtonPResid;

extern const Var kExtraProtonTruthP;

extern const Var kExtraPionTruthP;

extern const Var kExtraPionEndProcess;

extern const Var kSidebandPionStartProcess;

extern const Var kSidebandPionEndProcess;

extern const Var kSidebandPionParentPDG;

extern const Var kExtraPionTrueLength;

extern const Var kExtraPi0TruthP;

extern const Var kSidebandPionP;

extern const Var kSidebandPionTruthP;

extern const Var kSidebandPionPResid;

extern const Var kSidebandPionTrackLength;

extern const Var kSidebandPionTruePDG;

extern const Var kLeadingProtonPFrac;

extern const Var kLeadingProtonPFracOldTruth;

extern const Var kLeadingProtonPFracTruth;

extern const Var kLeadingProtonPFracFullTruth;

extern const Var kRecoMuonTruthP;

extern const Var kRecoMuonPResid;

extern const Var kRecoMuonPInverseResid;

extern const Var kRecoMuonTruthDirDotReco;

extern const Var kRecoMuonTruthDirDotRecoWide;

extern const Cut kBadMuon;

extern const Var kPrintBadMuon;

extern const Var kEHad_Proton;

extern const Var kEHad_CheatingProtonPs;

extern const Var kEHad_Pion;

extern const Var kEHad;

extern const Var kEHad_Truth;

extern const Var kEAvail_Truth;

extern const Var kHasKaon;

extern const Var kHasLambda;

extern const Var kHasSigma;

extern const Var kHasAr;

extern const Var kPrintExotics;

extern const Var kQ2;

extern const Var kQ2_CheatingMuonP;

extern const Var kQ2_CheatingAll;

extern const Var kQ2_Truth;

extern const Var kq3;

extern const Var kq3_Truth;

extern const Var kRecoENu;

extern const Var kW_Proton;

extern const Var kW_Pion;

extern const Var kW;

extern const Var kW_ProtonTruth;

extern const Var kW_PionTruth;

extern const Var kW_Truth;

extern const Var kW_Resid;

extern const Var kW_ProtonResid;

extern const Var kW_PionResid;

extern const Var kW_Exp;

extern const Var kW_Exp_CheatingMuonP;

extern const Var kW_Exp_CheatingAll;

extern const Var kW_Exp_Truth;

extern const Var kW_Exp2;

extern const Var kW_Exp2_CheatingMuonP;

extern const Var kW_Exp2_CheatingAll;

extern const Var kW_Exp2_Truth;

extern const Var kGENIEQ2;

extern const Var kGENIEW;

extern const Var kW_ExpResid;

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Nasty TKI section, each in proton mass, pion mass, and using truth variables

extern const Var kDeltaPT_Single;

extern const Var kDeltaPT_Proton;

extern const Var kDeltaPT_CheatingMuon;

extern const Var kDeltaPT_CheatingAngles;

extern const Var kDeltaPT_Pion;

extern const Var kDeltaPT_ProtonTruth;

extern const Var kDeltaPT_PionTruth;

extern const Var kDeltaPT;

extern const Var kDeltaPT_OldTruth;

extern const Var kDeltaPT_FullTruth;

extern const Var kDeltaPT_Truth;

extern const Var kDeltaPT_Resid;

extern const Var kDeltaPT_ProtonResid;

extern const Var kDeltaPT_CheatingMuonResid;

extern const Var kDeltaPT_CheatingAnglesResid;

extern const Var kDeltaAlphaT_Single;

extern const Var kDeltaAlphaT_Proton;

extern const Var kDeltaAlphaT_CheatingMuon;

extern const Var kDeltaAlphaT_CheatingAngles;

extern const Var kDeltaAlphaT_Pion;

extern const Var kDeltaAlphaT_ProtonTruth;

extern const Var kDeltaAlphaT_PionTruth;

extern const Var kDeltaAlphaT;

extern const Var kDeltaAlphaT_OldTruth;

extern const Var kDeltaAlphaT_FullTruth;

extern const Var kDeltaAlphaT_Truth;

extern const Var kDeltaAlphaT_Resid;

extern const Var kDeltaAlphaT_ProtonResid;

extern const Var kDeltaAlphaT_CheatingMuonResid;

extern const Var kDeltaAlphaT_CheatingAnglesResid;

extern const Var kDeltaPhiT_Single;

extern const Var kDeltaPhiT_Proton;

extern const Var kDeltaPhiT_CheatingMuon;

extern const Var kDeltaPhiT_CheatingAngles;

extern const Var kDeltaPhiT_Pion;

extern const Var kDeltaPhiT_ProtonTruth;

extern const Var kDeltaPhiT_PionTruth;

extern const Var kDeltaPhiT;

extern const Var kDeltaPhiT_OldTruth;

extern const Var kDeltaPhiT_FullTruth;

extern const Var kDeltaPhiT_Truth;

extern const Var kDeltaPhiT_Resid;

extern const Var kDeltaPhiT_ProtonResid;

extern const Var kDeltaPhiT_CheatingMuonResid;

extern const Var kDeltaPhiT_CheatingAnglesResid;

extern const Var kDeltaPTT_Single;

extern const Var kDeltaPTT_Proton;

extern const Var kDeltaPTT_CheatingMuon;

extern const Var kDeltaPTT_CheatingAngles;

extern const Var kDeltaPTT_Pion;

extern const Var kDeltaPTT_ProtonTruth;

extern const Var kDeltaPTT_PionTruth;

extern const Var kDeltaPTT;

extern const Var kDeltaPTT_OldTruth;

extern const Var kDeltaPTT_FullTruth;

extern const Var kDeltaPTT_Truth;

extern const Var kDeltaPTT_Resid;

extern const Var kDeltaPTT_ProtonResid;

extern const Var kDeltaPTT_CheatingMuonResid;

extern const Var kDeltaPTT_CheatingAnglesResid;

//! Truth Section
//! Here we use t<NAME> instead of k<NAME> convention to distinguish TruthCuts and TruthVars from Cuts and Vars

extern const TruthCut tIsSignal;

extern const TruthVar tTrueNeutrinoPDG;

extern const Var kTrueNeutrinoPDG;

extern const TruthVar tTrueNeutrinoE;

extern const Var kTrueNeutrinoE;

extern const TruthVar tTrueMuonP;

extern const Var kTrueMuonP;

extern const TruthVar tTrueLeadingProtonP;

extern const Var kTrueLeadingProtonP;

extern const TruthVar tTrueSecondProtonP;

extern const Var kTrueSecondProtonP;

extern const TruthVar tTrueThirdProtonP;

extern const TruthVar tTrueChargedPionP;

extern const TruthVar tTrueLeadingProtonPFrac;

extern const Var kTrueLeadingProtonPFrac;

extern const TruthVar tTrueHadronicOpeningAngle;

extern const Var kTrueHadronicOpeningAngle;

extern const TruthVar tTrueMuonHadronAngle;

extern const Var kTrueMuonHadronAngle;

extern const TruthVar tTrueDeltaPT;

extern const Var kTrueDeltaPT;

extern const TruthVar tTrueDeltaAlphaT;

extern const Var kTrueDeltaAlphaT;

extern const TruthVar tTrueDeltaPhiT;

extern const Var kTrueDeltaPhiT;

extern const TruthVar tTrueDeltaPTT;

extern const Var kTrueDeltaPTT;

extern const TruthVar tTrueEHad;

extern const Var kTrueEHad;

extern const TruthVar tTrueEAvail;

extern const Var kTrueEAvail;

extern const TruthVar tTrueQ2;

extern const Var kTrueQ2;

extern const TruthVar tTrueq3;

extern const Var kTrueq3;

extern const TruthVar tTrueNProtons;

extern const Var kTrueNProtons;

extern const TruthVar tTrueNeutronEnergy;

extern const Var kTrueNeutronEnergy;

extern const TruthVar tContained;

extern const TruthVar tStopping;

extern const TruthVar tAllStopping;








extern const Var kCategory;

extern const Var kClassLabel;

extern const Var kPrintExtraPrimaries;

extern const Var kContained;

extern const Var kGENIEMode;

std::vector<std::string> GetGENIEMultisigmaKnobNames(){

  return {
/* ZExpA_VariationResponse
"ZExpA1CCQE",
"ZExpA2CCQE",
"ZExpA3CCQE",
"ZExpA4CCQE",
*/
"RPA_CCQE",
"CoulombCCQE",
"NormCCMEC",
"NormNCMEC",
//"DecayAngMEC", --> MirrorSyst!
/* NCEL_VariationResponse
"MaNCEL",
"EtaNCEL",
*/
/* CCRES_VariationResponse
"MaCCRES",
"MvCCRES",
*/
/* NCRES_VariationResponse
"MaNCRES",
"MvNCRES",
*/
"NonRESBGvpCC1pi",
"NonRESBGvpCC2pi",
"NonRESBGvpNC1pi",
"NonRESBGvpNC2pi",
"NonRESBGvnCC1pi",
"NonRESBGvnCC2pi",
"NonRESBGvnNC1pi",
"NonRESBGvnNC2pi",
"NonRESBGvbarpCC1pi",
"NonRESBGvbarpCC2pi",
"NonRESBGvbarpNC1pi",
"NonRESBGvbarpNC2pi",
"NonRESBGvbarnCC1pi",
"NonRESBGvbarnCC2pi",
"NonRESBGvbarnNC1pi",
"NonRESBGvbarnNC2pi",
"RDecBR1gamma",
"RDecBR1eta",
"NormCCCOH",
"NormNCCOH",
/* DISBY_VariationResponse
"AhtBY",
"BhtBY",
"CV1uBY",
"CV2uBY",
*/
/* FIS_pi_VariationResponse
"MFP_pi",
"FrCEx_pi",
"FrInel_pi",
"FrAbs_pi",
"FrPiProd_pi",
*/
/* FSI_N_VariationResponse
"MFP_N",
"FrCEx_N",
"FrInel_N",
"FrAbs_N",
"FrPiProd_N",
*/
  };

}

std::vector<std::string> GetGENIEDependentKnobNames(){
  return {
"ZExpAVariationResponse",
"NCELVariationResponse",
"CCRESVariationResponse",
"NCRESVariationResponse",
"DISBYVariationResponse",
"FSI_pi_VariationResponse",
"FSI_N_VariationResponse",
"reinteractions_piminus_Geant4",
"reinteractions_piplus_Geant4",
"reinteractions_proton_Geant4",

  };
}

std::vector<std::string> GetFluxKnobNames(unsigned nPCs) {
  std::vector<std::string> knobs = {
"beam_div",
"beam_power",
"beam_shift_x",
"beam_spot",
"horn1_x",
"horn1_y",
"horn_current_plus",
"water_layer",
"beam_shift_y"
  };
  for ( unsigned i = 0; i < nPCs; i++ ) knobs.push_back( "flux_PCA_" + std::to_string(i) );
  knobs.push_back("flux_stat");

  return knobs;
}

}//END NAMESPACE
#endif
