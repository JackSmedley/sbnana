//###########  Vars and Cuts for the ICARUS NuMI 1mu>1p0pi XS analysis 

#include "sbnana/SBNAna/Cuts/MultiprotonXSVars.h"
#include "sbnana/SBNAna/Cuts/NuMIXSecSysts.h"

namespace ana{

double mMuon = 0.1056583755;
double mProton = 0.93827209;
double mNeutron = 0.93956542052;
double mPion = 0.13957039;


/////////////

// First attempts at NuMI nu direction vectors
TVector3 dFromNuMI(315120.380, 33644.912, 733632.532);
// also
//TVector3 NuDirection_NuMI(3.94583e-01, 4.26067e-02, 9.17677e-01);
//TVector3 NuDirection_NuMI(0.394654,0.042614,0.917841); //Slightly better normalized

/////////////////////////////////////////////////
// 10% UNBLINDED DATA
/*
const SpillCut k10PercentUnblinded([](const caf::SRSpillProxy* sr) {
  std::pair<int,int> thisRunEvt = {sr->hdr.run, sr->hdr.evt};

  bool unblindedRun1 = ( std::find(run1_10percent_RunEvtPairs.begin(), run1_10percent_RunEvtPairs.end(), thisRunEvt) != run1_10percent_RunEvtPairs.end() );
  bool unblindedRun2 = ( std::find(run2_10percent_RunEvtPairs.begin(), run2_10percent_RunEvtPairs.end(), thisRunEvt) != run2_10percent_RunEvtPairs.end() );

  return ( unblindedRun1 || unblindedRun2 );
  });
*/

/////////////////////////////////////////////////
// SLICE LEVEL

const Cut kIsNuSlice([](const caf::SRSliceProxy* slc) {
  return ( slc->truth.index >= 0 );
  });
 
const Cut kIsCosmic = ( !kIsNuSlice );

const Cut kNotClearCosmic([](const caf::SRSliceProxy* slc) {
    return !slc->is_clear_cosmic;
  });



//! Product of PPFX weight, single pion production reweight from MINERvA data, and track data-driven track split reweight
const Var kTotalCVWeight([](const caf::SRSliceProxy* slc) {
  return ( kGetNuMIFluxWeightG3Chase(slc) * kNuMISPPCVCorrection(slc) /** kNuMISplitTrackCVCorrection(slc)*/ );
  });

const Var kOne([](const caf::SRSliceProxy* slc) -> int {
  return 1;
  });

const Var kZero([](const caf::SRSliceProxy* slc) -> int {
  return 0;
  });

const Var kIsFHC([](const caf::SRSliceProxy* slc) -> int {
  return 1;
  });

const Var kIsRHC([](const caf::SRSliceProxy* slc) -> int {
  return 0;
  });

bool isInFV (double x, double y, double z)
{
  if ( std::isnan(x) || std::isnan(y) || std::isnan(z) ) return false;

  return (( ( x < -61.94 - 25 && x > -358.49 + 25 ) ||
	          ( x >  61.94 + 25 && x <  358.49 - 25 )) &&
	        (( y > -181.86 + 25 && y < 134.96 - 25 ) &&
	         ( z > -894.95 + 30 && z < 894.95 - 50 ) ));
}

const Cut isInFV_Vtx([](const caf::SRSliceProxy* sr)
		     {
		       const auto& vtx = sr->truth.position;

		       if ( std::isnan(vtx.x) || std::isnan(vtx.y) || std::isnan(vtx.z) ) return false;

		       return (( ( vtx.x < -61.94 - 25 && vtx.x > -358.49 + 25 ) ||
				             ( vtx.x >  61.94 + 25 && vtx.x <  358.49 - 25 )) &&
			             (( vtx.y > -181.86 + 25 && vtx.y < 134.96 - 25 ) &&
				            ( vtx.z > -894.95 + 30 && vtx.z < 894.95 - 50 ) ));
		     });

const Cut isInAV_Vtx([](const caf::SRSliceProxy* sr)
                     {
                       const auto& vtx = sr->truth.position;

                       if ( std::isnan(vtx.x) || std::isnan(vtx.y) || std::isnan(vtx.z) ) return false;

                       return (( ( vtx.x < -61.94 && vtx.x > -358.49 ) ||
                                 ( vtx.x >  61.94 && vtx.x <  358.49 )) &&
                               (( vtx.y > -181.86 && vtx.y < 134.96 ) &&
                                ( vtx.z > -894.95 && vtx.z < 894.95 ) ));
                     });

const Cut kFV([](const caf::SRSliceProxy* sr)
                     {
                       const auto& vtx = sr->vertex;

                       if ( std::isnan(vtx.x) || std::isnan(vtx.y) || std::isnan(vtx.z) ) return false;

                       return (( ( vtx.x < -61.94 - 25 && vtx.x > -358.49 + 25 ) ||
                                             ( vtx.x >  61.94 + 25 && vtx.x <  358.49 - 25 )) &&
                                     (( vtx.y > -181.86 + 25 && vtx.y < 134.96 - 25 ) &&
                                            ( vtx.z > -894.95 + 30 && vtx.z < 894.95 - 50 ) ));
                     });

/////////////////////////////////////////////////
// Slice Vars and Cuts (that aren't already in another header or in the fiducial volume stuff above)
/////////////////////////////////////////////////

/*
const Cut kIsNC([](const caf::SRSliceProxy* slc)
                 {
                   return slc->truth.isnc;
                 });


const Cut kIsCC([](const caf::SRSliceProxy* slc)
                 {
                   return slc->truth.iscc;
                 });
*/

// Redo the following two slice vars from NumuVars to use 10cm for containment (and to consider the other cryostat)

const Var kRecoMuonIdx([](const caf::SRSliceProxy* slc) -> int {
    // The (dis)qualification of a slice is based upon the track level information.
    float Longest(0);
    int PTrackInd(-1);
    for (std::size_t i(0); i < slc->reco.pfp.size(); ++i)
    {
        if( slc->reco.pfp.at(i).trackScore < 0.45 ) continue;
        auto const& trk = slc->reco.pfp.at(i).trk;
        if( std::isnan(trk.start.x) || trk.bestplane == -1 ) continue;

        // First we calculate the distance of each track to the slice vertex.
        const float Atslc = std::hypot(slc->vertex.x - trk.start.x,
                                       slc->vertex.y - trk.start.y,
                                       slc->vertex.z - trk.start.z);

        // We require that the distance of the track from the slice is less than
        // 10 cm and that the parent of the track has been marked as the primary.
        const bool AtSlice = ( Atslc < 10.0 && slc->reco.pfp.at(i).parent_is_primary);

        const float Chi2Proton = trk.chi2pid[2].chi2_proton;
        const float Chi2Muon = trk.chi2pid[2].chi2_muon;

        const bool Contained = (!isnan(trk.end.x) &&
                                ((trk.end.x < -61.94 - 10 && trk.end.x > -358.49 + 10) ||
                                 (trk.end.x >  61.94 + 10 && trk.end.x <  358.49 - 10)) &&
                                !isnan(trk.end.y) &&
                                ( trk.end.y > -181.86 + 10 && trk.end.y < 134.96 - 10 ) &&
                                !isnan(trk.end.z) &&
                                ( trk.end.z > -894.95 + 10 && trk.end.z < 894.95 - 10 ) );
        const bool MaybeMuonExiting = ( !Contained && trk.len > 50);
        const bool MaybeMuonContained = ( Contained && Chi2Proton > 60 && Chi2Muon < 30 && trk.len > 50 );
        if ( AtSlice && ( MaybeMuonExiting || MaybeMuonContained ) && trk.len > Longest )
    	  {
	        Longest = trk.len;
	        PTrackInd = i;
    	  }
    }
    return PTrackInd;
  });

const Cut kHasMuon([](const caf::SRSliceProxy* slc) {
    return ( kRecoMuonIdx(slc) >= 0 );
  });

const Var kRecoMuonPNew([](const caf::SRSliceProxy* slc) -> float {
    float p(-5.f);

    if ( kRecoMuonIdx(slc) >= 0 )
      {
        auto const& trk = slc->reco.pfp.at(kRecoMuonIdx(slc)).trk;
        const bool Contained = ( !isnan(trk.end.x) &&
				                         ((trk.end.x < -61.94 - 10 && trk.end.x > -358.49 + 10) ||
				                          (trk.end.x >  61.94 + 10 && trk.end.x <  358.49 - 10)) &&
                                 !isnan(trk.end.y) &&
                                 ( trk.end.y > -181.86 + 10 && trk.end.y < 134.96 - 10 ) &&
                                 !isnan(trk.end.z) &&
                                 ( trk.end.z > -894.95 + 10 && trk.end.z < 894.95 - 10 ) );
        if(Contained) p = trk.rangeP.p_muon;
        else p = trk.mcsP.fwdP_muon;
        if ( isnan(p) || p < 0.) std::cout << "NaN muon momentum! Contained?: " << Contained << ", track start: ("
            << trk.start.x << ", " << trk.start.y  << ", " << trk.start.z  << "), track end: (" <<  trk.end.x << ", " <<  trk.end.y << ", " <<  trk.end.z << ")" << std::endl;
      }
    return p;
  });

const Cut kRecoMuonContained([](const caf::SRSliceProxy* slc) {
    if ( kRecoMuonIdx(slc) >= 0 )
      {
        auto const& trk = slc->reco.pfp.at(kRecoMuonIdx(slc)).trk;
        const bool Contained = ( !isnan(trk.end.x) &&
                                 ((trk.end.x < -61.94 - 10 && trk.end.x > -358.49 + 10) ||
                                  (trk.end.x >  61.94 + 10 && trk.end.x <  358.49 - 10)) &&
                                 !isnan(trk.end.y) &&
                                 ( trk.end.y > -181.86 + 10 && trk.end.y < 134.96 - 10 ) &&
                                 !isnan(trk.end.z) &&
                                 ( trk.end.z > -894.95 + 10 && trk.end.z < 894.95 - 10 ) );
        if(Contained) return true;
        else return false;
      }
    return false;
  });

// PROTONs
// -- just picking the best proton track
const Var kRecoProtonIdx([](const caf::SRSliceProxy* slc) -> int {
    int primaryInd = kRecoMuonIdx(slc);

    int idxScdy = -1;
    float maxLength = -1;

    for ( unsigned int idxTrk = 0; idxTrk < slc->reco.pfp.size(); ++idxTrk ) {
      if ( (int)idxTrk == primaryInd ) continue; // need a different track...
  
//      if( slc->reco.pfp.at(idxTrk).trackScore < 0.5 ) continue; //TODO
      auto const& trk = slc->reco.pfp.at(idxTrk).trk;
      if( std::isnan(trk.start.x) || trk.bestplane == -1 ) continue;

      const float Atslc = std::hypot(slc->vertex.x - trk.start.x,
                                     slc->vertex.y - trk.start.y,
                                     slc->vertex.z - trk.start.z);
      const bool Contained = ( !isnan(trk.end.x) &&
                               ((trk.end.x < -61.94 - 10 && trk.end.x > -358.49 + 10) ||
                                (trk.end.x >  61.94 + 10 && trk.end.x <  358.49 - 10)) &&
                               !isnan(trk.end.y) &&
                               ( trk.end.y > -181.86 + 10 && trk.end.y < 134.96 - 10 ) &&
                               !isnan(trk.end.z) &&
                               ( trk.end.z > -894.95 + 10 && trk.end.z < 894.95 - 10 ) );

      const float Chi2Proton = trk.chi2pid[2].chi2_proton;
      const float Chi2Muon = trk.chi2pid[2].chi2_muon;


      float angle = 0.;
/*
      float angle = -5.0;
      if ( primaryInd >= 0 ) {
        const unsigned int idxPrim = (unsigned int)primaryInd;
        TVector3 muDir( slc->reco.pfp.at(idxPrim).trk.dir.x, slc->reco.pfp.at(idxPrim).trk.dir.y, slc->reco.pfp.at(idxPrim).trk.dir.z );
        TVector3 pDir( trk.dir.x, trk.dir.y, trk.dir.z );
        angle = TMath::Cos(muDir.Angle(pDir));
      }
*/

      // do we want to make the proton cut even tighter on PID
      if ( Atslc < 10.0 && slc->reco.pfp.at(idxTrk).parent_is_primary && Contained && /*Chi2Proton <= 100 && Chi2Muon >= 30*/  Chi2Proton < 50. && Chi2Muon != 0. && angle >= -0.9 && trk.len > maxLength ) {
        maxLength = trk.len;
        idxScdy = (int)idxTrk;
      }
    }

    return idxScdy;
  }); // kRecoProtonIdx

const Var kRecoProtonP([](const caf::SRSliceProxy* slc) -> float {
    float p(-5.f);

    if ( kRecoProtonIdx(slc) >= 0 )
    {
        auto const& trk = slc->reco.pfp.at(kRecoProtonIdx(slc)).trk;
        const bool Contained = ( !isnan(trk.end.x) &&
                                 ((trk.end.x < -61.94 - 10 && trk.end.x > -358.49 + 10) ||
                                  (trk.end.x >  61.94 + 10 && trk.end.x <  358.49 - 10)) &&
                                 !isnan(trk.end.y) &&
                                 ( trk.end.y > -181.86 + 10 && trk.end.y < 134.96 - 10 ) &&
                                 !isnan(trk.end.z) &&
                                 ( trk.end.z > -894.95 + 10 && trk.end.z < 894.95 - 10 ) );
        if(Contained) p = trk.rangeP.p_proton;
        else {
	        std::cout << "Currently kRecoProtonIdx requires a contained proton... Why am I trying to use MCS here??" << std::endl;
	        p = trk.mcsP.fwdP_proton;
      	}
    }
    return p;
  });

// Is selected proton (kRecoProtonIdx) actually a proton
const Cut kRecoProtonIsTrueProton([](const caf::SRSliceProxy* slc) {
    if ( kRecoProtonIdx(slc) >= 0 )
    {
      auto const& trk = slc->reco.pfp.at(kRecoProtonIdx(slc)).trk;

      if ( trk.truth.p.pdg == 2212 ) return true;
      else return false;
    }

    return false;
  });

// ---- reco w/ stubs
const Var kProtonMult_Reco_wStubs([](const caf::SRSliceProxy* slc) -> int {
    int primaryInd = kRecoMuonIdx(slc);
    int countP = 0;

    std::vector< int > proton_pfpids;

    for ( unsigned int idxTrk = 0; idxTrk < slc->reco.pfp.size(); ++idxTrk ) {
      if ( (int)idxTrk == primaryInd ) continue;

      if( slc->reco.pfp.at(idxTrk).trackScore < 0.5 ) continue;
      auto const& trk = slc->reco.pfp.at(idxTrk).trk;
      if( std::isnan(trk.start.x) || trk.bestplane == -1 ) continue;

      const float Atslc = std::hypot(slc->vertex.x - trk.start.x,
                                     slc->vertex.y - trk.start.y,
                                     slc->vertex.z - trk.start.z);
      const bool Contained = ( !isnan(trk.end.x) &&
                               ((trk.end.x < -61.94 - 10 && trk.end.x > -358.49 + 10) ||
                                (trk.end.x >  61.94 + 10 && trk.end.x <  358.49 - 10)) &&
                               !isnan(trk.end.y) &&
                               ( trk.end.y > -181.86 + 10 && trk.end.y < 134.96 - 10 ) &&
                               !isnan(trk.end.z) &&
                               ( trk.end.z > -894.95 + 10 && trk.end.z < 894.95 - 10 ) );

      const float Chi2Proton = trk.chi2pid[trk.bestplane].chi2_proton;
      const float Chi2Muon = trk.chi2pid[trk.bestplane].chi2_muon;

      float angle = -5.0;
      if ( primaryInd >= 0 ) {
        const unsigned int idxPrim = (unsigned int)primaryInd;
        TVector3 muDir( slc->reco.pfp.at(idxPrim).trk.dir.x, slc->reco.pfp.at(idxPrim).trk.dir.y, slc->reco.pfp.at(idxPrim).trk.dir.z );
        TVector3 pDir( trk.dir.x, trk.dir.y, trk.dir.z );
        angle = TMath::Cos(muDir.Angle(pDir));
      }

      if ( Atslc < 10.0 && slc->reco.pfp.at(idxTrk).parent_is_primary && Contained && Chi2Proton <= 100 && Chi2Muon >= 30 && angle >= -0.9 ) {
        countP += 1;
	      proton_pfpids.push_back( slc->reco.pfp.at(idxTrk).id );
      }
    }

    for ( unsigned int idxStub = 0; idxStub < slc->reco.stub.size(); ++idxStub ) {
      // Check if stub within 10cm of vertex, contained, and does NOT match one of our already
      // selected proton PFParticles. We'll assume it's a proton if so?

      auto const& stub = slc->reco.stub.at(idxStub);
      const float Atslc = std::hypot(slc->vertex.x - stub.vtx.x,
                                     slc->vertex.y - stub.vtx.y,
                                     slc->vertex.z - stub.vtx.z);
      const bool Contained = ( !isnan(stub.end.x) &&
                               ((stub.end.x < -61.94 - 10 && stub.end.x > -358.49 + 10) ||
                                (stub.end.x >  61.94 + 10 && stub.end.x <  358.49 - 10)) &&
                               !isnan(stub.end.y) &&
                               ( stub.end.y > -181.86 + 10 && stub.end.y < 134.96 - 10 ) &&
                               !isnan(stub.end.z) &&
                               ( stub.end.z > -894.95 + 10 && stub.end.z < 894.95 - 10 ) );

      bool matchesProtonPFP = false;
      if ( stub.pfpid >= 0 ) {
        for ( auto const& pfpid : proton_pfpids ) {
          if ( pfpid == stub.pfpid ) matchesProtonPFP = true;
        }
      }

      if ( Atslc < 10.0 && Contained && !matchesProtonPFP ) countP+=1;
    }

    return countP;
  });

// Modify cut looking for ONLY a secondary pion to instead be at LEAST a secondary pion
const Cut kHasScdyPionTrack([](const caf::SRSliceProxy* slc) {
    int primaryInd = kRecoMuonIdx(slc);
    if ( primaryInd < 0 ) return false;
    unsigned int idxPrim = (unsigned int)primaryInd;

    unsigned int idxScdy = 0;
    float maxLength = -1;

    for ( unsigned int idxTrk = 0; idxTrk < slc->reco.pfp.size(); ++idxTrk ) {
      if ( idxTrk == idxPrim ) continue; // need a different track...

      // Here I take the vertex determination and containment from kPTrackInd
      // and check this for meeting the qualificiations
      // NOTE: updated containment to match FV used here... Do I also want to change the padding?
      if( slc->reco.pfp.at(idxTrk).trackScore < 0.5 ) continue;
      auto const& trk = slc->reco.pfp.at(idxTrk).trk;
      if( std::isnan(trk.start.x) || trk.bestplane == -1 ) continue;

      const float Atslc = std::hypot(slc->vertex.x - trk.start.x,
                                     slc->vertex.y - trk.start.y,
                                     slc->vertex.z - trk.start.z);
      const bool Contained = ( !isnan(trk.end.x) &&
                               ((trk.end.x < -61.94 - 10 && trk.end.x > -358.49 + 10) ||
                                (trk.end.x >  61.94 + 10 && trk.end.x <  358.49 - 10)) &&
                               !isnan(trk.end.y) &&
                               ( trk.end.y > -181.86 + 10 && trk.end.y < 134.96 - 10 ) &&
                               !isnan(trk.end.z) &&
                               ( trk.end.z > -894.95 + 10 && trk.end.z < 894.95 - 10 ) );

      // Add Chi2:
      const float Chi2Proton = trk.chi2pid[trk.bestplane].chi2_proton;
      const float Chi2Pion = trk.chi2pid[trk.bestplane].chi2_pion;
      if ( Chi2Pion > 100. || Chi2Proton < 30. ) continue;
      ////////////

      if ( Atslc < 10.0 && slc->reco.pfp.at(idxTrk).parent_is_primary && Contained && trk.len > maxLength ) {
        maxLength = trk.len;
        idxScdy = idxTrk;
      }
    }

    if ( maxLength < 0. ) return false;

    TVector3 muDir( slc->reco.pfp.at(idxPrim).trk.dir.x, slc->reco.pfp.at(idxPrim).trk.dir.y, slc->reco.pfp.at(idxPrim).trk.dir.z );
    TVector3 piDir( slc->reco.pfp.at(idxScdy).trk.dir.x, slc->reco.pfp.at(idxScdy).trk.dir.y, slc->reco.pfp.at(idxScdy).trk.dir.z );
    if ( TMath::Cos(muDir.Angle(piDir)) < -0.9 ) return false;
    return true;
  }); // kHasScdyPionTrack

// ***************** *************** *******************
// ********** *************** *************** **********

/*
// Selects EXACTLY 3 tracks/showers
const Cut kThreePrimaryTracks([](const caf::SRSliceProxy* slc) {
    unsigned int nPrim=0;
    unsigned int nPrimTrks=0;

    for ( unsigned int idxTrk = 0; idxTrk < slc->reco.pfp.size(); ++idxTrk ) {
      if( slc->reco.pfp.at(idxTrk).parent_is_primary ) nPrim+=1;
//      if( slc->reco.pfp.at(idxTrk).trackScore < 0.5 ) continue; //TODO
      auto const& trk = slc->reco.pfp.at(idxTrk).trk;
      if( std::isnan(trk.start.x) || trk.bestplane == -1 ) continue;

      const float Atslc = std::hypot(slc->vertex.x - trk.start.x,
                                     slc->vertex.y - trk.start.y,
                                     slc->vertex.z - trk.start.z);
      if ( Atslc < 10.0 && slc->reco.pfp.at(idxTrk).parent_is_primary ) nPrimTrks+=1;
    }

//    return ( nPrimTrks==3 );
    return ( nPrim==4 && nPrimTrks==3 );
  });
*/

//Selects 3 + 1 below threshold
const Cut kThreePrimaryTracks([](const caf::SRSliceProxy* slc) {
  unsigned int nPrim=0;
  unsigned int nPrimTrks=0;
  double shwThreshold = .015;
  
  for ( unsigned int idxTrk = 0; idxTrk < slc->reco.pfp.size(); ++idxTrk ) {
    if ( !slc->reco.pfp.at(idxTrk).parent_is_primary ) continue;
    nPrim++;
    auto const& trk = slc->reco.pfp.at(idxTrk).trk;
    if ( std::isnan(trk.start.x) || trk.bestplane == -1 ) continue;
    const float Atslc = std::hypot(slc->vertex.x - trk.start.x,
                                     slc->vertex.y - trk.start.y,
                                     slc->vertex.z - trk.start.z);
    double shwEnergy = 0.;
    if ( !isnan(slc->reco.pfp.at(idxTrk).shw.bestplane_energy) ) shwEnergy = slc->reco.pfp.at(idxTrk).shw.bestplane_energy;
    if ( Atslc < 10.0 && shwEnergy > shwThreshold ) nPrimTrks+=1;
  }

  return ( nPrimTrks == 3 &&  nPrim < 6);
  });

const SpillCut kTriggeredEvent([](const caf::SRSpillProxy* sp) {
  return ( sp->hdr.triggerinfo.trigger_within_gate >= 0 );
  });


/*
const Cut k1mu2p0pi([](const caf::SRSliceProxy* slc) {
  if(slc->truth.index < 0) return false;

  int nVisProtons = 0;
  int nContained = 0;
  int nPions = 0;

  bool nuMuCCinFV = false;
  bool signalMuon = false;
  bool twoContainedProtons = false;
  bool pionConditon = false;

  if ( abs(slc->truth.pdg) == 14 &&
       slc->truth.iscc &&
       isInFV_Vtx(slc) ) {
         nuMuCCinFV = true;
         for ( auto const& prim : slc->truth.prim ) {
           if ( abs(prim.pdg) == 13 && sqrt( prim.startp.x*prim.startp.x + prim.startp.y*prim.startp.y + prim.startp.z*prim.startp.z ) >= 0.1 ) signalMuon = true;
           if ( prim.pdg == 2212 && sqrt( prim.startp.x*prim.startp.x + prim.startp.y*prim.startp.y + prim.startp.z*prim.startp.z ) >= 0.4 ) {
             nVisProtons++;
             if ( prim.contained ) nContained++;
           }
           if ( abs(prim.pdg) == 211  && sqrt( prim.startp.x*prim.startp.x + prim.startp.y*prim.startp.y + prim.startp.z*prim.startp.z ) >= 0.065 ) nPions++;
           if ( abs(prim.pdg) == 111 ) nPions++;
         } 
  }

  twoContainedProtons = ( nVisProtons == 2 && nContained == 2 );
  pionConditon = ( nPions == 0 );
  return ( nuMuCCinFV && signalMuon && twoContainedProtons && pionConditon );
  });
*/

const Cut k1mu2p0pi([](const caf::SRSliceProxy* slc) {
  if(slc->truth.index < 0) return false;

  int nVisProtons = 0;
  int nContained = 0;
  int nPions = 0;
  int protonsOver = 0;

  bool nuMuCCinFV = false;
  bool signalMuon = false;
  bool twoContainedProtons = false;
  bool pionConditon = false;

  if ( abs(slc->truth.pdg) == 14 &&
       slc->truth.iscc &&
       isInFV_Vtx(slc) ) {
         nuMuCCinFV = true;
         for ( auto const& prim : slc->truth.prim ) {
           if ( prim.start_process != 0 ) continue;
           if ( abs(prim.pdg) == 13 && sqrt( prim.startp.x*prim.startp.x + prim.startp.y*prim.startp.y + prim.startp.z*prim.startp.z ) >= 0.226 ) signalMuon = true;
           if ( prim.pdg == 2212 && sqrt( prim.startp.x*prim.startp.x + prim.startp.y*prim.startp.y + prim.startp.z*prim.startp.z ) >= 0.35  && sqrt( prim.startp.x*prim.startp.x + prim.startp.y*prim.startp.y + prim.startp.z*prim.startp.z ) < 2. ) {
             nVisProtons++;
             if ( prim.contained ) nContained++;
           }
           else if ( prim.pdg == 2212 && sqrt( prim.startp.x*prim.startp.x + prim.startp.y*prim.startp.y + prim.startp.z*prim.startp.z ) > 2. ) protonsOver++;
           if ( abs(prim.pdg) == 211 ) nPions++;
           if ( abs(prim.pdg) == 111 ) nPions++;
           if ( abs(prim.pdg) == 321 ) nPions++;
           if ( abs(prim.pdg) == 311 ) nPions++;
         } 
  }
  
  twoContainedProtons = ( nVisProtons == 2 /*&& nContained == 2*/ );
  pionConditon = ( nPions == 0 );
  return ( nuMuCCinFV && signalMuon && twoContainedProtons && pionConditon && protonsOver == 0);
  });

const SpillCut kSpill1mu2p0pi([](const caf::SRSpillProxy* sr) {
  bool hasSignal = false;
  for ( const auto &slc : sr->slc ) {
    if ( k1mu2p0pi(&slc) ) {
      hasSignal = true;
      break;
    }
  }

  return hasSignal;
  });

const Cut k1mu2p0piThreshold([](const caf::SRSliceProxy* slc) {
  if(slc->truth.index < 0) return false;

  int nVisProtons = 0;
  int nContained = 0;
  int nPions = 0;
  int nSoftPions = 0;

  bool nuMuCCinFV = false;
  bool signalMuon = false;
  bool twoContainedProtons = false;
  bool pionConditon = false;

  if ( abs(slc->truth.pdg) == 14 &&
       slc->truth.iscc &&
       isInFV_Vtx(slc) ) {
         nuMuCCinFV = true;
         for ( auto const& prim : slc->truth.prim ) {
           if ( prim.start_process != 0 ) continue;
           if ( abs(prim.pdg) == 13 && sqrt( prim.startp.x*prim.startp.x + prim.startp.y*prim.startp.y + prim.startp.z*prim.startp.z ) >= 0.226 ) signalMuon = true;
           if ( prim.pdg == 2212 && sqrt( prim.startp.x*prim.startp.x + prim.startp.y*prim.startp.y + prim.startp.z*prim.startp.z ) >= 0.35 ) {
             nVisProtons++;
             if ( prim.contained ) nContained++;
           }
           if ( abs(prim.pdg) == 211  ) {
             if ( sqrt( prim.startp.x*prim.startp.x + prim.startp.y*prim.startp.y + prim.startp.z*prim.startp.z ) >= 0.065 ) nPions++;
             else nSoftPions++;
           }
           if ( abs(prim.pdg) == 111 ) nPions++;
         } 
  }
  
  twoContainedProtons = ( nVisProtons == 2 /*&& nContained == 2*/ );
  pionConditon = ( nPions == 0 && nSoftPions != 0);
  return ( nuMuCCinFV && signalMuon && twoContainedProtons && pionConditon );
  });

const Cut k1mu2pNpi([](const caf::SRSliceProxy* slc) {
  if(slc->truth.index < 0) return false;

  int nVisProtons = 0;
  int nContained = 0;
  int nPions = 0;
  int protonsOver = 0;

  bool nuMuCCinFV = false;
  bool signalMuon = false;
  bool twoContainedProtons = false;
  bool pionConditon = false;

  if ( abs(slc->truth.pdg) == 14 &&
       slc->truth.iscc &&
       isInFV_Vtx(slc) ) {
         nuMuCCinFV = true;
         for ( auto const& prim : slc->truth.prim ) {
           if ( prim.start_process != 0 ) continue;
           if ( abs(prim.pdg) == 13 && sqrt( prim.startp.x*prim.startp.x + prim.startp.y*prim.startp.y + prim.startp.z*prim.startp.z ) >= 0.226 ) signalMuon = true;
           if ( prim.pdg == 2212 && sqrt( prim.startp.x*prim.startp.x + prim.startp.y*prim.startp.y + prim.startp.z*prim.startp.z ) >= 0.35 && sqrt( prim.startp.x*prim.startp.x + prim.startp.y*prim.startp.y + prim.startp.z*prim.startp.z ) < 2. ) {
             nVisProtons++;
             if ( prim.contained ) nContained++;
           }
           else if ( prim.pdg == 2212 && sqrt( prim.startp.x*prim.startp.x + prim.startp.y*prim.startp.y + prim.startp.z*prim.startp.z ) > 2. ) protonsOver++;
           if ( abs(prim.pdg) == 211 ) nPions++;
           if ( abs(prim.pdg) == 111 ) nPions++;
           if ( abs(prim.pdg) == 321 ) nPions++;
           if ( abs(prim.pdg) == 311 ) nPions++;
         } 
  }

  twoContainedProtons = ( nVisProtons == 2 /*&& nContained == 2*/ );
  pionConditon = ( nPions > 0 );
  return ( nuMuCCinFV && signalMuon && twoContainedProtons && pionConditon && protonsOver == 0);
  });

const Cut k1mu0p0pi([](const caf::SRSliceProxy* slc) {
  if(slc->truth.index < 0) return false;

  int nVisProtons = 0;
  int nContained = 0;
  int nPions = 0;

  bool nuMuCCinFV = false;
  bool signalMuon = false;
  bool protonCondition = false;
  bool pionConditon = false;

  if ( abs(slc->truth.pdg) == 14 &&
       slc->truth.iscc &&
       isInFV_Vtx(slc) ) {
         nuMuCCinFV = true;
         for ( auto const& prim : slc->truth.prim ) {
           if ( prim.start_process != 0 ) continue;
           if ( abs(prim.pdg) == 13 && sqrt( prim.startp.x*prim.startp.x + prim.startp.y*prim.startp.y + prim.startp.z*prim.startp.z ) >= 0.226 ) signalMuon = true;
           if ( prim.pdg == 2212 && sqrt( prim.startp.x*prim.startp.x + prim.startp.y*prim.startp.y + prim.startp.z*prim.startp.z ) >= 0.35 ) {
             nVisProtons++;
             if ( prim.contained ) nContained++;
           }
           if ( abs(prim.pdg) == 211 ) nPions++;
           if ( abs(prim.pdg) == 111 ) nPions++;
           if ( abs(prim.pdg) == 321 ) nPions++;
           if ( abs(prim.pdg) == 311 ) nPions++;
         }
  }

  protonCondition = ( nVisProtons == 0 );
  pionConditon = ( nPions == 0 );
  return ( nuMuCCinFV && signalMuon && protonCondition && pionConditon );
  });

const Cut k1mu0pNpi([](const caf::SRSliceProxy* slc) {
  if(slc->truth.index < 0) return false;
   
  int nVisProtons = 0;
  int nContained = 0;
  int nPions = 0;
   
  bool nuMuCCinFV = false;
  bool signalMuon = false;
  bool protonCondition = false;
  bool pionConditon = false;
    
  if ( abs(slc->truth.pdg) == 14 &&
       slc->truth.iscc &&
       isInFV_Vtx(slc) ) {
         nuMuCCinFV = true;
         for ( auto const& prim : slc->truth.prim ) {
           if ( prim.start_process != 0 ) continue;
           if ( abs(prim.pdg) == 13 && sqrt( prim.startp.x*prim.startp.x + prim.startp.y*prim.startp.y + prim.startp.z*prim.startp.z ) >= 0.226 ) signalMuon = true;
           if ( prim.pdg == 2212 && sqrt( prim.startp.x*prim.startp.x + prim.startp.y*prim.startp.y + prim.startp.z*prim.startp.z ) >= 0.35 ) {
             nVisProtons++;
             if ( prim.contained ) nContained++;
           }
           if ( abs(prim.pdg) == 211 ) nPions++;
           if ( abs(prim.pdg) == 111 ) nPions++;
           if ( abs(prim.pdg) == 321 ) nPions++;
           if ( abs(prim.pdg) == 311 ) nPions++;
         }
  }

  protonCondition = ( nVisProtons == 0 );
  pionConditon = ( nPions > 0 );
  return ( nuMuCCinFV && signalMuon && protonCondition && pionConditon );
  }); 

const Cut k1mu1p0pi([](const caf::SRSliceProxy* slc) {
  if(slc->truth.index < 0) return false;

  int nVisProtons = 0;
  int nContained = 0;
  int nPions = 0;

  bool nuMuCCinFV = false;
  bool signalMuon = false;
  bool protonCondition = false;
  bool pionConditon = false;

  if ( abs(slc->truth.pdg) == 14 &&
       slc->truth.iscc &&
       isInFV_Vtx(slc) ) {
         nuMuCCinFV = true;
         for ( auto const& prim : slc->truth.prim ) {
           if ( prim.start_process != 0 ) continue;
           if ( abs(prim.pdg) == 13 && sqrt( prim.startp.x*prim.startp.x + prim.startp.y*prim.startp.y + prim.startp.z*prim.startp.z ) >= 0.226 ) signalMuon = true;
           if ( prim.pdg == 2212 && sqrt( prim.startp.x*prim.startp.x + prim.startp.y*prim.startp.y + prim.startp.z*prim.startp.z ) >= 0.35 ) {
             nVisProtons++;
             if ( prim.contained ) nContained++;
           }
           if ( abs(prim.pdg) == 211 ) nPions++;
           if ( abs(prim.pdg) == 111 ) nPions++;
           if ( abs(prim.pdg) == 321 ) nPions++;
           if ( abs(prim.pdg) == 311 ) nPions++;
         }
  }

  protonCondition = ( nVisProtons == 1 );
  pionConditon = ( nPions == 0 );
  return ( nuMuCCinFV && signalMuon && protonCondition && pionConditon );
  });

const Cut k1mu1pNpi([](const caf::SRSliceProxy* slc) {
  if(slc->truth.index < 0) return false;
   
  int nVisProtons = 0;
  int nContained = 0;
  int nPions = 0;
   
  bool nuMuCCinFV = false;
  bool signalMuon = false;
  bool protonCondition = false;
  bool pionConditon = false;
    
  if ( abs(slc->truth.pdg) == 14 &&
       slc->truth.iscc &&
       isInFV_Vtx(slc) ) {
         nuMuCCinFV = true;
         for ( auto const& prim : slc->truth.prim ) {
           if ( prim.start_process != 0 ) continue;
           if ( abs(prim.pdg) == 13 && sqrt( prim.startp.x*prim.startp.x + prim.startp.y*prim.startp.y + prim.startp.z*prim.startp.z ) >= 0.226 ) signalMuon = true;
           if ( prim.pdg == 2212 && sqrt( prim.startp.x*prim.startp.x + prim.startp.y*prim.startp.y + prim.startp.z*prim.startp.z ) >= 0.35 ) {
             nVisProtons++;
             if ( prim.contained ) nContained++;
           }
           if ( abs(prim.pdg) == 211 ) nPions++;
           if ( abs(prim.pdg) == 111 ) nPions++;
           if ( abs(prim.pdg) == 321 ) nPions++;
           if ( abs(prim.pdg) == 311 ) nPions++;
         }
  }

  protonCondition = ( nVisProtons == 1 );
  pionConditon = ( nPions > 0 );
  return ( nuMuCCinFV && signalMuon && protonCondition && pionConditon );
  }); 

const Cut k1mu3p0pi([](const caf::SRSliceProxy* slc) {
  if(slc->truth.index < 0) return false;

  int nVisProtons = 0;
  int nContained = 0;
  int nPions = 0;

  bool nuMuCCinFV = false;
  bool signalMuon = false;
  bool protonCondition = false;
  bool pionConditon = false;

  if ( abs(slc->truth.pdg) == 14 &&
       slc->truth.iscc &&
       isInFV_Vtx(slc) ) {
         nuMuCCinFV = true;
         for ( auto const& prim : slc->truth.prim ) {
           if ( prim.start_process != 0 ) continue;
           if ( abs(prim.pdg) == 13 && sqrt( prim.startp.x*prim.startp.x + prim.startp.y*prim.startp.y + prim.startp.z*prim.startp.z ) >= 0.226 ) signalMuon = true;
           if ( prim.pdg == 2212 && sqrt( prim.startp.x*prim.startp.x + prim.startp.y*prim.startp.y + prim.startp.z*prim.startp.z ) >= 0.35 ) {
             nVisProtons++;
             if ( prim.contained ) nContained++;
           }
           if ( abs(prim.pdg) == 211 ) nPions++;
           if ( abs(prim.pdg) == 111 ) nPions++;
           if ( abs(prim.pdg) == 321 ) nPions++;
           if ( abs(prim.pdg) == 311 ) nPions++;
         }
  }

  protonCondition = ( nVisProtons > 2 );
  pionConditon = ( nPions == 0 );
  return ( nuMuCCinFV && signalMuon && protonCondition && pionConditon );
  });

const Cut k1mu3pNpi([](const caf::SRSliceProxy* slc) {
  if(slc->truth.index < 0) return false;
   
  int nVisProtons = 0;
  int nContained = 0;
  int nPions = 0;
   
  bool nuMuCCinFV = false;
  bool signalMuon = false;
  bool protonCondition = false;
  bool pionConditon = false;
    
  if ( abs(slc->truth.pdg) == 14 &&
       slc->truth.iscc &&
       isInFV_Vtx(slc) ) {
         nuMuCCinFV = true;
         for ( auto const& prim : slc->truth.prim ) {
           if ( prim.start_process != 0 ) continue;
           if ( abs(prim.pdg) == 13 && sqrt( prim.startp.x*prim.startp.x + prim.startp.y*prim.startp.y + prim.startp.z*prim.startp.z ) >= 0.226 ) signalMuon = true;
           if ( prim.pdg == 2212 && sqrt( prim.startp.x*prim.startp.x + prim.startp.y*prim.startp.y + prim.startp.z*prim.startp.z ) >= 0.35 ) {
             nVisProtons++;
             if ( prim.contained ) nContained++;
           }
           if ( abs(prim.pdg) == 211 ) nPions++;
           if ( abs(prim.pdg) == 111 ) nPions++;
           if ( abs(prim.pdg) == 321 ) nPions++;
           if ( abs(prim.pdg) == 311 ) nPions++;
         }
  }

  protonCondition = ( nVisProtons > 2 );
  pionConditon = ( nPions > 0 );
  return ( nuMuCCinFV && signalMuon && protonCondition && pionConditon );
  }); 

const Cut k1mu2p = ( k1mu2p0pi || k1mu2pNpi );

const Cut k1mu0p = ( k1mu0p0pi || k1mu0pNpi );

const Cut k1mu1p = ( k1mu1p0pi || k1mu1pNpi );

const Cut k1mu3p = ( k1mu3p0pi || k1mu3pNpi );

const MultiVar kSelectedCPiP([](const caf::SRSliceProxy* slc) {
  std::vector<double> pionPs;
  
  for ( auto const& prim : slc->truth.prim ) {
    if ( prim.start_process != 0 ) continue;
    if ( abs(prim.pdg) == 211 ) {
      pionPs.push_back( sqrt( prim.startp.x*prim.startp.x + prim.startp.y*prim.startp.y + prim.startp.z*prim.startp.z ) );
//      if ( sqrt( prim.startp.x*prim.startp.x + prim.startp.y*prim.startp.y + prim.startp.z*prim.startp.z ) <= 0.065 ) std::cout << std::endl << "~~~~~~~~~~~~~~~~~ Charged pion momentum: " <<  sqrt( prim.startp.x*prim.startp.x + prim.startp.y*prim.startp.y + prim.startp.z*prim.startp.z ) << " ~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
    }
  }

  return pionPs;
  });

const Cut kIsSignal = (k1mu2p0pi || k1mu3p0pi);

const MultiVar kSelectedPi0P([](const caf::SRSliceProxy* slc) {
  std::vector<double> pionPs;

  for ( auto const& prim : slc->truth.prim ) {
    if ( prim.start_process != 0 ) continue;
    if ( abs(prim.pdg) == 111 ) pionPs.push_back( sqrt( prim.startp.x*prim.startp.x + prim.startp.y*prim.startp.y + prim.startp.z*prim.startp.z ) );
  }

  return pionPs;
  });

const Var NonSignalPrinting([](const caf::SRSliceProxy* slc) -> float {
    int nProtons = 0;
    for ( const auto &prim : slc->truth.prim ) if ( prim.pdg == 2212 ) nProtons++;
    return nProtons;
  });

const Var kNPrimary([](const caf::SRSliceProxy* slc) -> float {
    return slc->truth.prim.size();
  });

const Var kNPrimaryProtons([](const caf::SRSliceProxy* slc) -> float {
    int nProtons = 0;
    for ( const auto &prim : slc->truth.prim ) if ( prim.pdg == 2212 ) nProtons++;
    return nProtons;
  });

const Var kNPrimaryPions([](const caf::SRSliceProxy* slc) -> float {
    int nPions = 0;
    for ( const auto &prim : slc->truth.prim ) if ( prim.pdg == 211 ) nPions++;
    return nPions;
  });

const Var kNPrimaryMuons([](const caf::SRSliceProxy* slc) -> float {
    int nMuons = 0;
    for ( const auto &prim : slc->truth.prim ) if ( prim.pdg == 13 ) nMuons++;
    return nMuons;
  });

const Var kNPrimaryNeutrons([](const caf::SRSliceProxy* slc) -> float {
    int nNeutrons = 0;
    for ( const auto &prim : slc->truth.prim ) if ( prim.pdg == 2112 ) nNeutrons++;
    return nNeutrons;
  });

const Var kNPrimaryPi0s([](const caf::SRSliceProxy* slc) -> float {
    int nPi0s = 0;
    for ( const auto &prim : slc->truth.prim ) if ( prim.pdg == 111 ) nPi0s++;
    return nPi0s;
  });

const Var kNPrimaryPiMins([](const caf::SRSliceProxy* slc) -> float {
    int nPiMins = 0;
    for ( const auto &prim : slc->truth.prim ) if ( prim.pdg == -211 ) nPiMins++;
    return nPiMins;
  });

const Var kNPrimaryOther([](const caf::SRSliceProxy* slc) -> float {
    int nOther = 0;
    for ( const auto &prim : slc->truth.prim ) { if ( prim.pdg != 13 && prim.pdg != 211 && prim.pdg != 2212 && prim.pdg != -211 &&prim.pdg != 111 && prim.pdg != 2112 ) nOther++; }
//      std::cout << "Other PDG: " << prim.pdg << std::endl;
    return nOther;
  });

const Var kNPrimaryPFPs([](const caf::SRSliceProxy* slc) -> float {
    int nPrimaries = 0;
    for ( const auto &pfp : slc->reco.pfp ) if ( pfp.parent_is_primary ) nPrimaries++;
    return nPrimaries;
  });

const Cut kFivePrimaryPFPs([](const caf::SRSliceProxy* slc) -> float {
    int nPrimaries = 0;
    for ( const auto &pfp : slc->reco.pfp ) if ( pfp.parent_is_primary ) nPrimaries++;
    return (nPrimaries >= 5);
  });

const MultiVar kPrimaryShwEnergies([](const caf::SRSliceProxy* slc) {
    std::vector<double> energies;
    for ( const auto &pfp : slc->reco.pfp ) {
      if ( pfp.parent_is_primary && pfp.trk.truth.p.pdg != -2147483648 && !isnan(pfp.shw.bestplane_energy) && pfp.shw.bestplane_energy >=0. ) {
        energies.push_back(pfp.shw.bestplane_energy);
      }
    }
    return energies;
  });

const Var kMinShwEnergy([](const caf::SRSliceProxy* slc) -> float {
    double min = 1.e9;
    for ( const auto &pfp : slc->reco.pfp ) if ( pfp.parent_is_primary && pfp.trk.truth.p.pdg != -2147483648 && !isnan(pfp.shw.bestplane_energy) && pfp.shw.bestplane_energy > 0. && pfp.shw.bestplane_energy < min ) min = pfp.shw.bestplane_energy;
    return min;
  });

const MultiVar kPrimaryTrkScores([](const caf::SRSliceProxy* slc) {
    std::vector<double> scores;
    for ( const auto &pfp : slc->reco.pfp ) {
      if ( pfp.parent_is_primary && pfp.trk.truth.p.pdg != -2147483648 && !isnan(pfp.trackScore) ) {
        scores.push_back(pfp.trackScore);
      }
    }
    return scores;
  });

const MultiVar kPrimaryCollectionSums([](const caf::SRSliceProxy* slc) {
    std::vector<double> sums;
    for ( const auto &pfp : slc->reco.pfp ) {
      if ( pfp.parent_is_primary && pfp.trk.truth.p.pdg != -2147483648 && !isnan(pfp.trk.calo[2].charge) ) {
        sums.push_back(pfp.trk.calo[2].charge);
      }
    }
    return sums;
  });

const Var kMinCollectionSum([](const caf::SRSliceProxy* slc) -> float {
    double min = 1.e9;
    for ( const auto &pfp : slc->reco.pfp ) if ( pfp.parent_is_primary && pfp.trk.truth.p.pdg != -2147483648 && !isnan(pfp.trk.calo[2].charge) && pfp.trk.calo[2].charge < min ) min = pfp.trk.calo[2].charge;
    return min;
  });

const MultiVar kPrimaryPDGs([](const caf::SRSliceProxy* slc) {
    std::vector<double> pdgs;
    int pdg;
    for ( const auto &pfp : slc->reco.pfp ) {
      if ( pfp.parent_is_primary && pfp.trk.truth.p.pdg != -2147483648 ) {
        if ( abs(pfp.trk.truth.p.pdg) == 11 )        pdg = 0.;
        else if ( abs(pfp.trk.truth.p.pdg) == 13 )   pdg = 1.;
        else if ( abs(pfp.trk.truth.p.pdg) == 22 )   pdg = 2.;
        else if ( abs(pfp.trk.truth.p.pdg) == 111 )  pdg = 3.;
        else if ( pfp.trk.truth.p.pdg == 211 )       pdg = 4.;
        else if ( pfp.trk.truth.p.pdg == -211 )      pdg = 5.;
        else if ( abs(pfp.trk.truth.p.pdg) == 2212 ) pdg = 6.;
        else                                     pdg = 7.;
        pdgs.push_back(pdg);
      }
    }
    return pdgs;
  });

const Cut kLeadingProtonThreshold([](const caf::SRSliceProxy* slc) {
  bool thresh = false;
  int protonInd = kRecoProtonIdx(slc);
  if ( protonInd >= 0 ) {
    const auto &trk = slc->reco.pfp.at(protonInd).trk;
    thresh = ( trk.rangeP.p_proton > 0.35 && trk.rangeP.p_proton < 2. );
  }

  return thresh;
  });

//Cuts for cut flow efficiencies

const Cut kPreselection = (kFV && kNotClearCosmic);

const Cut k1mu1p1X = ( kFV && kHasMuon && kLeadingProtonThreshold /*&& kThreePrimaryTracks*/ && kNotClearCosmic);

//For three track slices with one muon and one proton, the third track is the pion candidate
//Requiring hadronic containment
/*
const Var kPionCandidate([](const caf::SRSliceProxy* slc) -> float {
  int idxPion = -1;

  if ( !k1mu1p1X(slc) ) return idxPion;
  int idxMuon = kRecoMuonIdx(slc);
  int idxProton = kRecoProtonIdx(slc);

  int nPFPs = slc->reco.pfp.size();
  for ( int idxTrk = 0; idxTrk < nPFPs; ++idxTrk ) {
    if ( slc->reco.pfp.at(idxTrk).trackScore < 0.5 ) continue;
    const auto &trk = slc->reco.pfp.at(idxTrk).trk;
    const float Atslc = std::hypot(slc->vertex.x - trk.start.x,
                                   slc->vertex.y - trk.start.y,
                                   slc->vertex.z - trk.start.z);
    const bool AtSlice = ( Atslc < 10.0 && slc->reco.pfp.at(idxTrk).parent_is_primary);
    const bool Contained = (!isnan(trk.end.x) &&
                            ((trk.end.x < -61.94 - 10 && trk.end.x > -358.49 + 10) ||
                                (trk.end.x >  61.94 + 10 && trk.end.x <  358.49 - 10)) &&
                            !isnan(trk.end.y) &&
                            ( trk.end.y > -181.86 + 10 && trk.end.y < 134.96 - 10 ) &&
                            !isnan(trk.end.z) &&
                            ( trk.end.z > -894.95 + 10 && trk.end.z < 894.95 - 10 ) );

    
    const bool Chi2 = ( trk.chi2pid[2].chi2_proton > 60. && trk.chi2pid[2].chi2_muon < 30. );
    if ( idxTrk != idxMuon && idxTrk != idxProton && AtSlice && Contained && Chi2 ) idxPion = idxTrk;
  }

  return idxPion;
  });

const Var kScndProtonIdx([](const caf::SRSliceProxy* slc) -> float {
  int idxPion = -1;
  double shwThreshold = .015;

  if ( !k1mu1p1X(slc) ) return idxPion;
  int idxMuon = kRecoMuonIdx(slc);
  int idxProton = kRecoProtonIdx(slc);

  int nPFPs = slc->reco.pfp.size();
  for ( int idxTrk = 0; idxTrk < nPFPs; ++idxTrk ) {
//    if ( slc->reco.pfp.at(idxTrk).trackScore < 0.5 ) continue; //TODO
    double shwEnergy = 0.;
    if ( !isnan(slc->reco.pfp.at(idxTrk).shw.bestplane_energy) ) shwEnergy = slc->reco.pfp.at(idxTrk).shw.bestplane_energy;
    const auto &trk = slc->reco.pfp.at(idxTrk).trk;
    const float Atslc = std::hypot(slc->vertex.x - trk.start.x,
                                   slc->vertex.y - trk.start.y,
                                   slc->vertex.z - trk.start.z);
    const bool AtSlice = ( Atslc < 10.0 && slc->reco.pfp.at(idxTrk).parent_is_primary);
    if ( idxTrk != idxMuon && idxTrk != idxProton && AtSlice && shwEnergy > shwThreshold ) idxPion = idxTrk;
  }

  return idxPion;
  });
*/

const Var kScndProtonIdx([](const caf::SRSliceProxy* slc) -> float {
  int idxProton2 = -1;
  double longest = 0.;

  if ( !k1mu1p1X(slc) ) return idxProton2;
  int idxMuon = kRecoMuonIdx(slc);
  int idxProton = kRecoProtonIdx(slc);

  int nPFPs = slc->reco.pfp.size();
  for ( int idxTrk = 0; idxTrk < nPFPs; ++idxTrk ) {
    if ( idxTrk == idxMuon || idxTrk == idxProton ) continue;
    const auto &trk = slc->reco.pfp.at(idxTrk).trk;
    if( isnan(trk.start.x) || isnan(trk.chi2pid[2].chi2_proton) ) continue;
    const float Atslc = std::hypot(slc->vertex.x - trk.start.x,
                                   slc->vertex.y - trk.start.y,
                                   slc->vertex.z - trk.start.z);

    const bool Contained = (!isnan(trk.end.x) &&
                            ((trk.end.x < -61.94 - 10 && trk.end.x > -358.49 + 10) ||
                                (trk.end.x >  61.94 + 10 && trk.end.x <  358.49 - 10)) &&
                            !isnan(trk.end.y) &&
                            ( trk.end.y > -181.86 + 10 && trk.end.y < 134.96 - 10 ) &&
                            !isnan(trk.end.z) &&
                            ( trk.end.z > -894.95 + 10 && trk.end.z < 894.95 - 10 ) );

    const bool AtSlice = ( Atslc < 10.0 && slc->reco.pfp.at(idxTrk).parent_is_primary);
    const bool chi2PID = ( trk.chi2pid[2].chi2_proton < 50. && trk.chi2pid[2].chi2_muon != 0. );
    if ( AtSlice && chi2PID && trk.len > longest  && Contained ) {
      longest = trk.len;
      idxProton2 = idxTrk;
    }
  }

  return idxProton2;
  });

const Var kThirdProton([](const caf::SRSliceProxy* slc) -> int {
  int idx = -1;
  if ( !k1mu1p1X(slc) ) return idx;
  if ( kScndProtonIdx(slc) == -1 ) return idx;

  unsigned idxMuon = kRecoMuonIdx(slc);
  unsigned idxP1 = kRecoProtonIdx(slc);
  unsigned idxP2 = kScndProtonIdx(slc);

  double longest = -5.;
  for ( unsigned idxTrk = 0; idxTrk < slc->reco.pfp.size(); ++idxTrk ) {
    if ( !slc->reco.pfp.at(idxTrk).parent_is_primary ) continue;
    if ( idxTrk == idxMuon || idxTrk == idxP1 || idxTrk == idxP2 ) continue;
    const auto &pfp = slc->reco.pfp.at(idxTrk);
    if ( pfp.trk.chi2pid[2].chi2_proton <= 50. && pfp.trk.chi2pid[2].chi2_muon != 0. && pfp.trk.rangeP.p_proton >= .35 && pfp.trk.len > longest ) {
      longest = pfp.trk.rangeP.p_proton;
      idx = idxTrk;
    }
  }

  return idx;
  });

const MultiVar kRecoProtonIndices([](const caf::SRSliceProxy* slc) {
  std::vector<double> protons;

  for ( unsigned i = 0; i < slc->reco.pfp.size(); i++ ) {
    const auto &pfp = slc->reco.pfp.at(i);
    if ( !pfp.parent_is_primary ) continue;
    if ( isnan(pfp.trk.start.x) || isnan(pfp.trk.start.y) || isnan(pfp.trk.start.z) || isnan(pfp.trk.end.x) || isnan(pfp.trk.end.y) || isnan(pfp.trk.end.z) ) continue;

    bool contained = (
                       ((pfp.trk.end.x < -61.94 - 10 && pfp.trk.end.x > -358.49 + 10) ||
                         (pfp.trk.end.x >  61.94 + 10 && pfp.trk.end.x <  358.49 - 10)) &&
                       ( pfp.trk.end.y > -181.86 + 10 && pfp.trk.end.y < 134.96 - 10 ) &&
                       ( pfp.trk.end.z > -894.95 + 10 && pfp.trk.end.z < 894.95 - 10 )
                     );

    double atslc = std::hypot(slc->vertex.x - pfp.trk.start.x,
                                slc->vertex.y - pfp.trk.start.y,
                                slc->vertex.z - pfp.trk.start.z);

    if ( pfp.trk.chi2pid[2].chi2_proton < 50. && pfp.trk.chi2pid[2].chi2_muon != 0. && pfp.trk.rangeP.p_proton > .35 && contained && atslc < 10. ) protons.push_back(static_cast<double>(i));
  }

  return protons;
});

//For Minerba, just clean proton candidate vars
/*
const Var kLeadingProton([](const caf::SRSliceProxy* slc) -> int {
  int idxProton = -1;
  double longest = -5.;

  for ( unsigned i = 0; i < slc->reco.pfp.size(); i++ ) {
    const auto &pfp = slc->reco.pfp.at(i);
    if ( !pfp.parent_is_primary ) continue;
    if ( isnan(pfp.trk.start.x) || isnan(pfp.trk.start.y) || isnan(pfp.trk.start.z) || isnan(pfp.trk.end.x) || isnan(pfp.trk.end.y) || isnan(pfp.trk.end.z) ) continue;

    bool contained = (
                       ((pfp.trk.end.x < -61.94 - 10 && pfp.trk.end.x > -358.49 + 10) ||
                         (pfp.trk.end.x >  61.94 + 10 && pfp.trk.end.x <  358.49 - 10)) &&
                       ( pfp.trk.end.y > -181.86 + 10 && pfp.trk.end.y < 134.96 - 10 ) &&
                       ( pfp.trk.end.z > -894.95 + 10 && pfp.trk.end.z < 894.95 - 10 )
                     );

    double atslc = std::hypot(slc->vertex.x - pfp.trk.start.x,
                                slc->vertex.y - pfp.trk.start.y,
                                slc->vertex.z - pfp.trk.start.z);

    bool validProton = (pfp.trk.chi2pid[2].chi2_proton < 50. && pfp.trk.chi2pid[2].chi2_muon != 0. && pfp.trk.rangeP.p_proton > .35 && contained && atslc < 10.);
    if ( validProton && pfp.trk.len > longest ) {
      longest = pfp.trk.len;
      idxProton = i;
    }
  }

  return idxProton;
});

const Var kSecondProton([](const caf::SRSliceProxy* slc) -> int {
  int idxProton2 = -1;
  if ( kLeadingProton(slc) == -1 ) return idxProton2;
  unsigned idxProton1 = kLeadingProton(slc); 
  double longest = -5.;

  for ( unsigned i = 0; i < slc->reco.pfp.size(); i++ ) {
    if ( i == idxProton1 ) continue;
    const auto &pfp = slc->reco.pfp.at(i);
    if ( !pfp.parent_is_primary ) continue;
    if ( isnan(pfp.trk.start.x) || isnan(pfp.trk.start.y) || isnan(pfp.trk.start.z) || isnan(pfp.trk.end.x) || isnan(pfp.trk.end.y) || isnan(pfp.trk.end.z) ) continue;

    bool contained = (
                       ((pfp.trk.end.x < -61.94 - 10 && pfp.trk.end.x > -358.49 + 10) ||
                         (pfp.trk.end.x >  61.94 + 10 && pfp.trk.end.x <  358.49 - 10)) &&
                       ( pfp.trk.end.y > -181.86 + 10 && pfp.trk.end.y < 134.96 - 10 ) &&
                       ( pfp.trk.end.z > -894.95 + 10 && pfp.trk.end.z < 894.95 - 10 )
                     );

    double atslc = std::hypot(slc->vertex.x - pfp.trk.start.x,
                                slc->vertex.y - pfp.trk.start.y,
                                slc->vertex.z - pfp.trk.start.z);

    bool validProton = (pfp.trk.chi2pid[2].chi2_proton < 50. && pfp.trk.chi2pid[2].chi2_muon != 0. && pfp.trk.rangeP.p_proton > .35 && contained && atslc < 10.);
    if ( validProton && pfp.trk.len > longest ) { 
      longest = pfp.trk.len;
      idxProton2 = i;
    }
  }
  
  return idxProton2;
});
*/

//Copy/paste of kRecoMuonIdx used to search for a second MIP
const Var kSidebandPion([](const caf::SRSliceProxy* slc) -> int {
    // The (dis)qualification of a slice is based upon the track level information.
    float Longest(0);
    double lengthReq = 10.;
    int PTrackInd(-1);

    if ( !k1mu1p1X(slc) ) return PTrackInd;
    unsigned idxMuon = kRecoMuonIdx(slc);
//    unsigned idxP1 = kRecoProtonIdx(slc);
    std::vector<double> idcsProton = kRecoProtonIndices(slc);

    for (std::size_t i(0); i < slc->reco.pfp.size(); ++i)
    {   

        if ( i == idxMuon || std::find(idcsProton.begin(), idcsProton.end(), i) != idcsProton.end() ) continue;

        if( slc->reco.pfp.at(i).trackScore < .45 ) continue;
        auto const& trk = slc->reco.pfp.at(i).trk;
        if( std::isnan(trk.start.x) || trk.bestplane == -1 ) continue;
        
        // First we calculate the distance of each track to the slice vertex.
        const float Atslc = std::hypot(slc->vertex.x - trk.start.x,
                                       slc->vertex.y - trk.start.y,
                                       slc->vertex.z - trk.start.z);
        
        // We require that the distance of the track from the slice is less than
        // 10 cm and that the parent of the track has been marked as the primary.
        const bool AtSlice = ( Atslc < 10.0 && slc->reco.pfp.at(i).parent_is_primary);
        
        const float Chi2Proton = trk.chi2pid[2].chi2_proton;
        const float Chi2Muon = trk.chi2pid[2].chi2_muon;
        
        const bool Contained = (!isnan(trk.end.x) && 
                                ((trk.end.x < -61.94 - 10 && trk.end.x > -358.49 + 10) ||
                                 (trk.end.x >  61.94 + 10 && trk.end.x <  358.49 - 10)) &&
                                !isnan(trk.end.y) &&
                                ( trk.end.y > -181.86 + 10 && trk.end.y < 134.96 - 10 ) &&
                                !isnan(trk.end.z) &&
                                ( trk.end.z > -894.95 + 10 && trk.end.z < 894.95 - 10 ) );
        //const bool MaybeMuonExiting = ( !Contained && trk.len > 50);
        const bool MaybeMuonContained = ( Contained && Chi2Proton > 60 && Chi2Muon < 30 && trk.len > lengthReq );
        if ( AtSlice && MaybeMuonContained && trk.len > Longest )
          {     
                Longest = trk.len;
                PTrackInd = i;
          }
    }
    return PTrackInd;
  });

const Cut kThirdRENAMEIFNEEDEDPrimaryContainedTruth([](const caf::SRSliceProxy* slc) {
  bool Contained = false;
  if ( kScndProtonIdx(slc) >= 0 ) {
    const auto &trueParticle = slc->reco.pfp.at(kScndProtonIdx(slc)).trk.truth.p;
    Contained = (!isnan(trueParticle.end.x) &&
                            ((trueParticle.end.x < -61.94 - 10 && trueParticle.end.x > -358.49 + 10) ||
                                (trueParticle.end.x >  61.94 + 10 && trueParticle.end.x <  358.49 - 10)) &&
                            !isnan(trueParticle.end.y) &&
                            ( trueParticle.end.y > -181.86 + 10 && trueParticle.end.y < 134.96 - 10 ) &&
                            !isnan(trueParticle.end.z) &&
                            ( trueParticle.end.z > -894.95 + 10 && trueParticle.end.z < 894.95 - 10 ) );
    }
  return ( Contained );
  });

const Cut kMIPLikeTrack([](const caf::SRSliceProxy* slc) {
  bool Chi2 = false;
  int thirdPrim = kScndProtonIdx(slc);
  if ( thirdPrim >= 0 ) {
    const auto &trk = slc->reco.pfp.at(thirdPrim).trk;
    Chi2 = ( trk.chi2pid[2].chi2_proton > 60. && trk.chi2pid[2].chi2_muon < 30. );
  }

  return Chi2;
  });

const Cut kProtonLikeTrack([](const caf::SRSliceProxy* slc) {
  bool Chi2 = false;
  int thirdPrim = kScndProtonIdx(slc);
  if ( thirdPrim >= 0 ) {
    const auto &trk = slc->reco.pfp.at(thirdPrim).trk;
    Chi2 = ( /*trk.chi2pid[2].chi2_proton <= 100. && trk.chi2pid[2].chi2_muon >= 30.*/ trk.chi2pid[2].chi2_proton < 50. && trk.chi2pid[2].chi2_muon != 0. );
  }

  return Chi2;
  });

const Cut kScndProtonThreshold([](const caf::SRSliceProxy* slc) {
  bool thresh = false;
  int protonInd = kScndProtonIdx(slc);
  if ( protonInd >= 0 ) {
    const auto &trk = slc->reco.pfp.at(protonInd).trk;
    thresh = ( trk.rangeP.p_proton > 0.35 );
  }

  return thresh;
  });

const Cut kTruePion([](const caf::SRSliceProxy* slc) {
  bool PDG = false;
  if ( kScndProtonIdx(slc) >= 0 ) {
    const auto &trueParticle = slc->reco.pfp.at(kScndProtonIdx(slc)).trk.truth.p;
    PDG = ( abs(trueParticle.pdg) == 211 );
    }
  return ( PDG );
  });

const Cut kTrueMuon([](const caf::SRSliceProxy* slc) {
  bool PDG = false;
  if ( kScndProtonIdx(slc) >= 0 ) {
    const auto &trueParticle = slc->reco.pfp.at(kScndProtonIdx(slc)).trk.truth.p;
    PDG = ( abs(trueParticle.pdg) == 13 );
    }
  return ( PDG );
  });

const Cut kTrueProton([](const caf::SRSliceProxy* slc) {
  bool PDG = false;
  if ( kScndProtonIdx(slc) >= 0 ) {
    const auto &trueParticle = slc->reco.pfp.at(kScndProtonIdx(slc)).trk.truth.p;
    PDG = ( abs(trueParticle.pdg) == 2212 );
    }
  return ( PDG );
  });

const Cut kMuPiMixing([](const caf::SRSliceProxy* slc) {
  bool PDG = false;
  if ( kScndProtonIdx(slc) >= 0 ) {
    const auto &pionTruth = slc->reco.pfp.at(kScndProtonIdx(slc)).trk.truth.p;
    const auto &muonTruth = slc->reco.pfp.at(kRecoMuonIdx(slc)).trk.truth.p;
    PDG = ( abs(muonTruth.pdg) == 211 && abs(pionTruth.pdg) == 13 );
    }
  return ( PDG );
  });

const Var kRecoMuonThetaNuMI([](const caf::SRSliceProxy* slc) -> float {
    float costh = -9999.;

    if ( kRecoMuonIdx(slc) >= 0 ) {
      auto const& trk = slc->reco.pfp.at(kRecoMuonIdx(slc)).trk;
      TVector3 direction(trk.dir.x, trk.dir.y, trk.dir.z);
      const auto& vtx = slc->vertex;
      TVector3 vec_vtx(vtx.x, vtx.y, vtx.z);
      TVector3 vec_numi_to_vtx = (dFromNuMI + vec_vtx).Unit();

      costh = direction.Dot(vec_numi_to_vtx) / (direction.Mag() * vec_numi_to_vtx.Mag());
    }

    return costh;
  });

const Var kRecoMuonPNuMI([](const caf::SRSliceProxy* slc) -> float {
    float momentum = -9999.;

    if ( kRecoMuonIdx(slc) >= 0 ) {
      double pMu_mag = kRecoMuonPNew(slc);
      double costh = kRecoMuonThetaNuMI(slc);
      momentum = pMu_mag * costh;
    }

    return momentum;
  });

const Var kRecoMuonTrackLength([](const caf::SRSliceProxy* slc) -> float {
    float length = -9999.;

    if ( kRecoMuonIdx(slc) >= 0 ) {
      auto const& trk = slc->reco.pfp.at(kRecoMuonIdx(slc)).trk;
      length = trk.len;
    }

    return length;
  });

const Var kRecoMuonTrueLength([](const caf::SRSliceProxy* slc) -> float {
    float length = -9999.;

    if ( kRecoMuonIdx(slc) >= 0 ) {
      auto const& trk = slc->reco.pfp.at(kRecoMuonIdx(slc)).trk;
      length = trk.truth.p.length;
    }

    return length;
  });

const Cut kForwardMuon([](const caf::SRSliceProxy* slc) {
  double costh = kRecoMuonThetaNuMI(slc);
  return ( costh > -.5 );
  });

//const MultiVar kRecoProtonIndices -- moved up in file

const Var kNRecoProtons([](const caf::SRSliceProxy* slc) -> int {
  std::vector<double> protons = kRecoProtonIndices(slc);
  return protons.size();
});

const Var kNTrueProtons([](const caf::SRSliceProxy* slc) -> int {
  int nProtons = 0;

  for ( const auto &prim : slc->truth.prim ) {
    if ( prim.start_process != 0 ) continue;
    if ( prim.pdg == 2212 && std::hypot(prim.startp.x, prim.startp.y, prim.startp.z) > .35 ) nProtons++;
  }

  return nProtons;
});

const Var kNProtonsOverCap([](const caf::SRSliceProxy* slc) -> int {
  int nProtons = 0;

  for ( const auto &prim : slc->truth.prim ) {
    if ( prim.start_process != 0 ) continue;
    if ( prim.pdg == 2212 && std::hypot(prim.startp.x, prim.startp.y, prim.startp.z) > 2. ) nProtons++;
  }

  return nProtons;
});

const Var kNSoftProtons([](const caf::SRSliceProxy* slc) -> int {
  int nProtons = 0;

  for ( const auto &prim : slc->truth.prim ) {
    if ( prim.start_process != 0 ) continue;
    if ( prim.pdg == 2212 && std::hypot(prim.startp.x, prim.startp.y, prim.startp.z) < .35 ) nProtons++;
  }

  return nProtons;
});

const Var kNPions([](const caf::SRSliceProxy* slc) -> int {
  int nPions = 0;

  for ( const auto &prim : slc->truth.prim ) {
    if ( prim.start_process != 0 ) continue;
    if ( abs(prim.pdg) == 211 ) nPions++;
  }

  return nPions;
});

const Var kNTruePiPlus([](const caf::SRSliceProxy* slc) -> int {
  int nPions = 0;

  for ( const auto &prim : slc->truth.prim ) {
    if ( prim.start_process != 0 ) continue;
    if ( prim.pdg == 211 ) nPions++;
  }

  return nPions;
});

const Var kTruePiPlusP([](const caf::SRSliceProxy* slc) -> float {
  double max = -9999.;;

  for ( const auto &prim : slc->truth.prim ) {
    if ( prim.start_process != 0 ) continue;
    double thisP = std::hypot(prim.startp.x, prim.startp.y, prim.startp.z);
    if ( prim.pdg == 211 && thisP > max  && thisP < 10000.) max = thisP;
  }

  return max;
});

const Var kNTruePiMinus([](const caf::SRSliceProxy* slc) -> int {
  int nPions = 0;

  for ( const auto &prim : slc->truth.prim ) {
    if ( prim.start_process != 0 ) continue;
    if ( prim.pdg == -211 ) nPions++;
  }

  return nPions;
});

const Var kTruePiMinusP([](const caf::SRSliceProxy* slc) -> float {
  double max = -9999.;;

  for ( const auto &prim : slc->truth.prim ) {
    if ( prim.start_process != 0 ) continue; 
    double thisP = std::hypot(prim.startp.x, prim.startp.y, prim.startp.z);
    if ( prim.pdg == -211 && thisP > max  && thisP < 10000.) max = thisP;
  }

  return max;
});

const Var kNPi0s([](const caf::SRSliceProxy* slc) -> int {
  int nPions = 0;

  for ( const auto &prim : slc->truth.prim ) {
    if ( prim.start_process != 0 ) continue;
    if ( abs(prim.pdg) == 111 ) nPions++;
  }

  return nPions;
});

const Var kNTrueNeutrons([](const caf::SRSliceProxy* slc) -> int {
  int nNeutrons = 0;

  for ( const auto &prim : slc->truth.prim ) {
    if ( prim.start_process != 0 ) continue;
    if ( prim.pdg == 2112 ) nNeutrons++;
  }

  return nNeutrons;
});

const Var kTrueNeutronP([](const caf::SRSliceProxy* slc) -> float {
  double max = -9999.;;
    
  for ( const auto &prim : slc->truth.prim ) {
    if ( prim.start_process != 0 ) continue; 
    double thisP = std::hypot(prim.startp.x, prim.startp.y, prim.startp.z);
    if ( prim.pdg == 2112 && thisP > max  && thisP < 10000.) max = thisP;
  }
  
  return max;
});

const Cut kCountPrimaryPFPs([](const caf::SRSliceProxy* slc) {
  int nAbove = 0;
  double shwThreshold = 20.;

  unsigned int idxMuon = kRecoMuonIdx(slc);
  std::vector<double> protons = kRecoProtonIndices(slc);

  for ( unsigned idxTrk = 0; idxTrk < slc->reco.pfp.size(); ++idxTrk ) {
    if ( !slc->reco.pfp.at(idxTrk).parent_is_primary ) continue;
    if ( idxTrk == idxMuon || std::find(protons.begin(), protons.end(), static_cast<unsigned>(idxTrk)) != protons.end() ) continue;
    double shwLength = slc->reco.pfp.at(idxTrk).shw.len;
    if ( shwLength > shwThreshold ) nAbove+=1;
  }

  return 1;
//  return ( nAbove == 0 );
  });

// This should be moot
//const Cut kHasSidebandProton([](const caf::SRSliceProxy* slc) {
//  int idxP3 = kSidebandProton(slc);
//  return (idxP3 > -1);
//  });

const Cut kHasSidebandPion([](const caf::SRSliceProxy* slc) {
  int idxP3 = kSidebandPion(slc);
  return (idxP3 > -1);
  });

const Var kHasSidebandPionVar([](const caf::SRSliceProxy* slc) -> int {
  int idxP3 = kSidebandPion(slc);
  return (idxP3 > -1);
  });

/*
const Cut kThreeRecoProtonsONLY([](const caf::SRSliceProxy* slc) {
  int nPrim = 0;
  int nAbove = 0;
  double shwThreshold = .015;

  unsigned int idxMuon = kRecoMuonIdx(slc);
  unsigned int idxP1 = kRecoProtonIdx(slc);
  unsigned int idxP2 = kScndProtonIdx(slc);
  unsigned int idxP3 = kSidebandProton(slc);

  for ( unsigned int idxTrk = 0; idxTrk < slc->reco.pfp.size(); ++idxTrk ) {
    if ( !slc->reco.pfp.at(idxTrk).parent_is_primary ) continue;
    nPrim++;
    if ( idxTrk == idxMuon || idxTrk == idxP1 || idxTrk == idxP2 || idxTrk == idxP3 ) continue;
    double shwEnergy = slc->reco.pfp.at(idxTrk).shw.bestplane_energy;
    if ( shwEnergy > shwThreshold ) nAbove+=1;
  }

  return ( nAbove == 0 );
  });

const Cut kContainedSidebandHadron([](const caf::SRSliceProxy* slc) {
  bool Contained = false;
  int idx = -1;
  int idxProton = kSidebandProton(slc);
  int idxPion = kSidebandPion(slc);

  if ( (idxProton!=-1) == (idxPion!=-1) ) return Contained;

  if ( idxProton != -1 ) idx = idxProton;
  else if (idxPion != -1) idx = idxPion;
  const auto &trk = slc->reco.pfp.at(idx).trk;
  Contained = (!isnan(trk.end.x) &&
                            ((trk.end.x < -61.94 - 10 && trk.end.x > -358.49 + 10) ||
                                (trk.end.x >  61.94 + 10 && trk.end.x <  358.49 - 10)) &&
                            !isnan(trk.end.y) &&
                            ( trk.end.y > -181.86 + 10 && trk.end.y < 134.96 - 10 ) &&
                            !isnan(trk.end.z) &&
                            ( trk.end.z > -894.95 + 10 && trk.end.z < 894.95 - 10 ) );
  return ( Contained );
  });

const Cut kThreeRecoProtons([](const caf::SRSliceProxy* slc) {
  int nProtons = 0;
  for ( const auto &pfp : slc->reco.pfp ) {
    if ( !pfp.parent_is_primary ) continue;
    if ( pfp.trk.chi2pid[2].chi2_proton <= 50. && pfp.trk.chi2pid[2].chi2_muon != 0. && pfp.trk.rangeP.p_proton >= .35) nProtons++;
  }
  return (nProtons == 3);
});
*/

const Cut kPionSidebandPFPCut([](const caf::SRSliceProxy* slc) {
  int nPrim = 0;
  int nAbove = 0;
  double shwThreshold = .015;

  unsigned int idxMuon = kRecoMuonIdx(slc);
  unsigned int idxProton = kRecoProtonIdx(slc);
  unsigned int idxPion = kSidebandPion(slc);

  for ( unsigned int idxTrk = 0; idxTrk < slc->reco.pfp.size(); ++idxTrk ) {
    if ( !slc->reco.pfp.at(idxTrk).parent_is_primary ) continue;
    nPrim++;
    if ( idxTrk == idxMuon || idxTrk == idxProton || idxTrk == idxPion ) continue;
    double shwEnergy = slc->reco.pfp.at(idxTrk).shw.bestplane_energy;
    if ( shwEnergy > shwThreshold ) nAbove+=1;
  }

  return ( nAbove == 0 );
  });

/*
const Cut kThreeRecoProtonsONLY([](const caf::SRSliceProxy* slc) {
  unsigned idxMuon = kRecoMuonIdx(slc);
  int nProtons(0), nAbove(0), nPrim(0);
  double shwThreshold = .015;

  for ( unsigned int i = 0; i < slc->reco.pfp.size(); ++i ) {
    const auto &pfp = slc->reco.pfp.at(i);
    if ( !pfp.parent_is_primary ) continue;
    nPrim++;
    if ( i == idxMuon ) continue;
    if ( pfp.trk.chi2pid[2].chi2_proton <= 50. && pfp.trk.chi2pid[2].chi2_muon != 0. && pfp.trk.rangeP.p_proton >= .35) { nProtons++; continue; }
    double shwEnergy = pfp.shw.bestplane_energy;
    if ( shwEnergy > shwThreshold ) nAbove+=1;
  }
  return (nProtons == 3 && nAbove == 0);
});

const Cut kNoExtraMIP([](const caf::SRSliceProxy* slc) {
  int idxMIP = kRecoMuonIdxCOPY(slc);
  return ( idxMIP == -1 );
  });
*/

//>>>This was originally meant for signal region pion rejection, but I'm recycling it for multipions in the sidebands
const Var kNoExtraMIP([](const caf::SRSliceProxy* slc) -> float {
  bool noMIP = true;
  unsigned idxMuon = kRecoMuonIdx(slc);
  std::vector<double> idxProtons = kRecoProtonIndices(slc);
  unsigned idxPion = kSidebandPion(slc);//>>>
  //<<<const auto &muTrk = slc->reco.pfp.at(idxMuon).trk;

  double lengthReq = 10.;

  for ( unsigned i = 0; i < slc->reco.pfp.size(); i++ ) {
    if ( i == idxMuon || std::find(idxProtons.begin(), idxProtons.end(), i) != idxProtons.end() || i == idxPion ) continue;//>>>
    const auto &pfp = slc->reco.pfp.at(i);

//<<<
/*
    double trkStartDiff = 9999.;
    double shwStartDiff = 9999.;
    if ( !isnan(pfp.trk.start.x) && !isnan(pfp.trk.start.y) && !isnan(pfp.trk.start.z) )
      trkStartDiff = std::hypot( (pfp.trk.start.x - muTrk.end.x), (pfp.trk.start.y - muTrk.end.y), (pfp.trk.start.z - muTrk.end.z) );
    if ( !isnan(pfp.shw.start.x) && !isnan(pfp.shw.start.y) && !isnan(pfp.shw.start.z) )
      shwStartDiff = std::hypot( (pfp.shw.start.x - muTrk.end.x), (pfp.shw.start.y - muTrk.end.y), (pfp.shw.start.z - muTrk.end.z) );
    if ( std::min(trkStartDiff, shwStartDiff) < 5. ) continue;
*/

    const float Atslc = std::hypot(slc->vertex.x - pfp.trk.start.x,
                                       slc->vertex.y - pfp.trk.start.y,
                                       slc->vertex.z - pfp.trk.start.z);
    const bool AtSlice = ( Atslc < 10.0 && slc->reco.pfp.at(i).parent_is_primary);//>>>

    double trackLength = pfp.trk.len;
    double trackScore = pfp.trackScore;
    double chi2muon = ( !isnan(pfp.trk.chi2pid[2].chi2_muon) ) ? (double) pfp.trk.chi2pid[2].chi2_muon : -9999.;
    double chi2proton = ( !isnan(pfp.trk.chi2pid[2].chi2_proton) ) ? (double) pfp.trk.chi2pid[2].chi2_proton : -9999.;
 
    bool contained = (!isnan(pfp.trk.end.x) &&
                             ((pfp.trk.end.x < -61.94 - 10 && pfp.trk.end.x > -358.49 + 10) ||
                                 (pfp.trk.end.x >  61.94 + 10 && pfp.trk.end.x <  358.49 - 10)) &&
                             !isnan(pfp.trk.end.y) &&
                             ( pfp.trk.end.y > -181.86 + 10 && pfp.trk.end.y < 134.96 - 10 ) &&
                             !isnan(pfp.trk.end.z) &&
                             ( pfp.trk.end.z > -894.95 + 10 && pfp.trk.end.z < 894.95 - 10 ) );
   
    if ( contained && trackLength > lengthReq && chi2muon < 30 && chi2proton > 60 && trackScore > 0.45 && AtSlice ) {//>>>
      noMIP = false;
      break;
    }
  }

  return ( noMIP );
  });

const Cut kNoExtraMIPCut([](const caf::SRSliceProxy* slc) {
  bool noMIP = kNoExtraMIP(slc);
  return ( noMIP );
  });

const Var kExtraMIPIdx([](const caf::SRSliceProxy* slc) -> float {
  int idx = -5;
  unsigned idxMuon = kRecoMuonIdx(slc);
  std::vector<double> idxProtons = kRecoProtonIndices(slc);
  const auto &muTrk = slc->reco.pfp.at(idxMuon).trk;

  double lengthReq = 10.;

  for ( unsigned i = 0; i < slc->reco.pfp.size(); i++ ) {
    if ( i == idxMuon || std::find(idxProtons.begin(), idxProtons.end(), i) != idxProtons.end() ) continue;
    const auto &pfp = slc->reco.pfp.at(i);

    double trkStartDiff = 9999.;
    double shwStartDiff = 9999.;
    if ( !isnan(pfp.trk.start.x) && !isnan(pfp.trk.start.y) && !isnan(pfp.trk.start.z) )
      trkStartDiff = std::hypot( (pfp.trk.start.x - muTrk.end.x), (pfp.trk.start.y - muTrk.end.y), (pfp.trk.start.z - muTrk.end.z) );
    if ( !isnan(pfp.shw.start.x) && !isnan(pfp.shw.start.y) && !isnan(pfp.shw.start.z) )
      shwStartDiff = std::hypot( (pfp.shw.start.x - muTrk.end.x), (pfp.shw.start.y - muTrk.end.y), (pfp.shw.start.z - muTrk.end.z) );
    if ( std::min(trkStartDiff, shwStartDiff) < 5. ) continue;

    double trackLength = pfp.trk.len;
    double trackScore = pfp.trackScore;
    double chi2muon = ( !isnan(pfp.trk.chi2pid[2].chi2_muon) ) ? (double) pfp.trk.chi2pid[2].chi2_muon : -9999.;
    double chi2proton = ( !isnan(pfp.trk.chi2pid[2].chi2_proton) ) ? (double) pfp.trk.chi2pid[2].chi2_proton : -9999.;

    bool contained = (!isnan(pfp.trk.end.x) &&
                            ((pfp.trk.end.x < -61.94 - 10 && pfp.trk.end.x > -358.49 + 10) ||
                                (pfp.trk.end.x >  61.94 + 10 && pfp.trk.end.x <  358.49 - 10)) &&
                            !isnan(pfp.trk.end.y) &&
                            ( pfp.trk.end.y > -181.86 + 10 && pfp.trk.end.y < 134.96 - 10 ) &&
                            !isnan(pfp.trk.end.z) &&
                            ( pfp.trk.end.z > -894.95 + 10 && pfp.trk.end.z < 894.95 - 10 ) );

    if ( contained && trackLength > lengthReq && chi2muon < 30 && chi2proton > 60 && trackScore > 0.45 ) {
      idx = i;
      break;
    }
  }

  return ( idx );
  });

const Var kNoExtraMIPPrimariesOnly([](const caf::SRSliceProxy* slc) -> float {
  bool noMIP = true;
  unsigned idxMuon = kRecoMuonIdx(slc);
  std::vector<double> idxProtons = kRecoProtonIndices(slc);
  const auto &muTrk = slc->reco.pfp.at(idxMuon).trk;

  double lengthReq = 10.;

  for ( unsigned i = 0; i < slc->reco.pfp.size(); i++ ) {
    if ( i == idxMuon || std::find(idxProtons.begin(), idxProtons.end(), i) != idxProtons.end() ) continue;
    const auto &pfp = slc->reco.pfp.at(i);
    if ( !pfp.parent_is_primary ) continue;

    double trkStartDiff = 9999.;
    double shwStartDiff = 9999.;
    if ( !isnan(pfp.trk.start.x) && !isnan(pfp.trk.start.y) && !isnan(pfp.trk.start.z) )
      trkStartDiff = std::hypot( (pfp.trk.start.x - muTrk.end.x), (pfp.trk.start.y - muTrk.end.y), (pfp.trk.start.z - muTrk.end.z) );
    if ( !isnan(pfp.shw.start.x) && !isnan(pfp.shw.start.y) && !isnan(pfp.shw.start.z) )
      shwStartDiff = std::hypot( (pfp.shw.start.x - muTrk.end.x), (pfp.shw.start.y - muTrk.end.y), (pfp.shw.start.z - muTrk.end.z) );
    if ( std::min(trkStartDiff, shwStartDiff) < 5. ) continue;

    double trackLength = pfp.trk.len;
    double trackScore = pfp.trackScore;
    double chi2muon = ( !isnan(pfp.trk.chi2pid[2].chi2_muon) ) ? (double) pfp.trk.chi2pid[2].chi2_muon : -9999.;
    double chi2proton = ( !isnan(pfp.trk.chi2pid[2].chi2_proton) ) ? (double) pfp.trk.chi2pid[2].chi2_proton : -9999.;

    bool contained = (!isnan(pfp.trk.end.x) &&
                            ((pfp.trk.end.x < -61.94 - 10 && pfp.trk.end.x > -358.49 + 10) ||
                                (pfp.trk.end.x >  61.94 + 10 && pfp.trk.end.x <  358.49 - 10)) &&
                            !isnan(pfp.trk.end.y) &&
                            ( pfp.trk.end.y > -181.86 + 10 && pfp.trk.end.y < 134.96 - 10 ) &&
                            !isnan(pfp.trk.end.z) &&
                            ( pfp.trk.end.z > -894.95 + 10 && pfp.trk.end.z < 894.95 - 10 ) );

    if ( contained && trackLength > lengthReq && chi2muon < 30 && chi2proton > 60 && trackScore > 0.45 ) {
      noMIP = false;
      break;
    }
  }

  return ( noMIP );
  });

const Var kNoExtraShower([](const caf::SRSliceProxy* slc) -> float {
  bool noShower = true;
  unsigned idxMuon = kRecoMuonIdx(slc);
  std::vector<double> idxProtons = kRecoProtonIndices(slc);
  const auto &muTrk = slc->reco.pfp.at(idxMuon).trk;

  double lengthReq = 10.;

  for ( unsigned i = 0; i < slc->reco.pfp.size(); i++ ) {
    if ( i == idxMuon || std::find(idxProtons.begin(), idxProtons.end(), i) != idxProtons.end() ) continue;
    const auto &pfp = slc->reco.pfp.at(i);

    double trkStartDiff = 9999.;
    double shwStartDiff = 9999.;
    if ( !isnan(pfp.trk.start.x) && !isnan(pfp.trk.start.y) && !isnan(pfp.trk.start.z) )
      trkStartDiff = std::hypot( (pfp.trk.start.x - muTrk.end.x), (pfp.trk.start.y - muTrk.end.y), (pfp.trk.start.z - muTrk.end.z) );
    if ( !isnan(pfp.shw.start.x) && !isnan(pfp.shw.start.y) && !isnan(pfp.shw.start.z) )
      shwStartDiff = std::hypot( (pfp.shw.start.x - muTrk.end.x), (pfp.shw.start.y - muTrk.end.y), (pfp.shw.start.z - muTrk.end.z) );
    if ( std::min(trkStartDiff, shwStartDiff) < 5. ) continue;

    double showerLength = pfp.shw.len;
    double trackScore = pfp.trackScore;

    double atslc = -9999.;
    if ( !(isnan(pfp.shw.start.x) || isnan(pfp.shw.start.y) || isnan(pfp.shw.start.z)) ) {
      atslc = std::hypot(slc->vertex.x - pfp.shw.start.x,
                           slc->vertex.y - pfp.shw.start.y,
                           slc->vertex.z - pfp.shw.start.z);
    }

    if ( showerLength > lengthReq && atslc > 5. && trackScore < 0.45 ) {
      noShower = false;
      break;
    }
  }

  return ( noShower );
  });

const Cut kNoExtraShowerCut([](const caf::SRSliceProxy* slc) {
  bool noShower = kNoExtraShower(slc);
  return ( noShower );
  });

const Var kExtraShowerIdx([](const caf::SRSliceProxy* slc) -> float {
  int idx = -5;
  unsigned idxMuon = kRecoMuonIdx(slc);
  std::vector<double> idxProtons = kRecoProtonIndices(slc);
  const auto &muTrk = slc->reco.pfp.at(idxMuon).trk;

  double lengthReq = 10.;

  for ( unsigned i = 0; i < slc->reco.pfp.size(); i++ ) {
    if ( i == idxMuon || std::find(idxProtons.begin(), idxProtons.end(), i) != idxProtons.end() ) continue;
    const auto &pfp = slc->reco.pfp.at(i);

    double trkStartDiff = 9999.;
    double shwStartDiff = 9999.;
    if ( !isnan(pfp.trk.start.x) && !isnan(pfp.trk.start.y) && !isnan(pfp.trk.start.z) )
      trkStartDiff = std::hypot( (pfp.trk.start.x - muTrk.end.x), (pfp.trk.start.y - muTrk.end.y), (pfp.trk.start.z - muTrk.end.z) );
    if ( !isnan(pfp.shw.start.x) && !isnan(pfp.shw.start.y) && !isnan(pfp.shw.start.z) )
      shwStartDiff = std::hypot( (pfp.shw.start.x - muTrk.end.x), (pfp.shw.start.y - muTrk.end.y), (pfp.shw.start.z - muTrk.end.z) );
    if ( std::min(trkStartDiff, shwStartDiff) < 5. ) continue;

    double showerLength = pfp.shw.len;
    double trackScore = pfp.trackScore;
    
    double atslc = -9999.;
    if ( !(isnan(pfp.shw.start.x) || isnan(pfp.shw.start.y) || isnan(pfp.shw.start.z)) ) {
      atslc = std::hypot(slc->vertex.x - pfp.shw.start.x,
                           slc->vertex.y - pfp.shw.start.y,
                           slc->vertex.z - pfp.shw.start.z);
    }
    
    if ( showerLength > lengthReq && atslc > 5. && trackScore < 0.45 ) {
      idx = i;
      break;
    }
  }
  
  return ( idx );
  });

const Var kNoExtraShowerPrimariesOnly([](const caf::SRSliceProxy* slc) -> float {
  bool noShower = true;
  unsigned idxMuon = kRecoMuonIdx(slc);
  std::vector<double> idxProtons = kRecoProtonIndices(slc);
  const auto &muTrk = slc->reco.pfp.at(idxMuon).trk;

  double lengthReq = 10.;

  for ( unsigned i = 0; i < slc->reco.pfp.size(); i++ ) {
    if ( i == idxMuon || std::find(idxProtons.begin(), idxProtons.end(), i) != idxProtons.end() ) continue;
    const auto &pfp = slc->reco.pfp.at(i);
    if ( !pfp.parent_is_primary ) continue;

    double trkStartDiff = 9999.;
    double shwStartDiff = 9999.;
    if ( !isnan(pfp.trk.start.x) && !isnan(pfp.trk.start.y) && !isnan(pfp.trk.start.z) )
      trkStartDiff = std::hypot( (pfp.trk.start.x - muTrk.end.x), (pfp.trk.start.y - muTrk.end.y), (pfp.trk.start.z - muTrk.end.z) );
    if ( !isnan(pfp.shw.start.x) && !isnan(pfp.shw.start.y) && !isnan(pfp.shw.start.z) )
      shwStartDiff = std::hypot( (pfp.shw.start.x - muTrk.end.x), (pfp.shw.start.y - muTrk.end.y), (pfp.shw.start.z - muTrk.end.z) );
    if ( std::min(trkStartDiff, shwStartDiff) < 5. ) continue;

    double showerLength = pfp.shw.len;
    double trackScore = pfp.trackScore;

    double atslc = -9999.;
    if ( !(isnan(pfp.shw.start.x) || isnan(pfp.shw.start.y) || isnan(pfp.shw.start.z)) ) {
      atslc = std::hypot(slc->vertex.x - pfp.shw.start.x,
                           slc->vertex.y - pfp.shw.start.y,
                           slc->vertex.z - pfp.shw.start.z);
    }

    if ( showerLength > lengthReq && atslc > 5. && trackScore < 0.45 ) {
      noShower = false;
      break;
    }
  }

  return ( noShower );
  });

const Cut kHadronicContainment([](const caf::SRSliceProxy* slc) {
  bool contained = true;
  if ( kRecoMuonIdx(slc) < 0 )  return false;
  unsigned idxMuon = kRecoMuonIdx(slc);
  const auto &muTrk = slc->reco.pfp.at(idxMuon).trk;
  
  for ( unsigned i = 0; i < slc->reco.pfp.size(); i++ ) {
    if ( i == idxMuon ) continue;
    const auto &pfp = slc->reco.pfp.at(i);

    double trkStartDiff = 9999.;
    double shwStartDiff = 9999.;
    if ( !isnan(pfp.trk.start.x) && !isnan(pfp.trk.start.y) && !isnan(pfp.trk.start.z) )
      trkStartDiff = std::hypot( (pfp.trk.start.x - muTrk.end.x), (pfp.trk.start.y - muTrk.end.y), (pfp.trk.start.z - muTrk.end.z) );
    if ( !isnan(pfp.shw.start.x) && !isnan(pfp.shw.start.y) && !isnan(pfp.shw.start.z) )
      shwStartDiff = std::hypot( (pfp.shw.start.x - muTrk.end.x), (pfp.shw.start.y - muTrk.end.y), (pfp.shw.start.z - muTrk.end.z) );
    if ( std::min(trkStartDiff, shwStartDiff) < 5. ) continue;

    if ( !isnan(pfp.trk.end.x) && !isnan(pfp.trk.end.y) && !isnan(pfp.trk.end.z) ) {
      contained = (
                    ((pfp.trk.end.x < -61.94 - 10 && pfp.trk.end.x > -358.49 + 10) ||
                      (pfp.trk.end.x >  61.94 + 10 && pfp.trk.end.x <  358.49 - 10)) &&
                    ( pfp.trk.end.y > -181.86 + 10 && pfp.trk.end.y < 134.96 - 10 ) &&
                    ( pfp.trk.end.z > -894.95 + 10 && pfp.trk.end.z < 894.95 - 10 )
                  );
      if ( !contained ) break;
    }

    if ( !isnan(pfp.shw.end.x) && !isnan(pfp.shw.end.y) && !isnan(pfp.shw.end.z) ) {
      contained = (
                    ((pfp.shw.end.x < -61.94 - 10 && pfp.shw.end.x > -358.49 + 10) ||
                      (pfp.shw.end.x >  61.94 + 10 && pfp.shw.end.x <  358.49 - 10)) &&
                    ( pfp.shw.end.y > -181.86 + 10 && pfp.shw.end.y < 134.96 - 10 ) &&
                    ( pfp.shw.end.z > -894.95 + 10 && pfp.shw.end.z < 894.95 - 10 )
                  );
      if ( !contained ) break;
    }
  }
 
  return ( contained );
  });

const Var kHadronicContainmentPrimariesOnly([](const caf::SRSliceProxy* slc) -> float {
  bool contained = true;
  unsigned idxMuon = kRecoMuonIdx(slc);

  for ( unsigned i = 0; i < slc->reco.pfp.size(); i++ ) {
    if ( i == idxMuon ) continue;
    const auto &pfp = slc->reco.pfp.at(i);
    if ( !pfp.parent_is_primary ) continue;
    if ( !isnan(pfp.trk.end.x) && !isnan(pfp.trk.end.y) && !isnan(pfp.trk.end.x) ) {
      contained = (
                    ((pfp.trk.end.x < -61.94 - 10 && pfp.trk.end.x > -358.49 + 10) ||
                      (pfp.trk.end.x >  61.94 + 10 && pfp.trk.end.x <  358.49 - 10)) &&
                    ( pfp.trk.end.y > -181.86 + 10 && pfp.trk.end.y < 134.96 - 10 ) &&
                    ( pfp.trk.end.z > -894.95 + 10 && pfp.trk.end.z < 894.95 - 10 )
                  );
    }
    if ( !contained ) break;
  }

  return ( contained );
  });

const Var kFurthestFromVtx([](const caf::SRSliceProxy* slc) -> float {
  double dist = -1.;
  std::vector<double> idxProtons = kRecoProtonIndices(slc);

  for ( const auto &i : idxProtons ) {
    const auto &trk = slc->reco.pfp.at(i).trk;
    double atslc = std::hypot(slc->vertex.x - trk.start.x,
                                       slc->vertex.y - trk.start.y,
                                       slc->vertex.z - trk.start.z);
    if ( atslc > dist ) dist = atslc;
  }

  return ( dist );
  });

//PFP index of LONGEST primary that's not a reco muon/proton and not obviously a muon daughter PFP
const Var kExtraPrimaryIndex([](const caf::SRSliceProxy* slc) -> float {
  double longest = -9999.;
  int idx = -5;

  if ( kRecoMuonIdx(slc) < 0 ) return idx;
  unsigned idxMuon = kRecoMuonIdx(slc);
  const auto &muTrk = slc->reco.pfp.at(idxMuon).trk;
  std::vector<double> idxProtons = kRecoProtonIndices(slc);

  unsigned idxPion = 9999;
  if ( kHasSidebandPion(slc) ) idxPion = kSidebandPion(slc);

  for ( unsigned i = 1; i < slc->reco.pfp.size(); i++ ) {
    const auto &pfp = slc->reco.pfp.at(i);
    if ( !pfp.parent_is_primary || i == idxMuon || std::find(idxProtons.begin(), idxProtons.end(), i) != idxProtons.end() || i == idxPion ) continue;

    double trkStartDiff = 9999.;
    double shwStartDiff = 9999.;
    if ( !isnan(pfp.trk.start.x) && !isnan(pfp.trk.start.y) && !isnan(pfp.trk.start.z) )
      trkStartDiff = std::hypot( (pfp.trk.start.x - muTrk.end.x), (pfp.trk.start.y - muTrk.end.y), (pfp.trk.start.z - muTrk.end.z) );
    if ( !isnan(pfp.shw.start.x) && !isnan(pfp.shw.start.y) && !isnan(pfp.shw.start.z) )
      shwStartDiff = std::hypot( (pfp.shw.start.x - muTrk.end.x), (pfp.shw.start.y - muTrk.end.y), (pfp.shw.start.z - muTrk.end.z) );
    if ( std::min(trkStartDiff, shwStartDiff) < 5. ) continue;

    if ( !isnan(pfp.pfochar.linfitlen) && pfp.pfochar.linfitlen > longest ) {
      longest = pfp.pfochar.linfitlen;
      idx = i;
    }
  }

  return idx;
  });

const Var kExtraPrimaryLinFitLength([](const caf::SRSliceProxy* slc) -> float {
  double length = -9999.;
  int idx  = kExtraPrimaryIndex(slc);

  if ( idx > -1 ) {
    const auto &pfp = slc->reco.pfp.at(idx);
    length = pfp.pfochar.linfitlen;
  }

  return length;
  });

const Cut kExtraPrimaryLinFitLengthCut([](const caf::SRSliceProxy* slc) {
  double length = kExtraPrimaryLinFitLength(slc);
  return (length < 10.);
  });

/*
const Var kDynamicLengthVar([](const caf::SRSliceProxy* slc) -> int {
  double length = kExtraPrimaryLinFitLength(slc);
  double leadingProton = slc->reco.pfp.at(kRecoProtonIdx(slc)).pfochar.linfitlen;
  return (length < 10. && length < leadingProton);
  });
*/

const Cut kDynamicLengthCut([](const caf::SRSliceProxy* slc) {
  double length = kExtraPrimaryLinFitLength(slc);
  double leadingProton = (kRecoProtonIdx(slc) > -1 ) ? double(slc->reco.pfp.at(kRecoProtonIdx(slc)).pfochar.linfitlen) : 10.;
  return (length < 10. && length < leadingProton);
  });

const Cut kScndProtonCandidate = ( kScndProtonThreshold && kHadronicContainment && !kHasSidebandPion && kDynamicLengthCut );

//const Cut kThreePSideband = ( kThirdRENAMEIFNEEDEDPrimaryContained && kProtonLikeTrack && kScndProtonThreshold && kThreeRecoProtons && kHasSidebandProton && kThreeRecoProtonsONLY && kContainedSidebandHadron );

const Cut kPionSidebandBase = ( kLeadingProtonThreshold && kHasSidebandPion && kHadronicContainment && kExtraPrimaryLinFitLengthCut);

const Cut kPionSideband1p([](const caf::SRSliceProxy* slc) {
  return (kPionSidebandBase(slc) && kNRecoProtons(slc) == 1);
  });

const Cut kPionSidebandNp([](const caf::SRSliceProxy* slc) {
  return (kPionSidebandBase(slc) && kNRecoProtons(slc) > 1);
  });

const Cut kSignalOrSideband = ( kScndProtonCandidate || kPionSideband1p || kPionSidebandNp );

const Var kCutType([](const caf::SRSliceProxy* slc) -> int {
  if ( kScndProtonCandidate(slc) ) return 0;
  else if ( kPionSideband1p(slc) ) return 1;
  else if ( kPionSidebandNp(slc) ) return 2;
  else return -1;
  });

const Var kNPFPs([](const caf::SRSliceProxy* slc) -> float {
  return slc->reco.npfp;
  });

const Var kPrintNonprimaryPFPs([](const caf::SRSliceProxy* slc) -> float {
  if ( isnan(slc->vertex.x) || isnan(slc->vertex.y) || isnan(slc->vertex.z)  || slc->reco.npfp < 5) return 0.;

  std::cout << std::endl << "=============== Nonprimary PFPs";
  for ( const auto &pfp : slc->reco.pfp ) {
    if ( pfp.parent_is_primary ) continue;
    const float Atslc = std::hypot(slc->vertex.x - pfp.trk.start.x,
                                   slc->vertex.y - pfp.trk.start.y,
                                   slc->vertex.z - pfp.trk.start.z);
    std::cout << std::endl << "True particle: " << pfp.trk.truth.p.pdg << ", distance to vertex: " << Atslc << ", track score: " << pfp.trackScore << ", track length: " << pfp.trk.len << ", shower energy: " << pfp.shw.bestplane_energy;
  }
  return 0.;
  });

const Var kExtraEnergy([](const caf::SRSliceProxy* slc) -> float {
  double energy = 0.;
  std::vector<unsigned> trackedPrimaries;

  if ( kRecoMuonIdx(slc) != -1 ) trackedPrimaries.push_back(kRecoMuonIdx(slc));
  if ( kRecoProtonIdx(slc) != -1 ) trackedPrimaries.push_back(kRecoProtonIdx(slc));
  if ( kScndProtonIdx(slc) != -1 ) trackedPrimaries.push_back(kScndProtonIdx(slc));
  //if ( kSidebandProton(slc) != -1 ) trackedPrimaries.push_back(kSidebandProton(slc));
  if ( kSidebandPion(slc) != -1 ) trackedPrimaries.push_back(kSidebandPion(slc));

  for ( unsigned i = 0; i < slc->reco.pfp.size(); i++ ) {
    const auto &pfp = slc->reco.pfp.at(i);
    if ( !pfp.parent_is_primary || std::find(trackedPrimaries.begin(), trackedPrimaries.end(), i) != trackedPrimaries.end() ) continue;
    if ( !isnan(pfp.shw.bestplane_energy) && pfp.shw.bestplane_energy != -5. ) energy += pfp.shw.bestplane_energy;
  }

  return ( energy > 0. ? energy : -9999. );
  });

const Var kVertexX = SIMPLEVAR(vertex.x);

const Var kVertexY = SIMPLEVAR(vertex.y);

const Var kVertexZ = SIMPLEVAR(vertex.z);

const Var kDistToEdgeX([](const caf::SRSliceProxy* slc) -> float {
  double vertex = kVertexX(slc);
  int sign = (vertex < 0) ? -1 : 1;
  double upperBound(358.49*sign), lowerBound(61.94*sign);
  double dist = std::min( abs(vertex - lowerBound), abs(vertex - upperBound) );
  return dist;
  });

const Var kDistToEdgeY([](const caf::SRSliceProxy* slc) -> float {
  double vertex = kVertexY(slc);
  double upperBound(134.96), lowerBound(-181.86);
  double dist = std::min( abs(vertex - lowerBound), abs(vertex - upperBound) );
  return dist;
  });

const Var kDistToEdgeZ([](const caf::SRSliceProxy* slc) -> float {
  double vertex = kVertexZ(slc);
  double upperBound(894.95), lowerBound(-894.95);
  double dist = std::min( abs(vertex - lowerBound), abs(vertex - upperBound) );  
  return dist;
  });

const Var kRecoMuonTruePDG([](const caf::SRSliceProxy* slc) -> int {
  int pdg = -9999;;

  if ( kRecoMuonIdx(slc) >= 0 ) {
    const auto &trk = slc->reco.pfp.at(kRecoMuonIdx(slc)).trk;
    if ( !isnan(trk.truth.p.pdg) ) pdg = trk.truth.p.pdg;
  }

  return pdg;
  });

const Var kRecoProtonTruePDG([](const caf::SRSliceProxy* slc) -> int {
  int pdg = -9999;;

  if ( kRecoProtonIdx(slc) >= 0 ) {
    const auto &trk = slc->reco.pfp.at(kRecoProtonIdx(slc)).trk;
    if ( !isnan(trk.truth.p.pdg) ) pdg = trk.truth.p.pdg;
  }

  return pdg;
  });

const Var kScndProtonTruePDG([](const caf::SRSliceProxy* slc) -> int {
  int pdg = -9999;;

  if ( kScndProtonIdx(slc) >= 0 ) {
    const auto &trk = slc->reco.pfp.at(kScndProtonIdx(slc)).trk;
    if ( !isnan(trk.truth.p.pdg) ) pdg = trk.truth.p.pdg;
  }

  return pdg;
  });

const Var kRecoMuonChi2Muon([](const caf::SRSliceProxy* slc) -> float {
    float chi2muon(-5.f);

    if ( kRecoMuonIdx(slc) >= 0 )
    {
      auto const& trk = slc->reco.pfp.at(kRecoMuonIdx(slc)).trk;
      const bool Contained = ( !isnan(trk.end.x) &&
                               ((trk.end.x < -61.94 - 10 && trk.end.x > -358.49 + 10) ||
                                      (trk.end.x >  61.94 + 10 && trk.end.x <  358.49 - 10)) &&
                               !isnan(trk.end.y) &&
                               ( trk.end.y > -181.86 + 10 && trk.end.y < 134.96 - 10 ) &&
                               !isnan(trk.end.z) &&
                               ( trk.end.z > -894.95 + 10 && trk.end.z < 894.95 - 10 ) );
      if(Contained) chi2muon = trk.chi2pid[2].chi2_muon;
    }

    return chi2muon;
  });

const Var kRecoMuonChi2Proton([](const caf::SRSliceProxy* slc) -> float {
    float chi2proton(-5.f);

    if ( kRecoMuonIdx(slc) >= 0 )
    {
      auto const& trk = slc->reco.pfp.at(kRecoMuonIdx(slc)).trk;
      const bool Contained = ( !isnan(trk.end.x) &&
                               ((trk.end.x < -61.94 - 10 && trk.end.x > -358.49 + 10) ||
                                      (trk.end.x >  61.94 + 10 && trk.end.x <  358.49 - 10)) &&
                               !isnan(trk.end.y) &&
                               ( trk.end.y > -181.86 + 10 && trk.end.y < 134.96 - 10 ) &&
                               !isnan(trk.end.z) &&
                               ( trk.end.z > -894.95 + 10 && trk.end.z < 894.95 - 10 ) );
      if(Contained) chi2proton = trk.chi2pid[2].chi2_proton;
    }

    return chi2proton;
  });

const Var kRecoMuonChi2Pion([](const caf::SRSliceProxy* slc) -> float {
    float chi2pion(-5.f);

    if ( kRecoMuonIdx(slc) >= 0 )
    {
      auto const& trk = slc->reco.pfp.at(kRecoMuonIdx(slc)).trk;
      const bool Contained = ( !isnan(trk.end.x) &&
                               ((trk.end.x < -61.94 - 10 && trk.end.x > -358.49 + 10) ||
                                      (trk.end.x >  61.94 + 10 && trk.end.x <  358.49 - 10)) &&
                               !isnan(trk.end.y) &&
                               ( trk.end.y > -181.86 + 10 && trk.end.y < 134.96 - 10 ) &&
                               !isnan(trk.end.z) &&
                               ( trk.end.z > -894.95 + 10 && trk.end.z < 894.95 - 10 ) );
      if(Contained) chi2pion = trk.chi2pid[2].chi2_pion;
    }

    return chi2pion;
  });

const Var kRecoProtonChi2Muon([](const caf::SRSliceProxy* slc) -> float {
    float chi2muon(-5.f);

    if ( kRecoProtonIdx(slc) >= 0 )
    {
      auto const& trk = slc->reco.pfp.at(kRecoProtonIdx(slc)).trk;
      const bool Contained = ( !isnan(trk.end.x) &&
                               ((trk.end.x < -61.94 - 10 && trk.end.x > -358.49 + 10) ||
                                      (trk.end.x >  61.94 + 10 && trk.end.x <  358.49 - 10)) &&
                               !isnan(trk.end.y) &&
                               ( trk.end.y > -181.86 + 10 && trk.end.y < 134.96 - 10 ) &&
                               !isnan(trk.end.z) &&
                               ( trk.end.z > -894.95 + 10 && trk.end.z < 894.95 - 10 ) );
      if(Contained) chi2muon = trk.chi2pid[2].chi2_muon;
    }

    return chi2muon;
  });

const Var kRecoProtonChi2Proton([](const caf::SRSliceProxy* slc) -> float {
    float chi2proton(-5.f);

    if ( kRecoProtonIdx(slc) >= 0 )
    {
      auto const& trk = slc->reco.pfp.at(kRecoProtonIdx(slc)).trk;
      const bool Contained = ( !isnan(trk.end.x) &&
                               ((trk.end.x < -61.94 - 10 && trk.end.x > -358.49 + 10) ||
                                      (trk.end.x >  61.94 + 10 && trk.end.x <  358.49 - 10)) &&
                               !isnan(trk.end.y) &&
                               ( trk.end.y > -181.86 + 10 && trk.end.y < 134.96 - 10 ) &&
                               !isnan(trk.end.z) &&
                               ( trk.end.z > -894.95 + 10 && trk.end.z < 894.95 - 10 ) );
      if(Contained) chi2proton = trk.chi2pid[2].chi2_proton;
    }

    return chi2proton;
  });

const Var kRecoProtonChi2Pion([](const caf::SRSliceProxy* slc) -> float {
    float chi2pion(-5.f);

    if ( kRecoProtonIdx(slc) >= 0 )
    {
      auto const& trk = slc->reco.pfp.at(kRecoProtonIdx(slc)).trk;
      const bool Contained = ( !isnan(trk.end.x) &&
                               ((trk.end.x < -61.94 - 10 && trk.end.x > -358.49 + 10) ||
                                      (trk.end.x >  61.94 + 10 && trk.end.x <  358.49 - 10)) &&
                               !isnan(trk.end.y) &&
                               ( trk.end.y > -181.86 + 10 && trk.end.y < 134.96 - 10 ) &&
                               !isnan(trk.end.z) &&
                               ( trk.end.z > -894.95 + 10 && trk.end.z < 894.95 - 10 ) );
      if(Contained) chi2pion = trk.chi2pid[2].chi2_pion;
    }

    return chi2pion;
  });

const Var kScndProtonChi2Muon([](const caf::SRSliceProxy* slc) -> float {
    float chi2muon(-5.f);

    if ( kScndProtonIdx(slc) >= 0 )
    {
      auto const& trk = slc->reco.pfp.at(kScndProtonIdx(slc)).trk;
      const bool Contained = ( !isnan(trk.end.x) &&
                               ((trk.end.x < -61.94 - 10 && trk.end.x > -358.49 + 10) ||
                                      (trk.end.x >  61.94 + 10 && trk.end.x <  358.49 - 10)) &&
                               !isnan(trk.end.y) &&
                               ( trk.end.y > -181.86 + 10 && trk.end.y < 134.96 - 10 ) &&
                               !isnan(trk.end.z) &&
                               ( trk.end.z > -894.95 + 10 && trk.end.z < 894.95 - 10 ) );
      if( Contained && !isnan(trk.chi2pid[2].chi2_muon) && trk.chi2pid[2].chi2_muon != 0 ) chi2muon = trk.chi2pid[2].chi2_muon;
    }

    return chi2muon;
  });

const Var kScndProtonChi2Proton([](const caf::SRSliceProxy* slc) -> float {
    float chi2proton(-5.f);

    if ( kScndProtonIdx(slc) >= 0 )
    {
      auto const& trk = slc->reco.pfp.at(kScndProtonIdx(slc)).trk;
      const bool Contained = ( !isnan(trk.end.x) &&
                               ((trk.end.x < -61.94 - 10 && trk.end.x > -358.49 + 10) ||
                                      (trk.end.x >  61.94 + 10 && trk.end.x <  358.49 - 10)) &&
                               !isnan(trk.end.y) &&
                               ( trk.end.y > -181.86 + 10 && trk.end.y < 134.96 - 10 ) &&
                               !isnan(trk.end.z) &&
                               ( trk.end.z > -894.95 + 10 && trk.end.z < 894.95 - 10 ) );
      if( Contained && !isnan(trk.chi2pid[2].chi2_proton) && trk.chi2pid[2].chi2_proton != 0 ) chi2proton = trk.chi2pid[2].chi2_proton;
    }

    return chi2proton;
  });

const Var kScndProtonChi2Pion([](const caf::SRSliceProxy* slc) -> float {
    float chi2pion(-5.f);

    if ( kScndProtonIdx(slc) >= 0 )
    {
      auto const& trk = slc->reco.pfp.at(kScndProtonIdx(slc)).trk;
      const bool Contained = ( !isnan(trk.end.x) &&
                               ((trk.end.x < -61.94 - 10 && trk.end.x > -358.49 + 10) ||
                                      (trk.end.x >  61.94 + 10 && trk.end.x <  358.49 - 10)) &&
                               !isnan(trk.end.y) &&
                               ( trk.end.y > -181.86 + 10 && trk.end.y < 134.96 - 10 ) &&
                               !isnan(trk.end.z) &&
                               ( trk.end.z > -894.95 + 10 && trk.end.z < 894.95 - 10 ) );
      if( Contained && trk.chi2pid[2].chi2_pion != 0 ) chi2pion = trk.chi2pid[2].chi2_pion;
    }

    return chi2pion;
  });

const Var kScndProtonChi2ProtonNEGATIVE([](const caf::SRSliceProxy* slc) -> float {
    float chi2proton(-5.f);

    if ( kScndProtonIdx(slc) >= 0 )
    {
      auto const& trk = slc->reco.pfp.at(kScndProtonIdx(slc)).trk;
      const bool Contained = ( !isnan(trk.end.x) &&
                               ((trk.end.x < -61.94 - 10 && trk.end.x > -358.49 + 10) ||
                                      (trk.end.x >  61.94 + 10 && trk.end.x <  358.49 - 10)) &&
                               !isnan(trk.end.y) &&
                               ( trk.end.y > -181.86 + 10 && trk.end.y < 134.96 - 10 ) &&
                               !isnan(trk.end.z) &&
                               ( trk.end.z > -894.95 + 10 && trk.end.z < 894.95 - 10 ) );
      if( Contained && trk.chi2pid[2].chi2_proton != 0 ) chi2proton = trk.chi2pid[2].chi2_proton;
    }

    return -chi2proton;
  });

const Var kScndProtonNDaughters([](const caf::SRSliceProxy* slc) -> float {
    int nDaughters = -1;

    if ( kScndProtonIdx(slc) >= 0 )
    {
      auto const& pfp = slc->reco.pfp.at(kScndProtonIdx(slc));
      nDaughters = pfp.daughters.size();
    }

    return nDaughters;
  });

const Var kRecoProtonNDaughters([](const caf::SRSliceProxy* slc) -> float {
    int nDaughters = -1;

    if ( kRecoProtonIdx(slc) >= 0 )
    { 
      auto const& pfp = slc->reco.pfp.at(kRecoProtonIdx(slc));
      nDaughters = pfp.daughters.size();
    }

    return nDaughters;
  });

const Var kTotalNDaughters([](const caf::SRSliceProxy* slc) -> float {
    int nDaughters = -1;

    int muTrk = kRecoMuonIdx(slc);
    int pTrk1 = kRecoProtonIdx(slc);
    int pTrk2 = kScndProtonIdx(slc);

    if ( muTrk < 0 || pTrk1 < 0 ||  pTrk2 < 0 ) return nDaughters;

    nDaughters = slc->reco.pfp.at(muTrk).ndaughters;
    nDaughters += slc->reco.pfp.at(pTrk1).ndaughters;
    nDaughters += slc->reco.pfp.at(pTrk2).ndaughters;

    return nDaughters;
  });

const Var kScndProtonLength([](const caf::SRSliceProxy* slc) -> float {
    double length = -1.;

    if ( kScndProtonIdx(slc) >= 0 )
    {
      auto const& trk = slc->reco.pfp.at(kScndProtonIdx(slc)).trk;
      length = trk.len;
    }

    return length;
  });

const Var kScndProtonTrueLength([](const caf::SRSliceProxy* slc) -> float {
    double length = -1.;

    if ( kScndProtonIdx(slc) >= 0 )
    {
      auto const& trk = slc->reco.pfp.at(kScndProtonIdx(slc)).trk;
      length = trk.truth.p.length;
    }

    return length;
  });

const Var kScndProtonLengthResid([](const caf::SRSliceProxy* slc) -> float {
    double recoLength = kScndProtonLength(slc);
    double trueLength = kScndProtonTrueLength(slc);

    return (recoLength - trueLength) / trueLength;
  });

const Var kScndProtonPionP([](const caf::SRSliceProxy* slc) -> float {
    double momentum = -1.;

    if ( kScndProtonIdx(slc) >= 0 )
    {
      auto const& trk = slc->reco.pfp.at(kScndProtonIdx(slc)).trk;
      momentum = trk.rangeP.p_pion;
    }

    return momentum;
  });

const Var kScndProtonTrueP([](const caf::SRSliceProxy* slc) -> float {
    double momentum = -1.;

    if ( kScndProtonIdx(slc) >= 0 )
    {
      auto const& trk = slc->reco.pfp.at(kScndProtonIdx(slc)).trk;
      momentum = sqrt(std::pow( trk.truth.p.startp.x, 2 ) + std::pow( trk.truth.p.startp.y, 2 ) + std::pow( trk.truth.p.startp.z, 2 ));
    }

    return momentum;
  });

const Var kScndProtonPionPResid([](const caf::SRSliceProxy* slc) -> float {
    double recoMomentum = kScndProtonPionP(slc);
    double trueMomentum = kScndProtonTrueP(slc);
    
    return (recoMomentum - trueMomentum) / trueMomentum;
  });

/*
const Cut kGoodPionReco([](const caf::SRSliceProxy* slc) {
    if ( !kTrueContainedPion(slc) ) return false;
    double residualMomentum = kScndProtonPionPResid(slc);
    return ( abs(residualMomentum) <= .1 );
  });

const Cut kBadPionReco([](const caf::SRSliceProxy* slc) {
    if ( !kTrueContainedPion(slc) ) return false;
    double residualMomentum = kScndProtonPionPResid(slc);
    double residualLength = kScndProtonLengthResid(slc);
    return ( abs(residualMomentum) > .1 && abs(residualLength) <= .1);
  });
*/

const Cut kTruePiPlusUncontained([](const caf::SRSliceProxy* slc) {
  bool Contained = false;
  bool PDG = false;
  if ( kScndProtonIdx(slc) >= 0 ) {
    const auto &trueParticle = slc->reco.pfp.at(kScndProtonIdx(slc)).trk.truth.p;
    Contained = (!isnan(trueParticle.end.x) &&
                            ((trueParticle.end.x < -61.94 - 10 && trueParticle.end.x > -358.49 + 10) ||
                                (trueParticle.end.x >  61.94 + 10 && trueParticle.end.x <  358.49 - 10)) &&
                            !isnan(trueParticle.end.y) &&
                            ( trueParticle.end.y > -181.86 + 10 && trueParticle.end.y < 134.96 - 10 ) &&
                            !isnan(trueParticle.end.z) &&
                            ( trueParticle.end.z > -894.95 + 10 && trueParticle.end.z < 894.95 - 10 ) );
    PDG = ( trueParticle.pdg == 211 );
    }
  return ( !Contained && PDG );
  });

/*
const SpillVar kPrintFlashes([](const caf::SRSpillProxy* sr) {
  for ( const auto &slc : sr->slc ) {
    if ( kGoodPionReco(&slc) ) {
      const auto &trk = slc.reco.pfp.at(kScndProtonIdx(&slc)).trk;
      int cryo = slc.producer;
      std::cout << std::endl << "Pion Track End X: " << trk.end.x << ", Y: " << trk.end.y << ", Z:  " << trk.end.z;
      //std::cout << std::endl << "Flash times in spill:";
      for ( const auto &flash : sr->opflashes ) if (flash.cryo == cryo) std::cout << std::endl << "Flash Time: "<< flash.time << ", Center X: " << flash.center.x << ", Y:" << flash.center.y << ", Z: " << flash.center.z;
      break;
    }
  }

  return 1;
  });
*/

const Var kRecoProtonEnergyPerLength([](const caf::SRSliceProxy* slc) -> float {
    double ratio = -1.;
    double energy, length;

    if ( kRecoProtonIdx(slc) >= 0 ) {
      auto const& trk = slc->reco.pfp.at(kRecoProtonIdx(slc)).trk;
      energy = trk.calo[2].ke;
      length = trk.len;
      if (  !isnan(energy) && !isnan(length) && energy != 0 && length != 0 ) ratio = energy / length;
    }

    return ratio;
  });

const Var kScndProtonEnergyPerLength([](const caf::SRSliceProxy* slc) -> float {
    double ratio = -1.;
    double energy, length;

    if ( kScndProtonIdx(slc) >= 0 ) {
      auto const& trk = slc->reco.pfp.at(kScndProtonIdx(slc)).trk;
      energy = trk.calo[2].ke;
      length = trk.len;
      if (  !isnan(energy) && !isnan(length) && energy != 0 && length != 0 ) ratio = energy / length;
    }

    return ratio;
  });

const Var kRecoProtonChi2Comp([](const caf::SRSliceProxy* slc) -> float {
    double chi2Mu = kRecoProtonChi2Muon(slc);
    double chi2P = kRecoProtonChi2Proton(slc);

    return (chi2P - chi2Mu) / (chi2P + chi2Mu);
  });

const Var kScndProtonChi2Comp([](const caf::SRSliceProxy* slc) -> float {
    double chi2Mu = kScndProtonChi2Muon(slc);
    double chi2P = kScndProtonChi2Proton(slc);
    
    return (chi2P - chi2Mu) / (chi2P + chi2Mu);
  });

const Var kRecoMuonTrackScore([](const caf::SRSliceProxy* slc) -> float {
    double score = -1.;
    int idx =  kRecoMuonIdx(slc);
    if ( idx >= 0 ) {
    const auto &pfp = slc->reco.pfp.at(idx);
      if (  !isnan(pfp.trackScore) ) score = pfp.trackScore;
    }

    return score;
  });

const Var kRecoProtonTrackScore([](const caf::SRSliceProxy* slc) -> float {
    double score = -1.;
    int idx =  kRecoProtonIdx(slc);
    if ( idx >= 0 ) {
    const auto &pfp = slc->reco.pfp.at(idx);
      if (  !isnan(pfp.trackScore) ) score = pfp.trackScore;
    }

    return score;
  });

const Var kScndProtonTrackScore([](const caf::SRSliceProxy* slc) -> float {
    double score = -1.;
    int idx =  kScndProtonIdx(slc);
    if ( idx >= 0 ) {
    const auto &pfp = slc->reco.pfp.at(idx);
      if (  !isnan(pfp.trackScore) ) score = pfp.trackScore;
    }

    return score;
  });

const Var kRecoProtonTrackLength([](const caf::SRSliceProxy* slc) -> float {
    double length = -1.;
    int idx =  kRecoProtonIdx(slc);
    if ( idx >= 0 ) {
    const auto &pfp = slc->reco.pfp.at(idx);
      if (  !isnan(pfp.trk.len) ) length = pfp.trk.len;
    }

    return length;
  });

const Var kScndProtonTrackLength([](const caf::SRSliceProxy* slc) -> float {
    double length = -1.;
    int idx =  kScndProtonIdx(slc);
    if ( idx >= 0 ) {
    const auto &pfp = slc->reco.pfp.at(idx);
      if (  !isnan(pfp.trk.len) ) length = pfp.trk.len;
    }

    return length;
  });

const MultiVar kScndProtonDaughterPDG([](const caf::SRSliceProxy* slc) {
    std::vector<double> pdgs;
    if ( kScndProtonIdx(slc) >= 0 ) {
      auto const& piTrk = slc->reco.pfp.at(kScndProtonIdx(slc));
      for ( const auto &dtrID : piTrk.daughters ) {
        for ( const auto &pfp : slc->reco.pfp ) {
          if ( pfp.id != dtrID ) continue;
          if ( abs(pfp.trk.truth.p.pdg) == 11 )        pdgs.push_back(0.);
          else if ( abs(pfp.trk.truth.p.pdg) == 13 )   pdgs.push_back(1.);
          else if ( abs(pfp.trk.truth.p.pdg) == 22 )   pdgs.push_back(2.);
          else if ( abs(pfp.trk.truth.p.pdg) == 111 )  pdgs.push_back(3.);
          else if ( pfp.trk.truth.p.pdg == 211 )       pdgs.push_back(4.);
          else if ( pfp.trk.truth.p.pdg == -211 )      pdgs.push_back(5.);
          else if ( abs(pfp.trk.truth.p.pdg) == 2212 ) pdgs.push_back(6.);
          else                                         pdgs.push_back(7.);
        }
      }
    }

    return pdgs;
  });

const MultiVar kScndProtonDaughterTrkScore([](const caf::SRSliceProxy* slc) {
    std::vector<double> scores;
    if ( kScndProtonIdx(slc) >= 0 ) {
      auto const& piTrk = slc->reco.pfp.at(kScndProtonIdx(slc));
      for ( const auto &dtrID : piTrk.daughters ) {
        for ( const auto &pfp : slc->reco.pfp ) {
          if ( pfp.id != dtrID ) continue;
          scores.push_back(pfp.trackScore);
        }
      }
    }

    return scores;
  });

const Var kScndProtonDirChange([](const caf::SRSliceProxy* slc) -> float {
    double costh = -2.;

    if ( kScndProtonIdx(slc) >= 0 ) {
      auto const& trk = slc->reco.pfp.at(kScndProtonIdx(slc)).trk;
      costh = (trk.dir.x * trk.dir_end.x) + (trk.dir.y * trk.dir_end.y) + (trk.dir.z * trk.dir_end.z);
    }

    return (1. - costh);
  });

const Var kScndProtonP([](const caf::SRSliceProxy* slc) -> float {
    double momentum = -1.;

    if ( kScndProtonIdx(slc) >= 0 )
    { 
      auto const& trk = slc->reco.pfp.at(kScndProtonIdx(slc)).trk;
      momentum = trk.rangeP.p_proton;
    }

    return momentum;
  });

const Var kScndProtonPResid([](const caf::SRSliceProxy* slc) -> float {
    double recoMomentum = kScndProtonP(slc);
    double trueMomentum = kScndProtonTrueP(slc);

    return (recoMomentum - trueMomentum) / trueMomentum;
  });

const Var kRecoMuonEndProcess([](const caf::SRSliceProxy* slc) -> float {
    int proc = -5;
    int muonIdx =  kRecoMuonIdx(slc);
    if ( muonIdx >= 0 ) {
      auto const& truth = slc->reco.pfp.at(muonIdx).trk.truth.p;
      if ( !isnan(int(truth.end_process)) ) proc = truth.end_process;
    }
    return proc;
  });

const Var kRecoProtonEndProcess([](const caf::SRSliceProxy* slc) -> float {
    int proc = -5;
    int protonIdx =  kRecoProtonIdx(slc);
    if ( protonIdx >= 0 ) {
      auto const& truth = slc->reco.pfp.at(protonIdx).trk.truth.p;
      if ( !isnan(int(truth.end_process)) ) proc = truth.end_process;
    }
    return proc;
  });

const Var kScndProtonEndProcess([](const caf::SRSliceProxy* slc) -> float {
    int proc = -5;
    int protonIdx =  kScndProtonIdx(slc);
    if ( protonIdx >= 0 ) {
      auto const& truth = slc->reco.pfp.at(protonIdx).trk.truth.p;
      if ( !isnan(int(truth.end_process)) ) proc = truth.end_process;
    }
    return proc;
  });

const Var kRecoMuonStartProcess([](const caf::SRSliceProxy* slc) -> float {
    int proc = -5;
    int muonIdx =  kRecoMuonIdx(slc);
    if ( muonIdx >= 0 ) {
      auto const& truth = slc->reco.pfp.at(muonIdx).trk.truth.p;
      if ( !isnan(int(truth.start_process)) ) proc = truth.start_process;
    }
    return proc;
  });

const Var kRecoProtonStartProcess([](const caf::SRSliceProxy* slc) -> float {
    int proc = -5;
    int protonIdx =  kRecoProtonIdx(slc);
    if ( protonIdx >= 0 ) {
      auto const& truth = slc->reco.pfp.at(protonIdx).trk.truth.p;
      if ( !isnan(int(truth.start_process)) ) proc = truth.start_process;
    }
    return proc;
  });

const Var kScndProtonStartProcess([](const caf::SRSliceProxy* slc) -> float {
    int proc = -5;
    int protonIdx =  kScndProtonIdx(slc);
    if ( protonIdx >= 0 ) {
      auto const& truth = slc->reco.pfp.at(protonIdx).trk.truth.p;
      if ( !isnan(int(truth.start_process)) ) proc = truth.start_process;
    }
    return proc;
  });

const Var kRecoMuonParentPDG([](const caf::SRSliceProxy* slc) -> float {
    int parent = -5;
    int parentID = -5;
    int muonIdx =  kRecoMuonIdx(slc);
    if ( muonIdx >= 0 ) {
      if ( !isnan(slc->reco.pfp.at(muonIdx).trk.truth.p.parent) ) parentID = slc->reco.pfp.at(muonIdx).trk.truth.p.parent;
      if ( parentID >= 0 ) {
        for ( const auto &prim : slc->truth.prim ) {
          if ( prim.G4ID == parentID ) parent = prim.pdg;
        }
      }
    }
    return parent;
  });

const Var kRecoProtonParentPDG([](const caf::SRSliceProxy* slc) -> float {
    int parent = -5;
    int parentID = -5;
    int protonIdx =  kRecoProtonIdx(slc);
    if ( protonIdx >= 0 ) {
      if ( !isnan(slc->reco.pfp.at(protonIdx).trk.truth.p.parent) ) parentID = slc->reco.pfp.at(protonIdx).trk.truth.p.parent;
      if ( parentID >= 0 ) {
        for ( const auto &prim : slc->truth.prim ) {
          if ( prim.G4ID == parentID ) parent = prim.pdg;
        }
      }
    }
    return parent;
  });

const Var kScndProtonParentPDG([](const caf::SRSliceProxy* slc) -> float {
    int parent = -5;
    int parentID = -5;
    int protonIdx =  kScndProtonIdx(slc);
    if ( protonIdx >= 0 ) {
      if ( !isnan(slc->reco.pfp.at(protonIdx).trk.truth.p.parent) ) parentID = slc->reco.pfp.at(protonIdx).trk.truth.p.parent;
      if ( parentID >= 0 ) {
        for ( const auto &prim : slc->truth.prim ) {
          if ( prim.G4ID == parentID ) parent = prim.pdg;
        }
      }
    }
    return parent;
  });

const Var kRecoProtonChgEndFrac([](const caf::SRSliceProxy* slc) -> float {
    double mvaInput = -9999.;
    int protonIdx =  kRecoProtonIdx(slc);
    if ( protonIdx >= 0 ) {
      auto const& pfochar = slc->reco.pfp.at(protonIdx).pfochar;
      if ( !isnan(pfochar.chgendfrac) ) mvaInput = pfochar.chgendfrac;
    }
    return mvaInput;
  });

const Var kScndProtonChgEndFrac([](const caf::SRSliceProxy* slc) -> float {
    double mvaInput = -9999.;
    int protonIdx =  kScndProtonIdx(slc);
    if ( protonIdx >= 0 ) {
      auto const& pfochar = slc->reco.pfp.at(protonIdx).pfochar;
      if ( !isnan(pfochar.chgendfrac) ) mvaInput = pfochar.chgendfrac;
    }
    return mvaInput;
  });

const Var kRecoProtonChgFracSpread([](const caf::SRSliceProxy* slc) -> float {
    double mvaInput = -9999.;
    int protonIdx =  kRecoProtonIdx(slc);
    if ( protonIdx >= 0 ) {
      auto const& pfochar = slc->reco.pfp.at(protonIdx).pfochar;
      if ( !isnan(pfochar.chgfracspread) ) mvaInput = pfochar.chgfracspread;
    }
    return mvaInput;
  });

const Var kScndProtonChgFracSpread([](const caf::SRSliceProxy* slc) -> float {
    double mvaInput = -9999.;
    int protonIdx =  kScndProtonIdx(slc);
    if ( protonIdx >= 0 ) { 
      auto const& pfochar = slc->reco.pfp.at(protonIdx).pfochar;
      if ( !isnan(pfochar.chgfracspread) ) mvaInput = pfochar.chgfracspread;
    }
    return mvaInput;
  });

const Var kRecoProtonLinFitDiff([](const caf::SRSliceProxy* slc) -> float {
    double mvaInput = -9999.;
    int protonIdx =  kRecoProtonIdx(slc);
    if ( protonIdx >= 0 ) {
      auto const& pfochar = slc->reco.pfp.at(protonIdx).pfochar;
      if ( !isnan(pfochar.linfitdiff) ) mvaInput = pfochar.linfitdiff;
    }
    return mvaInput;
  });

const Var kScndProtonLinFitDiff([](const caf::SRSliceProxy* slc) -> float {
    double mvaInput = -9999.;
    int protonIdx =  kScndProtonIdx(slc);
    if ( protonIdx >= 0 ) { 
      auto const& pfochar = slc->reco.pfp.at(protonIdx).pfochar;
      if ( !isnan(pfochar.linfitdiff) ) mvaInput = pfochar.linfitdiff;
    }
    return mvaInput;
  });

const Var kRecoProtonLinFitLen([](const caf::SRSliceProxy* slc) -> float {
    double mvaInput = -9999.;
    int protonIdx =  kRecoProtonIdx(slc);
    if ( protonIdx >= 0 ) {
      auto const& pfochar = slc->reco.pfp.at(protonIdx).pfochar;
      if ( !isnan(pfochar.linfitlen) ) mvaInput = pfochar.linfitlen;
    }
    return mvaInput;
  });

const Var kScndProtonLinFitLen([](const caf::SRSliceProxy* slc) -> float {
    double mvaInput = -9999.;
    int protonIdx =  kScndProtonIdx(slc);
    if ( protonIdx >= 0 ) { 
      auto const& pfochar = slc->reco.pfp.at(protonIdx).pfochar;
      if ( !isnan(pfochar.linfitlen) ) mvaInput = pfochar.linfitlen;
    }
    return mvaInput;
  });

const Var kRecoProtonLinFitGapLen([](const caf::SRSliceProxy* slc) -> float {
    double mvaInput = -9999.;
    int protonIdx =  kRecoProtonIdx(slc);
    if ( protonIdx >= 0 ) {
      auto const& pfochar = slc->reco.pfp.at(protonIdx).pfochar;
      if ( !isnan(pfochar.linfitgaplen) ) mvaInput = pfochar.linfitgaplen;
    }
    return mvaInput;
  });

const Var kScndProtonLinFitGapLen([](const caf::SRSliceProxy* slc) -> float {
    double mvaInput = -9999.;
    int protonIdx =  kScndProtonIdx(slc);
    if ( protonIdx >= 0 ) { 
      auto const& pfochar = slc->reco.pfp.at(protonIdx).pfochar;
      if ( !isnan(pfochar.linfitgaplen) ) mvaInput = pfochar.linfitgaplen;
    }
    return mvaInput;
  });

const Var kRecoProtonLinFitRMS([](const caf::SRSliceProxy* slc) -> float {
    double mvaInput = -9999.;
    int protonIdx =  kRecoProtonIdx(slc);
    if ( protonIdx >= 0 ) {
      auto const& pfochar = slc->reco.pfp.at(protonIdx).pfochar;
      if ( !isnan(pfochar.linfitrms) ) mvaInput = pfochar.linfitrms;
    }
    return mvaInput;
  });

const Var kScndProtonLinFitRMS([](const caf::SRSliceProxy* slc) -> float {
    double mvaInput = -9999.;
    int protonIdx =  kScndProtonIdx(slc);
    if ( protonIdx >= 0 ) { 
      auto const& pfochar = slc->reco.pfp.at(protonIdx).pfochar;
      if ( !isnan(pfochar.linfitrms) ) mvaInput = pfochar.linfitrms;
    }
    return mvaInput;
  });

const Var kRecoProtonOpenAngleDiff([](const caf::SRSliceProxy* slc) -> float {
    double mvaInput = -9999.;
    int protonIdx =  kRecoProtonIdx(slc);
    if ( protonIdx >= 0 ) {
      auto const& pfochar = slc->reco.pfp.at(protonIdx).pfochar;
      if ( !isnan(pfochar.openanglediff) ) mvaInput = pfochar.openanglediff;
    }
    return mvaInput;
  });

const Var kScndProtonOpenAngleDiff([](const caf::SRSliceProxy* slc) -> float {
    double mvaInput = -9999.;
    int protonIdx =  kScndProtonIdx(slc);
    if ( protonIdx >= 0 ) { 
      auto const& pfochar = slc->reco.pfp.at(protonIdx).pfochar;
      if ( !isnan(pfochar.openanglediff) ) mvaInput = pfochar.openanglediff;
    }
    return mvaInput;
  });

const Var kRecoProtonPCA2Ratio([](const caf::SRSliceProxy* slc) -> float {
    double mvaInput = -9999.;
    int protonIdx =  kRecoProtonIdx(slc);
    if ( protonIdx >= 0 ) {
      auto const& pfochar = slc->reco.pfp.at(protonIdx).pfochar;
      if ( !isnan(pfochar.pca2ratio) ) mvaInput = pfochar.pca2ratio;
    }
    return mvaInput;
  });

const Var kScndProtonPCA2Ratio([](const caf::SRSliceProxy* slc) -> float {
    double mvaInput = -9999.;
    int protonIdx =  kScndProtonIdx(slc);
    if ( protonIdx >= 0 ) { 
      auto const& pfochar = slc->reco.pfp.at(protonIdx).pfochar;
      if ( !isnan(pfochar.pca2ratio) ) mvaInput = pfochar.pca2ratio;
    }
    return mvaInput;
  });

const Var kRecoProtonPCA3Ratio([](const caf::SRSliceProxy* slc) -> float {
    double mvaInput = -9999.;
    int protonIdx =  kRecoProtonIdx(slc);
    if ( protonIdx >= 0 ) {
      auto const& pfochar = slc->reco.pfp.at(protonIdx).pfochar;
      if ( !isnan(pfochar.pca3ratio) ) mvaInput = pfochar.pca3ratio;
    }
    return mvaInput;
  });

const Var kScndProtonPCA3Ratio([](const caf::SRSliceProxy* slc) -> float {
    double mvaInput = -9999.;
    int protonIdx =  kScndProtonIdx(slc);
    if ( protonIdx >= 0 ) { 
      auto const& pfochar = slc->reco.pfp.at(protonIdx).pfochar;
      if ( !isnan(pfochar.pca3ratio) ) mvaInput = pfochar.pca3ratio;
    }
    return mvaInput;
  });

const Var kRecoProtonVtxDist([](const caf::SRSliceProxy* slc) -> float {
    double mvaInput = -9999.;
    int protonIdx =  kRecoProtonIdx(slc);
    if ( protonIdx >= 0 ) {
      auto const& pfochar = slc->reco.pfp.at(protonIdx).pfochar;
      if ( !isnan(pfochar.vtxdist) ) mvaInput = pfochar.vtxdist;
    }
    return mvaInput;
  });

const Var kScndProtonVtxDist([](const caf::SRSliceProxy* slc) -> float {
    double mvaInput = -9999.;
    int protonIdx =  kScndProtonIdx(slc);
    if ( protonIdx >= 0 ) {
      auto const& pfochar = slc->reco.pfp.at(protonIdx).pfochar;
      if ( !isnan(pfochar.vtxdist) ) mvaInput = pfochar.vtxdist;
    }
    return mvaInput;
  });

const Var kRecoMuonShowerLength([](const caf::SRSliceProxy* slc) -> float {
    double showerLen = -9999.;
    int idx =  kRecoMuonIdx(slc);
    if ( idx >= 0 ) {
      auto const &shw = slc->reco.pfp.at(idx).shw;
      if ( !isnan(shw.len) ) showerLen = shw.len;
    }
    return showerLen;
  });

const Var kRecoProtonShowerLength([](const caf::SRSliceProxy* slc) -> float {
    double showerLen = -9999.;
    int idx =  kRecoProtonIdx(slc);
    if ( idx >= 0 ) {
      auto const &shw = slc->reco.pfp.at(idx).shw;
      if ( !isnan(shw.len) ) showerLen = shw.len;
    }
    return showerLen;
  });

const Var kScndProtonShowerLength([](const caf::SRSliceProxy* slc) -> float {
    double showerLen = -9999.;
    int idx =  kScndProtonIdx(slc);
    if ( idx >= 0 ) {
      auto const &shw = slc->reco.pfp.at(idx).shw;
      if ( !isnan(shw.len) ) showerLen = shw.len;
    }
    return showerLen;
  });

const Var kRecoProtonShowerEnergy([](const caf::SRSliceProxy* slc) -> float {
    double showerE = -9999.;
    int idx =  kRecoProtonIdx(slc);
    if ( idx >= 0 ) {
      auto const &shw = slc->reco.pfp.at(idx).shw;
      if ( !isnan(shw.bestplane_energy) ) showerE = shw.bestplane_energy;
    }
    return showerE;
  });

const Var kScndProtonShowerEnergy([](const caf::SRSliceProxy* slc) -> float {
    double showerE = -9999.;
    int idx =  kScndProtonIdx(slc);
    if ( idx >= 0 ) {
      auto const &shw = slc->reco.pfp.at(idx).shw;
      if ( !isnan(shw.bestplane_energy) ) showerE = shw.bestplane_energy;
    }
    return showerE;
  });

//ExtraPrimaryIndex and length cuts moved up in file

const Var kExtraPrimaryShowerLength([](const caf::SRSliceProxy* slc) -> float {
  double length = -9999.;
  int idx  = kExtraPrimaryIndex(slc);

  if ( idx > -1 ) {
    const auto &pfp = slc->reco.pfp.at(idx);
    length = pfp.shw.len;
  }

  return length;
  });

const Var kExtraPrimaryShowerEnergy([](const caf::SRSliceProxy* slc) -> float {
  double energy = -9999.;
  int idx  = kExtraPrimaryIndex(slc);

  if ( idx > -1 ) {
    const auto &pfp = slc->reco.pfp.at(idx);
    energy = pfp.shw.bestplane_energy;
  }

  return energy;
  });

const Var kExtraPrimaryShowerdEdx([](const caf::SRSliceProxy* slc) -> float {
  double dEdx = -9999.;
  int idx  = kExtraPrimaryIndex(slc);

  if ( idx > -1 ) {
    const auto &pfp = slc->reco.pfp.at(idx);
    dEdx = pfp.shw.bestplane_dEdx;
  }

  return dEdx;
  });

const Var kExtraPrimaryShowerConvGap([](const caf::SRSliceProxy* slc) -> float {
  double conv_gap = -9999.;
  int idx  = kExtraPrimaryIndex(slc);

  if ( idx > -1 ) {
    const auto &pfp = slc->reco.pfp.at(idx);
    conv_gap = pfp.shw.conversion_gap;
  }

  return conv_gap;
  });

const Var kExtraPrimaryShowerDensity([](const caf::SRSliceProxy* slc) -> float {
  double density = -9999.;
  int idx  = kExtraPrimaryIndex(slc);

  if ( idx > -1 ) {
    const auto &pfp = slc->reco.pfp.at(idx);
    density = pfp.shw.density;
  }

  return density;
  });

const Var kExtraPrimaryShowerOpenAngle([](const caf::SRSliceProxy* slc) -> float {
  double open_angle = -9999.;
  int idx  = kExtraPrimaryIndex(slc);

  if ( idx > -1 ) {
    const auto &pfp = slc->reco.pfp.at(idx);
    open_angle = pfp.shw.open_angle;
  }

  return open_angle;
  });

const Var kExtraPrimaryChi2Muon([](const caf::SRSliceProxy* slc) -> float {
  double chi2 = -9999.;
  int idx  = kExtraPrimaryIndex(slc);

  if ( idx > -1 ) {
    const auto &pfp = slc->reco.pfp.at(idx);
    if ( !isnan(pfp.trk.chi2pid[2].chi2_muon) ) chi2 = pfp.trk.chi2pid[2].chi2_muon;
  }

  return chi2;
  });

const Var kExtraPrimaryChi2Proton([](const caf::SRSliceProxy* slc) -> float {
  double chi2 = -9999.;
  int idx  = kExtraPrimaryIndex(slc);

  if ( idx > -1 ) {
    const auto &pfp = slc->reco.pfp.at(idx);
    if ( !isnan(pfp.trk.chi2pid[2].chi2_proton) )chi2 = pfp.trk.chi2pid[2].chi2_proton;
  }

  return chi2;
  });

const Var kExtraPrimaryTruePDG([](const caf::SRSliceProxy* slc) -> float {
  double pdg = -9999.;
  int idx  = kExtraPrimaryIndex(slc);

  if ( idx > -1 ) {
    const auto &pfp = slc->reco.pfp.at(idx);
    pdg = pfp.trk.truth.p.pdg;
  }

  return pdg;
  });


const Var kExtraPrimaryStartProcess([](const caf::SRSliceProxy* slc) -> float {
  double proc = -9999.;
  int idx  = kExtraPrimaryIndex(slc);

  if ( idx > -1 ) {
    const auto &pfp = slc->reco.pfp.at(idx);
    proc = pfp.trk.truth.p.start_process;
  }

  return proc;
  });

const Var kSplitMuonDistance([](const caf::SRSliceProxy* slc) -> float {
  double dist = -9999.;

  if ( abs(kExtraPrimaryTruePDG(slc)) == 11 || abs(kExtraPrimaryTruePDG(slc)) == 13 ) {
    const auto &muTrk = slc->reco.pfp.at(kRecoMuonIdx(slc)).trk;
    const auto &pfp = slc->reco.pfp.at(kExtraPrimaryIndex(slc));

    double trkStartDiff = 9999.;
    double shwStartDiff = 9999.;
    if ( !isnan(pfp.trk.start.x) && !isnan(pfp.trk.start.y) && !isnan(pfp.trk.start.z) )
      trkStartDiff = std::hypot( (pfp.trk.start.x - muTrk.end.x), (pfp.trk.start.y - muTrk.end.y), (pfp.trk.start.z - muTrk.end.z) );
    if ( !isnan(pfp.shw.start.x) && !isnan(pfp.shw.start.y) && !isnan(pfp.shw.start.z) )
      shwStartDiff = std::hypot( (pfp.shw.start.x - muTrk.end.x), (pfp.shw.start.y - muTrk.end.y), (pfp.shw.start.z - muTrk.end.z) );
    dist = std::min(trkStartDiff, shwStartDiff);
  }

  return dist;
  });

const Var kTrueSplitMuon([](const caf::SRSliceProxy* slc) -> int {
  bool correctMuon = false;
  int idx = kExtraPrimaryIndex(slc);
  if (idx < 0 ) return correctMuon;

  int muID = -9999;
  for (const auto &prim : slc->truth.prim) {
    if ( abs(prim.pdg) == 13 ) {
      muID = prim.G4ID;
      break;
    }
  }
  int pdg = kExtraPrimaryTruePDG(slc);
  int g4ID = slc->reco.pfp.at(idx).trk.truth.p.G4ID;
  int parentID = slc->reco.pfp.at(idx).trk.truth.p.parent;

  correctMuon = ( (abs(pdg) == 11 || abs(pdg) == 13) && (g4ID == muID || parentID == muID) );
  return correctMuon;
  });

/*
const Var kExtraPrimaryLinFitLength([](const caf::SRSliceProxy* slc) -> float {
  double longest = -9999.;

  unsigned idxMuon = kRecoMuonIdx(slc);
  const auto &muTrk = slc->reco.pfp.at(idxMuon).trk;
  std::vector<double> idxProtons = kRecoProtonIndices(slc);

  for ( unsigned i = 1; i < slc->reco.pfp.size(); i++ ) {
    const auto &pfp = slc->reco.pfp.at(i);
    if ( !pfp.parent_is_const auto &prim : slc->truth.primprimary || i == idxMuon || std::find(idxProtons.begin(), idxProtons.end(), i) != idxProtons.end() ) continue;

    double trkStartDiff = 9999.;
    double shwStartDiff = 9999.;
    if ( !isnan(pfp.trk.start.x) && !isnan(pfp.trk.start.y) && !isnan(pfp.trk.start.z) )
      trkStartDiff = std::hypot( (pfp.trk.start.x - muTrk.end.x), (pfp.trk.start.y - muTrk.end.y), (pfp.trk.start.z - muTrk.end.z) );
    if ( !isnan(pfp.shw.start.x) && !isnan(pfp.shw.start.y) && !isnan(pfp.shw.start.z) )
      shwStartDiff = std::hypot( (pfp.shw.start.x - muTrk.end.x), (pfp.shw.start.y - muTrk.end.y), (pfp.shw.start.z - muTrk.end.z) );
    if ( std::min(trkStartDiff, shwStartDiff) < 5. ) continue;

    if ( !isnan(pfp.pfochar.linfitlen) && pfp.pfochar.linfitlen > longest ) longest = pfp.pfochar.linfitlen;
  }

  return longest;
  });

const Cut kExtraPrimaryLinFitLengthCut([](const caf::SRSliceProxy* slc) {
  double length = kExtraPrimaryLinFitLength(slc);
  return (length < 15.);
  });

const Var kSplitMuonDistance([](const caf::SRSliceProxy* slc) -> int {
  double longest = -9999.;

  unsigned idxMuon = kRecoMuonIdx(slc);
  const auto &muTrk = slc->reco.pfp.at(idxMuon).trk;
  std::vector<double> idxProtons = kRecoProtonIndices(slc);

  for ( unsigned i = 1; i < slc->reco.pfp.size(); i++ ) {
    const auto &pfp = slc->reco.pfp.at(i);
    if ( !pfp.parent_is_primary || i == idxMuon || std::find(idxProtons.begin(), idxProtons.end(), i) != idxProtons.end() ) continue;

    double trkStartDiff = 9999.;
    double shwStartDiff = 9999.;
    if ( !isnan(pfp.trk.start.x) && !isnan(pfp.trk.start.y) && !isnan(pfp.trk.start.z) )
      trkStartDiff = std::hypot( (pfp.trk.start.x - muTrk.end.x), (pfp.trk.start.y - muTrk.end.y), (pfp.trk.start.z - muTrk.end.z) );
    if ( !isnan(pfp.shw.start.x) && !isnan(pfp.shw.start.y) && !isnan(pfp.shw.start.z) )
      shwStartDiff = std::hypot( (pfp.shw.start.x - muTrk.end.x), (pfp.shw.start.y - muTrk.end.y), (pfp.shw.start.z - muTrk.end.z) );
    if ( std::min(trkStartDiff, shwStartDiff) < 5. ) continue;

    if ( !isnan(pfp.pfochar.linfitlen) && pfp.pfochar.linfitlen > longest ) longest = pfp.pfochar.linfitlen;
  }

  return longest;
  });

const Var kDynamicLengthCut([](const caf::SRSliceProxy* slc) -> int {
  double longest = -9999.;

  unsigned idxMuon = kRecoMuonIdx(slc);
  const auto &muTrk = slc->reco.pfp.at(idxMuon).trk;
  std::vector<double> idxProtons = kRecoProtonIndices(slc);

  double leadingProton = slc->reco.pfp.at(kRecoProtonIdx(slc)).pfochar.linfitlen;

  for ( unsigned i = 1; i < slc->reco.pfp.size(); i++ ) {
    const auto &pfp = slc->reco.pfp.at(i);
    if ( !pfp.parent_is_primary || i == idxMuon || std::find(idxProtons.begin(), idxProtons.end(), i) != idxProtons.end() ) continue;

    double trkStartDiff = 9999.;
    double shwStartDiff = 9999.;
    if ( !isnan(pfp.trk.start.x) && !isnan(pfp.trk.start.y) && !isnan(pfp.trk.start.z) )
      trkStartDiff = std::hypot( (pfp.trk.start.x - muTrk.end.x), (pfp.trk.start.y - muTrk.end.y), (pfp.trk.start.z - muTrk.end.z) );
    if ( !isnan(pfp.shw.start.x) && !isnan(pfp.shw.start.y) && !isnan(pfp.shw.start.z) )
      shwStartDiff = std::hypot( (pfp.shw.start.x - muTrk.end.x), (pfp.shw.start.y - muTrk.end.y), (pfp.shw.start.z - muTrk.end.z) );
    if ( std::min(trkStartDiff, shwStartDiff) < 5. ) continue;

    if ( !isnan(pfp.pfochar.linfitlen) && pfp.pfochar.linfitlen > longest ) longest = pfp.pfochar.linfitlen;
  }

  return (longest < 15. && longest < leadingProton);
  });

const Var kExtraPrimaryShowerLength([](const caf::SRSliceProxy* slc) -> float {
    double longest = -9999.;

    std::vector<unsigned> trackedPrimaries;
    if ( kRecoMuonIdx(slc) != -1 ) trackedPrimaries.push_back(kRecoMuonIdx(slc));
    for ( const double i : kRecoProtonIndices(slc) ) trackedPrimaries.push_back(i);

    for ( unsigned i = 1; i < slc->reco.pfp.size(); i++ ) {
      const auto &pfp = slc->reco.pfp.at(i);
      if ( !pfp.parent_is_primary || std::find(trackedPrimaries.begin(), trackedPrimaries.end(), i) != trackedPrimaries.end() ) continue;
      if ( !isnan(pfp.shw.len) && pfp.shw.len > longest ) longest = pfp.shw.len;
    }

    return longest;
  });

const Var kExtraPrimaryShowerEnergy([](const caf::SRSliceProxy* slc) -> float {
    double longest = -9999.;
    double bestplane_energy = -9999.;

    std::vector<unsigned> trackedPrimaries;
    if ( kRecoMuonIdx(slc) != -1 ) trackedPrimaries.push_back(kRecoMuonIdx(slc));
    for ( const double i : kRecoProtonIndices(slc) ) trackedPrimaries.push_back(i);

    for ( unsigned i = 1; i < slc->reco.pfp.size(); i++ ) {
      const auto &pfp = slc->reco.pfp.at(i);
      if ( !pfp.parent_is_primary || std::find(trackedPrimaries.begin(), trackedPrimaries.end(), i) != trackedPrimaries.end() ) continue;
      if ( !isnan(pfp.shw.len) && pfp.shw.len > longest ) {
        longest = pfp.shw.len;
        bestplane_energy = pfp.shw.bestplane_energy;
      }
    }

    return bestplane_energy;
  });

const Var kExtraPrimaryShowerdEdx([](const caf::SRSliceProxy* slc) -> float {
    double longest = -9999.;
    double bestplane_dEdx = -9999.;

    std::vector<unsigned> trackedPrimaries;
    if ( kRecoMuonIdx(slc) != -1 ) trackedPrimaries.push_back(kRecoMuonIdx(slc));
    for ( const double i : kRecoProtonIndices(slc) ) trackedPrimaries.push_back(i);

    for ( unsigned i = 1; i < slc->reco.pfp.size(); i++ ) {
      const auto &pfp = slc->reco.pfp.at(i);
      if ( !pfp.parent_is_primary || std::find(trackedPrimaries.begin(), trackedPrimaries.end(), i) != trackedPrimaries.end() ) continue;
      if ( !isnan(pfp.shw.len) && pfp.shw.len > longest ) {
        longest = pfp.shw.len;
        bestplane_dEdx = pfp.shw.bestplane_dEdx;
      }
    }

    return bestplane_dEdx;
  });

const Var kExtraPrimaryShowerConvGap([](const caf::SRSliceProxy* slc) -> float {
    double longest = -9999.;
    double conversion_gap = -9999.;

    std::vector<unsigned> trackedPrimaries;
    if ( kRecoMuonIdx(slc) != -1 ) trackedPrimaries.push_back(kRecoMuonIdx(slc));
    for ( const double i : kRecoProtonIndices(slc) ) trackedPrimaries.push_back(i);

    for ( unsigned i = 1; i < slc->reco.pfp.size(); i++ ) {
      const auto &pfp = slc->reco.pfp.at(i);
      if ( !pfp.parent_is_primary || std::find(trackedPrimaries.begin(), trackedPrimaries.end(), i) != trackedPrimaries.end() ) continue;
      if ( !isnan(pfp.shw.len) && pfp.shw.len > longest ) {
        longest = pfp.shw.len;
        conversion_gap = pfp.shw.conversion_gap;
      }
    }

    return conversion_gap;
  });

const Var kExtraPrimaryShowerDensity([](const caf::SRSliceProxy* slc) -> float {
    double longest = -9999.;
    double density = -9999.;

    std::vector<unsigned> trackedPrimaries;
    if ( kRecoMuonIdx(slc) != -1 ) trackedPrimaries.push_back(kRecoMuonIdx(slc));
    for ( const double i : kRecoProtonIndices(slc) ) trackedPrimaries.push_back(i);

    for ( unsigned i = 1; i < slc->reco.pfp.size(); i++ ) {
      const auto &pfp = slc->reco.pfp.at(i);
      if ( !pfp.parent_is_primary || std::find(trackedPrimaries.begin(), trackedPrimaries.end(), i) != trackedPrimaries.end() ) continue;
      if ( !isnan(pfp.shw.len) && pfp.shw.len > longest ) {
        longest = pfp.shw.len;
        density = pfp.shw.density;
      }
    }

    return density;
  });

const Var kExtraPrimaryShowerOpenAngle([](const caf::SRSliceProxy* slc) -> float {
    double longest = -9999.;
    double open_angle = -9999.;

    std::vector<unsigned> trackedPrimaries;
    if ( kRecoMuonIdx(slc) != -1 ) trackedPrimaries.push_back(kRecoMuonIdx(slc));
    for ( const double i : kRecoProtonIndices(slc) ) trackedPrimaries.push_back(i);

    for ( unsigned i = 1; i < slc->reco.pfp.size(); i++ ) {
      const auto &pfp = slc->reco.pfp.at(i);
      if ( !pfp.parent_is_primary || std::find(trackedPrimaries.begin(), trackedPrimaries.end(), i) != trackedPrimaries.end() ) continue;
      if ( !isnan(pfp.shw.len) && pfp.shw.len > longest ) {
        longest = pfp.shw.len;
        open_angle = pfp.shw.open_angle;
      }
    }

    return open_angle;
  });

const Var kExtraPrimaryTruePDG([](const caf::SRSliceProxy* slc) -> float {
    double longest = -9999.;
    double pdg = -9999.;

    std::vector<unsigned> trackedPrimaries;
    if ( kRecoMuonIdx(slc) != -1 ) trackedPrimaries.push_back(kRecoMuonIdx(slc));
    for ( const double i : kRecoProtonIndices(slc) ) trackedPrimaries.push_back(i);

    for ( unsigned i = 1; i < slc->reco.pfp.size(); i++ ) {
      const auto &pfp = slc->reco.pfp.at(i);
      if ( !pfp.parent_is_primary || std::find(trackedPrimaries.begin(), trackedPrimaries.end(), i) != trackedPrimaries.end() ) continue;
      if ( !isnan(pfp.shw.len) && pfp.shw.len > longest ) {
        longest = pfp.shw.len;
        pdg = pfp.shw.truth.p.pdg;
      }
    }

    return pdg;
  });

const Var kExtraPrimaryChi2Proton([](const caf::SRSliceProxy* slc) -> float {
    double longest = -9999.;
    double chi2 = -9999.;

    std::vector<unsigned> trackedPrimaries;
    if ( kRecoMuonIdx(slc) != -1 ) trackedPrimaries.push_back(kRecoMuonIdx(slc));
    for ( const double i : kRecoProtonIndices(slc) ) trackedPrimaries.push_back(i);

    for ( unsigned i = 1; i < slc->reco.pfp.size(); i++ ) {
      const auto &pfp = slc->reco.pfp.at(i);
      if ( !pfp.parent_is_primary || std::find(trackedPrimaries.begin(), trackedPrimaries.end(), i) != trackedPrimaries.end() ) continue;
      if ( !isnan(pfp.shw.len) && pfp.shw.len > longest ) {
        longest = pfp.shw.len;
        if ( !isnan(pfp.trk.chi2pid[2].chi2_proton) ) chi2 = pfp.trk.chi2pid[2].chi2_proton;
      }
    }

    return chi2;
  });

const Var kExtraPrimaryChi2Muon([](const caf::SRSliceProxy* slc) -> float {
    double longest = -9999.;
    double chi2 = -9999.;

    std::vector<unsigned> trackedPrimaries;
    if ( kRecoMuonIdx(slc) != -1 ) trackedPrimaries.push_back(kRecoMuonIdx(slc));
    for ( const double i : kRecoProtonIndices(slc) ) trackedPrimaries.push_back(i);

    for ( unsigned i = 1; i < slc->reco.pfp.size(); i++ ) {
      const auto &pfp = slc->reco.pfp.at(i);
      if ( !pfp.parent_is_primary || std::find(trackedPrimaries.begin(), trackedPrimaries.end(), i) != trackedPrimaries.end() ) continue;
      if ( !isnan(pfp.shw.len) && pfp.shw.len > longest ) {
        longest = pfp.shw.len;
        if ( !isnan(pfp.trk.chi2pid[2].chi2_muon) ) chi2 = pfp.trk.chi2pid[2].chi2_muon;
      }
    }

    return chi2;
  });
*/

const Var kTrueNeutrinoEnergy([](const caf::SRSliceProxy* slc) -> float {
  if ( !isnan(slc->truth.E) ) return slc->truth.E;
  else return -9999.;
  });

const Cut kSignalQE([](const caf::SRSliceProxy* slc) {
  return ( k1mu2p0pi(slc) && slc->truth.genie_mode==caf::kQE );
  });

const Cut kSignalMEC([](const caf::SRSliceProxy* slc) {
  return ( k1mu2p0pi(slc) && slc->truth.genie_mode==caf::kMEC );
  });

const Cut kSignalRes([](const caf::SRSliceProxy* slc) {
  return ( k1mu2p0pi(slc) && slc->truth.genie_mode==caf::kRes );
  });

const Cut kSignalDIS([](const caf::SRSliceProxy* slc) {
  return ( k1mu2p0pi(slc) && slc->truth.genie_mode==caf::kDIS );
  });

//Nonsnse?
const Cut kSignalCoh([](const caf::SRSliceProxy* slc) {
  return ( k1mu2p(slc) && slc->truth.genie_mode==caf::kCoh );
  });

const Cut kCCBG = ( !k1mu2p && kIsCC );

const Var kCCBGPrintout([](const caf::SRSliceProxy* slc) -> float {
  if ( !kCCBG(slc) && kTrueProton(slc) ) return 0.;

  std::cout << std::endl << "================   Start Event";

  for ( auto const& prim : slc->truth.prim ) {
    int pdg = prim.pdg;
    double p = sqrt ( prim.startp.x*prim.startp.x + prim.startp.y*prim.startp.y + prim.startp.z*prim.startp.z );
    bool contained = prim.contained;
    std::cout << std::endl << "Primary particle: " << pdg << ", momentum (GeV): " <<  p <<  ", contained?: " << contained;
  }

  return 0.;
  });

const Var kPrintPrimaryPFPs([](const caf::SRSliceProxy* slc) -> float {
    bool selected = kScndProtonCandidate(slc);

    std::cout << std::endl << "~~~~~~~~~~~~~~~~~ G4 Primaries";
    for ( const auto &prim : slc->truth.prim ) {
      if ( prim.start_process != 0 ) continue;
      std::cout << std::endl << "G4ID: " << prim.G4ID << ", PDG: " << prim.pdg << ", momentum: " << sqrt( prim.startp.x*prim.startp.x + prim.startp.y*prim.startp.y + prim.startp.z*prim.startp.z );
    }

    std::cout << std::endl << "~~~~~~~~~~~~~~~~~ Pandora Primaries";
    for ( const auto &pfp : slc->reco.pfp ) {
      int pdg = -9999;
      int g4id = -9999;
      int nMatches = -9999;
      int pid_ndof = -9999;
      double trkLen = -9999.;
      double trkScore = -9999.;
      double shwEnergy = -9999.;
      double caloSum = -9999.;
      double chi2_proton = -9999.;
      double chi2_muon = -9999.;
      double protonP = -9999.;
//      if ( pfp.parent_is_primary && abs(pfp.trk.truth.p.pdg) != 12 && abs(pfp.trk.truth.p.pdg) != 14 ) {
        if ( !isnan(pfp.shw.truth.p.G4ID) ) g4id = pfp.shw.truth.p.G4ID;
        if ( !isnan(pfp.shw.truth.p.pdg) ) pdg = pfp.shw.truth.p.pdg;
        if ( !isnan(pfp.shw.truth.nmatches) ) nMatches = pfp.shw.truth.nmatches;
        if ( !isnan(pfp.trk.len) ) trkLen = pfp.trk.len;
        if ( !isnan(pfp.trackScore) ) trkScore = pfp.trackScore;
        if ( !isnan(pfp.trk.calo[2].charge) ) caloSum = pfp.trk.calo[2].charge;
        if ( !isnan(pfp.shw.bestplane_energy) ) shwEnergy = pfp.shw.bestplane_energy;
        if ( !isnan(pfp.trk.chi2pid[2].pid_ndof) ) pid_ndof = pfp.trk.chi2pid[2].pid_ndof;
        if ( !isnan(pfp.trk.chi2pid[2].chi2_proton) ) chi2_proton = pfp.trk.chi2pid[2].chi2_proton;
        if ( !isnan(pfp.trk.chi2pid[2].chi2_muon) ) chi2_muon = pfp.trk.chi2pid[2].chi2_muon;
        if ( !isnan(pfp.trk.rangeP.p_proton) ) protonP = pfp.trk.rangeP.p_proton;

        std::cout << std::endl << "Primary PFP?: " << pfp.parent_is_primary << ", Best match PDG: "<< pdg << ", Best match G4ID: " << g4id << ", nMatches: " << nMatches << ", track length: " << trkLen << ", track score: " << trkScore << ", collection charge: " << caloSum << ", shower energy: " << shwEnergy << ", pid_ndof: " << pid_ndof << ", chi2_p: " << chi2_proton << ", chi2_mu: " << chi2_muon;
        if ( pdg == 2212 ) std::cout << ", Proton range P: " << protonP;

//      }
    }
    std::cout << std::endl << "Selected slice?: " << selected;
    std::cout << std::endl;
    return 0.;
  });

const Var kCCBGReason([](const caf::SRSliceProxy* slc) -> int {

  if ( !kCCBG(slc) ) return -1;

  //Failure mode: nu_e
  if ( abs(slc->truth.pdg) == 12 ) return 0;

  //Fact finding...
  int nVisProtons = 0;
  int nContained = 0;

  bool nuMuCCinFV = ( isInFV_Vtx(slc) );
  bool signalMuon = false;
  bool badTruthInfo = false;
  std::vector<int> primG4ID;
  std::vector<int> primPDG;
  std::vector<double> primP;
  std::vector<bool> primContained;
  int thisG4ID;
  int thisPDG;
  double thisP;
  bool thisContained;
  for ( auto const& prim : slc->truth.prim ) {
    if ( prim.start_process != 0 ) continue;
    thisG4ID = prim.G4ID;
    thisPDG = prim.pdg;
    thisP = sqrt( prim.startp.x*prim.startp.x + prim.startp.y*prim.startp.y + prim.startp.z*prim.startp.z );
    if ( badTruthInfo == false && thisP > 1000. ) badTruthInfo = true;
    thisContained = prim.contained;
    primG4ID.push_back(thisG4ID);
    primPDG.push_back(thisPDG);
    primP.push_back(thisP);
    primContained.push_back(thisContained);
    if ( abs(thisPDG) == 13 && thisP > .226 ) signalMuon = true;
    if ( thisPDG == 2212 && thisP >= 0.35 ) {
      nVisProtons++;
      if ( thisContained ) nContained++;
    }
  }
  unsigned int idxMuon = (unsigned int) kRecoMuonIdx(slc);
  unsigned int idxProton = (unsigned int) kRecoProtonIdx(slc);
  unsigned int idxThird = (unsigned int) kScndProtonIdx(slc);
  const auto &muTrk = slc->reco.pfp.at(idxMuon).trk;
  const auto &proton1Trk = slc->reco.pfp.at(idxProton).trk;
  const auto &proton2Trk = slc->reco.pfp.at(idxThird).trk;

  bool truthMuon = ( abs(muTrk.truth.p.pdg) == 13 );
  bool truthProton1 = ( abs(proton1Trk.truth.p.pdg) == 2212 );
  bool truthProton2 = ( abs(proton2Trk.truth.p.pdg) == 2212 );
  double pProton1 = sqrt( proton1Trk.truth.p.startp.x*proton1Trk.truth.p.startp.x + proton1Trk.truth.p.startp.y*proton1Trk.truth.p.startp.y + proton1Trk.truth.p.startp.z*proton1Trk.truth.p.startp.z );
  double pProton2 = sqrt( proton2Trk.truth.p.startp.x*proton2Trk.truth.p.startp.x + proton2Trk.truth.p.startp.y*proton2Trk.truth.p.startp.y + proton2Trk.truth.p.startp.z*proton2Trk.truth.p.startp.z );
  int g4IDProton1 = proton1Trk.truth.p.G4ID;
  int g4IDProton2 = proton2Trk.truth.p.G4ID;
  bool isPrimaryProton1 = ( std::find(primP.begin(), primP.end(), pProton1) != primP.end() );
  bool isPrimaryProton2 = ( std::find(primP.begin(), primP.end(), pProton2) != primP.end() );


  //Failure mode: TeV particles?
  if ( badTruthInfo ) return 1;

  //Failure mode: misID leading proton
  if ( !truthProton1 ) return 2;

  //Failure mode: misID sub-leading proton
  if ( !truthProton2 ) return 3;

  //Failure mode: best match tracked proton below threshold
  if ( (pProton1 < .35) || (pProton2 < .35) ) return 4;

  //Failure mode: both protons matched to same G4ID
  if ( g4IDProton1 == g4IDProton2 ) return 5;

/*
  int nPrimaries = primG4ID.size();
  std::cout << std::endl << "================   Start Event";
  for ( int i = 0; i < nPrimaries; i++  ) std::cout << std::endl << "Primary particle: " << primPDG.at(i) << ", momentum (GeV): " << primP.at(i) <<  ", contained?: " << primContained.at(i) << ", G4ID: " << primG4ID.at(i);
  std::cout << std::endl << "Muon track momentum (GeV): " << pMuon << ", Best match G4ID: " << g4IDMuon;
  std::cout << std::endl << "Leading proton track momentum (GeV): " <<  pProton1 <<  ", Best match G4ID: " << g4IDProton1;
  for ( const auto &match : proton1Trk.truth.matches ) {
    if ( match.hit_completeness != 0 ) std::cout << std::endl << "Truth match G4ID: " << match.G4ID << ", Hit completeness: " << match.hit_completeness << ", Hit purity: " << match.hit_purity;
  }
  std::cout << std::endl << "Sub-leading proton track momentum (GeV): " << pProton2 <<  ", Best match G4ID: " << g4IDProton2;
  for ( const auto &match : proton2Trk.truth.matches ) {
    if ( match.hit_completeness != 0 ) std::cout << std::endl << "Truth match G4ID: " << match.G4ID << ", Hit completeness: " << match.hit_completeness << ", Hit purity: " << match.hit_purity; 
  }
*/

  //Failure mode: >2 protons above threshold
  if ( nVisProtons > 2 ) return 6;

  //Failure mode: true proton is not a G4 primary
  if ( !isPrimaryProton1 || !isPrimaryProton2 ) return 7;

  //Failure mode: just one true proton above threshold other?
  if ( nVisProtons == 1 || !nuMuCCinFV ) return 8;

  //Failure mode: misID muon
  if ( !truthMuon ) return 9;

  //Failure mode: truth uncontained particles
  if ( !nuMuCCinFV || !signalMuon ) return 10;

  //Not sure what else is wrong then...
  return  11;
  });

const Cut kCat0([](const caf::SRSliceProxy* slc) {
  if ( !kCCBG(slc) ) return false;
  int category = kCCBGReason(slc);
  return (category == 0);
  });

const Cut kCat1([](const caf::SRSliceProxy* slc) {
  if ( !kCCBG(slc) ) return false;
  int category = kCCBGReason(slc);
  return (category == 1);
  });

const Cut kCat2([](const caf::SRSliceProxy* slc) {
  if ( !kCCBG(slc) ) return false;
  int category = kCCBGReason(slc);
  return (category == 2);
  });

const Cut kCat3([](const caf::SRSliceProxy* slc) {
  if ( !kCCBG(slc) ) return false;
  int category = kCCBGReason(slc);
  return (category == 3);
  });

const Cut kCat4([](const caf::SRSliceProxy* slc) {
  if ( !kCCBG(slc) ) return false;
  int category = kCCBGReason(slc);
  return (category == 4);
  });

const Cut kCat5([](const caf::SRSliceProxy* slc) {
  if ( !kCCBG(slc) ) return false;
  int category = kCCBGReason(slc);
  return (category == 5);
  });

const Cut kCat6([](const caf::SRSliceProxy* slc) {
  if ( !kCCBG(slc) ) return false;
  int category = kCCBGReason(slc);
  return (category == 6);
  });

const Cut kCat7([](const caf::SRSliceProxy* slc) {
  if ( !kCCBG(slc) ) return false;
  int category = kCCBGReason(slc);
  return (category == 7);
  });

const Cut kCat8([](const caf::SRSliceProxy* slc) {
  if ( !kCCBG(slc) ) return false;
  int category = kCCBGReason(slc);
  return (category == 8);
  });

const Cut kCat9([](const caf::SRSliceProxy* slc) {
  if ( !kCCBG(slc) ) return false;
  int category = kCCBGReason(slc);
  return (category == 9);
  });

const Cut kCat10([](const caf::SRSliceProxy* slc) {
  if ( !kCCBG(slc) ) return false;
  int category = kCCBGReason(slc);
  return (category == 10);
  });

const Cut kCat11([](const caf::SRSliceProxy* slc) {
  if ( !kCCBG(slc) ) return false;
  int category = kCCBGReason(slc);
  return (category == 11);
  });

/*
const Cut kCCOther([](const caf::SRSliceProxy* slc) {
  if ( !kCCBG(slc) ) return false;
  int category = kCCBGReason(slc);
  return (category != 6);
  });
*/


const Cut kOOFV = ( !k1mu0p && !k1mu1p && !k1mu2p && !k1mu3p && kIsCC && !isInFV_Vtx);

const Cut kCCOther = ( !k1mu0p && !k1mu1p && !k1mu2p && !k1mu3p && kIsCC && isInFV_Vtx);

const Var kPrintCCBG([](const caf::SRSliceProxy* slc) -> float {
//    int category = kCCBGReason(slc);
//    if ( category == 6 ) return 0.;

    unsigned int muTrk = kRecoMuonIdx(slc);
    unsigned int pTrk1 = kRecoProtonIdx(slc);
    unsigned int pTrk2 = kScndProtonIdx(slc);


//    std::cout << std::endl << "~~~~~~~~~~~~~~~~~ CC Nonsignal Category: " << category;

    std::cout << std::endl << "~~~~~~~~~~~~~~~~~ True Primaries";
    for ( const auto &prim : slc->truth.prim ) {
      if ( prim.start_process != 0 ) continue;
      std::cout << std::endl << "G4ID: " << prim.G4ID << ", PDG: " << prim.pdg << ", momentum: " << std::hypot( prim.startp.x, prim.startp.y, prim.startp.z ) << ", end process: " << prim.end_process;
      if ( std::hypot( prim.startp.x, prim.startp.y, prim.startp.z ) > 1000. || std::hypot( prim.startp.x, prim.startp.y, prim.startp.z ) == 0. )
          std::cout << ", momentum vector: (" << prim.startp.x << ", " << prim.startp.y << ", " << prim.startp.z << "), energy: " << prim.startE << ", gen position: (" << 
          prim.gen.x << ", " << prim.gen.y << ", " << prim.gen.z << "), start position: (" << prim.start.x << ", " << prim.start.y << ", " << prim.start.z << ")";
    }
    std::cout << std::endl << "~~~~~~~~~~~~~~~~~ True Secondaries";
    for ( const auto &prim : slc->truth.prim ) {
      if ( prim.start_process == 0 ) continue;
      std::cout << std::endl << "G4ID: " << prim.G4ID << ", PDG: " << prim.pdg << ", momentum: " << std::hypot( prim.startp.x, prim.startp.y, prim.startp.z ) << ", parent ID: " << 
          prim.parent << ", start process: " << prim.start_process << ", end process: " << prim.end_process;
    }

    std::cout << std::endl << "~~~~~~~~~~~~~~~~~ Pandora PFPs";
    unsigned int nPFPs = slc->reco.npfp;
    for ( unsigned int i = 0; i < nPFPs; i++ ) {
      const auto &pfp = slc->reco.pfp.at(i);
      int pfpid = -9999;
      int pdg = -9999;
      int g4id = -9999;
      int nMatches = -9999;
      int parent = -9999;
      double trkLen = -9999.;
      double trkScore = -9999.;
      double shwEnergy = -9999.;
      double caloSum = -9999.;
      double chi2_proton = -9999.;
      double chi2_muon = -9999.;
        if ( !isnan(pfp.id) ) pfpid = pfp.id;
        if ( !isnan(pfp.shw.truth.p.G4ID) ) g4id = pfp.shw.truth.p.G4ID;
        if ( !isnan(pfp.shw.truth.p.pdg) ) pdg = pfp.shw.truth.p.pdg;
        if ( !isnan(pfp.shw.truth.nmatches) ) nMatches = pfp.shw.truth.nmatches;
        if ( !isnan(pfp.parent) ) parent = pfp.parent;
        if ( !isnan(pfp.trk.len) ) trkLen = pfp.trk.len;
        if ( !isnan(pfp.trackScore) ) trkScore = pfp.trackScore;
        if ( !isnan(pfp.trk.calo[2].charge) ) caloSum = pfp.trk.calo[2].charge;
        if ( !isnan(pfp.shw.bestplane_energy) ) shwEnergy = pfp.shw.bestplane_energy;
        if ( !isnan(pfp.trk.chi2pid[2].chi2_proton) ) chi2_proton = pfp.trk.chi2pid[2].chi2_proton;
        if ( !isnan(pfp.trk.chi2pid[2].chi2_muon) ) chi2_muon = pfp.trk.chi2pid[2].chi2_muon;

        std::cout << std::endl << "PFP index: " << i << ", PFPID: " << pfpid << ", Best match PDG: "<< pdg << ", Best match G4ID: " << g4id << ", nMatches: " << nMatches <<  ", parent PFPID: " << 
            parent <<  ", track length: " << trkLen << ", track score: " << trkScore << ", collection charge: " << caloSum << ", shower energy: " << shwEnergy << ", chi2_p: " << 
            chi2_proton << ", chi2_mu: " << chi2_muon;

    }
    std::cout << std::endl << "Muon index: " << muTrk << ", Leading proton index: " << pTrk1 << ", Sub-leading proton index: " << pTrk2;

    std::cout << std::endl << "Reco vertex: (" << slc->vertex.x << ", " << slc->vertex.y << ", " << slc->vertex.z << "), True vertex: (" << 
        slc->truth.position.x << ", " << slc->truth.position.y << ", " << slc->truth.position.z << "), Difference: " << 
        std::hypot(slc->truth.position.x-slc->vertex.x, slc->truth.position.y-slc->vertex.y, slc->truth.position.z-slc->vertex.z) << ", In FV?: " <<  isInFV_Vtx(slc);
        if ( !isInFV_Vtx(slc) ) { std::cout << ", In AV?: " <<  isInAV_Vtx(slc); }

//    if ( category == 1 ) std::cout << std::endl << "Neutrino energy: " << slc->truth.E << ", GENIE mode: " <<  slc->truth.genie_mode << ", Nuclear target: " << slc->truth.targetPDG;

    std::cout << std::endl;
    return 0.;
  });

/*
const SpillVar kPrintDoubleSlice([](const caf::SRSPillProxy* spill) -> float {
    std::cout << std::endl << "~~~~~~~~~~~~~~~~~ True Primaries";
    for ( const auto &prim : spill->truth.nu.prim ) {
      if ( prim.start_process != 0 ) continue;
      std::cout << std::endl << "G4ID: " << prim.G4ID << ", PDG: " << prim.pdg << ", momentum: " << std::hypot( prim.startp.x, prim.startp.y, prim.startp.z ) << ", end process: " << prim.end_process;
      if ( std::hypot( prim.startp.x, prim.startp.y, prim.startp.z ) > 1000. || std::hypot( prim.startp.x, prim.startp.y, prim.startp.z ) == 0. )
          std::cout << ", momentum vector: (" << prim.startp.x << ", " << prim.startp.y << ", " << prim.startp.z << "), energy: " << prim.startE << ", gen position: (" <<
          prim.gen.x << ", " << prim.gen.y << ", " << prim.gen.z << "), start position: (" << prim.start.x << ", " << prim.start.y << ", " << prim.start.z << ")";
    }
    std::cout << std::endl << "~~~~~~~~~~~~~~~~~ True Secondaries";
    for ( const auto &prim : spill->truth.nu.prim ) {
      if ( prim.start_process == 0 ) continue;
      std::cout << std::endl << "G4ID: " << prim.G4ID << ", PDG: " << prim.pdg << ", momentum: " << std::hypot( prim.startp.x, prim.startp.y, prim.startp.z ) << ", parent ID: " <<
          prim.parent << ", start process: " << prim.start_process << ", end process: " << prim.end_process;
    }

    

 
  });
*/

const Var kRecoMuonIsPrimary([](const caf::SRSliceProxy* slc) -> int {
  bool isPrimary = false;

  int g4ID = slc->reco.pfp.at(kRecoMuonIdx(slc)).trk.truth.p.G4ID;
  
  for ( const auto &prim : slc->truth.prim ) {
    if ( prim.start_process != 0 ) continue;
    if ( prim.G4ID == g4ID ) { isPrimary = true; break; }
  }

  return isPrimary;
  });

const Var kRecoProtonIsPrimary([](const caf::SRSliceProxy* slc) -> int {
  bool isPrimary = false;

  int g4ID = slc->reco.pfp.at(kRecoProtonIdx(slc)).trk.truth.p.G4ID;

  for ( const auto &prim : slc->truth.prim ) {
    if ( prim.start_process != 0 ) continue;
    if ( prim.G4ID == g4ID ) { isPrimary = true; break; }
  }

  return isPrimary;
  });

const Var kScndProtonIsPrimary([](const caf::SRSliceProxy* slc) -> int {
  bool isPrimary = false;
  if ( kScndProtonIdx(slc) < 0 ) return isPrimary;

  int g4ID = slc->reco.pfp.at(kScndProtonIdx(slc)).trk.truth.p.G4ID;
  
  for ( const auto &prim : slc->truth.prim ) {
    if ( prim.start_process != 0 ) continue;
    if ( prim.G4ID == g4ID ) { isPrimary = true; break; }
  }

  return isPrimary;
  });

const Var kCountTeVParticles([](const caf::SRSliceProxy* slc) -> int {
  int nTeV = 0;

  for ( auto const& prim : slc->truth.prim ) {
    if ( prim.start_process != 0 ) continue;
    double thisP = sqrt( prim.startp.x*prim.startp.x + prim.startp.y*prim.startp.y + prim.startp.z*prim.startp.z );
    if ( thisP > 1000. ) nTeV++;
  }

  return nTeV;
  });

const Var kPrintTeVParticles([](const caf::SRSliceProxy* slc) -> int {
  const auto &muTrk = slc->reco.pfp.at(kRecoMuonIdx(slc)).trk;
  const auto &pTrk1 = slc->reco.pfp.at(kRecoProtonIdx(slc)).trk;
  const auto &pTrk2 = slc->reco.pfp.at(kScndProtonIdx(slc)).trk;

  if ( std::hypot(muTrk.truth.p.startp.x, muTrk.truth.p.startp.y, muTrk.truth.p.startp.z) > 1000. )
      std::cout << std::endl << "Muon candidate with TeV momentum! Length: " << muTrk.len << ", Reco momentum: " << muTrk.rangeP.p_muon;
  if ( std::hypot(pTrk1.truth.p.startp.x, pTrk1.truth.p.startp.y, pTrk1.truth.p.startp.z) > 1000. )
      std::cout << std::endl << "Leading proton candidate with TeV momentum! Length: " << pTrk1.len << ", Reco momentum: " << pTrk1.rangeP.p_proton;
  if ( std::hypot(pTrk2.truth.p.startp.x, pTrk2.truth.p.startp.y, pTrk2.truth.p.startp.z) > 1000. )
      std::cout << std::endl << "Second candidate with TeV momentum! Length: " << pTrk2.len << ", Reco momentum: " << pTrk2.rangeP.p_proton;

  return 0;
  });

const Var kNeutronEnergy([](const caf::SRSliceProxy* slc) -> float {
  double sum = 0.;

  for ( const auto &prim : slc->truth.prim ) {
    if ( prim.start_process != 0 ) continue;
    if ( prim.startE == -9999. ) return -9999.;
    if ( prim.pdg == 2112 ) sum += prim.startE - mNeutron;
  }

  return sum;
  });

const Var kRecoMuonDistToVtx([](const caf::SRSliceProxy* slc) -> float {
  const auto &trk = slc->reco.pfp.at(kRecoProtonIdx(slc)).trk;
  const float Atslc = std::hypot(slc->vertex.x - trk.start.x,
                                 slc->vertex.y - trk.start.y,
                                 slc->vertex.z - trk.start.z);
  return Atslc;
  });

const Var kRecoProtonDistToVtx([](const caf::SRSliceProxy* slc) -> float {
  const auto &trk = slc->reco.pfp.at(kRecoProtonIdx(slc)).trk;
  const float Atslc = std::hypot(slc->vertex.x - trk.start.x,
                                 slc->vertex.y - trk.start.y,
                                 slc->vertex.z - trk.start.z);
  return Atslc;
  });

const Var kScndProtonDistToVtx([](const caf::SRSliceProxy* slc) -> float {
  if ( kScndProtonIdx(slc) < 0 ) return -9999.;
  const auto &trk = slc->reco.pfp.at(kScndProtonIdx(slc)).trk;
  const float Atslc = std::hypot(slc->vertex.x - trk.start.x,
                                 slc->vertex.y - trk.start.y,
                                 slc->vertex.z - trk.start.z);
  return Atslc;
  });

const MultiVar kScndProtondEdx([](const caf::SRSliceProxy* slc) {
  std::vector<double> dedx;
  const auto &trk = slc->reco.pfp.at(kScndProtonIdx(slc)).trk;
  for ( const auto &point : trk.calo[2].points ) dedx.push_back(point.dedx);
  return dedx;
  });

const MultiVar kScndProtonRR([](const caf::SRSliceProxy* slc) {
  std::vector<double> rr;
  const auto &trk = slc->reco.pfp.at(kScndProtonIdx(slc)).trk;
  for ( const auto &point : trk.calo[2].points ) rr.push_back(point.rr);
  return rr;
  });

const MultiVar kRecoProtondEdx([](const caf::SRSliceProxy* slc) {
  std::vector<double> dedx;
  const auto &trk = slc->reco.pfp.at(kRecoProtonIdx(slc)).trk;
  for ( const auto &point : trk.calo[2].points ) dedx.push_back(point.dedx);
  return dedx;
  });

const MultiVar kRecoProtonRR([](const caf::SRSliceProxy* slc) {
  std::vector<double> rr;
  const auto &trk = slc->reco.pfp.at(kRecoProtonIdx(slc)).trk;
  for ( const auto &point : trk.calo[2].points ) rr.push_back(point.rr);
  return rr;
  });

const MultiVar kBothProtonsdEdx([](const caf::SRSliceProxy* slc) {
  std::vector<double> dedx;
  const auto &trk1 = slc->reco.pfp.at(kRecoProtonIdx(slc)).trk;
  const auto &trk2 = slc->reco.pfp.at(kScndProtonIdx(slc)).trk;
  for ( const auto &point : trk1.calo[2].points ) dedx.push_back(point.dedx);
  for ( const auto &point : trk2.calo[2].points ) dedx.push_back(point.dedx);
  return dedx;
  });

const MultiVar kBothProtonsRR([](const caf::SRSliceProxy* slc) {
  std::vector<double> rr;
  const auto &trk1 = slc->reco.pfp.at(kRecoProtonIdx(slc)).trk;
  const auto &trk2 = slc->reco.pfp.at(kScndProtonIdx(slc)).trk;
  for ( const auto &point : trk1.calo[2].points ) rr.push_back(point.rr);
  for ( const auto &point : trk2.calo[2].points ) rr.push_back(point.rr);
  return rr;
  });

const Var kTrueNeutrinoDirection([](const caf::SRSliceProxy* slc) -> float {
    float costh = -9999.;

    if ( !kIsNuSlice(slc) ) return costh;

    TVector3 nuP(slc->truth.momentum.x, slc->truth.momentum.y, slc->truth.momentum.z);
    const auto& vtx = slc->vertex;
    TVector3 vec_vtx(vtx.x, vtx.y, vtx.z);
    TVector3 vec_numi_to_vtx = (dFromNuMI + vec_vtx).Unit();

    costh = nuP.Dot(vec_numi_to_vtx) / nuP.Mag();
    return costh;
  });

const Var kThirdProtonP([](const caf::SRSliceProxy* slc) -> float {
  double momentum = -9999.;
  int idxP3 = kThirdProton(slc);

  if ( idxP3 >= 0 ) momentum = slc->reco.pfp.at(idxP3).trk.rangeP.p_proton;

  return momentum;
  });


const Var kThirdProtonTruthP([](const caf::SRSliceProxy* slc) -> float {
  double momentum = -9999.;
  int idxP3 = kThirdProton(slc);

  if ( idxP3 >= 0 ) {
    const auto &vec = slc->reco.pfp.at(idxP3).trk.truth.p.startp;
    momentum = std::hypot(vec.x, vec.y, vec.z);
  }

  return momentum;
  });

const Var kThirdProtonPResid([](const caf::SRSliceProxy* slc) -> float {
  double trueP = kThirdProtonTruthP(slc);
  double recoP = kThirdProtonP(slc);

  return (recoP-trueP) / trueP;
  });

const Var kHadronicOpeningAngle([](const caf::SRSliceProxy* slc) -> float {
    float costh = -9999.;

    if ( kNRecoProtons(slc) < 2 ) return costh;

    unsigned int idxProton = (unsigned int) kRecoProtonIdx(slc); //Leading proton
    unsigned int idxThird = (unsigned int) kScndProtonIdx(slc); //Charged pion or sub-leading proton

    TVector3 dirP(slc->reco.pfp.at(idxProton).trk.dir.x, slc->reco.pfp.at(idxProton).trk.dir.y, slc->reco.pfp.at(idxProton).trk.dir.z);
    TVector3 dirThird(slc->reco.pfp.at(idxThird).trk.dir.x, slc->reco.pfp.at(idxThird).trk.dir.y, slc->reco.pfp.at(idxThird).trk.dir.z);

    costh = dirP.Dot(dirThird);
    return costh;
  });

const Var kHadronicOpeningAngleTruth([](const caf::SRSliceProxy* slc) -> float {
    float costh = -9999.;
    if ( slc->truth.index < 0 ) return costh;
    const auto &nu = slc->truth;

    TVector3 thisP;
    TVector3 p1 = {0., 0., 0.};
    TVector3 p2 = {0., 0., 0.};

    for ( const auto &prim : nu.prim ) {
      if ( prim.start_process != 0 ) continue;
      thisP = {prim.startp.x, prim.startp.y, prim.startp.z};
      if ( thisP.Mag() == std::hypot(-9999., -9999., -9999.) ) return costh;
      if ( prim.pdg == 2212 && thisP.Mag() > .35 ) {
        if ( thisP.Mag() > p1.Mag() ) {
          p2 = p1;
          p1 = thisP;
        }
        else if ( thisP.Mag() > p2.Mag() ) p2 = thisP;
      }
    }

    if ( p1.Mag() == 0. || p2.Mag() == 0. ) return costh; 
    costh = p1.Dot(p2) / (  p1.Mag() * p2.Mag() );

    return costh;
  });

const Var kHadronicOpeningAngleOldTruth([](const caf::SRSliceProxy* slc) -> float {
    float costh = -9999.;

    if ( kNRecoProtons(slc) < 2 ) return costh;

    unsigned int idxProton = (unsigned int) kRecoProtonIdx(slc); //Leading proton
    unsigned int idxThird = (unsigned int) kScndProtonIdx(slc); //Charged pion or sub-leading proton

    TVector3 dirP(slc->reco.pfp.at(idxProton).trk.truth.p.startp.x, slc->reco.pfp.at(idxProton).trk.truth.p.startp.y, slc->reco.pfp.at(idxProton).trk.truth.p.startp.z);
    TVector3 dirThird(slc->reco.pfp.at(idxThird).trk.truth.p.startp.x, slc->reco.pfp.at(idxThird).trk.truth.p.startp.y, slc->reco.pfp.at(idxThird).trk.truth.p.startp.z);

    costh = dirP.Dot(dirThird) / (dirP.Mag() * dirThird.Mag());
    return costh;
  });

const Var kHadronicOpeningAngleResid([](const caf::SRSliceProxy* slc) -> float {
    double recoAngle = kHadronicOpeningAngle(slc);
    double trueAngle = kHadronicOpeningAngleTruth(slc);

    return (recoAngle - trueAngle) / trueAngle;
  });

const Var kMuonHadronAngle([](const caf::SRSliceProxy* slc) -> float {
    float costh = -9999.;

    if ( kNRecoProtons(slc) < 2 ) return costh;

    unsigned int idxMuon = (unsigned int) kRecoMuonIdx(slc);
//    std::vector<double> idcsProton = kRecoProtonIndices(slc);
    std::vector<double> idcsProton = {kRecoProtonIdx(slc), kScndProtonIdx(slc)};

    double pMu_mag = kRecoMuonPNew(slc);
    TVector3 pMu(pMu_mag*slc->reco.pfp.at(idxMuon).trk.dir.x, pMu_mag*slc->reco.pfp.at(idxMuon).trk.dir.y, pMu_mag*slc->reco.pfp.at(idxMuon).trk.dir.z);

    TVector3 thisP;
    TVector3 pHad(0., 0., 0.);
    for ( const auto &i : idcsProton ) {
      const auto &trk = slc->reco.pfp.at(i).trk;
      thisP = {trk.dir.x, trk.dir.y, trk.dir.z};
      thisP *= trk.rangeP.p_proton;
      pHad += thisP;
    }

    costh = pHad.Dot(pMu) / pMu.Mag() / pHad.Mag();
    return costh;
  });

/*
For Minerba
const Var kMuonHadronAngle([](const caf::SRSliceProxy* slc) -> float {
    float costh = -9999.;

    if ( kNRecoProtons(slc) < 2 ) return costh;

    unsigned int idxMuon = (unsigned int) kRecoMuonIdx(slc);
    unsigned int idxProton1 = (unsigned int) kLeadingProton(slc);
    unsigned int idxProton2 = (unsigned int) kSecondProton(slc);

    //Get the momentum vector for each particle
    double pMu_mag = slc->reco.pfp.at(idxMuon).trk.rangeP.p_muon;
    TVector3 pMu(pMu_mag*slc->reco.pfp.at(idxMuon).trk.dir.x, pMu_mag*slc->reco.pfp.at(idxMuon).trk.dir.y, pMu_mag*slc->reco.pfp.at(idxMuon).trk.dir.z);

    double pP1_mag = slc->reco.pfp.at(idxProton1).trk.rangeP.p_muon;
    TVector3 pP1(pP1_mag*slc->reco.pfp.at(idxProton1).trk.dir.x, pP1_mag*slc->reco.pfp.at(idxProton1).trk.dir.y, pP1_mag*slc->reco.pfp.at(idxProton1).trk.dir.z);

    double pP2_mag = slc->reco.pfp.at(idxProton2).trk.rangeP.p_muon;
    TVector3 pP2(pP2_mag*slc->reco.pfp.at(idxProton2).trk.dir.x, pP2_mag*slc->reco.pfp.at(idxProton2).trk.dir.y, pP2_mag*slc->reco.pfp.at(idxProton2).trk.dir.z);

    //Reutrn cosine of the angle between the muon momentum vector and hadronic (proton pair) momentum vector
    TVector3 pHad = pP1 + pP2;
    costh = pHad.Dot(pMu) / pMu.Mag() / pHad.Mag();
    return costh;
  });
*/

const Var kMuonHadronAngleTruth([](const caf::SRSliceProxy* slc) -> float {
    float costh = -9999.;
    if ( slc->truth.index < 0 ) return costh;
    const auto &nu = slc->truth;

    TVector3 thisP;
    TVector3 p1 = {0., 0., 0.};
    TVector3 p2 = {0., 0., 0.};
    TVector3 pMu = {0., 0., 0.};

    for ( const auto &prim : nu.prim ) {
      if ( prim.start_process != 0 ) continue;
      thisP = {prim.startp.x, prim.startp.y, prim.startp.z};
      if ( thisP.Mag() == std::hypot(-9999., -9999., -9999.) ) return costh;
      if ( abs(prim.pdg) == 13 ) pMu = thisP;
      else if ( prim.pdg == 2212 && thisP.Mag() > .35 ) {
        if ( thisP.Mag() > p1.Mag() ) {
          p2 = p1;
          p1 = thisP;
        }
        else if ( thisP.Mag() > p2.Mag() ) p2 = thisP;
      }
    }

    if ( pMu.Mag() == 0. || p2.Mag() == 0. ) return costh;

    TVector3 pHad = p1 + p2;
    costh = pMu.Dot(pHad) / (  pMu.Mag() * pHad.Mag() );

    return costh;
  });

const Var kMuonHadronAngleFullTruth([](const caf::SRSliceProxy* slc) -> float {
    float costh = -9999.;
    if ( slc->truth.index < 0 ) return costh;
    const auto &nu = slc->truth;

    TVector3 p3;
    TVector3 pMu(0., 0., 0.);
    TVector3 pHad(0., 0., 0.);

    for ( unsigned int i = 0; i < nu.prim.size(); i++ ) {
      const auto &prim = nu.prim.at(i);
      if ( prim.start_process != 0 ) continue;
      p3 = {prim.startp.x, prim.startp.y, prim.startp.z};
      if ( p3.Mag() == std::hypot(-9999., -9999., -9999.) ) return costh;
      if ( abs(prim.pdg) == 13 ) pMu = p3;
      else if ( abs(prim.pdg) == 2212 && p3.Mag() > .35 ) pHad += p3;
      else if ( abs(prim.pdg) == 211 || abs(prim.pdg) == 321 ) pHad += p3;
    }

    if ( pMu.Mag() == 0. || pHad.Mag() == 0. ) return costh;
    costh = pMu.Dot(pHad) / (  pMu.Mag() * pHad.Mag() );

    return costh;
  });

const Var kMuonHadronAngleOldTruth([](const caf::SRSliceProxy* slc) -> float {
    float costh = -9999.;

    //if ( !kThirdRENAMEIFNEEDEDPrimaryContained(slc) ) return costh;

    unsigned int idxMuon = (unsigned int) kRecoMuonIdx(slc);
    unsigned int idxProton = (unsigned int) kRecoProtonIdx(slc); //Leading proton
    unsigned int idxThird = (unsigned int) kScndProtonIdx(slc); //Charged pion or sub-leading proton

    TVector3 pMu(slc->reco.pfp.at(idxMuon).trk.truth.p.startp.x, slc->reco.pfp.at(idxMuon).trk.truth.p.startp.y, slc->reco.pfp.at(idxMuon).trk.truth.p.startp.z);
    TVector3 pP(slc->reco.pfp.at(idxProton).trk.truth.p.startp.x, slc->reco.pfp.at(idxProton).trk.truth.p.startp.y, slc->reco.pfp.at(idxProton).trk.truth.p.startp.z);
    TVector3 pThird(slc->reco.pfp.at(idxThird).trk.truth.p.startp.x, slc->reco.pfp.at(idxThird).trk.truth.p.startp.y, slc->reco.pfp.at(idxThird).trk.truth.p.startp.z);
    TVector3 pHad = pP + pThird;

    costh = pHad.Dot(pMu) / pMu.Mag() / pHad.Mag();
    return costh;
  });

const Var kMuonHadronAngleResid([](const caf::SRSliceProxy* slc) -> float {
    double recoAngle = kMuonHadronAngle(slc);
    double trueAngle = kMuonHadronAngleTruth(slc);

    return (recoAngle - trueAngle) / trueAngle;
  });

const Var kMuonPionAngle([](const caf::SRSliceProxy* slc) -> float {
    float costh = -9999.;

    if ( !kHasSidebandPion(slc) ) return costh;

    unsigned int idxMuon = (unsigned int) kRecoMuonIdx(slc);
    unsigned int idxPion = (unsigned int) kSidebandPion(slc);

    TVector3 pMu(slc->reco.pfp.at(idxMuon).trk.dir.x, slc->reco.pfp.at(idxMuon).trk.dir.y, slc->reco.pfp.at(idxMuon).trk.dir.z);
    TVector3 pPi(slc->reco.pfp.at(idxPion).trk.dir.x, slc->reco.pfp.at(idxPion).trk.dir.y, slc->reco.pfp.at(idxPion).trk.dir.z);

    costh = pPi.Dot(pMu) / pMu.Mag() / pPi.Mag();
    return costh;
  });

const Var kMuonPionAngleTruth([](const caf::SRSliceProxy* slc) -> float {
    float costh = -9999.;
    if ( slc->truth.index < 0 ) return costh;
    const auto &nu = slc->truth;

    TVector3 thisP;
    TVector3 pPi = {0., 0., 0.};
    TVector3 pMu = {0., 0., 0.};

    for ( const auto &prim : nu.prim ) {
      if ( prim.start_process != 0 ) continue;
      thisP = {prim.startp.x, prim.startp.y, prim.startp.z};
      if ( thisP.Mag() == std::hypot(-9999., -9999., -9999.) ) return costh;
      if ( abs(prim.pdg) == 13 ) pMu = thisP;
      else if ( abs(prim.pdg) == 211 ) {
        if ( thisP.Mag() > pPi.Mag() )  pPi = thisP;
      }
    }

    if ( pMu.Mag() == 0. || pPi.Mag() == 0. ) return costh;
    costh = pMu.Dot(pPi) / (  pMu.Mag() * pPi.Mag() );

    return costh;
  });

const Var kSidebandPionThetaNuMI([](const caf::SRSliceProxy* slc) -> float {
    float theta = -9999.;

    if ( kSidebandPion(slc) >= 0 ) {
      auto const& trk = slc->reco.pfp.at(kSidebandPion(slc)).trk;
      TVector3 direction(trk.dir.x, trk.dir.y, trk.dir.z);
      const auto& vtx = slc->vertex;
      TVector3 vec_vtx(vtx.x, vtx.y, vtx.z);
      TVector3 vec_numi_to_vtx = (dFromNuMI + vec_vtx).Unit();
      theta = direction.Dot(vec_numi_to_vtx) / (direction.Mag() * vec_numi_to_vtx.Mag());
    }

    return theta;
  });

const Var kSidebandPionTruthThetaNuMI([](const caf::SRSliceProxy* slc) -> float {
    float costh = -9999.;
    if ( slc->truth.index < 0 ) return costh;
    const auto &nu = slc->truth;
    TVector3 NuDirection_Truth = (TVector3(nu.momentum.x, nu.momentum.y, nu.momentum.z)).Unit().Unit();

    TVector3 thisP;
    TVector3 pPi = {0., 0., 0.};

    for ( const auto &prim : nu.prim ) {
      if ( prim.start_process != 0 ) continue;
      thisP = {prim.startp.x, prim.startp.y, prim.startp.z};
      if ( thisP.Mag() == std::hypot(-9999., -9999., -9999.) ) return costh;
      if ( abs(prim.pdg) == 211 ) {
        if ( thisP.Mag() > pPi.Mag() )  pPi = thisP;
      }
    }

    if ( pPi.Mag() == 0. ) return costh;
    costh = NuDirection_Truth.Dot(pPi) / (  NuDirection_Truth.Mag() * pPi.Mag() );

    return costh;
  });

const Var kRecoMuonTruthThetaNuMI([](const caf::SRSliceProxy* slc) -> float {
    float costh = -9999.;

    if ( slc->truth.index < 0 ) return costh;
    const auto &nu = slc->truth;
    TVector3 NuDirection_Truth = (TVector3(nu.momentum.x, nu.momentum.y, nu.momentum.z)).Unit().Unit();

    if ( kRecoMuonIdx(slc) >= 0 ) {
      auto const& trk = slc->reco.pfp.at(kRecoMuonIdx(slc)).trk;
      TVector3 direction(trk.truth.p.startp.x, trk.truth.p.startp.y, trk.truth.p.startp.z);
      costh = direction.Dot(NuDirection_Truth) / (direction.Mag() * NuDirection_Truth.Mag());
    }

    return costh;
  });

const Var kRecoMuonThetaNuMIResid([](const caf::SRSliceProxy* slc) -> float {
    double recoAngle = kRecoMuonThetaNuMI(slc);
    double trueAngle = kRecoMuonTruthThetaNuMI(slc);

    return (recoAngle - trueAngle) / trueAngle;
  });

const Var kRecoMuonStartDotEnd([](const caf::SRSliceProxy* slc) -> float {
    float costheta = -9999.;

    if ( kRecoMuonIdx(slc) >= 0 ) {
      auto const& trk = slc->reco.pfp.at(kRecoProtonIdx(slc)).trk;
      TVector3 start(trk.dir.x, trk.dir.y, trk.dir.z);
      TVector3 end(trk.dir_end.x, trk.dir.y, trk.dir.z);
      start = start.Unit();
      end = end.Unit();
      costheta = start.Dot(end);
    }

    return costheta;
  });

const Var kRecoProtonThetaNuMI([](const caf::SRSliceProxy* slc) -> float {
    float theta = -9999.;

    if ( kRecoProtonIdx(slc) >= 0 ) {
      auto const& trk = slc->reco.pfp.at(kRecoProtonIdx(slc)).trk;
      TVector3 direction(trk.dir.x, trk.dir.y, trk.dir.z);
      const auto& vtx = slc->vertex;
      TVector3 vec_vtx(vtx.x, vtx.y, vtx.z);
      TVector3 vec_numi_to_vtx = (dFromNuMI + vec_vtx).Unit();
      theta = direction.Dot(vec_numi_to_vtx) / (direction.Mag() * vec_numi_to_vtx.Mag());
    }

    return theta;
  });

const Var kRecoProtonTruthThetaNuMI([](const caf::SRSliceProxy* slc) -> float {
    float costh = -9999.;

    if ( slc->truth.index < 0 ) return costh;
    const auto &nu = slc->truth;
    TVector3 NuDirection_Truth = (TVector3(nu.momentum.x, nu.momentum.y, nu.momentum.z)).Unit().Unit();

    if ( kRecoProtonIdx(slc) >= 0 ) {
      auto const& trk = slc->reco.pfp.at(kRecoProtonIdx(slc)).trk;
      TVector3 direction(trk.truth.p.startp.x, trk.truth.p.startp.y, trk.truth.p.startp.z);
      costh = direction.Dot(NuDirection_Truth) / (direction.Mag() * NuDirection_Truth.Mag());
    }

    return costh;
  });

const Var kRecoProtonThetaNuMIResid([](const caf::SRSliceProxy* slc) -> float {
    double recoAngle = kRecoProtonThetaNuMI(slc);
    double trueAngle = kRecoProtonTruthThetaNuMI(slc);

    return (recoAngle - trueAngle) / trueAngle;
  });

const Var kRecoProtonTruthDirDotReco([](const caf::SRSliceProxy* slc) -> float {
    float costheta = -9999.;

    if ( kRecoProtonIdx(slc) >= 0 ) {
      auto const& trk = slc->reco.pfp.at(kRecoProtonIdx(slc)).trk;
      TVector3 truth(trk.truth.p.startp.x, trk.truth.p.startp.y, trk.truth.p.startp.z);
      truth = truth.Unit();
      TVector3 reco(trk.dir.x, trk.dir.y, trk.dir.z);
      costheta = truth.Dot(reco);
    }

    if ( costheta < 0.8 ) return .801;
    else return costheta;
  });

const Var kRecoProtonTruthDirDotRecoWide([](const caf::SRSliceProxy* slc) -> float {
    float costheta = -9999.;

    if ( kRecoProtonIdx(slc) >= 0 ) {
      auto const& trk = slc->reco.pfp.at(kRecoProtonIdx(slc)).trk;
      TVector3 truth(trk.truth.p.startp.x, trk.truth.p.startp.y, trk.truth.p.startp.z);
      truth = truth.Unit();
      TVector3 reco(trk.dir.x, trk.dir.y, trk.dir.z);
      costheta = truth.Dot(reco);
    }

    return costheta;
  });

const Var kRecoProtonStartDotEnd([](const caf::SRSliceProxy* slc) -> float {
    float costheta = -9999.;

    if ( kRecoProtonIdx(slc) >= 0 ) {
      auto const& trk = slc->reco.pfp.at(kRecoProtonIdx(slc)).trk;
      TVector3 start(trk.dir.x, trk.dir.y, trk.dir.z);
      TVector3 end(trk.dir_end.x, trk.dir.y, trk.dir.z);
      start = start.Unit();
      end = end.Unit();
      costheta = start.Dot(end);
    }

    return costheta;
  });

const Var kScndProtonThetaNuMI([](const caf::SRSliceProxy* slc) -> float {
    float costh = -9999.;

    if ( kScndProtonIdx(slc) >= 0 ) {
      auto const& trk = slc->reco.pfp.at(kScndProtonIdx(slc)).trk;
      TVector3 direction(trk.dir.x, trk.dir.y, trk.dir.z);
      const auto& vtx = slc->vertex;
      TVector3 vec_vtx(vtx.x, vtx.y, vtx.z);
      TVector3 vec_numi_to_vtx = (dFromNuMI + vec_vtx).Unit();
      costh = direction.Dot(vec_numi_to_vtx) / direction.Mag() / vec_numi_to_vtx.Mag();
    }

    return costh;
  });

const Var kScndProtonTruthThetaNuMI([](const caf::SRSliceProxy* slc) -> float {
    float costh = -9999.;

    if ( slc->truth.index < 0 ) return costh;
    const auto &nu = slc->truth;
    TVector3 NuDirection_Truth = (TVector3(nu.momentum.x, nu.momentum.y, nu.momentum.z)).Unit().Unit();

    if ( kScndProtonIdx(slc) >= 0 ) {
      auto const& trk = slc->reco.pfp.at(kScndProtonIdx(slc)).trk;
      TVector3 direction(trk.truth.p.startp.x, trk.truth.p.startp.y, trk.truth.p.startp.z);
      costh = direction.Angle(NuDirection_Truth);
    }

    return costh;
  });
    
const Var kScndProtonThetaNuMIResid([](const caf::SRSliceProxy* slc) -> float {
    double recoAngle = kScndProtonThetaNuMI(slc);
    double trueAngle = kScndProtonTruthThetaNuMI(slc);

    return (recoAngle - trueAngle) / trueAngle;
  });

const Var kScndProtonTruthDirDotReco([](const caf::SRSliceProxy* slc) -> float {
    float costheta = -9999.;

    if ( kScndProtonIdx(slc) >= 0 ) {
      auto const& trk = slc->reco.pfp.at(kScndProtonIdx(slc)).trk;
      TVector3 truth(trk.truth.p.startp.x, trk.truth.p.startp.y, trk.truth.p.startp.z);
      truth = truth.Unit();
      TVector3 reco(trk.dir.x, trk.dir.y, trk.dir.z);
      costheta = truth.Dot(reco);
    }

    if ( costheta < 0.8 ) return .801;
    else return costheta;
  });

const Var kScndProtonTruthDirDotRecoWide([](const caf::SRSliceProxy* slc) -> float {
    float costheta = -9999.;

    if ( kScndProtonIdx(slc) >= 0 ) {
      auto const& trk = slc->reco.pfp.at(kScndProtonIdx(slc)).trk;
      TVector3 truth(trk.truth.p.startp.x, trk.truth.p.startp.y, trk.truth.p.startp.z);
      truth = truth.Unit();
      TVector3 reco(trk.dir.x, trk.dir.y, trk.dir.z);
      costheta = truth.Dot(reco);
    }

    return costheta;
  });

const Var kScndProtonStartDotEnd([](const caf::SRSliceProxy* slc) -> float {
    float costheta = -9999.;

    if ( kScndProtonIdx(slc) >= 0 ) {
      auto const& trk = slc->reco.pfp.at(kScndProtonIdx(slc)).trk;
      TVector3 start(trk.dir.x, trk.dir.y, trk.dir.z);
      TVector3 end(trk.dir_end.x, trk.dir.y, trk.dir.z);
      start = start.Unit();
      end = end.Unit();
      costheta = start.Dot(end);
    }

    return costheta;
  });

const Var kRecoProtonTruthP([](const caf::SRSliceProxy* slc) -> float {
    float momentum = -9999.;

    if ( kRecoProtonIdx(slc) >= 0 ) {
      auto const& trk = slc->reco.pfp.at(kRecoProtonIdx(slc)).trk;
      momentum = sqrt(std::pow( trk.truth.p.startp.x, 2 ) + std::pow( trk.truth.p.startp.y, 2 ) + std::pow( trk.truth.p.startp.z, 2 ));
    }

    return momentum;
  }); 

const Var kRecoProtonPResid([](const caf::SRSliceProxy* slc) -> float {
    double recoMomentum = kRecoProtonP(slc);
    double trueMomentum = kRecoProtonTruthP(slc);

    return (recoMomentum - trueMomentum) / trueMomentum;
  });

const Var kExtraProtonTruthP([](const caf::SRSliceProxy* slc) -> float {
  double max = -9999.;
  if ( !k1mu3p(slc) ) return max;

  std::vector<double> pMomenta;
  for ( auto const& prim : slc->truth.prim ) {
    if ( prim.start_process != 0 || prim.pdg != 2212 ) continue;
    double momentum = sqrt( prim.startp.x*prim.startp.x + prim.startp.y*prim.startp.y + prim.startp.z*prim.startp.z );
    if ( momentum >= 0.35 ) pMomenta.push_back(momentum);
  }

  double p1 = kRecoProtonTruthP(slc);
  double p2 = kScndProtonTrueP(slc);

  //Could add something to exlude the ~10% of cases where our tracked protons are technically scatters
  for ( const auto &p : pMomenta ) if ( !(p == p1 || p == p2) && p > max ) max = p;
  return max;
  });

const Var kExtraPionTruthP([](const caf::SRSliceProxy* slc) -> float {
  double max = -9999.;
//  if ( !(k1mu3pNpi(slc) || k1mu2pNpi(slc) || k1mu1pNpi(slc) || k1mu0pNpi(slc)) ) return max;

  for ( auto const& prim : slc->truth.prim ) {
    if ( prim.start_process != 0 || abs(prim.pdg) != 211 ) continue;
    double momentum = sqrt( prim.startp.x*prim.startp.x + prim.startp.y*prim.startp.y + prim.startp.z*prim.startp.z );
    if ( momentum > max ) max = momentum;
  }

  return max;
  });

const Var kExtraPionEndProcess([](const caf::SRSliceProxy* slc) -> float {
  double max = -9999.;
  double proc = -9999.;

  for ( auto const& prim : slc->truth.prim ) {
    if ( prim.start_process != 0 || abs(prim.pdg) != 211 ) continue;
    double momentum = sqrt( prim.startp.x*prim.startp.x + prim.startp.y*prim.startp.y + prim.startp.z*prim.startp.z );
    if ( momentum > max ) {
      max = momentum;
      proc = prim.end_process;
    }
  }

  return proc;
  });

const Var kSidebandPionStartProcess([](const caf::SRSliceProxy* slc) -> float {
  double proc = -9999.;
  int idxPion = kSidebandPion(slc);
  if ( idxPion > -1 ) proc = slc->reco.pfp.at(idxPion).trk.truth.p.start_process;
  return proc;
  });

const Var kSidebandPionEndProcess([](const caf::SRSliceProxy* slc) -> float {
  double proc = -9999.;
  int idxPion = kSidebandPion(slc);
  if ( idxPion > -1 ) proc = slc->reco.pfp.at(idxPion).trk.truth.p.end_process;
  return proc;
  });

const Var kSidebandPionParentPDG([](const caf::SRSliceProxy* slc) -> float {
    int parent = -5;
    int parentID = -5;
    int pionIdx =  kSidebandPion(slc);
    if ( pionIdx >= 0 ) {
      if ( !isnan(slc->reco.pfp.at(pionIdx).trk.truth.p.parent) ) parentID = slc->reco.pfp.at(pionIdx).trk.truth.p.parent;
      if ( parentID >= 0 ) {
        for ( const auto &prim : slc->truth.prim ) {
          if ( prim.G4ID == parentID ) parent = prim.pdg;
        }
      }
    }
    return parent;
  });

const Var kExtraPionTrueLength([](const caf::SRSliceProxy* slc) -> float {
  double max = -9999.;
  double length = -9999.;

  for ( auto const& prim : slc->truth.prim ) {
    if ( prim.start_process != 0 || abs(prim.pdg) != 211 ) continue;
    double momentum = sqrt( prim.startp.x*prim.startp.x + prim.startp.y*prim.startp.y + prim.startp.z*prim.startp.z );
    if ( momentum > max ) {
      max = momentum;
      length = prim.length;
    }
  }

  return length;
  });

const Var kExtraPi0TruthP([](const caf::SRSliceProxy* slc) -> float {
  double max = -9999.;
//  if ( !(k1mu3pNpi(slc) || k1mu2pNpi(slc) || k1mu1pNpi(slc) || k1mu0pNpi(slc)) ) return max;

  for ( auto const& prim : slc->truth.prim ) {
    if ( prim.start_process != 0 || abs(prim.pdg) != 111 ) continue;
    double momentum = sqrt( prim.startp.x*prim.startp.x + prim.startp.y*prim.startp.y + prim.startp.z*prim.startp.z );
    if ( momentum > max ) max = momentum;
  }

  return max;
  });

/*
const Var kSidebandProtonP([](const caf::SRSliceProxy* slc) -> float {
  double momentum = -9999.;
  int idxP3 = kSidebandProton(slc);

  if ( idxP3 >= 0 ) momentum = slc->reco.pfp.at(idxP3).trk.rangeP.p_proton;

  return momentum;
  });


const Var kSidebandProtonTruthP([](const caf::SRSliceProxy* slc) -> float {
  double momentum = -9999.;
  int idxP3 = kSidebandProton(slc);

  if ( idxP3 >= 0 ) {
    const auto &vec = slc->reco.pfp.at(idxP3).trk.truth.p.startp;
    momentum = std::hypot(vec.x, vec.y, vec.z);
  }

  return momentum;
  });

const Var kSidebandProtonPResid([](const caf::SRSliceProxy* slc) -> float {
  double trueP = kSidebandProtonTruthP(slc);
  double recoP = kSidebandProtonP(slc);

  return (recoP-trueP) / trueP;
  });
*/

const Var kSidebandPionP([](const caf::SRSliceProxy* slc) -> float {
  double momentum = -9999.;
  int idxP3 = kSidebandPion(slc);

  if ( idxP3 >= 0 ) momentum = slc->reco.pfp.at(idxP3).trk.rangeP.p_pion;

  return momentum;
  });


const Var kSidebandPionTruthP([](const caf::SRSliceProxy* slc) -> float {
  double momentum = -9999.;
  int idxP3 = kSidebandPion(slc);

  if ( idxP3 >= 0 ) {
    const auto &vec = slc->reco.pfp.at(idxP3).trk.truth.p.startp;
    momentum = std::hypot(vec.x, vec.y, vec.z);
  }

  return momentum;
  });

const Var kSidebandPionPResid([](const caf::SRSliceProxy* slc) -> float {
  double trueP = kSidebandPionTruthP(slc);
  double recoP = kSidebandPionP(slc);

  return (recoP-trueP) / trueP;
  });

const Var kSidebandPionTrackLength([](const caf::SRSliceProxy* slc) -> float {
  double length = -9999.;
  int idxP3 = kSidebandPion(slc);

  if ( idxP3 >= 0 ) length = slc->reco.pfp.at(idxP3).trk.len;

  return length;
  });

const Var kSidebandPionTruePDG([](const caf::SRSliceProxy* slc) -> float {
  int pdg = -9999;
  int idxP3 = kSidebandPion(slc);

  if ( idxP3 >= 0 ) pdg = slc->reco.pfp.at(idxP3).trk.truth.p.pdg;

  return pdg;
  });

const Var kLeadingProtonPFrac([](const caf::SRSliceProxy* slc) -> float {
    float frac = -9999.;
    double sum = 0.;
    double leading = -1.;

    if ( kNRecoProtons(slc) < 2 ) return frac;

//    std::vector<double> idcsProton = kRecoProtonIndices(slc);
    std::vector<double> idcsProton = {kRecoProtonIdx(slc), kScndProtonIdx(slc)};
    for ( const auto &i : idcsProton ) {
      const auto &trk = slc->reco.pfp.at(i).trk;
      sum += trk.rangeP.p_proton;
      if ( trk.rangeP.p_proton > leading ) leading  = trk.rangeP.p_proton;
    }

    frac = leading / sum;

    return frac;
  });

const Var kLeadingProtonPFracOldTruth([](const caf::SRSliceProxy* slc) -> float {
    float frac = -9999.;

    //if ( !kScndProtonCandidate(slc) ) return frac;

    const auto & leadingProton = slc->reco.pfp.at(kRecoProtonIdx(slc)).trk.truth.p;
    const auto & scndProton = slc->reco.pfp.at(kScndProtonIdx(slc)).trk.truth.p;
    if ( leadingProton.pdg != 2212 || scndProton.pdg != 2212 ) return frac;

    double leadingP = std::hypot(leadingProton.startp.x, leadingProton.startp.y, leadingProton.startp.z);
    double scndP = std::hypot(scndProton.startp.x, scndProton.startp.y, scndProton.startp.z);
    frac = std::max(leadingP, scndP) / (leadingP+scndP);
    return frac;
  });

const Var kLeadingProtonPFracTruth([](const caf::SRSliceProxy* slc) -> float {
    float frac = -9999.;

    double leading = -1.;
    double second = -1.;
    TVector3 p3;

    for ( unsigned int i = 0; i < slc->truth.prim.size(); i++ ) {
      const auto &prim = slc->truth.prim.at(i);
      if ( prim.start_process != 0 ) continue;
      p3 = {prim.startp.x, prim.startp.y, prim.startp.z};
      if ( p3.Mag() == std::hypot(-9999., -9999., -9999.) ) return frac;
      if ( abs(prim.pdg) == 2212 && p3.Mag() > .35 ) {
        if ( p3.Mag() > leading ) {
          second = leading;
          leading = p3.Mag();
        }
        else if ( p3.Mag() > second ) {
          second = p3.Mag();
        }
      }
    }
    if ( (leading+second) == 0. || second == -1. ) return frac;

    frac = leading / (leading+second);
    return frac;
  });

const Var kLeadingProtonPFracFullTruth([](const caf::SRSliceProxy* slc) -> float {
    float frac = -9999.;

    double sum = 0.;
    double leading = -1.;
    TVector3 p3;

    for ( unsigned int i = 0; i < slc->truth.prim.size(); i++ ) {
      const auto &prim = slc->truth.prim.at(i);
      if ( prim.start_process != 0 ) continue;
      p3 = {prim.startp.x, prim.startp.y, prim.startp.z};
      if ( p3.Mag() == std::hypot(-9999., -9999., -9999.) ) return frac;
      if ( abs(prim.pdg) == 2212 && p3.Mag() > .35 ) {
        sum += p3.Mag();
        if ( p3.Mag() > leading ) leading = p3.Mag();
      }
    }
    if ( sum == 0. || leading == -1. ) return frac;

    frac = leading / sum;
    return frac;
  });

const Var kRecoMuonTruthP([](const caf::SRSliceProxy* slc) -> float {
    float momentum = -9999.;

    if ( kRecoMuonIdx(slc) >= 0 ) {
      auto const& trk = slc->reco.pfp.at(kRecoMuonIdx(slc)).trk;
      momentum = sqrt(std::pow( trk.truth.p.startp.x, 2 ) + std::pow( trk.truth.p.startp.y, 2 ) + std::pow( trk.truth.p.startp.z, 2 ));
    }

    return momentum;
  });

const Var kRecoMuonPResid([](const caf::SRSliceProxy* slc) -> float {
    double recoMomentum = kRecoMuonPNew(slc);
    double trueMomentum = kRecoMuonTruthP(slc);

    return (recoMomentum - trueMomentum) / trueMomentum;
  });

const Var kRecoMuonPInverseResid([](const caf::SRSliceProxy* slc) -> float {
    double recoMomentum = 1. / kRecoMuonPNew(slc);
    double trueMomentum = 1. / kRecoMuonTruthP(slc);

    return (recoMomentum - trueMomentum) / trueMomentum;
  });

const Var kRecoMuonTruthDirDotReco([](const caf::SRSliceProxy* slc) -> float {
    float costheta = -9999.;

    if ( kRecoMuonIdx(slc) >= 0 ) {
      auto const& trk = slc->reco.pfp.at(kRecoMuonIdx(slc)).trk;
      TVector3 truth(trk.truth.p.startp.x, trk.truth.p.startp.y, trk.truth.p.startp.z);
      truth = truth.Unit();
      TVector3 reco(trk.dir.x, trk.dir.y, trk.dir.z);
      costheta = truth.Dot(reco);
    }

    if ( costheta < 0.8 ) return .801;
    else return costheta;
  });

const Var kRecoMuonTruthDirDotRecoWide([](const caf::SRSliceProxy* slc) -> float {
    float costheta = -9999.;

    if ( kRecoMuonIdx(slc) >= 0 ) {
      auto const& trk = slc->reco.pfp.at(kRecoMuonIdx(slc)).trk;
      TVector3 truth(trk.truth.p.startp.x, trk.truth.p.startp.y, trk.truth.p.startp.z);
      truth = truth.Unit();
      TVector3 reco(trk.dir.x, trk.dir.y, trk.dir.z);
      costheta = truth.Dot(reco);
    }

    return costheta;
  });

const Cut kBadMuon([](const caf::SRSliceProxy* slc) {
  double momentum = kRecoMuonPNew(slc);
  return( isnan(momentum) || momentum < 0 );
  });

const Var kPrintBadMuon([](const caf::SRSliceProxy* slc) -> float {
  double momentum = kRecoMuonPNew(slc);
  bool contained = kRecoMuonContained(slc);
  auto const& trk = slc->reco.pfp.at(kRecoMuonIdx(slc)).trk;

  std::cout << std::endl << "~~~~~~~~~~~~~~~~~~~Bad Muon!";
  std::cout << std::endl << "Var: " << momentum << ", Contained?: " << contained << ", Range P: " << trk.rangeP.p_muon << ", MCS P: " << trk.mcsP.fwdP_muon;

  std::cout << std::endl << "cos(theta): " << kRecoMuonThetaNuMI(slc) << ", Track End - Start: " << std::hypot((trk.start.x-trk.end.x), (trk.start.y-trk.end.y), (trk.start.z-trk.end.z)) << ", trk.len" << trk.len << ", Muon track true species: " << trk.truth.p.pdg;
//  if ( slc->truth.isnc ) std::cout << ", CC?/NC? : Neutral current";
//  else if ( slc->truth.iscc ) std::cout << ", CC?/NC? : Charged current, CC Other category: " << kCCBGReason(slc);
  std::cout << std::endl << "Reco track start: (" << trk.start.x << ", " << trk.start.y << ", " << trk.start.z << "), Reco track end: (" << trk.end.x << ", " << trk.end.y << ", " << trk.end.z << "), True track start: ("
    << trk.truth.p.start.x << ", " << trk.truth.p.start.y << ", " << trk.truth.p.start.z << "), True track end: (" << trk.truth.p.end.x << ", " << trk.truth.p.end.y << ", " << trk.truth.p.end.z << ")";
  std::cout << std::endl << "Reco vertex: (" << slc->vertex.x << ", " << slc->vertex.y << ", " << slc->vertex.z << "), True vertex: (" << slc->truth.position.x << ", " << slc->truth.position.y << ", " << slc->truth.position.z << ")";

/*
  auto const& pTrk1 = slc->reco.pfp.at(kRecoProtonIdx(slc)).trk;
  auto const& pTrk2 = slc->reco.pfp.at(kScndProtonIdx(slc)).trk;
  std::cout << std::endl << "Reco proton track starts: (" << pTrk1.start.x << ", "  << pTrk1.start.y << ", "  << pTrk1.start.z << "), (" << pTrk2.start.x << ", "  << pTrk2.start.y << ", "  << pTrk2.start.z << ")";
*/

  std::cout << std::endl;

  return 0;
  });

/* Come back to this!
const Var kTrueSecondary([](const caf::SRSliceProxy* slc) -> float {
    float dpTMag = -9999.;
    if ( !kThirdRENAMEIFNEEDEDPrimaryContained(slc) ) return dpTMag;

    unsigned int idxMuon = (unsigned int) kRecoMuonIdx(slc);
    unsigned int idxProton = (unsigned int) kRecoProtonIdx(slc); //Leading proton
    unsigned int idxThird = (unsigned int) kScndProtonIdx(slc); //Charged pion or sub-leading proton

  });
*/

const Var kEHad_Proton([](const caf::SRSliceProxy* slc) -> float {
    double eHad = 0.;

    std::vector<double> idcsProton = kRecoProtonIndices(slc);
    for ( const auto &i : idcsProton ) {
      const auto &trk = slc->reco.pfp.at(i).trk;
      double thisP = trk.rangeP.p_proton;
      eHad += std::hypot(thisP, mProton) - mProton;
    }

    return eHad;
  });

const Var kEHad_CheatingProtonPs([](const caf::SRSliceProxy* slc) -> float {
    double eHad = -9999.;
    if ( kScndProtonIdx(slc) < 0 ) return eHad;

    double pP_mag = kRecoProtonTruthP(slc);
    double pThird_mag = kScndProtonTrueP(slc);

    double leadingT = std::hypot(pP_mag, mProton) - mProton;
    double secondT = std::hypot(pThird_mag, mProton) - mProton;

    eHad = leadingT + secondT;
    return eHad;
  });

/*
const Var kEHad_ThreeP([](const caf::SRSliceProxy* slc) -> float {
    double eHad = -9999.;

    unsigned int idxP1 = (unsigned int) kRecoProtonIdx(slc);
    unsigned int idxP2 = (unsigned int) kScndProtonIdx(slc);
    unsigned int idxP3 = (unsigned int) kSidebandProton(slc);

    double pP1_mag = kRecoProtonP(slc);
    double pP2_mag = kScndProtonP(slc);
    double pP3_mag = kSidebandProtonP(slc);

    double leadingT = std::hypot(pP1_mag, mProton) - mProton;
    double secondT = std::hypot(pP2_mag, mProton) - mProton;
    double thirdT = std::hypot(pP3_mag, mProton) - mProton;

    eHad = leadingT + secondT + thirdT;
    return eHad;
  });
*/

const Var kEHad_Pion([](const caf::SRSliceProxy* slc) -> float {
    double eHad = -9999.;

    double pP_mag = kRecoProtonP(slc);
    double pThird_mag = kSidebandPionP(slc);

    double protonT = std::hypot(pP_mag, mProton) - mProton;
    double pionE = std::hypot(pThird_mag, mPion);

    eHad = protonT + pionE;
    return eHad;
  });

const Var kEHad([](const caf::SRSliceProxy* slc) -> float {
    double eHad = -9999.;
    if ( kScndProtonCandidate(slc) ) eHad = kEHad_Proton(slc);
    //else if ( kThreePSideband(slc) ) eHad = kEHad_ThreeP(slc);
    else if ( kPionSidebandBase(slc) ) eHad = kEHad_Pion(slc);
    assert(((void)"No valid selection for EHad value", eHad != -9999.));

    return eHad;
  });

const Var kEHad_Truth([](const caf::SRSliceProxy* slc) -> float {
    double eHad = 0.;

    for ( auto const& prim : slc->truth.prim ) {
      if ( prim.start_process != 0 ) continue;
      if ( prim.startE == -9999. )  return -9999.;
      if ( (abs(prim.pdg) >= 11 && abs(prim.pdg) <= 14) || prim.pdg == 22 || prim.pdg == 1000180400 ) continue;
      else if ( abs(prim.pdg) == 2212 ) eHad += prim.startE - mProton;
      else if ( abs(prim.pdg) == 2112 ) eHad += prim.startE - mNeutron;
      //Just subtract off proton mass for strange hadrons. It's close enough to right, and they're < 0.5% of the sample anyway
      else if ( abs(prim.pdg) > 3000 && abs(prim.pdg) < 4000) eHad += prim.startE - mProton;
      else eHad += prim.startE;
    }

    if ( eHad < 0. ) return -9999.;
    else return eHad;
  });

//Same as above, minus neutrons
const Var kEAvail_Truth([](const caf::SRSliceProxy* slc) -> float {
    double eAvail = 0.;

    for ( auto const& prim : slc->truth.prim ) {
      if ( prim.start_process != 0 ) continue;
      if ( prim.startE == -9999. )  return -9999.;
      if ( (abs(prim.pdg) >= 11 && abs(prim.pdg) <= 14) || prim.pdg == 22 || prim.pdg == 1000180400 || prim.pdg == 2112 ) continue;
      else if ( abs(prim.pdg) == 2212 ) {
        if ( std::hypot(prim.startp.x,prim.startp.y,prim.startp.z) > .35 ) eAvail += prim.startE - mProton;
        else continue;
      }
      else if ( abs(prim.pdg) > 3000 && abs(prim.pdg) < 4000) eAvail += prim.startE - mProton;
      else eAvail += prim.startE;
    }

    if ( eAvail < 0. ) return -9999.;
    else return eAvail;
  });

const Var kHasKaon([](const caf::SRSliceProxy* slc) -> int {
  bool hasKaon = false;
  for ( auto const& prim : slc->truth.prim ) {
    if ( prim.start_process != 0 ) continue;
    if ( abs(prim.pdg) == 321 ||  abs(prim.pdg) == 311 ) {
      hasKaon = true;
      break;
    }
  }

  return hasKaon;
  });

const Var kHasLambda([](const caf::SRSliceProxy* slc) -> int {
  bool hasLambda = false;

  for ( auto const& prim : slc->truth.prim ) {
    if ( prim.start_process != 0 ) continue;
    if ( abs(prim.pdg) == 3122 ) {
      hasLambda = true;
      break;
    }
  }

  return hasLambda;
  });

const Var kHasSigma([](const caf::SRSliceProxy* slc) -> int {
  bool hasSigma = false;

  for ( auto const& prim : slc->truth.prim ) {
    if ( prim.start_process != 0 ) continue;
    if ( abs(prim.pdg) == 3222 ||  abs(prim.pdg) == 3212||  abs(prim.pdg) == 3112 ) {
      hasSigma = true;
      break;
    }
  }

  return hasSigma;
  });

const Var kHasAr([](const caf::SRSliceProxy* slc) -> int {
  bool hasAr = false;

  for ( auto const& prim : slc->truth.prim ) {
    if ( prim.start_process != 0 ) continue;
    if ( abs(prim.pdg) == 1000180400 ) {
      hasAr = true;
      break;
    }
  }

  return hasAr;
  });

const Var kPrintExotics([](const caf::SRSliceProxy* slc) -> float {
  if ( !(kHasKaon(slc) || kHasLambda(slc) || kHasSigma(slc) || kHasAr(slc)) ) return 0.;

  int cat = -1;
  if ( k1mu2p0pi(slc) ) cat = 0;
  else if ( k1mu2pNpi(slc) ) cat = 1;
  else if ( k1mu3p0pi(slc) ) cat = 2;
  else if ( k1mu3pNpi(slc) ) cat = 3;
  else if ( k1mu1p0pi(slc) ) cat = 4;
  else if ( k1mu1pNpi(slc) ) cat = 5;
  else if ( k1mu0p0pi(slc) ) cat = 6;
  else if ( k1mu0pNpi(slc) ) cat = 7;
  else if ( kCCOther(slc) ) cat = 8;
  else if ( kIsNC(slc) ) cat = 9;
  else if ( kIsCosmic(slc) ) cat = 10;


  int reason = kCCBGReason(slc);
  std::cout << std::endl << "~~~~~~~~~~~~~~~~~ Event category: " << cat << ", reason for selection: " << reason;
  kPrintCCBG(slc);
  return 0.;
  });

const Var kQ2([](const caf::SRSliceProxy* slc) -> float {
    double Q2 = -9999.;
    double eHad = kEHad(slc);

    unsigned int idxMuon = (unsigned int) kRecoMuonIdx(slc); 
    double pMu_mag = kRecoMuonPNew(slc);
    TVector3 pMu(pMu_mag*slc->reco.pfp.at(idxMuon).trk.dir.x, pMu_mag*slc->reco.pfp.at(idxMuon).trk.dir.y, pMu_mag*slc->reco.pfp.at(idxMuon).trk.dir.z);

    double eMu = std::hypot(pMu_mag, mMuon);
    double eNu = eMu + eHad;

    const auto& vtx = slc->vertex;
    TVector3 vec_vtx(vtx.x, vtx.y, vtx.z);
    TVector3 vec_numi_to_vtx = (dFromNuMI + vec_vtx).Unit();

    Q2 = 2.*eNu*(eMu - pMu.Dot(vec_numi_to_vtx)) - mMuon*mMuon;
    return Q2;
  });

const Var kQ2_CheatingMuonP([](const caf::SRSliceProxy* slc) -> float {
    double Q2 = -9999.;
    double eHad = kEHad(slc);

    unsigned int idxMuon = (unsigned int) kRecoMuonIdx(slc);
    double pMu_mag = kRecoMuonTruthP(slc);
    TVector3 pMu(pMu_mag*slc->reco.pfp.at(idxMuon).trk.dir.x, pMu_mag*slc->reco.pfp.at(idxMuon).trk.dir.y, pMu_mag*slc->reco.pfp.at(idxMuon).trk.dir.z);

    double eMu = std::hypot(pMu_mag, mMuon);
    double eNu = eMu + eHad;

    const auto& vtx = slc->vertex;
    TVector3 vec_vtx(vtx.x, vtx.y, vtx.z);
    TVector3 vec_numi_to_vtx = (dFromNuMI + vec_vtx).Unit();

    Q2 = 2.*eNu*(eMu - pMu.Dot(vec_numi_to_vtx)) - mMuon*mMuon;
    return Q2;
  });

const Var kQ2_CheatingAll([](const caf::SRSliceProxy* slc) -> float {
    double Q2 = -9999.;
    double eHad = kEHad_CheatingProtonPs(slc);

    unsigned int idxMuon = (unsigned int) kRecoMuonIdx(slc);
    double pMu_mag = kRecoMuonTruthP(slc);
    TVector3 pMu(pMu_mag*slc->reco.pfp.at(idxMuon).trk.dir.x, pMu_mag*slc->reco.pfp.at(idxMuon).trk.dir.y, pMu_mag*slc->reco.pfp.at(idxMuon).trk.dir.z);

    double eMu = std::hypot(pMu_mag, mMuon);
    double eNu = eMu + eHad;

    const auto& vtx = slc->vertex;
    TVector3 vec_vtx(vtx.x, vtx.y, vtx.z);
    TVector3 vec_numi_to_vtx = (dFromNuMI + vec_vtx).Unit();

    Q2 = 2.*eNu*(eMu - pMu.Dot(vec_numi_to_vtx)) - mMuon*mMuon;
    return Q2;
  });

const Var kQ2_Truth([](const caf::SRSliceProxy* slc) -> float {
    double Q2 = -9999.;
    TVector3 NuDirection_Truth = (TVector3(slc->truth.momentum.x, slc->truth.momentum.y, slc->truth.momentum.z)).Unit().Unit();

    double eHad = kEHad_Truth(slc);
    if ( eHad == -9999.  || slc->truth.index < 0 ) return Q2;

    TVector3 pMu(0., 0., 0.);
    for ( unsigned int i = 0; i < slc->truth.prim.size(); i++ ) {
      const auto &prim = slc->truth.prim.at(i);
      if ( prim.start_process != 0 ) continue;
      if ( abs(prim.pdg) == 13 ) {
        pMu = {prim.startp.x, prim.startp.y, prim.startp.z};
        break;
      }
    }
    if ( pMu.Mag() == 0. || pMu.Mag() == std::hypot(-9999., -9999., -9999.) ) return Q2;

    double eMu = std::hypot(pMu.Mag(), mMuon);
    double eNu = eMu + eHad;

    Q2 = 2.*eNu*(eMu - pMu.Dot(NuDirection_Truth)) - mMuon*mMuon;
    return Q2;
  });

const Var kq3([](const caf::SRSliceProxy* slc) -> float {
    double q3 = -9999.;
    double eHad = kEHad(slc);
    q3 = std::sqrt(kQ2(slc) + eHad*eHad);
    return q3;
  });

const Var kq3_Truth([](const caf::SRSliceProxy* slc) -> float {
    double q3 = -9999.;
    double eHad = kEHad_Truth(slc);
    double Q2 = kQ2_Truth(slc);
    if ( Q2 == -9999.  || slc->truth.index < 0 ) return q3;
    q3 = std::sqrt(Q2 + eHad*eHad);
    return q3;
  });

const Var kRecoENu([](const caf::SRSliceProxy* slc) -> float {
    double eNu = -9999.;
    double eHad = kEHad(slc);

    unsigned int idxMuon = (unsigned int) kRecoMuonIdx(slc);
    double pMu_mag = kRecoMuonPNew(slc);
    TVector3 pMu(pMu_mag*slc->reco.pfp.at(idxMuon).trk.dir.x, pMu_mag*slc->reco.pfp.at(idxMuon).trk.dir.y, pMu_mag*slc->reco.pfp.at(idxMuon).trk.dir.z);

    double eMu = std::hypot(pMu_mag, mMuon);
    eNu = eMu + eHad;

    return eNu;
  });

const Var kW_Proton([](const caf::SRSliceProxy* slc) -> float {
    float W = -9999.;

    unsigned int idxProton = (unsigned int) kRecoProtonIdx(slc); //Leading proton
    unsigned int idxThird = (unsigned int) kScndProtonIdx(slc); //Charged pion or sub-leading proton

    double pP_mag = kRecoProtonP(slc);
    double pThird_mag = kScndProtonP(slc);
    double EP = sqrt(mProton*mProton + pP_mag*pP_mag);
    double EThird = sqrt(mProton*mProton + pThird_mag*pThird_mag);

    TVector3 pP(pP_mag*slc->reco.pfp.at(idxProton).trk.dir.x, pP_mag*slc->reco.pfp.at(idxProton).trk.dir.y, pP_mag*slc->reco.pfp.at(idxProton).trk.dir.z);
    TVector3 pThird(pThird_mag*slc->reco.pfp.at(idxThird).trk.dir.x, pThird_mag*slc->reco.pfp.at(idxThird).trk.dir.y, pThird_mag*slc->reco.pfp.at(idxThird).trk.dir.z);
    TVector3 pTotal = pP + pThird;

    W = sqrt( (EP + EThird)*(EP + EThird) - pTotal.Mag2() );

    return W;
  });

const Var kW_Pion([](const caf::SRSliceProxy* slc) -> float {
    float W = -9999.;

    unsigned int idxProton = (unsigned int) kRecoProtonIdx(slc); //Leading proton
    unsigned int idxThird = (unsigned int) kSidebandPion(slc); //Charged pion or sub-leading proton

    double pP_mag = kRecoProtonP(slc);
    double pThird_mag = kSidebandPionP(slc);
    double EP = sqrt(mProton*mProton + pP_mag*pP_mag);
    double EThird = sqrt(mPion*mPion + pThird_mag*pThird_mag);

    TVector3 pP(pP_mag*slc->reco.pfp.at(idxProton).trk.dir.x, pP_mag*slc->reco.pfp.at(idxProton).trk.dir.y, pP_mag*slc->reco.pfp.at(idxProton).trk.dir.z);
    TVector3 pThird(pThird_mag*slc->reco.pfp.at(idxThird).trk.dir.x, pThird_mag*slc->reco.pfp.at(idxThird).trk.dir.y, pThird_mag*slc->reco.pfp.at(idxThird).trk.dir.z);
    TVector3 pTotal = pP + pThird;

    W = sqrt( (EP + EThird)*(EP + EThird) - pTotal.Mag2() );

    return W;
  });

/*
const Var kW_ThreeP([](const caf::SRSliceProxy* slc) -> float {
    float W = -9999.;

    unsigned int idxP1 = (unsigned int) kRecoProtonIdx(slc); //Leading proton
    unsigned int idxP2 = (unsigned int) kScndProtonIdx(slc); //Charged pion or sub-leading proton
    unsigned int idxP3 = (unsigned int) kSidebandProton(slc);

    double pP1_mag = kRecoProtonP(slc);
    double pP2_mag = kScndProtonP(slc);
    double pP3_mag = kSidebandProtonP(slc);
    double EP1 = sqrt(mProton*mProton + pP1_mag*pP1_mag);
    double EP2 = sqrt(mProton*mProton + pP2_mag*pP2_mag);
    double EP3 = sqrt(mProton*mProton + pP3_mag*pP3_mag);

    TVector3 pP1(pP1_mag*slc->reco.pfp.at(idxP1).trk.dir.x, pP1_mag*slc->reco.pfp.at(idxP1).trk.dir.y, pP1_mag*slc->reco.pfp.at(idxP1).trk.dir.z);
    TVector3 pP2(pP2_mag*slc->reco.pfp.at(idxP2).trk.dir.x, pP2_mag*slc->reco.pfp.at(idxP2).trk.dir.y, pP2_mag*slc->reco.pfp.at(idxP2).trk.dir.z);
    TVector3 pP3(pP3_mag*slc->reco.pfp.at(idxP3).trk.dir.x, pP3_mag*slc->reco.pfp.at(idxP3).trk.dir.y, pP3_mag*slc->reco.pfp.at(idxP3).trk.dir.z);
    TVector3 pTotal = pP1 + pP2 + pP3;

    W = sqrt( (EP1+EP2+EP3)*(EP1+EP2+EP3) - pTotal.Mag2() );

    return W;
  });
*/

const Var kW([](const caf::SRSliceProxy* slc) -> float {
    double w = -9999.;
    if ( kScndProtonCandidate(slc) ) w = kW_Proton(slc);
    //else if ( kThreePSideband(slc) ) w = kW_ThreeP(slc);
    else if ( kPionSidebandBase(slc) ) w = kW_Pion(slc);
    assert(((void)"No valid selection for W value", w != -9999.));

    return w;
  });

const Var kW_ProtonTruth([](const caf::SRSliceProxy* slc) -> float {
    float W = -9999.;

    unsigned int idxProton = (unsigned int) kRecoProtonIdx(slc); //Leading proton
    unsigned int idxThird = (unsigned int) kScndProtonIdx(slc); //Charged pion or sub-leading proton

    double EP = slc->reco.pfp.at(idxProton).trk.truth.p.startE;
    double EThird = slc->reco.pfp.at(idxThird).trk.truth.p.startE;

    TVector3 pP(slc->reco.pfp.at(idxProton).trk.truth.p.startp.x, slc->reco.pfp.at(idxProton).trk.truth.p.startp.y, slc->reco.pfp.at(idxProton).trk.truth.p.startp.z);
    TVector3 pThird(slc->reco.pfp.at(idxThird).trk.truth.p.startp.x, slc->reco.pfp.at(idxThird).trk.truth.p.startp.y, slc->reco.pfp.at(idxThird).trk.truth.p.startp.z);
    TVector3 pTotal = pP + pThird;

    W = sqrt( (EP + EThird)*(EP + EThird) - pTotal.Mag2() );

    return W;
  });

const Var kW_PionTruth([](const caf::SRSliceProxy* slc) -> float {
    float W = -9999.;

    unsigned int idxProton = (unsigned int) kRecoProtonIdx(slc); //Leading proton
    unsigned int idxThird = (unsigned int) kSidebandPion(slc); //Charged pion or sub-leading proton

    double EP = slc->reco.pfp.at(idxProton).trk.truth.p.startE;
    double EThird = slc->reco.pfp.at(idxThird).trk.truth.p.startE;

    TVector3 pP(slc->reco.pfp.at(idxProton).trk.truth.p.startp.x, slc->reco.pfp.at(idxProton).trk.truth.p.startp.y, slc->reco.pfp.at(idxProton).trk.truth.p.startp.z);
    TVector3 pThird(slc->reco.pfp.at(idxThird).trk.truth.p.startp.x, slc->reco.pfp.at(idxThird).trk.truth.p.startp.y, slc->reco.pfp.at(idxThird).trk.truth.p.startp.z);
    TVector3 pTotal = pP + pThird;

    W = sqrt( (EP + EThird)*(EP + EThird) - pTotal.Mag2() );

    return W;
  });

/*
const Var kW_ThreePTruth([](const caf::SRSliceProxy* slc) -> float {
    float W = -9999.;

    unsigned int idxP1 = (unsigned int) kRecoProtonIdx(slc); //Leading proton
    unsigned int idxP2 = (unsigned int) kScndProtonIdx(slc); //Charged pion or sub-leading proton
    unsigned int idxP3 = (unsigned int) kSidebandProton(slc);

    double EP1 = slc->reco.pfp.at(idxP1).trk.truth.p.startE;
    double EP2 = slc->reco.pfp.at(idxP2).trk.truth.p.startE;
    double EP3 = slc->reco.pfp.at(idxP3).trk.truth.p.startE;

    TVector3 pP1(slc->reco.pfp.at(idxP1).trk.truth.p.startp.x, slc->reco.pfp.at(idxP1).trk.truth.p.startp.y, slc->reco.pfp.at(idxP1).trk.truth.p.startp.z);
    TVector3 pP2(slc->reco.pfp.at(idxP2).trk.truth.p.startp.x, slc->reco.pfp.at(idxP2).trk.truth.p.startp.y, slc->reco.pfp.at(idxP2).trk.truth.p.startp.z);
    TVector3 pP3(slc->reco.pfp.at(idxP3).trk.truth.p.startp.x, slc->reco.pfp.at(idxP3).trk.truth.p.startp.y, slc->reco.pfp.at(idxP3).trk.truth.p.startp.z);
    TVector3 pTotal = pP1 + pP2 + pP3;

    W = sqrt( (EP1+EP2+EP3)*(EP1+EP2+EP3) - pTotal.Mag2() );

    return W;
  });
*/

const Var kW_Truth([](const caf::SRSliceProxy* slc) -> float {
    double w = -9999.;
    if ( kScndProtonCandidate(slc) ) w = kW_ProtonTruth(slc);
    //else if ( kThreePSideband(slc) ) w = kW_ThreePTruth(slc);
    else if ( kPionSidebandBase(slc) ) w = kW_PionTruth(slc);
    assert(((void)"No valid selection for W_Truth value", w != -9999.));

    return w;
  });

const Var kW_Resid([](const caf::SRSliceProxy* slc) -> float {
    double recoW = kW(slc);
    double trueW = kW_Truth(slc);

    return (recoW - trueW) / trueW;
  });

const Var kW_ProtonResid([](const caf::SRSliceProxy* slc) -> float {
    double recoW = kW_Proton(slc);
    double trueW = kW_Truth(slc);

    return (recoW - trueW) / trueW;
  });

const Var kW_PionResid([](const caf::SRSliceProxy* slc) -> float {
    double recoW = kW_Pion(slc);
    double trueW = kW_Truth(slc);

    return (recoW - trueW) / trueW;
  });

const Var kW_Exp([](const caf::SRSliceProxy* slc) -> float {
    double wExp = -9999.;
    double eHad = kEHad(slc);
    double q2 = kQ2(slc);

    double wExp2 = mNeutron*mNeutron + 2*mNeutron*eHad - q2;
    if ( wExp2 >= 0. ) wExp = sqrt(wExp2);

    return wExp;
  });

const Var kW_Exp_CheatingMuonP([](const caf::SRSliceProxy* slc) -> float {
    double wExp = -9999.;
    double eHad = kEHad(slc);
    double q2 = kQ2_CheatingMuonP(slc);

    double wExp2 = mNeutron*mNeutron + 2*mNeutron*eHad - q2;
    if ( wExp2 >= 0. ) wExp = sqrt(wExp2);

    return wExp;
  });

const Var kW_Exp_CheatingAll([](const caf::SRSliceProxy* slc) -> float {
    double wExp = -9999.;
    double eHad = kEHad_CheatingProtonPs(slc);
    double q2 = kQ2_CheatingAll(slc);

    double wExp2 = mNeutron*mNeutron + 2*mNeutron*eHad - q2;
    if ( wExp2 >= 0. ) wExp = sqrt(wExp2);

    return wExp;
  });

const Var kW_Exp_Truth([](const caf::SRSliceProxy* slc) -> float {
    double wExp = -9999.;
    double eHad = kEHad_Truth(slc);
    if ( eHad == -9999. ) return wExp;
    double q2 = kQ2_Truth(slc);

    double wExp2 = mNeutron*mNeutron + 2*mNeutron*eHad - q2;
    if ( wExp2 >= 0. ) wExp = sqrt(wExp2);

    return wExp;
  });

const Var kW_Exp2([](const caf::SRSliceProxy* slc) -> float {
    double eHad = kEHad(slc);
    double q2 = kQ2(slc);

    double wExp2 = mNeutron*mNeutron + 2*mNeutron*eHad - q2;

    return wExp2;
  });

const Var kW_Exp2_CheatingMuonP([](const caf::SRSliceProxy* slc) -> float {
    double eHad = kEHad(slc);
    double q2 = kQ2_CheatingMuonP(slc);

    double wExp2 = mNeutron*mNeutron + 2*mNeutron*eHad - q2;

    return wExp2;
  });

const Var kW_Exp2_CheatingAll([](const caf::SRSliceProxy* slc) -> float {
    double eHad = kEHad_CheatingProtonPs(slc);
    double q2 = kQ2_CheatingAll(slc);

    double wExp2 = mNeutron*mNeutron + 2*mNeutron*eHad - q2;

    return wExp2;
  });

const Var kW_Exp2_Truth([](const caf::SRSliceProxy* slc) -> float {
    double wExp2 = -9999.;
    if ( slc->truth.index < 0 ) return wExp2;
    double eHad = kEHad_Truth(slc);
    double q2 = kQ2_Truth(slc);

    wExp2 = mNeutron*mNeutron + 2*mNeutron*eHad - q2;

    return wExp2;
  });

const Var kGENIEQ2([](const caf::SRSliceProxy* slc) -> float {
    double q2 = -9999.;

    if ( !isnan(slc->truth.Q2) ) q2 = slc->truth.Q2;
    return q2;
  });

const Var kGENIEW([](const caf::SRSliceProxy* slc) -> float {
    double w = -9999.;

    if ( !isnan(slc->truth.w) ) w = slc->truth.w;
    return w;
  });

const Var kW_ExpResid([](const caf::SRSliceProxy* slc) -> float {
    double recoW = kW_Exp(slc);
    double trueW = kGENIEW(slc);

    if ( trueW == -9999. ) return trueW;
    else return (recoW - trueW) / trueW;
  });

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Nasty TKI section, each in proton mass, pion mass, and using truth variables

//Single-transverse momentum imbalance
const Var kDeltaPT_Single([](const caf::SRSliceProxy* slc) -> float {
    float dpTMag = -9999.;

    unsigned int idxMuon = (unsigned int) kRecoMuonIdx(slc);
    unsigned int idxProton = (unsigned int) kRecoProtonIdx(slc); //Leading proton

    double pMu_mag = kRecoMuonPNew(slc);
    double pP_mag = kRecoProtonP(slc);

    TVector3 pMu(pMu_mag*slc->reco.pfp.at(idxMuon).trk.dir.x, pMu_mag*slc->reco.pfp.at(idxMuon).trk.dir.y, pMu_mag*slc->reco.pfp.at(idxMuon).trk.dir.z);
    TVector3 pP(pP_mag*slc->reco.pfp.at(idxProton).trk.dir.x, pP_mag*slc->reco.pfp.at(idxProton).trk.dir.y, pP_mag*slc->reco.pfp.at(idxProton).trk.dir.z);

    const auto& vtx = slc->vertex;
    TVector3 vec_vtx(vtx.x, vtx.y, vtx.z);
    TVector3 vec_numi_to_vtx = (dFromNuMI + vec_vtx).Unit();

    //pT = p - pL
    TVector3 pTMu = pMu - (pMu.Dot(vec_numi_to_vtx) * vec_numi_to_vtx);
    TVector3 pTP = pP - (pP.Dot(vec_numi_to_vtx) * vec_numi_to_vtx);
    TVector3 dpT = pTMu + pTP;

    dpTMag = dpT.Mag();

    return dpTMag;
  });

const Var kDeltaPT_Proton([](const caf::SRSliceProxy* slc) -> float {
    float dpTMag = -9999.;

    unsigned int idxMuon = (unsigned int) kRecoMuonIdx(slc);
    unsigned int idxProton = (unsigned int) kRecoProtonIdx(slc); //Leading proton
    unsigned int idxThird = (unsigned int) kScndProtonIdx(slc); //Charged pion or sub-leading proton

    double pMu_mag = kRecoMuonPNew(slc);
    double pP_mag = kRecoProtonP(slc);
    double pThird_mag = kScndProtonP(slc);

    TVector3 pMu(pMu_mag*slc->reco.pfp.at(idxMuon).trk.dir.x, pMu_mag*slc->reco.pfp.at(idxMuon).trk.dir.y, pMu_mag*slc->reco.pfp.at(idxMuon).trk.dir.z);
    TVector3 pP(pP_mag*slc->reco.pfp.at(idxProton).trk.dir.x, pP_mag*slc->reco.pfp.at(idxProton).trk.dir.y, pP_mag*slc->reco.pfp.at(idxProton).trk.dir.z);
    TVector3 pThird(pThird_mag*slc->reco.pfp.at(idxThird).trk.dir.x, pThird_mag*slc->reco.pfp.at(idxThird).trk.dir.y, pThird_mag*slc->reco.pfp.at(idxThird).trk.dir.z);

    const auto& vtx = slc->vertex;
    TVector3 vec_vtx(vtx.x, vtx.y, vtx.z);
    TVector3 vec_numi_to_vtx = (dFromNuMI + vec_vtx).Unit();

    //pT = p - pL
    TVector3 pTMu = pMu - (pMu.Dot(vec_numi_to_vtx) * vec_numi_to_vtx);
    TVector3 pTP = pP - (pP.Dot(vec_numi_to_vtx) * vec_numi_to_vtx);
    TVector3 pTThird = pThird - (pThird.Dot(vec_numi_to_vtx) * vec_numi_to_vtx);
    TVector3 dpT = pTMu + pTP + pTThird;

    dpTMag = dpT.Mag();

    if ( isnan(dpTMag) ) std::cout << std::endl << "~~~~~~~~~~~~~~~~~~~~~~~~~ NaN deltapT!";
    return dpTMag;
  });

/*
For Minerba
const Var kDeltaPT([](const caf::SRSliceProxy* slc) -> float {
    float dpTMag = -9999.;

    //Is this the beam direction in BNB detector coordinates? 
    TVector3 neutrinoDirection(0., 0., 1.);

    if ( kNRecoProtons(slc) < 2 ) return dpTMag;
    
    unsigned int idxMuon = (unsigned int) kRecoMuonIdx(slc);
    unsigned int idxProton1 = (unsigned int) kLeadingProton(slc);
    unsigned int idxProton2 = (unsigned int) kSecondProton(slc);
    
    double pMu_mag = slc->reco.pfp.at(idxMuon).trk.rangeP.p_muon;
    TVector3 pMu(pMu_mag*slc->reco.pfp.at(idxMuon).trk.dir.x, pMu_mag*slc->reco.pfp.at(idxMuon).trk.dir.y, pMu_mag*slc->reco.pfp.at(idxMuon).trk.dir.z);
    
    double pP1_mag = slc->reco.pfp.at(idxProton1).trk.rangeP.p_muon;
    TVector3 pP1(pP1_mag*slc->reco.pfp.at(idxProton1).trk.dir.x, pP1_mag*slc->reco.pfp.at(idxProton1).trk.dir.y, pP1_mag*slc->reco.pfp.at(idxProton1).trk.dir.z);
    
    double pP2_mag = slc->reco.pfp.at(idxProton2).trk.rangeP.p_muon;
    TVector3 pP2(pP2_mag*slc->reco.pfp.at(idxProton2).trk.dir.x, pP2_mag*slc->reco.pfp.at(idxProton2).trk.dir.y, pP2_mag*slc->reco.pfp.at(idxProton2).trk.dir.z);

    //Find the transver component (pT) by subtructiong the longitudinal component (pL) from the total vector, pT = p - pL
    TVector3 pTMu = pMu - (pMu.Dot(NuDirection_NuMI) * NuDirection_NuMI);
    TVector3 pTP1 = pP1 - (pP1.Dot(NuDirection_NuMI) * NuDirection_NuMI);
    TVector3 pTP2 = pP2 - (pP2.Dot(NuDirection_NuMI) * NuDirection_NuMI);

    //Return the total transverse momentum
    TVector3 dpT = pTMu + pTP1 + pTP2;
    dpTMag = dpT.Mag();
    return dptMag;
  });
*/

const Var kDeltaPT_CheatingMuon([](const caf::SRSliceProxy* slc) -> float {
    float dpTMag = -9999.;

    unsigned int idxMuon = (unsigned int) kRecoMuonIdx(slc);
    unsigned int idxProton = (unsigned int) kRecoProtonIdx(slc); //Leading proton
    if ( kScndProtonIdx(slc) < 0 ) return dpTMag;
    unsigned int idxThird = (unsigned int) kScndProtonIdx(slc); //Charged pion or sub-leading proton

    double pP_mag = kRecoProtonP(slc);
    double pThird_mag = kScndProtonP(slc);

    TVector3 pMu(slc->reco.pfp.at(idxMuon).trk.truth.p.startp.x, slc->reco.pfp.at(idxMuon).trk.truth.p.startp.y, slc->reco.pfp.at(idxMuon).trk.truth.p.startp.z);
    TVector3 pP(pP_mag*slc->reco.pfp.at(idxProton).trk.dir.x, pP_mag*slc->reco.pfp.at(idxProton).trk.dir.y, pP_mag*slc->reco.pfp.at(idxProton).trk.dir.z);
    TVector3 pThird(pThird_mag*slc->reco.pfp.at(idxThird).trk.dir.x, pThird_mag*slc->reco.pfp.at(idxThird).trk.dir.y, pThird_mag*slc->reco.pfp.at(idxThird).trk.dir.z);

    const auto& vtx = slc->vertex;
    TVector3 vec_vtx(vtx.x, vtx.y, vtx.z);
    TVector3 vec_numi_to_vtx = (dFromNuMI + vec_vtx).Unit();

    //pT = p - pL
    TVector3 pTMu = pMu - (pMu.Dot(vec_numi_to_vtx) * vec_numi_to_vtx);
    TVector3 pTP = pP - (pP.Dot(vec_numi_to_vtx) * vec_numi_to_vtx);
    TVector3 pTThird = pThird - (pThird.Dot(vec_numi_to_vtx) * vec_numi_to_vtx);
    TVector3 dpT = pTMu + pTP + pTThird;

    dpTMag = dpT.Mag();

    return dpTMag;
  });

const Var kDeltaPT_CheatingAngles([](const caf::SRSliceProxy* slc) -> float {
    float dpTMag = -9999.;

    unsigned int idxMuon = (unsigned int) kRecoMuonIdx(slc);
    unsigned int idxProton = (unsigned int) kRecoProtonIdx(slc); //Leading proton
    if ( kScndProtonIdx(slc) < 0 ) return dpTMag;
    unsigned int idxThird = (unsigned int) kScndProtonIdx(slc); //Charged pion or sub-leading proton

    double pMu_mag = kRecoMuonPNew(slc);
    double pP_mag = kRecoProtonP(slc);
    double pThird_mag = kScndProtonP(slc);

    TVector3 pMu(slc->reco.pfp.at(idxMuon).trk.truth.p.startp.x, slc->reco.pfp.at(idxMuon).trk.truth.p.startp.y, slc->reco.pfp.at(idxMuon).trk.truth.p.startp.z);
    TVector3 pP(slc->reco.pfp.at(idxProton).trk.truth.p.startp.x, slc->reco.pfp.at(idxProton).trk.truth.p.startp.y, slc->reco.pfp.at(idxProton).trk.truth.p.startp.z);
    TVector3 pThird(slc->reco.pfp.at(idxThird).trk.truth.p.startp.x, slc->reco.pfp.at(idxThird).trk.truth.p.startp.y, slc->reco.pfp.at(idxThird).trk.truth.p.startp.z);

    pMu *= pMu_mag / pMu.Mag();
    pP *= pP_mag / pP.Mag();
    pThird *= pThird_mag / pThird.Mag();

    const auto& vtx = slc->vertex;
    TVector3 vec_vtx(vtx.x, vtx.y, vtx.z);
    TVector3 vec_numi_to_vtx = (dFromNuMI + vec_vtx).Unit();

    //pT = p - pL
    TVector3 pTMu = pMu - (pMu.Dot(vec_numi_to_vtx) * vec_numi_to_vtx);
    TVector3 pTP = pP - (pP.Dot(vec_numi_to_vtx) * vec_numi_to_vtx);
    TVector3 pTThird = pThird - (pThird.Dot(vec_numi_to_vtx) * vec_numi_to_vtx);
    TVector3 dpT = pTMu + pTP + pTThird;

    dpTMag = dpT.Mag();

    return dpTMag;
  });

const Var kDeltaPT_Pion([](const caf::SRSliceProxy* slc) -> float {
    float dpTMag = -9999.;

    unsigned int idxMuon = (unsigned int) kRecoMuonIdx(slc);
    unsigned int idxProton = (unsigned int) kRecoProtonIdx(slc); //Leading proton
    unsigned int idxThird = (unsigned int) kSidebandPion(slc); //Charged pion or sub-leading proton

    double pMu_mag = kRecoMuonPNew(slc);
    double pP_mag = kRecoProtonP(slc);
    double pThird_mag = kSidebandPionP(slc);

    TVector3 pMu(pMu_mag*slc->reco.pfp.at(idxMuon).trk.dir.x, pMu_mag*slc->reco.pfp.at(idxMuon).trk.dir.y, pMu_mag*slc->reco.pfp.at(idxMuon).trk.dir.z);
    TVector3 pP(pP_mag*slc->reco.pfp.at(idxProton).trk.dir.x, pP_mag*slc->reco.pfp.at(idxProton).trk.dir.y, pP_mag*slc->reco.pfp.at(idxProton).trk.dir.z);
    TVector3 pThird(pThird_mag*slc->reco.pfp.at(idxThird).trk.dir.x, pThird_mag*slc->reco.pfp.at(idxThird).trk.dir.y, pThird_mag*slc->reco.pfp.at(idxThird).trk.dir.z);

    const auto& vtx = slc->vertex;
    TVector3 vec_vtx(vtx.x, vtx.y, vtx.z);
    TVector3 vec_numi_to_vtx = (dFromNuMI + vec_vtx).Unit();

    //pT = p - pL
    TVector3 pTMu = pMu - (pMu.Dot(vec_numi_to_vtx) * vec_numi_to_vtx);
    TVector3 pTP = pP - (pP.Dot(vec_numi_to_vtx) * vec_numi_to_vtx);
    TVector3 pTThird = pThird - (pThird.Dot(vec_numi_to_vtx) * vec_numi_to_vtx);
    TVector3 dpT = pTMu + pTP + pTThird;

    dpTMag = dpT.Mag();

    return dpTMag;
  });

/*
const Var kDeltaPT_ThreeP([](const caf::SRSliceProxy* slc) -> float {
    float dpTMag = -9999.;

    unsigned int idxMuon = (unsigned int) kRecoMuonIdx(slc);
    unsigned int idxP1 = (unsigned int) kRecoProtonIdx(slc); //Leading proton
    unsigned int idxP2 = (unsigned int) kScndProtonIdx(slc); //Charged pion or sub-leading proton
    unsigned int idxP3 = (unsigned int) kSidebandProton(slc);

    double pMu_mag = kRecoMuonPNew(slc);
    double pP1_mag = kRecoProtonP(slc);
    double pP2_mag = kScndProtonP(slc);
    double pP3_mag = kSidebandProtonP(slc);

    TVector3 pMu(pMu_mag*slc->reco.pfp.at(idxMuon).trk.dir.x, pMu_mag*slc->reco.pfp.at(idxMuon).trk.dir.y, pMu_mag*slc->reco.pfp.at(idxMuon).trk.dir.z);
    TVector3 pP1(pP1_mag*slc->reco.pfp.at(idxP1).trk.dir.x, pP1_mag*slc->reco.pfp.at(idxP1).trk.dir.y, pP1_mag*slc->reco.pfp.at(idxP1).trk.dir.z);
    TVector3 pP2(pP2_mag*slc->reco.pfp.at(idxP2).trk.dir.x, pP2_mag*slc->reco.pfp.at(idxP2).trk.dir.y, pP2_mag*slc->reco.pfp.at(idxP2).trk.dir.z);
    TVector3 pP3(pP3_mag*slc->reco.pfp.at(idxP3).trk.dir.x, pP3_mag*slc->reco.pfp.at(idxP3).trk.dir.y, pP3_mag*slc->reco.pfp.at(idxP3).trk.dir.z);

    const auto& vtx = slc->vertex;
    TVector3 vec_vtx(vtx.x, vtx.y, vtx.z);
    TVector3 vec_numi_to_vtx = (dFromNuMI + vec_vtx).Unit();

    //pT = p - pL
    TVector3 pTMu = pMu - (pMu.Dot(vec_numi_to_vtx) * vec_numi_to_vtx);
    TVector3 pTP1 = pP1 - (pP1.Dot(vec_numi_to_vtx) * vec_numi_to_vtx);
    TVector3 pTP2 = pP2 - (pP2.Dot(vec_numi_to_vtx) * vec_numi_to_vtx);
    TVector3 pTP3 = pP3 - (pP3.Dot(vec_numi_to_vtx) * vec_numi_to_vtx);
    TVector3 dpT = pTMu + pTP1 + pTP2 + pTP3;

    dpTMag = dpT.Mag();

    return dpTMag;
  });
*/

const Var kDeltaPT_ProtonTruth([](const caf::SRSliceProxy* slc) -> float {
    float dpTMag = -9999.;

    if ( slc->truth.index < 0 ) return dpTMag;
    TVector3 NuDirection_Truth = (TVector3(slc->truth.momentum.x, slc->truth.momentum.y, slc->truth.momentum.z)).Unit().Unit();

    unsigned int idxMuon = (unsigned int) kRecoMuonIdx(slc);
    unsigned int idxProton = (unsigned int) kRecoProtonIdx(slc); //Leading proton
    unsigned int idxThird = (unsigned int) kScndProtonIdx(slc); //Charged pion or sub-leading proton

    TVector3 pMu(slc->reco.pfp.at(idxMuon).trk.truth.p.startp.x, slc->reco.pfp.at(idxMuon).trk.truth.p.startp.y, slc->reco.pfp.at(idxMuon).trk.truth.p.startp.z);
    TVector3 pP(slc->reco.pfp.at(idxProton).trk.truth.p.startp.x, slc->reco.pfp.at(idxProton).trk.truth.p.startp.y, slc->reco.pfp.at(idxProton).trk.truth.p.startp.z);
    TVector3 pThird(slc->reco.pfp.at(idxThird).trk.truth.p.startp.x, slc->reco.pfp.at(idxThird).trk.truth.p.startp.y, slc->reco.pfp.at(idxThird).trk.truth.p.startp.z);

    //pT = p - pL
    TVector3 pTMu = pMu - (pMu.Dot(NuDirection_Truth) * NuDirection_Truth);
    TVector3 pTP = pP - (pP.Dot(NuDirection_Truth) * NuDirection_Truth);
    TVector3 pTThird = pThird - (pThird.Dot(NuDirection_Truth) * NuDirection_Truth);
    TVector3 dpT = pTMu + pTP + pTThird;

    dpTMag = dpT.Mag();

    return dpTMag;
  });

const Var kDeltaPT_PionTruth([](const caf::SRSliceProxy* slc) -> float {
    float dpTMag = -9999.;

    if ( slc->truth.index < 0 ) return dpTMag;
    TVector3 NuDirection_Truth = (TVector3(slc->truth.momentum.x, slc->truth.momentum.y, slc->truth.momentum.z)).Unit().Unit();

    unsigned int idxMuon = (unsigned int) kRecoMuonIdx(slc);
    unsigned int idxProton = (unsigned int) kRecoProtonIdx(slc); //Leading proton
    unsigned int idxThird = (unsigned int) kSidebandPion(slc); //Charged pion or sub-leading proton

    TVector3 pMu(slc->reco.pfp.at(idxMuon).trk.truth.p.startp.x, slc->reco.pfp.at(idxMuon).trk.truth.p.startp.y, slc->reco.pfp.at(idxMuon).trk.truth.p.startp.z);
    TVector3 pP(slc->reco.pfp.at(idxProton).trk.truth.p.startp.x, slc->reco.pfp.at(idxProton).trk.truth.p.startp.y, slc->reco.pfp.at(idxProton).trk.truth.p.startp.z);
    TVector3 pThird(slc->reco.pfp.at(idxThird).trk.truth.p.startp.x, slc->reco.pfp.at(idxThird).trk.truth.p.startp.y, slc->reco.pfp.at(idxThird).trk.truth.p.startp.z);

    //pT = p - pL
    TVector3 pTMu = pMu - (pMu.Dot(NuDirection_Truth) * NuDirection_Truth);
    TVector3 pTP = pP - (pP.Dot(NuDirection_Truth) * NuDirection_Truth);
    TVector3 pTThird = pThird - (pThird.Dot(NuDirection_Truth) * NuDirection_Truth);
    TVector3 dpT = pTMu + pTP + pTThird;

    dpTMag = dpT.Mag();

    return dpTMag;
  });

/*
const Var kDeltaPT_ThreePTruth([](const caf::SRSliceProxy* slc) -> float {
    float dpTMag = -9999.;

    if ( slc->truth.index < 0 ) return dpTMag;
    TVector3 NuDirection_Truth = (TVector3(slc->truth.momentum.x, slc->truth.momentum.y, slc->truth.momentum.z)).Unit().Unit();

    unsigned int idxMuon = (unsigned int) kRecoMuonIdx(slc);
    unsigned int idxP1 = (unsigned int) kRecoProtonIdx(slc); //Leading proton
    unsigned int idxP2 = (unsigned int) kScndProtonIdx(slc); //Charged pion or sub-leading proton
    unsigned int idxP3 = (unsigned int) kSidebandProton(slc);

    TVector3 pMu(slc->reco.pfp.at(idxMuon).trk.truth.p.startp.x, slc->reco.pfp.at(idxMuon).trk.truth.p.startp.y, slc->reco.pfp.at(idxMuon).trk.truth.p.startp.z);
    TVector3 pP1(slc->reco.pfp.at(idxP1).trk.truth.p.startp.x, slc->reco.pfp.at(idxP1).trk.truth.p.startp.y, slc->reco.pfp.at(idxP1).trk.truth.p.startp.z);
    TVector3 pP2(slc->reco.pfp.at(idxP2).trk.truth.p.startp.x, slc->reco.pfp.at(idxP2).trk.truth.p.startp.y, slc->reco.pfp.at(idxP2).trk.truth.p.startp.z);
    TVector3 pP3(slc->reco.pfp.at(idxP3).trk.truth.p.startp.x, slc->reco.pfp.at(idxP3).trk.truth.p.startp.y, slc->reco.pfp.at(idxP3).trk.truth.p.startp.z);

    //pT = p - pL
    TVector3 pTMu = pMu - (pMu.Dot(NuDirection_Truth) * NuDirection_Truth);
    TVector3 pTP1 = pP1 - (pP1.Dot(NuDirection_Truth) * NuDirection_Truth);
    TVector3 pTP2 = pP2 - (pP2.Dot(NuDirection_Truth) * NuDirection_Truth);
    TVector3 pTP3 = pP3 - (pP3.Dot(NuDirection_Truth) * NuDirection_Truth);
    TVector3 dpT = pTMu + pTP1 + pTP2 + pTP3;

    dpTMag = dpT.Mag();

    return dpTMag;
  });
*/

const Var kDeltaPT([](const caf::SRSliceProxy* slc) -> float {
    double dPT = -9999.;
    if ( kScndProtonCandidate(slc) ) dPT = kDeltaPT_Proton(slc);
    //else if ( kThreePSideband(slc) ) dPT = kDeltaPT_ThreeP(slc);
    else if ( kPionSidebandBase(slc) ) dPT = kDeltaPT_Pion(slc);
    assert(((void)"No valid selection for DeltaPT value", dPT != -9999.));

    return dPT;
  });

const Var kDeltaPT_OldTruth([](const caf::SRSliceProxy* slc) -> float {
    double dPT = -99999.; //Note the extra 9, as -9999. is a valid response for selected cosmic slices
    if ( kScndProtonCandidate(slc) ) dPT = kDeltaPT_ProtonTruth(slc);
    //else if ( kThreePSideband(slc) ) dPT = kDeltaPT_ThreePTruth(slc);
    else if ( kPionSidebandBase(slc) ) dPT = kDeltaPT_PionTruth(slc);
    assert(((void)"No valid selection for DeltaPT value", dPT != -99999.));

    return dPT;
  });

//Built out of SRTrueInteraction rather than truth info matched to reco tracks
//Considers all charged pions and protons over threshold, so should be correct for all selections
const Var kDeltaPT_FullTruth([](const caf::SRSliceProxy* slc) -> float {
    float dpTMag = -9999.;
    if ( slc->truth.index < 0 ) return dpTMag;
    const auto &nu = slc->truth;

    TVector3 NuDirection_Truth = (TVector3(nu.momentum.x, nu.momentum.y, nu.momentum.z)).Unit().Unit();
    TVector3 p3;
    TVector3 pTMu(0., 0., 0.);
    TVector3 dpT(0., 0., 0.);

    for ( unsigned int i = 0; i < nu.prim.size(); i++ ) {
      const auto &prim = nu.prim.at(i);
      if ( prim.start_process != 0 ) continue;
      p3 = {prim.startp.x, prim.startp.y, prim.startp.z};
      if ( p3.Mag() == std::hypot(-9999., -9999., -9999.) ) return dpTMag;
      if ( abs(prim.pdg) == 13 ) {
        pTMu = p3 - (p3.Dot(NuDirection_Truth) * NuDirection_Truth);
        dpT += pTMu;
      }
      else if ( abs(prim.pdg) == 2212 && p3.Mag() > .35 ) dpT += p3 - (p3.Dot(NuDirection_Truth) * NuDirection_Truth);
      else if ( abs(prim.pdg) == 211 || abs(prim.pdg) == 321 ) dpT += p3 - (p3.Dot(NuDirection_Truth) * NuDirection_Truth);
    }

    if ( pTMu.Mag() == 0. ) return dpTMag;
    dpTMag = dpT.Mag();

    return dpTMag;
  });

const Var kDeltaPT_Truth([](const caf::SRSliceProxy* slc) -> float {
    float dpTMag = -9999.;
    if ( slc->truth.index < 0 ) return dpTMag;
    const auto &nu = slc->truth;

    TVector3 NuDirection_Truth = (TVector3(nu.momentum.x, nu.momentum.y, nu.momentum.z)).Unit().Unit();
    TVector3 p3;
    TVector3 pTMu(0., 0., 0.);
    TVector3 dpT(0., 0., 0.);

    int idxP1 = -1;
    int idxP2 = -1;
    double leading = -9999.;
    double second = -9999.;

    for ( unsigned int i = 0; i < nu.prim.size(); i++ ) {
      const auto &prim = nu.prim.at(i);
      if ( prim.start_process != 0 ) continue;
      p3 = {prim.startp.x, prim.startp.y, prim.startp.z};
      if ( p3.Mag() == std::hypot(-9999., -9999., -9999.) ) return dpTMag;
      if ( abs(prim.pdg) == 13 ) {
        pTMu = p3 - (p3.Dot(NuDirection_Truth) * NuDirection_Truth);
        dpT += pTMu;
      }
      else if ( abs(prim.pdg) == 2212 && p3.Mag() > .35 ) {
        if ( p3.Mag() > leading ) {
          second = leading;
          leading = p3.Mag();
          idxP2 = idxP1;
          idxP1 = i;
        }
        else if ( p3.Mag() > second ) {
          second = p3.Mag();
          idxP2 = i;
        }
      }
    }

    if ( idxP2 == -1 ) return dpTMag;
    TVector3 P1 = {nu.prim.at(idxP1).startp.x, nu.prim.at(idxP1).startp.y, nu.prim.at(idxP1).startp.z};
    TVector3 P2 = {nu.prim.at(idxP2).startp.x, nu.prim.at(idxP2).startp.y, nu.prim.at(idxP2).startp.z};
    dpT += P1 - (P1.Dot(NuDirection_Truth) * NuDirection_Truth);
    dpT += P2 - (P2.Dot(NuDirection_Truth) * NuDirection_Truth);

    if ( pTMu.Mag() == 0. ) return dpTMag;
    dpTMag = dpT.Mag();

    return dpTMag;
  });

const Var kDeltaPT_Resid([](const caf::SRSliceProxy* slc) -> float {
    double recoDeltaPT = kDeltaPT(slc);
    double trueDeltaPT = kDeltaPT_Truth(slc);

    if ( trueDeltaPT == -9999. || trueDeltaPT == 0. ) return -9999.;
    return (recoDeltaPT - trueDeltaPT) / trueDeltaPT;
  });

const Var kDeltaPT_ProtonResid([](const caf::SRSliceProxy* slc) -> float {
    double recoDeltaPT = kDeltaPT_Proton(slc);
    double trueDeltaPT = kDeltaPT_Truth(slc);

    if ( trueDeltaPT == -9999. || trueDeltaPT == 0. ) return -9999.;
    return (recoDeltaPT - trueDeltaPT) / trueDeltaPT;
  });

const Var kDeltaPT_CheatingMuonResid([](const caf::SRSliceProxy* slc) -> float {
    double recoDeltaPT = kDeltaPT_CheatingMuon(slc);
    double trueDeltaPT = kDeltaPT_Truth(slc);

    if ( trueDeltaPT == -9999. || trueDeltaPT == 0. ) return -9999.;
    return (recoDeltaPT - trueDeltaPT) / trueDeltaPT;
  });

const Var kDeltaPT_CheatingAnglesResid([](const caf::SRSliceProxy* slc) -> float {
    double recoDeltaPT = kDeltaPT_CheatingAngles(slc);
    double trueDeltaPT = kDeltaPT_Truth(slc);

    if ( trueDeltaPT == -9999. || trueDeltaPT == 0. ) return -9999.;
    return (recoDeltaPT - trueDeltaPT) / trueDeltaPT;
  });

//Transverse boosting angle, "This observable quantifies whether the hadronic system is accelerated or decelerated by nuclear effects."
const Var kDeltaAlphaT_Single([](const caf::SRSliceProxy* slc) -> float {
    float daT = -9999.;

    unsigned int idxMuon = (unsigned int) kRecoMuonIdx(slc);
    unsigned int idxProton = (unsigned int) kRecoProtonIdx(slc); //Leading proton

    double pMu_mag = kRecoMuonPNew(slc);
    double pP_mag = kRecoProtonP(slc);

    TVector3 pMu(pMu_mag*slc->reco.pfp.at(idxMuon).trk.dir.x, pMu_mag*slc->reco.pfp.at(idxMuon).trk.dir.y, pMu_mag*slc->reco.pfp.at(idxMuon).trk.dir.z);
    TVector3 pP(pP_mag*slc->reco.pfp.at(idxProton).trk.dir.x, pP_mag*slc->reco.pfp.at(idxProton).trk.dir.y, pP_mag*slc->reco.pfp.at(idxProton).trk.dir.z);

    const auto& vtx = slc->vertex;
    TVector3 vec_vtx(vtx.x, vtx.y, vtx.z);
    TVector3 vec_numi_to_vtx = (dFromNuMI + vec_vtx).Unit();

    //pT = p - pL
    TVector3 pTMu = pMu - (pMu.Dot(vec_numi_to_vtx) * vec_numi_to_vtx);
    TVector3 pTP = pP - (pP.Dot(vec_numi_to_vtx) * vec_numi_to_vtx);
    TVector3 dpT = pTMu + pTP;

    daT = TMath::ACos( -1. * (pTMu.Dot(dpT)) / (pTMu.Mag() * dpT.Mag()) );

    return daT;
  });

const Var kDeltaAlphaT_Proton([](const caf::SRSliceProxy* slc) -> float {
    float daT = -9999.;

    unsigned int idxMuon = (unsigned int) kRecoMuonIdx(slc);
    unsigned int idxProton = (unsigned int) kRecoProtonIdx(slc); //Leading proton
    unsigned int idxThird = (unsigned int) kScndProtonIdx(slc); //Charged pion or sub-leading proton

    double pMu_mag = kRecoMuonPNew(slc);
    double pP_mag = kRecoProtonP(slc);
    double pThird_mag = kScndProtonP(slc);

    TVector3 pMu(pMu_mag*slc->reco.pfp.at(idxMuon).trk.dir.x, pMu_mag*slc->reco.pfp.at(idxMuon).trk.dir.y, pMu_mag*slc->reco.pfp.at(idxMuon).trk.dir.z);
    TVector3 pP(pP_mag*slc->reco.pfp.at(idxProton).trk.dir.x, pP_mag*slc->reco.pfp.at(idxProton).trk.dir.y, pP_mag*slc->reco.pfp.at(idxProton).trk.dir.z);
    TVector3 pThird(pThird_mag*slc->reco.pfp.at(idxThird).trk.dir.x, pThird_mag*slc->reco.pfp.at(idxThird).trk.dir.y, pThird_mag*slc->reco.pfp.at(idxThird).trk.dir.z);

    const auto& vtx = slc->vertex;
    TVector3 vec_vtx(vtx.x, vtx.y, vtx.z);
    TVector3 vec_numi_to_vtx = (dFromNuMI + vec_vtx).Unit();

    //pT = p - pL
    TVector3 pTMu = pMu - (pMu.Dot(vec_numi_to_vtx) * vec_numi_to_vtx);
    TVector3 pTP = pP - (pP.Dot(vec_numi_to_vtx) * vec_numi_to_vtx);
    TVector3 pTThird = pThird - (pThird.Dot(vec_numi_to_vtx) * vec_numi_to_vtx);
    TVector3 dpT = pTMu + pTP + pTThird;

    daT = TMath::ACos( -1. * (pTMu.Dot(dpT)) / (pTMu.Mag() * dpT.Mag()) );

    return daT;
  });

//Temporary, use true muon momentum to characterize the impact of getting it wrong
const Var kDeltaAlphaT_CheatingMuon([](const caf::SRSliceProxy* slc) -> float {
    float daT = -9999.;

    unsigned int idxMuon = (unsigned int) kRecoMuonIdx(slc);
    unsigned int idxProton = (unsigned int) kRecoProtonIdx(slc); //Leading proton
    if ( kScndProtonIdx(slc) < 0 ) return daT;
    unsigned int idxThird = (unsigned int) kScndProtonIdx(slc); //Charged pion or sub-leading proton

    double pP_mag = kRecoProtonP(slc);
    double pThird_mag = kScndProtonP(slc);

    TVector3 pMu(slc->reco.pfp.at(idxMuon).trk.truth.p.startp.x, slc->reco.pfp.at(idxMuon).trk.truth.p.startp.y, slc->reco.pfp.at(idxMuon).trk.truth.p.startp.z);
    TVector3 pP(pP_mag*slc->reco.pfp.at(idxProton).trk.dir.x, pP_mag*slc->reco.pfp.at(idxProton).trk.dir.y, pP_mag*slc->reco.pfp.at(idxProton).trk.dir.z);
    TVector3 pThird(pThird_mag*slc->reco.pfp.at(idxThird).trk.dir.x, pThird_mag*slc->reco.pfp.at(idxThird).trk.dir.y, pThird_mag*slc->reco.pfp.at(idxThird).trk.dir.z);

    const auto& vtx = slc->vertex;
    TVector3 vec_vtx(vtx.x, vtx.y, vtx.z);
    TVector3 vec_numi_to_vtx = (dFromNuMI + vec_vtx).Unit();

    //pT = p - pL
    TVector3 pTMu = pMu - (pMu.Dot(vec_numi_to_vtx) * vec_numi_to_vtx);
    TVector3 pTP = pP - (pP.Dot(vec_numi_to_vtx) * vec_numi_to_vtx);
    TVector3 pTThird = pThird - (pThird.Dot(vec_numi_to_vtx) * vec_numi_to_vtx);
    TVector3 dpT = pTMu + pTP + pTThird;

    daT = TMath::ACos( -1. * (pTMu.Dot(dpT)) / (pTMu.Mag() * dpT.Mag()) );

    return daT;
  });

const Var kDeltaAlphaT_CheatingAngles([](const caf::SRSliceProxy* slc) -> float {
    float daT = -9999.;

    unsigned int idxMuon = (unsigned int) kRecoMuonIdx(slc);
    unsigned int idxProton = (unsigned int) kRecoProtonIdx(slc); //Leading proton
    if ( kScndProtonIdx(slc) < 0 ) return daT;
    unsigned int idxThird = (unsigned int) kScndProtonIdx(slc); //Charged pion or sub-leading proton

    double pMu_mag = kRecoMuonPNew(slc);
    double pP_mag = kRecoProtonP(slc);
    double pThird_mag = kScndProtonP(slc);

    TVector3 pMu(slc->reco.pfp.at(idxMuon).trk.truth.p.startp.x, slc->reco.pfp.at(idxMuon).trk.truth.p.startp.y, slc->reco.pfp.at(idxMuon).trk.truth.p.startp.z);
    TVector3 pP(slc->reco.pfp.at(idxProton).trk.truth.p.startp.x, slc->reco.pfp.at(idxProton).trk.truth.p.startp.y, slc->reco.pfp.at(idxProton).trk.truth.p.startp.z);
    TVector3 pThird(slc->reco.pfp.at(idxThird).trk.truth.p.startp.x, slc->reco.pfp.at(idxThird).trk.truth.p.startp.y, slc->reco.pfp.at(idxThird).trk.truth.p.startp.z);

    pMu *= pMu_mag / pMu.Mag();
    pP *= pP_mag / pP.Mag();
    pThird *= pThird_mag / pThird.Mag();

    const auto& vtx = slc->vertex;
    TVector3 vec_vtx(vtx.x, vtx.y, vtx.z);
    TVector3 vec_numi_to_vtx = (dFromNuMI + vec_vtx).Unit();

    //pT = p - pL
    TVector3 pTMu = pMu - (pMu.Dot(vec_numi_to_vtx) * vec_numi_to_vtx);
    TVector3 pTP = pP - (pP.Dot(vec_numi_to_vtx) * vec_numi_to_vtx);
    TVector3 pTThird = pThird - (pThird.Dot(vec_numi_to_vtx) * vec_numi_to_vtx);
    TVector3 dpT = pTMu + pTP + pTThird;

    daT = TMath::ACos( -1. * (pTMu.Dot(dpT)) / (pTMu.Mag() * dpT.Mag()) );

    if ( isnan(daT) ) daT = -9999.;
    return daT;
  });

const Var kDeltaAlphaT_Pion([](const caf::SRSliceProxy* slc) -> float {
    float daT = -9999.;

    unsigned int idxMuon = (unsigned int) kRecoMuonIdx(slc);
    unsigned int idxProton = (unsigned int) kRecoProtonIdx(slc); //Leading proton
    unsigned int idxThird = (unsigned int) kSidebandPion(slc); //Charged pion or sub-leading proton

    double pMu_mag = kRecoMuonPNew(slc);
    double pP_mag = kRecoProtonP(slc);
    double pThird_mag = kSidebandPionP(slc);

    TVector3 pMu(pMu_mag*slc->reco.pfp.at(idxMuon).trk.dir.x, pMu_mag*slc->reco.pfp.at(idxMuon).trk.dir.y, pMu_mag*slc->reco.pfp.at(idxMuon).trk.dir.z);
    TVector3 pP(pP_mag*slc->reco.pfp.at(idxProton).trk.dir.x, pP_mag*slc->reco.pfp.at(idxProton).trk.dir.y, pP_mag*slc->reco.pfp.at(idxProton).trk.dir.z);
    TVector3 pThird(pThird_mag*slc->reco.pfp.at(idxThird).trk.dir.x, pThird_mag*slc->reco.pfp.at(idxThird).trk.dir.y, pThird_mag*slc->reco.pfp.at(idxThird).trk.dir.z);

    const auto& vtx = slc->vertex;
    TVector3 vec_vtx(vtx.x, vtx.y, vtx.z);
    TVector3 vec_numi_to_vtx = (dFromNuMI + vec_vtx).Unit();

    //pT = p - pL
    TVector3 pTMu = pMu - (pMu.Dot(vec_numi_to_vtx) * vec_numi_to_vtx);
    TVector3 pTP = pP - (pP.Dot(vec_numi_to_vtx) * vec_numi_to_vtx);
    TVector3 pTThird = pThird - (pThird.Dot(vec_numi_to_vtx) * vec_numi_to_vtx);
    TVector3 dpT = pTMu + pTP + pTThird;

    daT = TMath::ACos( -1. * (pTMu.Dot(dpT)) / (pTMu.Mag() * dpT.Mag()) );

    return daT;
  });

/*
const Var kDeltaAlphaT_ThreeP([](const caf::SRSliceProxy* slc) -> float {
    float daT = -9999.;

    unsigned int idxMuon = (unsigned int) kRecoMuonIdx(slc);
    unsigned int idxP1 = (unsigned int) kRecoProtonIdx(slc); //Leading proton
    unsigned int idxP2 = (unsigned int) kScndProtonIdx(slc); //Charged pion or sub-leading proton
    unsigned int idxP3 = (unsigned int) kSidebandProton(slc);

    double pMu_mag = kRecoMuonPNew(slc);
    double pP1_mag = kRecoProtonP(slc);
    double pP2_mag = kScndProtonP(slc);
    double pP3_mag = kSidebandProtonP(slc);

    TVector3 pMu(pMu_mag*slc->reco.pfp.at(idxMuon).trk.dir.x, pMu_mag*slc->reco.pfp.at(idxMuon).trk.dir.y, pMu_mag*slc->reco.pfp.at(idxMuon).trk.dir.z);
    TVector3 pP1(pP1_mag*slc->reco.pfp.at(idxP1).trk.dir.x, pP1_mag*slc->reco.pfp.at(idxP1).trk.dir.y, pP1_mag*slc->reco.pfp.at(idxP1).trk.dir.z);
    TVector3 pP2(pP2_mag*slc->reco.pfp.at(idxP2).trk.dir.x, pP2_mag*slc->reco.pfp.at(idxP2).trk.dir.y, pP2_mag*slc->reco.pfp.at(idxP2).trk.dir.z);
    TVector3 pP3(pP3_mag*slc->reco.pfp.at(idxP3).trk.dir.x, pP3_mag*slc->reco.pfp.at(idxP3).trk.dir.y, pP3_mag*slc->reco.pfp.at(idxP3).trk.dir.z);

    const auto& vtx = slc->vertex;
    TVector3 vec_vtx(vtx.x, vtx.y, vtx.z);
    TVector3 vec_numi_to_vtx = (dFromNuMI + vec_vtx).Unit();

    //pT = p - pL
    TVector3 pTMu = pMu - (pMu.Dot(vec_numi_to_vtx) * vec_numi_to_vtx);
    TVector3 pTP1 = pP1 - (pP1.Dot(vec_numi_to_vtx) * vec_numi_to_vtx);
    TVector3 pTP2 = pP2 - (pP2.Dot(vec_numi_to_vtx) * vec_numi_to_vtx);
    TVector3 pTP3 = pP3 - (pP3.Dot(vec_numi_to_vtx) * vec_numi_to_vtx);
    TVector3 dpT = pTMu + pTP1 + pTP2 + pTP3;

    daT = TMath::ACos( -1. * (pTMu.Dot(dpT)) / (pTMu.Mag() * dpT.Mag()) );

    return daT;
  });
*/

const Var kDeltaAlphaT_ProtonTruth([](const caf::SRSliceProxy* slc) -> float {
    float daT = -9999.;

    if ( slc->truth.index < 0 ) return daT;
    TVector3 NuDirection_Truth = (TVector3(slc->truth.momentum.x, slc->truth.momentum.y, slc->truth.momentum.z)).Unit().Unit();

    unsigned int idxMuon = (unsigned int) kRecoMuonIdx(slc);
    unsigned int idxProton = (unsigned int) kRecoProtonIdx(slc); //Leading proton
    unsigned int idxThird = (unsigned int) kScndProtonIdx(slc); //Charged pion or sub-leading proton

    TVector3 pMu(slc->reco.pfp.at(idxMuon).trk.truth.p.startp.x, slc->reco.pfp.at(idxMuon).trk.truth.p.startp.y, slc->reco.pfp.at(idxMuon).trk.truth.p.startp.z);
    TVector3 pP(slc->reco.pfp.at(idxProton).trk.truth.p.startp.x, slc->reco.pfp.at(idxProton).trk.truth.p.startp.y, slc->reco.pfp.at(idxProton).trk.truth.p.startp.z);
    TVector3 pThird(slc->reco.pfp.at(idxThird).trk.truth.p.startp.x, slc->reco.pfp.at(idxThird).trk.truth.p.startp.y, slc->reco.pfp.at(idxThird).trk.truth.p.startp.z);

    //pT = p - pL
    TVector3 pTMu = pMu - (pMu.Dot(NuDirection_Truth) * NuDirection_Truth);
    TVector3 pTP = pP - (pP.Dot(NuDirection_Truth) * NuDirection_Truth);
    TVector3 pTThird = pThird - (pThird.Dot(NuDirection_Truth) * NuDirection_Truth);
    TVector3 dpT = pTMu + pTP + pTThird;

    daT = TMath::ACos( -1. * (pTMu.Dot(dpT)) / (pTMu.Mag() * dpT.Mag()) );

    return daT;
  });

const Var kDeltaAlphaT_PionTruth([](const caf::SRSliceProxy* slc) -> float {
    float daT = -9999.;

    if ( slc->truth.index < 0 ) return daT;
    TVector3 NuDirection_Truth = (TVector3(slc->truth.momentum.x, slc->truth.momentum.y, slc->truth.momentum.z)).Unit().Unit();

    unsigned int idxMuon = (unsigned int) kRecoMuonIdx(slc);
    unsigned int idxProton = (unsigned int) kRecoProtonIdx(slc); //Leading proton
    unsigned int idxThird = (unsigned int) kSidebandPion(slc); //Charged pion or sub-leading proton

    TVector3 pMu(slc->reco.pfp.at(idxMuon).trk.truth.p.startp.x, slc->reco.pfp.at(idxMuon).trk.truth.p.startp.y, slc->reco.pfp.at(idxMuon).trk.truth.p.startp.z);
    TVector3 pP(slc->reco.pfp.at(idxProton).trk.truth.p.startp.x, slc->reco.pfp.at(idxProton).trk.truth.p.startp.y, slc->reco.pfp.at(idxProton).trk.truth.p.startp.z);
    TVector3 pThird(slc->reco.pfp.at(idxThird).trk.truth.p.startp.x, slc->reco.pfp.at(idxThird).trk.truth.p.startp.y, slc->reco.pfp.at(idxThird).trk.truth.p.startp.z);

    //pT = p - pL
    TVector3 pTMu = pMu - (pMu.Dot(NuDirection_Truth) * NuDirection_Truth);
    TVector3 pTP = pP - (pP.Dot(NuDirection_Truth) * NuDirection_Truth);
    TVector3 pTThird = pThird - (pThird.Dot(NuDirection_Truth) * NuDirection_Truth);
    TVector3 dpT = pTMu + pTP + pTThird;

    daT = TMath::ACos( -1. * (pTMu.Dot(dpT)) / (pTMu.Mag() * dpT.Mag()) );

    return daT;
  });

/*
const Var kDeltaAlphaT_ThreePTruth([](const caf::SRSliceProxy* slc) -> float {
    float daT = -9999.;

    if ( slc->truth.index < 0 ) return daT;
    TVector3 NuDirection_Truth = (TVector3(slc->truth.momentum.x, slc->truth.momentum.y, slc->truth.momentum.z)).Unit().Unit();

    unsigned int idxMuon = (unsigned int) kRecoMuonIdx(slc);
    unsigned int idxP1 = (unsigned int) kRecoProtonIdx(slc); //Leading proton
    unsigned int idxP2 = (unsigned int) kScndProtonIdx(slc); //Charged pion or sub-leading proton
    unsigned int idxP3 = (unsigned int) kSidebandProton(slc);

    TVector3 pMu(slc->reco.pfp.at(idxMuon).trk.truth.p.startp.x, slc->reco.pfp.at(idxMuon).trk.truth.p.startp.y, slc->reco.pfp.at(idxMuon).trk.truth.p.startp.z);
    TVector3 pP1(slc->reco.pfp.at(idxP1).trk.truth.p.startp.x, slc->reco.pfp.at(idxP1).trk.truth.p.startp.y, slc->reco.pfp.at(idxP1).trk.truth.p.startp.z);
    TVector3 pP2(slc->reco.pfp.at(idxP2).trk.truth.p.startp.x, slc->reco.pfp.at(idxP2).trk.truth.p.startp.y, slc->reco.pfp.at(idxP2).trk.truth.p.startp.z);
    TVector3 pP3(slc->reco.pfp.at(idxP3).trk.truth.p.startp.x, slc->reco.pfp.at(idxP3).trk.truth.p.startp.y, slc->reco.pfp.at(idxP3).trk.truth.p.startp.z);

    //pT = p - pL
    TVector3 pTMu = pMu - (pMu.Dot(NuDirection_Truth) * NuDirection_Truth);
    TVector3 pTP1 = pP1 - (pP1.Dot(NuDirection_Truth) * NuDirection_Truth);
    TVector3 pTP2 = pP2 - (pP2.Dot(NuDirection_Truth) * NuDirection_Truth);
    TVector3 pTP3 = pP3 - (pP3.Dot(NuDirection_Truth) * NuDirection_Truth);
    TVector3 dpT = pTMu + pTP1 + pTP2 + pTP3;

    daT = TMath::ACos( -1. * (pTMu.Dot(dpT)) / (pTMu.Mag() * dpT.Mag()) );

    return daT;
  });
*/

const Var kDeltaAlphaT([](const caf::SRSliceProxy* slc) -> float {
    double daT = -9999.;
    if ( kScndProtonCandidate(slc) ) daT = kDeltaAlphaT_Proton(slc);
    //else if ( kThreePSideband(slc) ) daT = kDeltaAlphaT_ThreeP(slc);
    else if ( kPionSidebandBase(slc) ) daT = kDeltaAlphaT_Pion(slc);
    assert(((void)"No valid selection for DeltaAlphaT value", daT != -9999.));

    return daT;
  });

const Var kDeltaAlphaT_OldTruth([](const caf::SRSliceProxy* slc) -> float {
    double daT = -99999.; //Note the extra 9, as -9999. is a valid response for selected cosmic slices
    if ( kScndProtonCandidate(slc) ) daT = kDeltaAlphaT_ProtonTruth(slc);
    //else if ( kThreePSideband(slc) ) daT = kDeltaAlphaT_ThreePTruth(slc);
    else if ( kPionSidebandBase(slc) ) daT = kDeltaAlphaT_PionTruth(slc);
    assert(((void)"No valid selection for DeltaAlphaT value", daT != -99999.));

    return daT;
  });

//Built out of SRTrueInteraction rather than truth info matched to reco tracks
//Considers all charged pions and protons over threshold, so should be correct for all selections
const Var kDeltaAlphaT_FullTruth([](const caf::SRSliceProxy* slc) -> float {
    float daT = -9999.;
    if ( slc->truth.index < 0 ) return daT;
    const auto &nu = slc->truth;

    TVector3 NuDirection_Truth = (TVector3(nu.momentum.x, nu.momentum.y, nu.momentum.z)).Unit().Unit();
    TVector3 p3;
    TVector3 pTMu(0., 0., 0.);
    TVector3 dpT(0., 0., 0.);
 
    for ( unsigned int i = 0; i < nu.prim.size(); i++ ) {
      const auto &prim = nu.prim.at(i);
      if ( prim.start_process != 0 ) continue;
      p3 = {prim.startp.x, prim.startp.y, prim.startp.z};
      if ( p3.Mag() == std::hypot(-9999., -9999., -9999.) ) return daT;
      if ( abs(prim.pdg) == 13 ) {
        pTMu = p3 - (p3.Dot(NuDirection_Truth) * NuDirection_Truth);
        dpT += pTMu;
      }
      else if ( abs(prim.pdg) == 2212 && p3.Mag() > .35 ) dpT += p3 - (p3.Dot(NuDirection_Truth) * NuDirection_Truth);
      else if ( abs(prim.pdg) == 211 || abs(prim.pdg) == 321 ) dpT += p3 - (p3.Dot(NuDirection_Truth) * NuDirection_Truth);
    }

    if ( pTMu.Mag() == 0. ) return daT;
    daT = TMath::ACos( -1. * (pTMu.Dot(dpT)) / (pTMu.Mag() * dpT.Mag()) );

    return daT;
  });

const Var kDeltaAlphaT_Truth([](const caf::SRSliceProxy* slc) -> float {
    float daT = -9999.;
    if ( slc->truth.index < 0 ) return daT;
    const auto &nu = slc->truth;

    TVector3 NuDirection_Truth = (TVector3(nu.momentum.x, nu.momentum.y, nu.momentum.z)).Unit().Unit();
    TVector3 p3;
    TVector3 pTMu(0., 0., 0.);
    TVector3 dpT(0., 0., 0.);

    int idxP1 = -1;
    int idxP2 = -1;
    double leading = -9999.;
    double second = -9999.;

    for ( unsigned int i = 0; i < nu.prim.size(); i++ ) {
      const auto &prim = nu.prim.at(i);
      if ( prim.start_process != 0 ) continue;
      p3 = {prim.startp.x, prim.startp.y, prim.startp.z};
      if ( p3.Mag() == std::hypot(-9999., -9999., -9999.) ) return daT;
      if ( abs(prim.pdg) == 13 ) {
        pTMu = p3 - (p3.Dot(NuDirection_Truth) * NuDirection_Truth);
        dpT += pTMu;
      }
      else if ( abs(prim.pdg) == 2212 && p3.Mag() > .35 ) {
        if ( p3.Mag() > leading ) {
          second = leading;
          leading = p3.Mag();
          idxP2 = idxP1;
          idxP1 = i;
        }
        else if ( p3.Mag() > second ) {
          second = p3.Mag();
          idxP2 = i;
        }
      }
    }

    if ( idxP2 == -1 ) return daT;
    TVector3 P1 = {nu.prim.at(idxP1).startp.x, nu.prim.at(idxP1).startp.y, nu.prim.at(idxP1).startp.z};
    TVector3 P2 = {nu.prim.at(idxP2).startp.x, nu.prim.at(idxP2).startp.y, nu.prim.at(idxP2).startp.z};
    dpT += P1 - (P1.Dot(NuDirection_Truth) * NuDirection_Truth);
    dpT += P2 - (P2.Dot(NuDirection_Truth) * NuDirection_Truth);

    if ( pTMu.Mag() == 0. ) return daT;
    daT = TMath::ACos( -1. * (pTMu.Dot(dpT)) / (pTMu.Mag() * dpT.Mag()) );

    return daT;
  });

const Var kDeltaAlphaT_Resid([](const caf::SRSliceProxy* slc) -> float {
    double recoDeltaAlphaT = kDeltaAlphaT(slc);
    double trueDeltaAlphaT = kDeltaAlphaT_Truth(slc);

    if ( trueDeltaAlphaT == -9999. || trueDeltaAlphaT == 0. ) return -9999.;
    return (recoDeltaAlphaT - trueDeltaAlphaT) / trueDeltaAlphaT;
  });

const Var kDeltaAlphaT_ProtonResid([](const caf::SRSliceProxy* slc) -> float {
    double recoDeltaAlphaT = kDeltaAlphaT_Proton(slc);
    double trueDeltaAlphaT = kDeltaAlphaT_Truth(slc);

    if ( trueDeltaAlphaT == -9999. || trueDeltaAlphaT == 0. ) return -9999.;
    return (recoDeltaAlphaT - trueDeltaAlphaT) / trueDeltaAlphaT;
  });

//Temporary, use true muon momentum to characterize the impact of getting it wrong
const Var kDeltaAlphaT_CheatingMuonResid([](const caf::SRSliceProxy* slc) -> float {
    double recoDeltaAlphaT = kDeltaAlphaT_CheatingMuon(slc);
    double trueDeltaAlphaT = kDeltaAlphaT_Truth(slc);

    if ( trueDeltaAlphaT == -9999. || trueDeltaAlphaT == 0. ) return -9999.;
    return (recoDeltaAlphaT - trueDeltaAlphaT) / trueDeltaAlphaT;
  });

const Var kDeltaAlphaT_CheatingAnglesResid([](const caf::SRSliceProxy* slc) -> float {
    double recoDeltaAlphaT = kDeltaAlphaT_CheatingAngles(slc);
    double trueDeltaAlphaT = kDeltaAlphaT_Truth(slc);

    if ( trueDeltaAlphaT == -9999. || trueDeltaAlphaT == 0. ) return -9999.;
    return (recoDeltaAlphaT - trueDeltaAlphaT) / trueDeltaAlphaT;
  });

//Transverse opening angle, "measures the deflection of N [hadronic system] with respect to q in the transverse plane ... depends on the lepton kinematics which are sensitive to the neutrino energy"
const Var kDeltaPhiT_Single([](const caf::SRSliceProxy* slc) -> float {
    float dphiT = -9999.;

    unsigned int idxMuon = (unsigned int) kRecoMuonIdx(slc);
    unsigned int idxProton = (unsigned int) kRecoProtonIdx(slc); //Leading proton

    double pMu_mag = kRecoMuonPNew(slc);
    double pP_mag = kRecoProtonP(slc);

    TVector3 pMu(pMu_mag*slc->reco.pfp.at(idxMuon).trk.dir.x, pMu_mag*slc->reco.pfp.at(idxMuon).trk.dir.y, pMu_mag*slc->reco.pfp.at(idxMuon).trk.dir.z);
    TVector3 pP(pP_mag*slc->reco.pfp.at(idxProton).trk.dir.x, pP_mag*slc->reco.pfp.at(idxProton).trk.dir.y, pP_mag*slc->reco.pfp.at(idxProton).trk.dir.z);

    const auto& vtx = slc->vertex;
    TVector3 vec_vtx(vtx.x, vtx.y, vtx.z);
    TVector3 vec_numi_to_vtx = (dFromNuMI + vec_vtx).Unit();

    //pT = p - pL
    TVector3 pTMu = pMu - (pMu.Dot(vec_numi_to_vtx) * vec_numi_to_vtx);
    TVector3 pTP = pP - (pP.Dot(vec_numi_to_vtx) * vec_numi_to_vtx);

    dphiT = TMath::ACos( -1. * (pTMu.Dot(pTP)) / (pTMu.Mag() * pTP.Mag()) );

    return dphiT;
  });

const Var kDeltaPhiT_Proton([](const caf::SRSliceProxy* slc) -> float {
    float dphiT = -9999.;

    unsigned int idxMuon = (unsigned int) kRecoMuonIdx(slc);
    unsigned int idxProton = (unsigned int) kRecoProtonIdx(slc); //Leading proton
    unsigned int idxThird = (unsigned int) kScndProtonIdx(slc); //Charged pion or sub-leading proton

    double pMu_mag = kRecoMuonPNew(slc);
    double pP_mag = kRecoProtonP(slc);
    double pThird_mag = kScndProtonP(slc);

    TVector3 pMu(pMu_mag*slc->reco.pfp.at(idxMuon).trk.dir.x, pMu_mag*slc->reco.pfp.at(idxMuon).trk.dir.y, pMu_mag*slc->reco.pfp.at(idxMuon).trk.dir.z);
    TVector3 pP(pP_mag*slc->reco.pfp.at(idxProton).trk.dir.x, pP_mag*slc->reco.pfp.at(idxProton).trk.dir.y, pP_mag*slc->reco.pfp.at(idxProton).trk.dir.z);
    TVector3 pThird(pThird_mag*slc->reco.pfp.at(idxThird).trk.dir.x, pThird_mag*slc->reco.pfp.at(idxThird).trk.dir.y, pThird_mag*slc->reco.pfp.at(idxThird).trk.dir.z);

    const auto& vtx = slc->vertex;
    TVector3 vec_vtx(vtx.x, vtx.y, vtx.z);
    TVector3 vec_numi_to_vtx = (dFromNuMI + vec_vtx).Unit();

    //pT = p - pL
    TVector3 pTMu = pMu - (pMu.Dot(vec_numi_to_vtx) * vec_numi_to_vtx);
    TVector3 pTP = pP - (pP.Dot(vec_numi_to_vtx) * vec_numi_to_vtx);
    TVector3 pTThird = pThird - (pThird.Dot(vec_numi_to_vtx) * vec_numi_to_vtx);
    TVector3 pTH = pTP + pTThird;

    dphiT = TMath::ACos( -1. * (pTMu.Dot(pTH)) / (pTMu.Mag() * pTH.Mag()) );

    return dphiT;
  });

const Var kDeltaPhiT_CheatingMuon([](const caf::SRSliceProxy* slc) -> float {
    float dphiT = -9999.;

    unsigned int idxMuon = (unsigned int) kRecoMuonIdx(slc);
    unsigned int idxProton = (unsigned int) kRecoProtonIdx(slc); //Leading proton
    if ( kScndProtonIdx(slc) < 0 ) return dphiT;
    unsigned int idxThird = (unsigned int) kScndProtonIdx(slc); //Charged pion or sub-leading proton

    double pP_mag = kRecoProtonP(slc);
    double pThird_mag = kScndProtonP(slc);

    TVector3 pMu(slc->reco.pfp.at(idxMuon).trk.truth.p.startp.x, slc->reco.pfp.at(idxMuon).trk.truth.p.startp.y, slc->reco.pfp.at(idxMuon).trk.truth.p.startp.z);
    TVector3 pP(pP_mag*slc->reco.pfp.at(idxProton).trk.dir.x, pP_mag*slc->reco.pfp.at(idxProton).trk.dir.y, pP_mag*slc->reco.pfp.at(idxProton).trk.dir.z);
    TVector3 pThird(pThird_mag*slc->reco.pfp.at(idxThird).trk.dir.x, pThird_mag*slc->reco.pfp.at(idxThird).trk.dir.y, pThird_mag*slc->reco.pfp.at(idxThird).trk.dir.z);

    const auto& vtx = slc->vertex;
    TVector3 vec_vtx(vtx.x, vtx.y, vtx.z);
    TVector3 vec_numi_to_vtx = (dFromNuMI + vec_vtx).Unit();

    //pT = p - pL
    TVector3 pTMu = pMu - (pMu.Dot(vec_numi_to_vtx) * vec_numi_to_vtx);
    TVector3 pTP = pP - (pP.Dot(vec_numi_to_vtx) * vec_numi_to_vtx);
    TVector3 pTThird = pThird - (pThird.Dot(vec_numi_to_vtx) * vec_numi_to_vtx);
    TVector3 pTH = pTP + pTThird;

    dphiT = TMath::ACos( -1. * (pTMu.Dot(pTH)) / (pTMu.Mag() * pTH.Mag()) );

    return dphiT;
  });

const Var kDeltaPhiT_CheatingAngles([](const caf::SRSliceProxy* slc) -> float {
    float dphiT = -9999.;

    unsigned int idxMuon = (unsigned int) kRecoMuonIdx(slc);
    unsigned int idxProton = (unsigned int) kRecoProtonIdx(slc); //Leading proton
    if ( kScndProtonIdx(slc) < 0 ) return dphiT;
    unsigned int idxThird = (unsigned int) kScndProtonIdx(slc); //Charged pion or sub-leading proton

    double pMu_mag = kRecoMuonPNew(slc);
    double pP_mag = kRecoProtonP(slc);
    double pThird_mag = kScndProtonP(slc);

    TVector3 pMu(slc->reco.pfp.at(idxMuon).trk.truth.p.startp.x, slc->reco.pfp.at(idxMuon).trk.truth.p.startp.y, slc->reco.pfp.at(idxMuon).trk.truth.p.startp.z);
    TVector3 pP(slc->reco.pfp.at(idxProton).trk.truth.p.startp.x, slc->reco.pfp.at(idxProton).trk.truth.p.startp.y, slc->reco.pfp.at(idxProton).trk.truth.p.startp.z);
    TVector3 pThird(slc->reco.pfp.at(idxThird).trk.truth.p.startp.x, slc->reco.pfp.at(idxThird).trk.truth.p.startp.y, slc->reco.pfp.at(idxThird).trk.truth.p.startp.z);

    pMu *= pMu_mag / pMu.Mag();
    pP *= pP_mag / pP.Mag();
    pThird *= pThird_mag / pThird.Mag();

    const auto& vtx = slc->vertex;
    TVector3 vec_vtx(vtx.x, vtx.y, vtx.z);
    TVector3 vec_numi_to_vtx = (dFromNuMI + vec_vtx).Unit();

    //pT = p - pL
    TVector3 pTMu = pMu - (pMu.Dot(vec_numi_to_vtx) * vec_numi_to_vtx);
    TVector3 pTP = pP - (pP.Dot(vec_numi_to_vtx) * vec_numi_to_vtx);
    TVector3 pTThird = pThird - (pThird.Dot(vec_numi_to_vtx) * vec_numi_to_vtx);
    TVector3 pTH = pTP + pTThird;

    dphiT = TMath::ACos( -1. * (pTMu.Dot(pTH)) / (pTMu.Mag() * pTH.Mag()) );

    if ( isnan(dphiT) ) dphiT = -9999.;
    return dphiT;
  });

const Var kDeltaPhiT_Pion([](const caf::SRSliceProxy* slc) -> float {
    float dphiT = -9999.;

    unsigned int idxMuon = (unsigned int) kRecoMuonIdx(slc);
    unsigned int idxProton = (unsigned int) kRecoProtonIdx(slc); //Leading proton
    unsigned int idxThird = (unsigned int) kSidebandPion(slc); //Charged pion or sub-leading proton

    double pMu_mag = kRecoMuonPNew(slc);
    double pP_mag = kRecoProtonP(slc);
    double pThird_mag = kSidebandPionP(slc);

    TVector3 pMu(pMu_mag*slc->reco.pfp.at(idxMuon).trk.dir.x, pMu_mag*slc->reco.pfp.at(idxMuon).trk.dir.y, pMu_mag*slc->reco.pfp.at(idxMuon).trk.dir.z);
    TVector3 pP(pP_mag*slc->reco.pfp.at(idxProton).trk.dir.x, pP_mag*slc->reco.pfp.at(idxProton).trk.dir.y, pP_mag*slc->reco.pfp.at(idxProton).trk.dir.z);
    TVector3 pThird(pThird_mag*slc->reco.pfp.at(idxThird).trk.dir.x, pThird_mag*slc->reco.pfp.at(idxThird).trk.dir.y, pThird_mag*slc->reco.pfp.at(idxThird).trk.dir.z);

    const auto& vtx = slc->vertex;
    TVector3 vec_vtx(vtx.x, vtx.y, vtx.z);
    TVector3 vec_numi_to_vtx = (dFromNuMI + vec_vtx).Unit();

    //pT = p - pL
    TVector3 pTMu = pMu - (pMu.Dot(vec_numi_to_vtx) * vec_numi_to_vtx);
    TVector3 pTP = pP - (pP.Dot(vec_numi_to_vtx) * vec_numi_to_vtx);
    TVector3 pTThird = pThird - (pThird.Dot(vec_numi_to_vtx) * vec_numi_to_vtx);
    TVector3 pTH = pTP + pTThird;

    dphiT = TMath::ACos( -1. * (pTMu.Dot(pTH)) / (pTMu.Mag() * pTH.Mag()) );

    return dphiT;
  });

/*
const Var kDeltaPhiT_ThreeP([](const caf::SRSliceProxy* slc) -> float {
    float dphiT = -9999.;

    unsigned int idxMuon = (unsigned int) kRecoMuonIdx(slc);
    unsigned int idxP1 = (unsigned int) kRecoProtonIdx(slc); //Leading proton
    unsigned int idxP2 = (unsigned int) kScndProtonIdx(slc); //Charged pion or sub-leading proton
    unsigned int idxP3 = (unsigned int) kSidebandProton(slc);

    double pMu_mag = kRecoMuonPNew(slc);
    double pP1_mag = kRecoProtonP(slc);
    double pP2_mag = kScndProtonP(slc);
    double pP3_mag = kSidebandProtonP(slc);

    TVector3 pMu(pMu_mag*slc->reco.pfp.at(idxMuon).trk.dir.x, pMu_mag*slc->reco.pfp.at(idxMuon).trk.dir.y, pMu_mag*slc->reco.pfp.at(idxMuon).trk.dir.z);
    TVector3 pP1(pP1_mag*slc->reco.pfp.at(idxP1).trk.dir.x, pP1_mag*slc->reco.pfp.at(idxP1).trk.dir.y, pP1_mag*slc->reco.pfp.at(idxP1).trk.dir.z);
    TVector3 pP2(pP2_mag*slc->reco.pfp.at(idxP2).trk.dir.x, pP2_mag*slc->reco.pfp.at(idxP2).trk.dir.y, pP2_mag*slc->reco.pfp.at(idxP2).trk.dir.z);
    TVector3 pP3(pP3_mag*slc->reco.pfp.at(idxP3).trk.dir.x, pP3_mag*slc->reco.pfp.at(idxP3).trk.dir.y, pP3_mag*slc->reco.pfp.at(idxP3).trk.dir.z);

    const auto& vtx = slc->vertex;
    TVector3 vec_vtx(vtx.x, vtx.y, vtx.z);
    TVector3 vec_numi_to_vtx = (dFromNuMI + vec_vtx).Unit();

    //pT = p - pL
    TVector3 pTMu = pMu - (pMu.Dot(vec_numi_to_vtx) * vec_numi_to_vtx);
    TVector3 pTP1 = pP1 - (pP1.Dot(vec_numi_to_vtx) * vec_numi_to_vtx);
    TVector3 pTP2 = pP2 - (pP2.Dot(vec_numi_to_vtx) * vec_numi_to_vtx);
    TVector3 pTP3 = pP3 - (pP3.Dot(vec_numi_to_vtx) * vec_numi_to_vtx);
    TVector3 pTH = pTP1 + pTP2 + pTP3;

    dphiT = TMath::ACos( -1. * (pTMu.Dot(pTH)) / (pTMu.Mag() * pTH.Mag()) );

    return dphiT;
  });
*/

const Var kDeltaPhiT_ProtonTruth([](const caf::SRSliceProxy* slc) -> float {
    float dphiT = -9999.;

    if ( slc->truth.index < 0 ) return dphiT;
    TVector3 NuDirection_Truth = (TVector3(slc->truth.momentum.x, slc->truth.momentum.y, slc->truth.momentum.z)).Unit().Unit();

    unsigned int idxMuon = (unsigned int) kRecoMuonIdx(slc);
    unsigned int idxProton = (unsigned int) kRecoProtonIdx(slc); //Leading proton
    unsigned int idxThird = (unsigned int) kScndProtonIdx(slc); //Charged pion or sub-leading proton

    TVector3 pMu(slc->reco.pfp.at(idxMuon).trk.truth.p.startp.x, slc->reco.pfp.at(idxMuon).trk.truth.p.startp.y, slc->reco.pfp.at(idxMuon).trk.truth.p.startp.z);
    TVector3 pP(slc->reco.pfp.at(idxProton).trk.truth.p.startp.x, slc->reco.pfp.at(idxProton).trk.truth.p.startp.y, slc->reco.pfp.at(idxProton).trk.truth.p.startp.z);
    TVector3 pThird(slc->reco.pfp.at(idxThird).trk.truth.p.startp.x, slc->reco.pfp.at(idxThird).trk.truth.p.startp.y, slc->reco.pfp.at(idxThird).trk.truth.p.startp.z);

    //pT = p - pL
    TVector3 pTMu = pMu - (pMu.Dot(NuDirection_Truth) * NuDirection_Truth);
    TVector3 pTP = pP - (pP.Dot(NuDirection_Truth) * NuDirection_Truth);
    TVector3 pTThird = pThird - (pThird.Dot(NuDirection_Truth) * NuDirection_Truth);
    TVector3 pTH = pTP + pTThird;

    dphiT = TMath::ACos( -1. * (pTMu.Dot(pTH)) / (pTMu.Mag() * pTH.Mag()) );

    return dphiT;
  });

const Var kDeltaPhiT_PionTruth([](const caf::SRSliceProxy* slc) -> float {
    float dphiT = -9999.;

    if ( slc->truth.index < 0 ) return dphiT;
    TVector3 NuDirection_Truth = (TVector3(slc->truth.momentum.x, slc->truth.momentum.y, slc->truth.momentum.z)).Unit().Unit();

    unsigned int idxMuon = (unsigned int) kRecoMuonIdx(slc);
    unsigned int idxProton = (unsigned int) kRecoProtonIdx(slc); //Leading proton
    unsigned int idxThird = (unsigned int) kSidebandPion(slc); //Charged pion or sub-leading proton

    TVector3 pMu(slc->reco.pfp.at(idxMuon).trk.truth.p.startp.x, slc->reco.pfp.at(idxMuon).trk.truth.p.startp.y, slc->reco.pfp.at(idxMuon).trk.truth.p.startp.z);
    TVector3 pP(slc->reco.pfp.at(idxProton).trk.truth.p.startp.x, slc->reco.pfp.at(idxProton).trk.truth.p.startp.y, slc->reco.pfp.at(idxProton).trk.truth.p.startp.z);
    TVector3 pThird(slc->reco.pfp.at(idxThird).trk.truth.p.startp.x, slc->reco.pfp.at(idxThird).trk.truth.p.startp.y, slc->reco.pfp.at(idxThird).trk.truth.p.startp.z);

    //pT = p - pL
    TVector3 pTMu = pMu - (pMu.Dot(NuDirection_Truth) * NuDirection_Truth);
    TVector3 pTP = pP - (pP.Dot(NuDirection_Truth) * NuDirection_Truth);
    TVector3 pTThird = pThird - (pThird.Dot(NuDirection_Truth) * NuDirection_Truth);
    TVector3 pTH = pTP + pTThird;

    dphiT = TMath::ACos( -1. * (pTMu.Dot(pTH)) / (pTMu.Mag() * pTH.Mag()) );

    return dphiT;
  });

/*
const Var kDeltaPhiT_ThreePTruth([](const caf::SRSliceProxy* slc) -> float {
    float dphiT = -9999.;

    if ( slc->truth.index < 0 ) return dphiT;
    TVector3 NuDirection_Truth = (TVector3(slc->truth.momentum.x, slc->truth.momentum.y, slc->truth.momentum.z)).Unit().Unit();

    unsigned int idxMuon = (unsigned int) kRecoMuonIdx(slc);
    unsigned int idxP1 = (unsigned int) kRecoProtonIdx(slc); //Leading proton
    unsigned int idxP2 = (unsigned int) kScndProtonIdx(slc); //Charged pion or sub-leading proton
    unsigned int idxP3 = (unsigned int) kSidebandProton(slc);

    TVector3 pMu(slc->reco.pfp.at(idxMuon).trk.truth.p.startp.x, slc->reco.pfp.at(idxMuon).trk.truth.p.startp.y, slc->reco.pfp.at(idxMuon).trk.truth.p.startp.z);
    TVector3 pP1(slc->reco.pfp.at(idxP1).trk.truth.p.startp.x, slc->reco.pfp.at(idxP1).trk.truth.p.startp.y, slc->reco.pfp.at(idxP1).trk.truth.p.startp.z);
    TVector3 pP2(slc->reco.pfp.at(idxP2).trk.truth.p.startp.x, slc->reco.pfp.at(idxP2).trk.truth.p.startp.y, slc->reco.pfp.at(idxP2).trk.truth.p.startp.z);
    TVector3 pP3(slc->reco.pfp.at(idxP3).trk.truth.p.startp.x, slc->reco.pfp.at(idxP3).trk.truth.p.startp.y, slc->reco.pfp.at(idxP3).trk.truth.p.startp.z);

    //pT = p - pL
    TVector3 pTMu = pMu - (pMu.Dot(NuDirection_Truth) * NuDirection_Truth);
    TVector3 pTP1 = pP1 - (pP1.Dot(NuDirection_Truth) * NuDirection_Truth);
    TVector3 pTP2 = pP2 - (pP2.Dot(NuDirection_Truth) * NuDirection_Truth);
    TVector3 pTP3 = pP3 - (pP3.Dot(NuDirection_Truth) * NuDirection_Truth);
    TVector3 pTH = pTP1 + pTP2 + pTP3;

    dphiT = TMath::ACos( -1. * (pTMu.Dot(pTH)) / (pTMu.Mag() * pTH.Mag()) );

    return dphiT;
  });
*/

const Var kDeltaPhiT([](const caf::SRSliceProxy* slc) -> float {
    double dphiT = -9999.;
    if ( kScndProtonCandidate(slc) ) dphiT = kDeltaPhiT_Proton(slc);
    //else if ( kThreePSideband(slc) ) dphiT = kDeltaPhiT_ThreeP(slc);
    else if ( kPionSidebandBase(slc) ) dphiT = kDeltaPhiT_Pion(slc);
    assert(((void)"No valid selection for DeltaPhiT value", dphiT != -9999.));

    return dphiT;
  });

const Var kDeltaPhiT_OldTruth([](const caf::SRSliceProxy* slc) -> float {
    double dphiT = -99999.; //Note the extra 9, as -9999. is a valid response for selected cosmic slices
    if ( kScndProtonCandidate(slc) ) dphiT = kDeltaPhiT_ProtonTruth(slc);
    //else if ( kThreePSideband(slc) ) dphiT = kDeltaPhiT_ThreePTruth(slc);
    else if ( kPionSidebandBase(slc) ) dphiT = kDeltaPhiT_PionTruth(slc);
    assert(((void)"No valid selection for DeltaPhiT value", dphiT != -99999.));

    return dphiT;
  });

//Built out of SRTrueInteraction rather than truth info matched to reco tracks
//Considers all charged pions and protons over threshold, so should be correct for all selections
const Var kDeltaPhiT_FullTruth([](const caf::SRSliceProxy* slc) -> float {
    float dphiT = -9999.;
    if ( slc->truth.index < 0 ) return dphiT;
    const auto &nu = slc->truth;

    TVector3 NuDirection_Truth = (TVector3(nu.momentum.x, nu.momentum.y, nu.momentum.z)).Unit().Unit();
    TVector3 p3;
    TVector3 pTMu(0., 0., 0.);
    TVector3 pTH(0., 0., 0.);
    
    for ( unsigned int i = 0; i < nu.prim.size(); i++ ) {
      const auto &prim = nu.prim.at(i);
      if ( prim.start_process != 0 ) continue;
      p3 = {prim.startp.x, prim.startp.y, prim.startp.z};
      if ( p3.Mag() == std::hypot(-9999., -9999., -9999.) ) return dphiT;
      if ( abs(prim.pdg) == 13 ) {
        pTMu = p3 - (p3.Dot(NuDirection_Truth) * NuDirection_Truth);
      }
      else if ( abs(prim.pdg) == 2212 && p3.Mag() > .35 ) pTH += p3 - (p3.Dot(NuDirection_Truth) * NuDirection_Truth);
      else if ( abs(prim.pdg) == 211 || abs(prim.pdg) == 321 ) pTH += p3 - (p3.Dot(NuDirection_Truth) * NuDirection_Truth);
    }

    if ( pTMu.Mag() == 0. ) return dphiT;
    if ( pTH.Mag() == 0. ) { /*std::cout << "pTH == 0. !" << std::endl;*/ return dphiT; }
    dphiT = TMath::ACos( -1. * (pTMu.Dot(pTH)) / (pTMu.Mag() * pTH.Mag()) );

    return dphiT;
  });

const Var kDeltaPhiT_Truth([](const caf::SRSliceProxy* slc) -> float {
    float dphiT = -9999.;
    if ( slc->truth.index < 0 ) return dphiT;
    const auto &nu = slc->truth;

    TVector3 NuDirection_Truth = (TVector3(nu.momentum.x, nu.momentum.y, nu.momentum.z)).Unit().Unit();
    TVector3 p3;
    TVector3 pTMu(0., 0., 0.);
    TVector3 pTH(0., 0., 0.);

    int idxP1 = -1;
    int idxP2 = -1;
    double leading = -9999.;
    double second = -9999.;

    for ( unsigned int i = 0; i < nu.prim.size(); i++ ) {
      const auto &prim = nu.prim.at(i);
      if ( prim.start_process != 0 ) continue;
      p3 = {prim.startp.x, prim.startp.y, prim.startp.z};
      if ( p3.Mag() == std::hypot(-9999., -9999., -9999.) ) return dphiT;
      if ( abs(prim.pdg) == 13 ) {
        pTMu = p3 - (p3.Dot(NuDirection_Truth) * NuDirection_Truth);
      }
      else if ( abs(prim.pdg) == 2212 && p3.Mag() > .35 ) {
        if ( p3.Mag() > leading ) {
          second = leading;
          leading = p3.Mag();
          idxP2 = idxP1;
          idxP1 = i;
        }
        else if ( p3.Mag() > second ) {
          second = p3.Mag();
          idxP2 = i;
        }
      }
    }

    if ( idxP2 == -1 ) return dphiT;
    TVector3 P1 = {nu.prim.at(idxP1).startp.x, nu.prim.at(idxP1).startp.y, nu.prim.at(idxP1).startp.z};
    TVector3 P2 = {nu.prim.at(idxP2).startp.x, nu.prim.at(idxP2).startp.y, nu.prim.at(idxP2).startp.z};
    pTH += P1 - (P1.Dot(NuDirection_Truth) * NuDirection_Truth);
    pTH += P2 - (P2.Dot(NuDirection_Truth) * NuDirection_Truth);

    if ( pTMu.Mag() == 0. ||  pTH.Mag() == 0. ) return dphiT;
    dphiT = TMath::ACos( -1. * (pTMu.Dot(pTH)) / (pTMu.Mag() * pTH.Mag()) );

    return dphiT;
  });

const Var kDeltaPhiT_Resid([](const caf::SRSliceProxy* slc) -> float {
    double recoDeltaPhiT = kDeltaPhiT(slc);
    double trueDeltaPhiT = kDeltaPhiT_Truth(slc);

    if ( trueDeltaPhiT == -9999. || trueDeltaPhiT == 0. ) return -9999.;
    return (recoDeltaPhiT - trueDeltaPhiT) / trueDeltaPhiT;
  });

const Var kDeltaPhiT_ProtonResid([](const caf::SRSliceProxy* slc) -> float {
    double recoDeltaPhiT = kDeltaPhiT_Proton(slc);
    double trueDeltaPhiT = kDeltaPhiT_Truth(slc);

    if ( trueDeltaPhiT == -9999. || trueDeltaPhiT == 0. ) return -9999.;
    return (recoDeltaPhiT - trueDeltaPhiT) / trueDeltaPhiT;
  });

const Var kDeltaPhiT_CheatingMuonResid([](const caf::SRSliceProxy* slc) -> float {
    double recoDeltaPhiT = kDeltaPhiT_CheatingMuon(slc);
    double trueDeltaPhiT = kDeltaPhiT_Truth(slc);

    if ( trueDeltaPhiT == -9999. || trueDeltaPhiT == 0. ) return -9999.;
    return (recoDeltaPhiT - trueDeltaPhiT) / trueDeltaPhiT;
  });

const Var kDeltaPhiT_CheatingAnglesResid([](const caf::SRSliceProxy* slc) -> float {
    double recoDeltaPhiT = kDeltaPhiT_CheatingAngles(slc);
    double trueDeltaPhiT = kDeltaPhiT_Truth(slc);

    if ( trueDeltaPhiT == -9999. || trueDeltaPhiT == 0. ) return -9999.;
    return (recoDeltaPhiT - trueDeltaPhiT) / trueDeltaPhiT;
  });

//Double-transverse momentum imbalance
const Var kDeltaPTT_Single([](const caf::SRSliceProxy* slc) -> float {
    float dpTT = -9999.;

    unsigned int idxMuon = (unsigned int) kRecoMuonIdx(slc);
    unsigned int idxProton = (unsigned int) kRecoProtonIdx(slc); //Leading proton

    double pMu_mag = kRecoMuonPNew(slc);
    double pP_mag = kRecoProtonP(slc);

    TVector3 pMu(pMu_mag*slc->reco.pfp.at(idxMuon).trk.dir.x, pMu_mag*slc->reco.pfp.at(idxMuon).trk.dir.y, pMu_mag*slc->reco.pfp.at(idxMuon).trk.dir.z);
    TVector3 pP(pP_mag*slc->reco.pfp.at(idxProton).trk.dir.x, pP_mag*slc->reco.pfp.at(idxProton).trk.dir.y, pP_mag*slc->reco.pfp.at(idxProton).trk.dir.z);

    const auto& vtx = slc->vertex;
    TVector3 vec_vtx(vtx.x, vtx.y, vtx.z);
    TVector3 vec_numi_to_vtx = (dFromNuMI + vec_vtx).Unit();

    TVector3 doubleTransverse = ( vec_numi_to_vtx.Cross(pMu) ).Unit();

    double pTTP = pP.Dot(doubleTransverse);
    dpTT = pTTP;

    return dpTT;
  });

const Var kDeltaPTT_Proton([](const caf::SRSliceProxy* slc) -> float {
    float dpTT = -9999.;

    unsigned int idxMuon = (unsigned int) kRecoMuonIdx(slc);
    unsigned int idxProton = (unsigned int) kRecoProtonIdx(slc); //Leading proton
    unsigned int idxThird = (unsigned int) kScndProtonIdx(slc); //Charged pion or sub-leading proton

    double pMu_mag = kRecoMuonPNew(slc);
    double pP_mag = kRecoProtonP(slc);
    double pThird_mag = kScndProtonP(slc);

    TVector3 pMu(pMu_mag*slc->reco.pfp.at(idxMuon).trk.dir.x, pMu_mag*slc->reco.pfp.at(idxMuon).trk.dir.y, pMu_mag*slc->reco.pfp.at(idxMuon).trk.dir.z);
    TVector3 pP(pP_mag*slc->reco.pfp.at(idxProton).trk.dir.x, pP_mag*slc->reco.pfp.at(idxProton).trk.dir.y, pP_mag*slc->reco.pfp.at(idxProton).trk.dir.z);
    TVector3 pThird(pThird_mag*slc->reco.pfp.at(idxThird).trk.dir.x, pThird_mag*slc->reco.pfp.at(idxThird).trk.dir.y, pThird_mag*slc->reco.pfp.at(idxThird).trk.dir.z);

    const auto& vtx = slc->vertex;
    TVector3 vec_vtx(vtx.x, vtx.y, vtx.z);
    TVector3 vec_numi_to_vtx = (dFromNuMI + vec_vtx).Unit();

    TVector3 doubleTransverse = ( vec_numi_to_vtx.Cross(pMu) ).Unit();

    double pTTP = pP.Dot(doubleTransverse);
    double pTTThird = pThird.Dot(doubleTransverse);
    dpTT = pTTP + pTTThird;

    return dpTT;
  });

const Var kDeltaPTT_CheatingMuon([](const caf::SRSliceProxy* slc) -> float {
    float dpTT = -9999.;

    unsigned int idxMuon = (unsigned int) kRecoMuonIdx(slc);
    unsigned int idxProton = (unsigned int) kRecoProtonIdx(slc); //Leading proton
    if ( kScndProtonIdx(slc) < 0 ) return dpTT;
    unsigned int idxThird = (unsigned int) kScndProtonIdx(slc); //Charged pion or sub-leading proton

    double pP_mag = kRecoProtonP(slc);
    double pThird_mag = kScndProtonP(slc);

    TVector3 pMu(slc->reco.pfp.at(idxMuon).trk.truth.p.startp.x, slc->reco.pfp.at(idxMuon).trk.truth.p.startp.y, slc->reco.pfp.at(idxMuon).trk.truth.p.startp.z);
    TVector3 pP(pP_mag*slc->reco.pfp.at(idxProton).trk.dir.x, pP_mag*slc->reco.pfp.at(idxProton).trk.dir.y, pP_mag*slc->reco.pfp.at(idxProton).trk.dir.z);
    TVector3 pThird(pThird_mag*slc->reco.pfp.at(idxThird).trk.dir.x, pThird_mag*slc->reco.pfp.at(idxThird).trk.dir.y, pThird_mag*slc->reco.pfp.at(idxThird).trk.dir.z);

    const auto& vtx = slc->vertex;
    TVector3 vec_vtx(vtx.x, vtx.y, vtx.z);
    TVector3 vec_numi_to_vtx = (dFromNuMI + vec_vtx).Unit();

    TVector3 doubleTransverse = ( vec_numi_to_vtx.Cross(pMu) ).Unit();

    double pTTP = pP.Dot(doubleTransverse);
    double pTTThird = pThird.Dot(doubleTransverse);
    dpTT = pTTP + pTTThird;

    return dpTT;
  });

const Var kDeltaPTT_CheatingAngles([](const caf::SRSliceProxy* slc) -> float {
    float dpTT = -9999.;

    unsigned int idxMuon = (unsigned int) kRecoMuonIdx(slc);
    unsigned int idxProton = (unsigned int) kRecoProtonIdx(slc); //Leading proton
    if ( kScndProtonIdx(slc) < 0 ) return dpTT;
    unsigned int idxThird = (unsigned int) kScndProtonIdx(slc); //Charged pion or sub-leading proton

    double pMu_mag = kRecoMuonPNew(slc);
    double pP_mag = kRecoProtonP(slc);
    double pThird_mag = kScndProtonP(slc);

    TVector3 pMu(slc->reco.pfp.at(idxMuon).trk.truth.p.startp.x, slc->reco.pfp.at(idxMuon).trk.truth.p.startp.y, slc->reco.pfp.at(idxMuon).trk.truth.p.startp.z);
    TVector3 pP(slc->reco.pfp.at(idxProton).trk.truth.p.startp.x, slc->reco.pfp.at(idxProton).trk.truth.p.startp.y, slc->reco.pfp.at(idxProton).trk.truth.p.startp.z);
    TVector3 pThird(slc->reco.pfp.at(idxThird).trk.truth.p.startp.x, slc->reco.pfp.at(idxThird).trk.truth.p.startp.y, slc->reco.pfp.at(idxThird).trk.truth.p.startp.z);

    pMu *= pMu_mag / pMu.Mag();
    pP *= pP_mag / pP.Mag();
    pThird *= pThird_mag / pThird.Mag();

    const auto& vtx = slc->vertex;
    TVector3 vec_vtx(vtx.x, vtx.y, vtx.z);
    TVector3 vec_numi_to_vtx = (dFromNuMI + vec_vtx).Unit();

    TVector3 doubleTransverse = ( vec_numi_to_vtx.Cross(pMu) ).Unit();

    double pTTP = pP.Dot(doubleTransverse);
    double pTTThird = pThird.Dot(doubleTransverse);
    dpTT = pTTP + pTTThird;

    return dpTT;
  });

const Var kDeltaPTT_Pion([](const caf::SRSliceProxy* slc) -> float {
    float dpTT = -9999.;

    unsigned int idxMuon = (unsigned int) kRecoMuonIdx(slc);
    unsigned int idxProton = (unsigned int) kRecoProtonIdx(slc); //Leading proton
    unsigned int idxThird = (unsigned int) kSidebandPion(slc); //Charged pion or sub-leading proton

    double pMu_mag = kRecoMuonPNew(slc);
    double pP_mag = kRecoProtonP(slc);
    double pThird_mag = kSidebandPionP(slc);

    TVector3 pMu(pMu_mag*slc->reco.pfp.at(idxMuon).trk.dir.x, pMu_mag*slc->reco.pfp.at(idxMuon).trk.dir.y, pMu_mag*slc->reco.pfp.at(idxMuon).trk.dir.z);
    TVector3 pP(pP_mag*slc->reco.pfp.at(idxProton).trk.dir.x, pP_mag*slc->reco.pfp.at(idxProton).trk.dir.y, pP_mag*slc->reco.pfp.at(idxProton).trk.dir.z);
    TVector3 pThird(pThird_mag*slc->reco.pfp.at(idxThird).trk.dir.x, pThird_mag*slc->reco.pfp.at(idxThird).trk.dir.y, pThird_mag*slc->reco.pfp.at(idxThird).trk.dir.z);
    
    const auto& vtx = slc->vertex;
    TVector3 vec_vtx(vtx.x, vtx.y, vtx.z);
    TVector3 vec_numi_to_vtx = (dFromNuMI + vec_vtx).Unit();

    TVector3 doubleTransverse = ( vec_numi_to_vtx.Cross(pMu) ).Unit();

    double pTTP = pP.Dot(doubleTransverse);
    double pTTThird = pThird.Dot(doubleTransverse);
    dpTT = pTTP + pTTThird;

    return dpTT;
  });

/*
const Var kDeltaPTT_ThreeP([](const caf::SRSliceProxy* slc) -> float {
    float dpTT = -9999.;

    unsigned int idxMuon = (unsigned int) kRecoMuonIdx(slc);
    unsigned int idxP1 = (unsigned int) kRecoProtonIdx(slc); //Leading proton
    unsigned int idxP2 = (unsigned int) kScndProtonIdx(slc); //Charged pion or sub-leading proton
    unsigned int idxP3 = (unsigned int) kSidebandProton(slc);

    double pMu_mag = kRecoMuonPNew(slc);
    double pP1_mag = kRecoProtonP(slc);
    double pP2_mag = kScndProtonP(slc);
    double pP3_mag = kSidebandProtonP(slc);

    TVector3 pMu(pMu_mag*slc->reco.pfp.at(idxMuon).trk.dir.x, pMu_mag*slc->reco.pfp.at(idxMuon).trk.dir.y, pMu_mag*slc->reco.pfp.at(idxMuon).trk.dir.z);
    TVector3 pP1(pP1_mag*slc->reco.pfp.at(idxP1).trk.dir.x, pP1_mag*slc->reco.pfp.at(idxP1).trk.dir.y, pP1_mag*slc->reco.pfp.at(idxP1).trk.dir.z);
    TVector3 pP2(pP2_mag*slc->reco.pfp.at(idxP2).trk.dir.x, pP2_mag*slc->reco.pfp.at(idxP2).trk.dir.y, pP2_mag*slc->reco.pfp.at(idxP2).trk.dir.z);
    TVector3 pP3(pP3_mag*slc->reco.pfp.at(idxP3).trk.dir.x, pP3_mag*slc->reco.pfp.at(idxP3).trk.dir.y, pP3_mag*slc->reco.pfp.at(idxP3).trk.dir.z);

    const auto& vtx = slc->vertex;
    TVector3 vec_vtx(vtx.x, vtx.y, vtx.z);
    TVector3 vec_numi_to_vtx = (dFromNuMI + vec_vtx).Unit();

    TVector3 doubleTransverse = ( vec_numi_to_vtx.Cross(pMu) ).Unit();

    double pTTP1 = pP1.Dot(doubleTransverse);
    double pTTP2 = pP2.Dot(doubleTransverse);
    double pTTP3 = pP3.Dot(doubleTransverse);
    dpTT = pTTP1 + pTTP2 + pTTP3;

    return dpTT;
  });
*/

const Var kDeltaPTT_ProtonTruth([](const caf::SRSliceProxy* slc) -> float {
    float dpTT = -9999.;

    if ( slc->truth.index < 0 ) return dpTT;
    TVector3 NuDirection_Truth = (TVector3(slc->truth.momentum.x, slc->truth.momentum.y, slc->truth.momentum.z)).Unit().Unit();

    unsigned int idxMuon = (unsigned int) kRecoMuonIdx(slc);
    unsigned int idxProton = (unsigned int) kRecoProtonIdx(slc); //Leading proton
    unsigned int idxThird = (unsigned int) kScndProtonIdx(slc); //Charged pion or sub-leading proton

    TVector3 pMu(slc->reco.pfp.at(idxMuon).trk.truth.p.startp.x, slc->reco.pfp.at(idxMuon).trk.truth.p.startp.y, slc->reco.pfp.at(idxMuon).trk.truth.p.startp.z);
    TVector3 pP(slc->reco.pfp.at(idxProton).trk.truth.p.startp.x, slc->reco.pfp.at(idxProton).trk.truth.p.startp.y, slc->reco.pfp.at(idxProton).trk.truth.p.startp.z);
    TVector3 pThird(slc->reco.pfp.at(idxThird).trk.truth.p.startp.x, slc->reco.pfp.at(idxThird).trk.truth.p.startp.y, slc->reco.pfp.at(idxThird).trk.truth.p.startp.z);

    TVector3 doubleTransverse = ( NuDirection_Truth.Cross(pMu) ).Unit();

    double pTTP = pP.Dot(doubleTransverse);
    double pTTThird = pThird.Dot(doubleTransverse);
    dpTT = pTTP + pTTThird;

    return dpTT;
  });

const Var kDeltaPTT_PionTruth([](const caf::SRSliceProxy* slc) -> float {
    float dpTT = -9999.;

    if ( slc->truth.index < 0 ) return dpTT;
    TVector3 NuDirection_Truth = (TVector3(slc->truth.momentum.x, slc->truth.momentum.y, slc->truth.momentum.z)).Unit().Unit();

    unsigned int idxMuon = (unsigned int) kRecoMuonIdx(slc);
    unsigned int idxProton = (unsigned int) kRecoProtonIdx(slc); //Leading proton
    unsigned int idxThird = (unsigned int) kSidebandPion(slc); //Charged pion or sub-leading proton

    TVector3 pMu(slc->reco.pfp.at(idxMuon).trk.truth.p.startp.x, slc->reco.pfp.at(idxMuon).trk.truth.p.startp.y, slc->reco.pfp.at(idxMuon).trk.truth.p.startp.z);
    TVector3 pP(slc->reco.pfp.at(idxProton).trk.truth.p.startp.x, slc->reco.pfp.at(idxProton).trk.truth.p.startp.y, slc->reco.pfp.at(idxProton).trk.truth.p.startp.z);
    TVector3 pThird(slc->reco.pfp.at(idxThird).trk.truth.p.startp.x, slc->reco.pfp.at(idxThird).trk.truth.p.startp.y, slc->reco.pfp.at(idxThird).trk.truth.p.startp.z);

    TVector3 doubleTransverse = ( NuDirection_Truth.Cross(pMu) ).Unit();

    double pTTP = pP.Dot(doubleTransverse);
    double pTTThird = pThird.Dot(doubleTransverse);
    dpTT = pTTP + pTTThird;

    return dpTT;
  });

/*
const Var kDeltaPTT_ThreePTruth([](const caf::SRSliceProxy* slc) -> float {
    float dpTT = -9999.;

    if ( slc->truth.index < 0 ) return dpTT;
    TVector3 NuDirection_Truth = (TVector3(slc->truth.momentum.x, slc->truth.momentum.y, slc->truth.momentum.z)).Unit().Unit();

    unsigned int idxMuon = (unsigned int) kRecoMuonIdx(slc);
    unsigned int idxP1 = (unsigned int) kRecoProtonIdx(slc); //Leading proton
    unsigned int idxP2 = (unsigned int) kScndProtonIdx(slc); //Charged pion or sub-leading proton
    unsigned int idxP3 = (unsigned int) kSidebandProton(slc);

    TVector3 pMu(slc->reco.pfp.at(idxMuon).trk.truth.p.startp.x, slc->reco.pfp.at(idxMuon).trk.truth.p.startp.y, slc->reco.pfp.at(idxMuon).trk.truth.p.startp.z);
    TVector3 pP1(slc->reco.pfp.at(idxP1).trk.truth.p.startp.x, slc->reco.pfp.at(idxP1).trk.truth.p.startp.y, slc->reco.pfp.at(idxP1).trk.truth.p.startp.z);
    TVector3 pP2(slc->reco.pfp.at(idxP2).trk.truth.p.startp.x, slc->reco.pfp.at(idxP2).trk.truth.p.startp.y, slc->reco.pfp.at(idxP2).trk.truth.p.startp.z);
    TVector3 pP3(slc->reco.pfp.at(idxP3).trk.truth.p.startp.x, slc->reco.pfp.at(idxP3).trk.truth.p.startp.y, slc->reco.pfp.at(idxP3).trk.truth.p.startp.z);

    const auto& vtx = slc->vertex;
    TVector3 vec_vtx(vtx.x, vtx.y, vtx.z);
    TVector3 vec_numi_to_vtx = (dFromNuMI + vec_vtx).Unit();

    TVector3 doubleTransverse = ( vec_numi_to_vtx.Cross(pMu) ).Unit();

    double pTTP1 = pP1.Dot(doubleTransverse);
    double pTTP2 = pP2.Dot(doubleTransverse);
    double pTTP3 = pP3.Dot(doubleTransverse);
    dpTT = pTTP1 + pTTP2 + pTTP3;

    return dpTT;
  });
*/

const Var kDeltaPTT([](const caf::SRSliceProxy* slc) -> float {
    double dPTT = -9999.;
    if ( kScndProtonCandidate(slc) ) dPTT = kDeltaPTT_Proton(slc);
    //else if ( kThreePSideband(slc) ) dPTT = kDeltaPTT_ThreeP(slc);
    else if ( kPionSidebandBase(slc) ) dPTT = kDeltaPTT_Pion(slc);
    assert(((void)"No valid selection for DeltaPTT value", dPTT != -9999.));

    return dPTT;
  });

const Var kDeltaPTT_OldTruth([](const caf::SRSliceProxy* slc) -> float {
    double dPTT = -99999.; //Note the extra 9, as -9999. is a valid response for selected cosmic slices
    if ( kScndProtonCandidate(slc) ) dPTT = kDeltaPTT_ProtonTruth(slc);
    //else if ( kThreePSideband(slc) ) dPTT = kDeltaPTT_ThreePTruth(slc);
    else if ( kPionSidebandBase(slc) ) dPTT = kDeltaPTT_PionTruth(slc);
    assert(((void)"No valid selection for DeltaPTT value", dPTT != -99999.));

    return dPTT;
  });

//Built out of SRTrueInteraction rather than truth info matched to reco tracks
//Considers all charged pions and protons over threshold, so should be correct for all selections
const Var kDeltaPTT_FullTruth([](const caf::SRSliceProxy* slc) -> float {
    float dpTT = -9999.;
    if ( slc->truth.index < 0 ) return dpTT;
    const auto &nu = slc->truth;

    TVector3 NuDirection_Truth = (TVector3(nu.momentum.x, nu.momentum.y, nu.momentum.z)).Unit().Unit();
    TVector3 p3;
    TVector3 pMu(0., 0., 0.);
    TVector3 doubleTransverse(0., 0., 0.);

    for ( unsigned int i = 0; i < nu.prim.size(); i++ ) {
      const auto &prim = nu.prim.at(i);
      if ( prim.start_process != 0 ) continue;
      p3 = {prim.startp.x, prim.startp.y, prim.startp.z};
      if ( p3.Mag() == std::hypot(-9999., -9999., -9999.) ) return dpTT;
      if ( abs(prim.pdg) == 13 ) {
        pMu = p3;
        break;
      }
    }
    if ( pMu.Mag() == 0. ) return dpTT;
    doubleTransverse = ( NuDirection_Truth.Cross(pMu) ).Unit();

    dpTT = 0.;
    for ( unsigned int i = 0; i < nu.prim.size(); i++ ) {
      const auto &prim = nu.prim.at(i);
      if ( prim.start_process != 0 ) continue;
      p3 = {prim.startp.x, prim.startp.y, prim.startp.z};
      if ( abs(prim.pdg) == 13 ) continue;
      else if ( abs(prim.pdg) == 2212 && p3.Mag() > .35 ) dpTT += p3.Dot(doubleTransverse);
      else if ( abs(prim.pdg) == 211 || abs(prim.pdg) == 321 ) dpTT += p3.Dot(doubleTransverse);
    }

    return dpTT;
  });

const Var kDeltaPTT_Truth([](const caf::SRSliceProxy* slc) -> float {
    float dpTT = -9999.;
    if ( slc->truth.index < 0 ) return dpTT;
    const auto &nu = slc->truth;

    TVector3 NuDirection_Truth = (TVector3(nu.momentum.x, nu.momentum.y, nu.momentum.z)).Unit().Unit();
    TVector3 p3;
    TVector3 pMu(0., 0., 0.);
    TVector3 doubleTransverse(0., 0., 0.);

    int idxP1 = -1;
    int idxP2 = -1;
    double leading = -9999.;
    double second = -9999.;

    for ( unsigned int i = 0; i < nu.prim.size(); i++ ) {
      const auto &prim = nu.prim.at(i);
      if ( prim.start_process != 0 ) continue;
      p3 = {prim.startp.x, prim.startp.y, prim.startp.z};
      if ( p3.Mag() == std::hypot(-9999., -9999., -9999.) ) return dpTT;
      if ( abs(prim.pdg) == 13 ) {
        pMu = p3;
      }
      else if ( abs(prim.pdg) == 2212 && p3.Mag() > .35 ) {
        if ( p3.Mag() > leading ) {
          second = leading;
          leading = p3.Mag();
          idxP2 = idxP1;
          idxP1 = i;
        }
        else if ( p3.Mag() > second ) {
          second = p3.Mag();
          idxP2 = i;
        }
      }
    }

    if ( pMu.Mag() == 0. ) return dpTT;
    doubleTransverse = ( NuDirection_Truth.Cross(pMu) ).Unit();
    if ( idxP2 == -1 ) return dpTT;

    dpTT = 0.;
    TVector3 P1 = {nu.prim.at(idxP1).startp.x, nu.prim.at(idxP1).startp.y, nu.prim.at(idxP1).startp.z};
    TVector3 P2 = {nu.prim.at(idxP2).startp.x, nu.prim.at(idxP2).startp.y, nu.prim.at(idxP2).startp.z};
    dpTT += P1.Dot(doubleTransverse);
    dpTT += P2.Dot(doubleTransverse);

    return dpTT;
  });

const Var kDeltaPTT_Resid([](const caf::SRSliceProxy* slc) -> float {
    double recoDeltaPTT = kDeltaPTT(slc);
    double trueDeltaPTT = kDeltaPTT_Truth(slc);

    if ( trueDeltaPTT == -9999. || trueDeltaPTT == 0. ) return -9999.;
    return (recoDeltaPTT - trueDeltaPTT) / trueDeltaPTT;
  });

const Var kDeltaPTT_ProtonResid([](const caf::SRSliceProxy* slc) -> float {
    double recoDeltaPTT = kDeltaPTT_Proton(slc);
    double trueDeltaPTT = kDeltaPTT_Truth(slc);

    if ( trueDeltaPTT == -9999. || trueDeltaPTT == 0. ) return -9999.;
    return (recoDeltaPTT - trueDeltaPTT) / trueDeltaPTT;
  });

const Var kDeltaPTT_CheatingMuonResid([](const caf::SRSliceProxy* slc) -> float {
    double recoDeltaPTT = kDeltaPTT_CheatingMuon(slc);
    double trueDeltaPTT = kDeltaPTT_Truth(slc);

    if ( trueDeltaPTT == -9999. || trueDeltaPTT == 0. ) return -9999.;
    return (recoDeltaPTT - trueDeltaPTT) / trueDeltaPTT;
  });

const Var kDeltaPTT_CheatingAnglesResid([](const caf::SRSliceProxy* slc) -> float {
    double recoDeltaPTT = kDeltaPTT_CheatingAngles(slc);
    double trueDeltaPTT = kDeltaPTT_Truth(slc);

    if ( trueDeltaPTT == -9999. || trueDeltaPTT == 0. ) return -9999.;
    return (recoDeltaPTT - trueDeltaPTT) / trueDeltaPTT;
  });


//! Truth Section
//! Here we sue t<NAME> instead of k<NAME> convention to distinguish TruthCuts and TruthVars from Cuts and Vars

bool containedMuon = false;
bool containedProton = false;
double muonThresh = 0.226;
double protonThresh = 0.35;
double protonMax = 2.;
double pionThresh=-1.;
double pi0Thresh=-1.;


const TruthCut tIsSignal([](const caf::SRTrueInteractionProxy* nu) {
  if ( abs(nu->pdg) != 14 ||
       !nu->iscc ||
       std::isnan(nu->position.x) || std::isnan(nu->position.y) || std::isnan(nu->position.z) ||
       !isInFV(nu->position.x,nu->position.y,nu->position.z) )
    return false;

  unsigned int nMu(0), nP(0), nPMax(0), nPi0(0), nChgPi(0);

  for ( auto const& prim : nu->prim ) {
    if ( prim.start_process != 0 ) continue;
    TVector3 p3(prim.startp.x, prim.startp.y, prim.startp.z);
    double momentum = p3.Mag();
    if ( abs(prim.pdg) == 13 ) {
      if ( containedMuon && !prim.contained ) continue;
      if ( momentum < muonThresh ) continue;
      nMu+=1;
    }
    if ( abs(prim.pdg) == 2212 ){
      if ( containedProton && !prim.contained ) continue;
      if ( momentum < protonThresh ) continue;
      if ( momentum > protonMax ) nPMax+=1;
      else nP+=1;
    }
    if ( abs(prim.pdg) == 111 ){
      if ( momentum < pi0Thresh ) continue;
      nPi0+=1;
    }
    if ( abs(prim.pdg) == 211 ){
      if ( momentum < pionThresh ) continue;
      nChgPi+=1;
    }
  }

  return ( nMu == 1 && nP >= 2 && nChgPi == 0 && nPi0 == 0 && nPMax == 0 );
  });


int trueMuon_idx(const caf::SRTrueInteractionProxy* nu) {
  int idx = -1;

  for ( unsigned int i = 0; i < nu->prim.size(); i++ ) {
    const auto &prim = nu->prim.at(i);
    if ( prim.start_process != 0 ) continue;
    if ( abs(prim.pdg) != 13 ) continue;
    TVector3 p3(prim.startp.x, prim.startp.y, prim.startp.z);
    double momentum = p3.Mag();
    if ( momentum < muonThresh ) continue;
    idx = i;
  }
  return idx;
}

int trueLeadingProton_idx(const caf::SRTrueInteractionProxy* nu) {
  int idx = -1;
  double leadingP(0.);

  for ( unsigned int i = 0; i < nu->prim.size(); i++ ) {
    const auto &prim = nu->prim.at(i);
    if ( prim.start_process != 0 ) continue;
    if ( abs(prim.pdg) != 2212 ) continue;
    TVector3 p3(prim.startp.x, prim.startp.y, prim.startp.z);
    double momentum = p3.Mag();
    if ( momentum < protonThresh ) continue;
    if ( momentum > leadingP ) {
      leadingP = momentum;
      idx = i;
    }
  }
  return idx;
}

int trueScndProton_idx(const caf::SRTrueInteractionProxy* nu) {
  int idx = -1;
  double scndP(0.);
  unsigned int leadingProton = trueLeadingProton_idx(nu);

  for ( unsigned int i = 0; i < nu->prim.size(); i++ ) {
    if ( i == leadingProton ) continue;
    const auto &prim = nu->prim.at(i);
    if ( prim.start_process != 0 ) continue;
    if ( abs(prim.pdg) != 2212 ) continue;
    TVector3 p3(prim.startp.x, prim.startp.y, prim.startp.z);
    double momentum = p3.Mag();
    if ( momentum < protonThresh ) continue;
    if ( momentum > scndP ) {
      scndP = momentum;
      idx = i;
    }
  }
  return idx;
}

int trueThirdProton_idx(const caf::SRTrueInteractionProxy* nu) {
  int idx = -1;
  double thirdP(0.);
  unsigned int leadingProton = trueLeadingProton_idx(nu);
  unsigned int scndProton = trueScndProton_idx(nu);

  for ( unsigned int i = 0; i < nu->prim.size(); i++ ) {
    if ( i == leadingProton || i == scndProton ) continue;
    const auto &prim = nu->prim.at(i);
    if ( prim.start_process != 0 ) continue;
    if ( abs(prim.pdg) != 2212 ) continue;
    TVector3 p3(prim.startp.x, prim.startp.y, prim.startp.z);
    double momentum = p3.Mag();
    if ( momentum < protonThresh ) continue;
    if ( momentum > thirdP ) {
      thirdP = momentum;
      idx = i;
    }
  }
  return idx;
}

int trueChargedPion_idx(const caf::SRTrueInteractionProxy* nu) {
  int idx = -1;
  double leadingP(0.);

  for ( unsigned int i = 0; i < nu->prim.size(); i++ ) {
    const auto &prim = nu->prim.at(i);
    if ( prim.start_process != 0 ) continue;
    if ( abs(prim.pdg) != 211 ) continue;
    TVector3 p3(prim.startp.x, prim.startp.y, prim.startp.z);
    double momentum = p3.Mag();
    if ( momentum > leadingP ) {
      leadingP = momentum;
      idx = i;
    }
  }
  return idx;
}

const TruthVar tTrueNeutrinoPDG([](const caf::SRTrueInteractionProxy* nu) -> double {
  return ( nu->pdg );
  });

const Var kTrueNeutrinoPDG ([](const caf::SRSliceProxy* slc) -> float {
    if ( slc->truth.index < 0 ) return -9999.;
    else return tTrueNeutrinoPDG(&slc->truth);
  });

const TruthVar tTrueNeutrinoE([](const caf::SRTrueInteractionProxy* nu) -> double {
  return ( nu->E );
  });

const Var kTrueNeutrinoE ([](const caf::SRSliceProxy* slc) -> float {
    if ( !(k1mu2p0pi(slc) || k1mu3p0pi(slc)) ) return -9999.;
    else return tTrueNeutrinoE(&slc->truth);
  });

const TruthVar tTrueMuonP([](const caf::SRTrueInteractionProxy* nu) -> double {
  int muon = trueMuon_idx(nu);

  TVector3 muP3 (nu->prim.at(muon).startp.x, nu->prim.at(muon).startp.y, nu->prim.at(muon).startp.z);

  return ( muP3.Mag() );
  });

const Var kTrueMuonP ([](const caf::SRSliceProxy* slc) -> float {
    if ( !(tIsSignal(&slc->truth)) ) return -9999.;
    else return tTrueMuonP(&slc->truth);
  });

const TruthVar tTrueLeadingProtonP([](const caf::SRTrueInteractionProxy* nu) -> double {
  int leadingProton = trueLeadingProton_idx(nu);

  TVector3 leadingP3 (nu->prim.at(leadingProton).startp.x, nu->prim.at(leadingProton).startp.y, nu->prim.at(leadingProton).startp.z);

  return ( leadingP3.Mag() );
  });

const Var kTrueLeadingProtonP ([](const caf::SRSliceProxy* slc) -> float {
    if ( trueLeadingProton_idx(&slc->truth) < 0 ) return -9999.;
    else return tTrueLeadingProtonP(&slc->truth);
  });

const TruthVar tTrueSecondProtonP([](const caf::SRTrueInteractionProxy* nu) -> double {
  int secondProton = trueScndProton_idx(nu);

  TVector3 secondP3 (nu->prim.at(secondProton).startp.x, nu->prim.at(secondProton).startp.y, nu->prim.at(secondProton).startp.z);

  return ( secondP3.Mag() );
  });

const Var kTrueSecondProtonP ([](const caf::SRSliceProxy* slc) -> float {
    if ( trueScndProton_idx(&slc->truth) < 0 ) return -9999.;
    else return tTrueSecondProtonP(&slc->truth);
  });

const TruthVar tTrueThirdProtonP([](const caf::SRTrueInteractionProxy* nu) -> double {
  int thirdProton = trueThirdProton_idx(nu);

  TVector3 thirdP3 (nu->prim.at(thirdProton).startp.x, nu->prim.at(thirdProton).startp.y, nu->prim.at(thirdProton).startp.z);

  return ( thirdP3.Mag() );
  });

const TruthVar tTrueChargedPionP([](const caf::SRTrueInteractionProxy* nu) -> double {
  int pion = trueChargedPion_idx(nu);

  TVector3 pionP3 (nu->prim.at(pion).startp.x, nu->prim.at(pion).startp.y, nu->prim.at(pion).startp.z);

  return ( pionP3.Mag() );
  });

const TruthVar tTrueLeadingProtonPFrac([](const caf::SRTrueInteractionProxy* nu) -> double {
  int leadingProton = trueLeadingProton_idx(nu);
  int secondProton = trueScndProton_idx(nu);

  TVector3 leadingP3 (nu->prim.at(leadingProton).startp.x, nu->prim.at(leadingProton).startp.y, nu->prim.at(leadingProton).startp.z);
  TVector3 secondP3 (nu->prim.at(secondProton).startp.x, nu->prim.at(secondProton).startp.y, nu->prim.at(secondProton).startp.z);

  return ( leadingP3.Mag() / (leadingP3.Mag() + secondP3.Mag()) );
  });

const Var kTrueLeadingProtonPFrac ([](const caf::SRSliceProxy* slc) -> float {
    if ( !(tIsSignal(&slc->truth)) ) return -9999.;
    else return tTrueLeadingProtonPFrac(&slc->truth);
  });

const TruthVar tTrueHadronicOpeningAngle([](const caf::SRTrueInteractionProxy* nu) -> double {
  int leadingProton = trueLeadingProton_idx(nu);
  int secondProton = trueScndProton_idx(nu);

  TVector3 leadingP3 (nu->prim.at(leadingProton).startp.x, nu->prim.at(leadingProton).startp.y, nu->prim.at(leadingProton).startp.z);
  TVector3 secondP3 (nu->prim.at(secondProton).startp.x, nu->prim.at(secondProton).startp.y, nu->prim.at(secondProton).startp.z);

  return ( leadingP3.Dot(secondP3) / (leadingP3.Mag() * secondP3.Mag()) );
  });

const Var kTrueHadronicOpeningAngle ([](const caf::SRSliceProxy* slc) -> float {
    if ( !(tIsSignal(&slc->truth)) ) return -9999.;
    else return tTrueHadronicOpeningAngle(&slc->truth);
  });

const TruthVar tTrueMuonHadronAngle([](const caf::SRTrueInteractionProxy* nu) -> double {
  int muon = trueMuon_idx(nu);
  int leadingProton = trueLeadingProton_idx(nu);
  int secondProton = trueScndProton_idx(nu);

  TVector3 muP3 (nu->prim.at(muon).startp.x, nu->prim.at(muon).startp.y, nu->prim.at(muon).startp.z);
  TVector3 leadingP3 (nu->prim.at(leadingProton).startp.x, nu->prim.at(leadingProton).startp.y, nu->prim.at(leadingProton).startp.z);
  TVector3 secondP3 (nu->prim.at(secondProton).startp.x, nu->prim.at(secondProton).startp.y, nu->prim.at(secondProton).startp.z);
  

  return ( muP3.Dot( leadingP3+secondP3 ) / (muP3.Mag() * (leadingP3+secondP3).Mag()) );
  });

const Var kTrueMuonHadronAngle ([](const caf::SRSliceProxy* slc) -> float {
    if ( !(tIsSignal(&slc->truth)) ) return -9999.;
    else return tTrueMuonHadronAngle(&slc->truth);
  });

const TruthVar tTrueDeltaPT([](const caf::SRTrueInteractionProxy* nu) -> double {
  int muon = trueMuon_idx(nu);
  int leadingProton = trueLeadingProton_idx(nu);
  int secondProton = trueScndProton_idx(nu);

  TVector3 muP3 (nu->prim.at(muon).startp.x, nu->prim.at(muon).startp.y, nu->prim.at(muon).startp.z);
  TVector3 leadingP3 (nu->prim.at(leadingProton).startp.x, nu->prim.at(leadingProton).startp.y, nu->prim.at(leadingProton).startp.z);
  TVector3 secondP3 (nu->prim.at(secondProton).startp.x, nu->prim.at(secondProton).startp.y, nu->prim.at(secondProton).startp.z);
  TVector3 NuDirection_Truth = (TVector3(nu->momentum.x, nu->momentum.y, nu->momentum.z)).Unit().Unit();  

  TVector3 muPT = muP3 - (muP3.Dot(NuDirection_Truth) * NuDirection_Truth);
  TVector3 leadingPT = leadingP3 - (leadingP3.Dot(NuDirection_Truth) * NuDirection_Truth);
  TVector3 secondPT = secondP3 - (secondP3.Dot(NuDirection_Truth) * NuDirection_Truth);
  TVector3 dpT = muPT + leadingPT + secondPT;

  return dpT.Mag();
  });

const Var kTrueDeltaPT ([](const caf::SRSliceProxy* slc) -> float {
    if ( !(tIsSignal(&slc->truth)) ) return -9999.;
    else return tTrueDeltaPT(&slc->truth);
  });

const TruthVar tTrueDeltaAlphaT([](const caf::SRTrueInteractionProxy* nu) -> double {
  int muon = trueMuon_idx(nu);
  int leadingProton = trueLeadingProton_idx(nu);
  int secondProton = trueScndProton_idx(nu);

  TVector3 muP3 (nu->prim.at(muon).startp.x, nu->prim.at(muon).startp.y, nu->prim.at(muon).startp.z);
  TVector3 leadingP3 (nu->prim.at(leadingProton).startp.x, nu->prim.at(leadingProton).startp.y, nu->prim.at(leadingProton).startp.z);
  TVector3 secondP3 (nu->prim.at(secondProton).startp.x, nu->prim.at(secondProton).startp.y, nu->prim.at(secondProton).startp.z);
  TVector3 NuDirection_Truth = (TVector3(nu->momentum.x, nu->momentum.y, nu->momentum.z)).Unit().Unit();

  TVector3 muPT = muP3 - (muP3.Dot(NuDirection_Truth) * NuDirection_Truth);
  TVector3 leadingPT = leadingP3 - (leadingP3.Dot(NuDirection_Truth) * NuDirection_Truth);
  TVector3 secondPT = secondP3 - (secondP3.Dot(NuDirection_Truth) * NuDirection_Truth);
  TVector3 dpT = muPT + leadingPT + secondPT;

  double dalphaT = TMath::ACos( -1. * (muPT.Dot(dpT)) / (muPT.Mag() * dpT.Mag()) );

  return dalphaT;
  });

const Var kTrueDeltaAlphaT ([](const caf::SRSliceProxy* slc) -> float {
    if ( !(tIsSignal(&slc->truth)) ) return -9999.;
    else return tTrueDeltaAlphaT(&slc->truth);
  });

const TruthVar tTrueDeltaPhiT([](const caf::SRTrueInteractionProxy* nu) -> double {
  int muon = trueMuon_idx(nu);
  int leadingProton = trueLeadingProton_idx(nu);
  int secondProton = trueScndProton_idx(nu);
 
  TVector3 muP3 (nu->prim.at(muon).startp.x, nu->prim.at(muon).startp.y, nu->prim.at(muon).startp.z);
  TVector3 leadingP3 (nu->prim.at(leadingProton).startp.x, nu->prim.at(leadingProton).startp.y, nu->prim.at(leadingProton).startp.z);
  TVector3 secondP3 (nu->prim.at(secondProton).startp.x, nu->prim.at(secondProton).startp.y, nu->prim.at(secondProton).startp.z);
  TVector3 NuDirection_Truth = (TVector3(nu->momentum.x, nu->momentum.y, nu->momentum.z)).Unit().Unit();
 
  TVector3 muPT = muP3 - (muP3.Dot(NuDirection_Truth) * NuDirection_Truth);
  TVector3 leadingPT = leadingP3 - (leadingP3.Dot(NuDirection_Truth) * NuDirection_Truth);
  TVector3 secondPT = secondP3 - (secondP3.Dot(NuDirection_Truth) * NuDirection_Truth);
  TVector3 hadronPT = leadingPT + secondPT;
 
  double dphiT = TMath::ACos( -1. * (muPT.Dot(hadronPT)) / (muPT.Mag() * hadronPT.Mag()) );

  return dphiT;
  });

const Var kTrueDeltaPhiT ([](const caf::SRSliceProxy* slc) -> float {
    if ( !(tIsSignal(&slc->truth)) ) return -9999.;
    else return tTrueDeltaPhiT(&slc->truth);
  });

const TruthVar tTrueDeltaPTT([](const caf::SRTrueInteractionProxy* nu) -> double {
  int muon = trueMuon_idx(nu);
  int leadingProton = trueLeadingProton_idx(nu);
  int secondProton = trueScndProton_idx(nu);

  TVector3 muP3 (nu->prim.at(muon).startp.x, nu->prim.at(muon).startp.y, nu->prim.at(muon).startp.z);
  TVector3 leadingP3 (nu->prim.at(leadingProton).startp.x, nu->prim.at(leadingProton).startp.y, nu->prim.at(leadingProton).startp.z);
  TVector3 secondP3 (nu->prim.at(secondProton).startp.x, nu->prim.at(secondProton).startp.y, nu->prim.at(secondProton).startp.z);
  TVector3 NuDirection_Truth = (TVector3(nu->momentum.x, nu->momentum.y, nu->momentum.z)).Unit().Unit();
  TVector3 doubleTransverse = ( NuDirection_Truth.Cross(muP3) ).Unit(); 

  double leadingPTT = leadingP3.Dot(doubleTransverse);
  double secondPTT = secondP3.Dot(doubleTransverse);
  double dpTT = leadingPTT + secondPTT;

  return dpTT;
  });

const Var kTrueDeltaPTT ([](const caf::SRSliceProxy* slc) -> float {
    if ( !(tIsSignal(&slc->truth)) ) return -9999.;
    else return tTrueDeltaPTT(&slc->truth);
  });

const TruthVar tTrueEHad([](const caf::SRTrueInteractionProxy* nu) -> double {
    double eHad = 0.;

    for ( auto const& prim : nu->prim ) {
      if ( prim.start_process != 0 ) continue;
      if ( prim.startE == -9999. )  return -9999.;
      if ( (abs(prim.pdg) >= 11 && abs(prim.pdg) <= 14) || prim.pdg == 22 || prim.pdg == 1000180400 ) continue;
      else if ( abs(prim.pdg) == 2212 ) eHad += prim.startE - mProton;
      else if ( abs(prim.pdg) == 2112 ) eHad += prim.startE - mNeutron;
      //Just subtract off proton mass for strange hadrons. It's close enough to right, and they're < 0.5% of the sample anyway
      else if ( abs(prim.pdg) > 3000 && abs(prim.pdg) < 4000) eHad += prim.startE - mProton;
      else eHad += prim.startE;
    }

    if ( eHad < 0. ) return -9999.;
    else return eHad;
  });

const Var kTrueEHad ([](const caf::SRSliceProxy* slc) -> float {
    if ( !(tIsSignal(&slc->truth)) ) return -9999.;
    else return tTrueEHad(&slc->truth);
  });

const TruthVar tTrueEAvail([](const caf::SRTrueInteractionProxy* nu) -> double {
    double eAvail = 0.;

    for ( auto const& prim : nu->prim ) {
      if ( prim.start_process != 0 ) continue;
      if ( prim.startE == -9999. )  return -9999.;
      if ( (abs(prim.pdg) >= 11 && abs(prim.pdg) <= 14) || prim.pdg == 22 || prim.pdg == 1000180400 || prim.pdg == 2112 ) continue;
      else if ( abs(prim.pdg) == 2212 ) {
        if ( std::hypot(prim.startp.x,prim.startp.y,prim.startp.z) > .35 ) eAvail += prim.startE - mProton;
        else continue;
      }
      else if ( abs(prim.pdg) > 3000 && abs(prim.pdg) < 4000) eAvail += prim.startE - mProton;
      else eAvail += prim.startE;
    }

    if ( eAvail < 0. ) return -9999.;
    else return eAvail;
  });

const Var kTrueEAvail ([](const caf::SRSliceProxy* slc) -> float {
    if ( !(tIsSignal(&slc->truth)) ) return -9999.;
    else return tTrueEAvail(&slc->truth);
  });

const TruthVar tTrueQ2([](const caf::SRTrueInteractionProxy* nu) -> double {
    double Q2 = -9999.;
    TVector3 NuDirection_Truth = (TVector3(nu->momentum.x, nu->momentum.y, nu->momentum.z)).Unit().Unit();

    double eHad = tTrueEHad(nu);
    if ( eHad == -9999. ) return Q2;

    TVector3 pMu(0., 0., 0.);
    for ( unsigned int i = 0; i < nu->prim.size(); i++ ) {
      const auto &prim = nu->prim.at(i);
      if ( prim.start_process != 0 ) continue;
      if ( abs(prim.pdg) == 13 ) {
        pMu = {prim.startp.x, prim.startp.y, prim.startp.z};
        break;
      }
    }
    if ( pMu.Mag() == 0. || pMu.Mag() == std::hypot(-9999., -9999., -9999.) ) return Q2;

    double eMu = std::hypot(pMu.Mag(), mMuon);
    double eNu = eMu + eHad;

    Q2 = 2.*eNu*(eMu - pMu.Dot(NuDirection_Truth)) - mMuon*mMuon;
    return Q2;
  });

const Var kTrueQ2 ([](const caf::SRSliceProxy* slc) -> float {
    if ( !(tIsSignal(&slc->truth)) ) return -9999.;
    else return tTrueQ2(&slc->truth);
  });

const TruthVar tTrueq3 ([](const caf::SRTrueInteractionProxy* nu) -> double {
    double q3 = -9999.;
    double eHad = tTrueEHad(nu);
    double Q2 = tTrueQ2(nu);
    if ( Q2 == -9999. ) return q3;
    q3 = std::sqrt(Q2 + eHad*eHad);
    return q3;
  });

const Var kTrueq3 ([](const caf::SRSliceProxy* slc) -> float {
    if ( !(tIsSignal(&slc->truth)) ) return -9999.;
    else return tTrueq3(&slc->truth);
  });

const TruthVar tTrueNProtons([](const caf::SRTrueInteractionProxy* nu) -> double {
    int nProtons = 0;

    for ( auto const& prim : nu->prim ) {
      if ( prim.start_process != 0 ) continue;
      if ( prim.startE == -9999. )  return -9999.;
      if ( prim.pdg == 2212 && std::hypot(prim.startp.x,prim.startp.y,prim.startp.z) > .35 ) nProtons++;
    }

    return nProtons;
  });

const Var kTrueNProtons ([](const caf::SRSliceProxy* slc) -> float {
    if ( !(tIsSignal(&slc->truth)) ) return -9999.;
    else return tTrueNProtons(&slc->truth);
  });

const TruthVar tTrueNeutronEnergy([](const caf::SRTrueInteractionProxy* nu) -> double {
    double sum = 0.;

    for ( const auto &prim : nu->prim ) {
      if ( prim.start_process != 0 ) continue;
      if ( prim.startE == -9999. ) return -9999.;
      if ( prim.pdg == 2112 ) sum += prim.startE - mNeutron;
    }

    return sum;
  });

const Var kTrueNeutronEnergy ([](const caf::SRSliceProxy* slc) -> float {
    if ( !(tIsSignal(&slc->truth)) ) return -9999.;
    else return tTrueNeutronEnergy(&slc->truth);
  });

const TruthVar tContained([](const caf::SRTrueInteractionProxy* nu) -> int {
  int leadingProton = trueLeadingProton_idx(nu);
  int secondProton = trueScndProton_idx(nu);

  return ( nu->prim.at(leadingProton).contained && nu->prim.at(secondProton).contained );
  });

const TruthVar tStopping([](const caf::SRTrueInteractionProxy* nu) -> int {
  int leadingProton = trueLeadingProton_idx(nu);
  int secondProton = trueScndProton_idx(nu);

  return ( nu->prim.at(leadingProton).end_process == 2 && nu->prim.at(secondProton).contained );
  });

const TruthVar tAllStopping([](const caf::SRTrueInteractionProxy* nu) -> int {
  bool allStopping = true;

  for ( auto const& prim : nu->prim ) {
    if ( prim.start_process != 0 ) continue;
    if ( prim.pdg == 2212 && std::hypot(prim.startp.x,prim.startp.y,prim.startp.z) > .35 ) {
      if ( !prim.contained || prim.end_process != 2 ) {
        allStopping = false;
        break;
      }
    }
  }

  return allStopping;
  });

/*
bool cutPreselection(const caf::SRSliceProxy &slc) { return (kPreselection(&slc)); }
bool cutMuon(const caf::SRSliceProxy &slc) { return (kPreselection(&slc) && kHasMuon(&slc)); }
bool cutHadronicContainment(const caf::SRSliceProxy &slc) { return (kPreselection(&slc) && kHasMuon(&slc) && kHadronicContainment(&slc)); }
bool cutLeadingProton(const caf::SRSliceProxy &slc) { return (kPreselection(&slc) && kHasMuon(&slc) && kHadronicContainment(&slc) && kLeadingProtonThreshold(&slc)); }
bool cutScndProton(const caf::SRSliceProxy &slc) { return (kPreselection(&slc) && kHasMuon(&slc) && kHadronicContainment(&slc) && kLeadingProtonThreshold(&slc) && kScndProtonThreshold(&slc)); }
bool cutCountPrimaryPFPs(const caf::SRSliceProxy &slc) { return (kPreselection(&slc) && kHasMuon(&slc) && kHadronicContainment(&slc) && kLeadingProtonThreshold(&slc) && kScndProtonThreshold(&slc) && kCountPrimaryPFPs(&slc)/ && !kThreeRecoProtons(&slc) && !kHasSidebandPion(&slc)/); }
bool cutNoExtraMIP(const caf::SRSliceProxy &slc) { return (kPreselection(&slc) && kHasMuon(&slc) && kHadronicContainment(&slc) && kLeadingProtonThreshold(&slc) && kScndProtonThreshold(&slc) && kNoExtraMIP(&slc)); }
bool cutNoExtraShower(const caf::SRSliceProxy &slc) { return (kPreselection(&slc) && kHasMuon(&slc) && kHadronicContainment(&slc) && kLeadingProtonThreshold(&slc) && kScndProtonThreshold(&slc) && kNoExtraMIPCut(&slc) && kNoExtraShowerCut(&slc)); }
bool cutExtraPrimaryLinFitLength(const caf::SRSliceProxy &slc) { return (kPreselection(&slc) && kHasMuon(&slc) && kHadronicContainment(&slc) && kLeadingProtonThreshold(&slc) && kScndProtonThreshold(&slc) && kExtraPrimaryLinFitLengthCut(&slc)); }
bool cutDynamicLength(const caf::SRSliceProxy &slc) { return (kPreselection(&slc) && kHasMuon(&slc) && kHadronicContainment(&slc) && kLeadingProtonThreshold(&slc) && kScndProtonThreshold(&slc) && kDynamicLengthCut(&slc)); }
bool cutNoPion(const caf::SRSliceProxy &slc) { return (kPreselection(&slc) && kHasMuon(&slc) && kHadronicContainment(&slc) && kLeadingProtonThreshold(&slc) && kScndProtonThreshold(&slc) && kDynamicLengthCut(&slc) && !kHasSidebandPion(&slc)); }

const std::vector<std::pair< TString, std::function<bool(const caf::SRTrueInteractionProxy&)> >> probabilityCuts = {
{"NoCut", cutNoCut},
{"Contained", cutContained},
{"Stopping", cutStopping},
{"AllStopping", cutAllStopping},

};

const std::vector<std::pair< TString, std::function<bool(const caf::SRSliceProxy&)> >> efficiencyCuts = {
{"Preselection", cutPreselection},
{"Muon", cutMuon},
{"HadronicContainment", cutHadronicContainment},
{"LeadingProton", cutLeadingProton},
{"ScndProton", cutScndProton},
{"Selected", cutCountPrimaryPFPs},
{"NoExtraMIP", cutNoExtraMIP},
{"NoExtraShower", cutNoExtraShower},
{"ExtraPrimaryLinFitLength", cutExtraPrimaryLinFitLength},
{"DynamicLength", cutDynamicLength},
{"NoPion", cutNoPion},

};

std::vector<int> getSignalIndices ( const caf::SRSpillProxy* sr, bool containedLepton=containedMuon, bool containedProton=containedProtons,
    double muThresh=muonThresh, double pThresh=protonThresh, double pMaxThresh=protonMax, double piThresh=pionThresh, double pi0Thresh=nPionThresh ) {
  std::vector<int> signalNus;

  for ( auto const& nu : sr->mc.nu ) {
    if ( abs(nu.pdg) != 14 ||
                                 !nu.iscc ||
                                 std::isnan(nu.position.x) || std::isnan(nu.position.y) || std::isnan(nu.position.z) ||
                                 !isInFV(nu.position.x,nu.position.y,nu.position.z) )
      continue;

    unsigned int nMu(0), nSigMu(0), nP(0), nPMax(0), nPi0(0), nChgPi(0);

    for ( auto const& prim : nu.prim ) {
      if ( prim.start_process != 0 ) continue;
      TVector3 p3(prim.startp.x, prim.startp.y, prim.startp.z);
      double momentum = p3.Mag();
      if ( abs(prim.pdg) == 13 ) {
        if ( containedLepton && !prim.contained ) continue;
        if ( momentum < muThresh ) continue;
        nMu+=1;
      }
      if ( abs(prim.pdg) == 2212 ){
        if ( containedProton && !prim.contained ) continue;
        if ( momentum < pThresh ) continue;
        if ( momentum > pMaxThresh ) nPMax+=1;
        else nP+=1;
      }
      if ( abs(prim.pdg) == 111 ){
        if ( momentum < pi0Thresh ) continue;
        nPi0+=1;
      }
      if ( abs(prim.pdg) == 211 ){
        if ( momentum < piThresh ) continue;
        nChgPi+=1;
      }
    }
    if ( nMu == 1 && nP >= 2 && nChgPi == 0 && nPi0 == 0 && nPMax == 0) {
      signalNus.push_back(nu.index);
    }
  }

  return signalNus;
}

std::vector<double> getTrueVarVectorOld(const caf::SRSpillProxy* sr, std::function<double(const caf::SRTrueInteractionProxy&)> trueVar,
    std::function<bool(const caf::SRTrueInteractionProxy&)> qualified) {
  std::vector<double> vals;
  std::vector<int> signalNus = getSignalIndices(sr);
  if ( signalNus.size() == 0 ) return vals;

  for ( auto const& nu : sr->mc.nu ) {
    if ( std::find(signalNus.begin(), signalNus.end(), nu.index) == signalNus.end() ) continue;
    if ( qualified(nu) ) vals.push_back( trueVar(nu) );
    else vals.push_back( -9999. ); //Ensure one entry per signal nu
  }

  if ( vals.size() != signalNus.size() ) std::cout << std::endl << "Size mismatch in True Interaction var! vals.size(): " << vals.size() << ", signalNus.size(): " << signalNus.size() << std::endl;
  assert (vals.size() == signalNus.size() );
  return vals;
  }

std::vector<double> getTrueVarVectorSelected(const caf::SRSpillProxy* sr, std::function<double(const caf::SRTrueInteractionProxy&) > trueVar, 
    std::function<bool(const caf::SRSliceProxy&)> selected) {
  std::vector<double> vals;
  std::vector<int> signalNus = getSignalIndices(sr);
  if ( signalNus.size() == 0 ) return vals;

  unsigned totalSignalNus = signalNus.size();
  unsigned seenSignalNus = 0;
  for ( auto const& nu : sr->mc.nu ) {
    if ( std::find(signalNus.begin(), signalNus.end(), nu.index) == signalNus.end() ) continue;
    seenSignalNus++;
    for ( unsigned i = 0; i < sr->slc.size(); i++ ) {
      auto const &slc = sr->slc.at(i);
      if ( slc.truth.index != nu.index ) continue;
      if ( selected(slc) ) { vals.push_back( trueVar(slc.truth) ); break; }
    }
    if ( vals.size() < seenSignalNus ) vals.push_back( -9999. ); //Ensure one entry per signal nu
  }

  if ( vals.size() != signalNus.size() ) std::cout << std::endl << "Size mismatch in Slice var! vals.size(): " << vals.size() << ", signalNus.size(): " << signalNus.size() << std::endl;
  assert (vals.size() == signalNus.size() );
  return vals;
}

std::vector<double> getTrueVarVector(const caf::SRSpillProxy* sr, std::function<double(const caf::SRTrueInteractionProxy&)> trueVar) {
  std::vector<double> vals;
  std::vector<int> signalNus = getSignalIndices(sr);
  if ( signalNus.size() == 0 ) return vals;

  for ( auto const& nu : sr->mc.nu ) {
    if ( std::find(signalNus.begin(), signalNus.end(), nu.index) == signalNus.end() ) continue;
    vals.push_back( trueVar(nu) );
  }

  return vals;
  }

std::vector<double> applyTruthCut(const caf::SRSpillProxy* sr, std::function<bool(const caf::SRTrueInteractionProxy&)> qualified) {
  std::vector<double> vals;
  std::vector<int> signalNus = getSignalIndices(sr);
  if ( signalNus.size() == 0 ) return vals;

  for ( auto const& nu : sr->mc.nu ) {
    if ( std::find(signalNus.begin(), signalNus.end(), nu.index) == signalNus.end() ) continue;
    if ( qualified(nu) ) vals.push_back( 1. );
    else vals.push_back( 0. );
  }

  return vals;
  }

std::vector<double> applyRecoCut(const caf::SRSpillProxy* sr, std::function<bool(const caf::SRSliceProxy&)> selected) {
  std::vector<double> vals;
  std::vector<int> signalNus = getSignalIndices(sr);
  if ( signalNus.size() == 0 ) return vals;

  for ( auto const& nu : sr->mc.nu ) {
    if ( std::find(signalNus.begin(), signalNus.end(), nu.index) == signalNus.end() ) continue;
    bool selectedSlc = false;
    for ( const auto &slc : sr->slc ) {
      if ( slc.truth.index != nu.index ) continue;
      if ( selected(slc) ) { selectedSlc = true; break; }
    }
    vals.push_back(selectedSlc);
  }

  return vals;
}

SpillMultiVar getTrueSpillMultiVarOld(std::function<double(const caf::SRTrueInteractionProxy&) > trueVar, std::function<bool(const caf::SRTrueInteractionProxy&)> qualified) {
  return SpillMultiVar( [=](const caf::SRSpillProxy *sr) -> std::vector<double> {
    return getTrueVarVectorOld(sr, trueVar, qualified);
  });
}

SpillMultiVar getTrueSpillMultiVarSelected(std::function<double(const caf::SRTrueInteractionProxy&) > trueVar, std::function<bool(const caf::SRSliceProxy&)> selected) {
  return SpillMultiVar( [=](const caf::SRSpillProxy *sr) -> std::vector<double> {
    return getTrueVarVectorSelected(sr, trueVar, selected);
  });
}

SpillMultiVar getTrueSpillMultiVar(std::function<double(const caf::SRTrueInteractionProxy&) > trueVar) {
  return SpillMultiVar( [=](const caf::SRSpillProxy *sr) -> std::vector<double> {
    return getTrueVarVector(sr, trueVar);
  });
}

SpillMultiVar applyTruthCutSpillMultiVar(std::function<bool(const caf::SRTrueInteractionProxy&)> qualified) {
  return SpillMultiVar( [=](const caf::SRSpillProxy *sr) -> std::vector<double> {
    return applyTruthCut(sr, qualified);
  });
}

SpillMultiVar applyRecoCutSpillMultiVar(std::function<bool(const caf::SRSliceProxy&)> selected) {
  return SpillMultiVar( [=](const caf::SRSpillProxy *sr) -> std::vector<double> {
    return applyRecoCut(sr, selected);
  });
}

*/

const Var kCategory([](const caf::SRSliceProxy* slc) -> int {
  int cat = -1;

  if ( k1mu2p0pi(slc) ) cat = 0;
  else if ( k1mu2pNpi(slc) ) cat = 1;
  else if ( k1mu3p0pi(slc) ) cat = 2;
  else if ( k1mu3pNpi(slc) ) cat = 3;
  else if ( k1mu1p0pi(slc) ) cat = 4;
  else if ( k1mu1pNpi(slc) ) cat = 5;
  else if ( k1mu0p0pi(slc) ) cat = 6;
  else if ( k1mu0pNpi(slc) ) cat = 7;
  else if ( kCCOther(slc) ) cat = 8;
  else if ( kIsNC(slc) ) cat = 9;
  else if ( kOOFV(slc) ) cat = 10;
  else if ( kIsCosmic(slc) ) cat = 11;
  else std::cout << std::endl << "!!!!!!!!!! NO CATEGORY FOR EVENT! !!!!!!!!!!!!!";
  return cat;
  });

//Concatenation of kCategory for easier plotting, listed in stack order
const Var kClassLabel([](const caf::SRSliceProxy* slc) -> int {
  int cat = kCategory(slc);

  if ( cat == 2 ) return 0;                  //CC>2p0pi
  else if ( cat == 0 ) return 1;             //CC2p0pi
  else if ( cat == 1 || cat == 3 ) return 2; //CC>1pNpi
  else if ( cat < 8 ) return 3;              //CC<2p
  else if ( cat == 8 ) return 4;             //CC Other
  else if ( cat == 10 ) return 5;            //OOFV
  else if ( cat == 9 ) return 6;             //NC
  else return 7;                             //Cosmic
  });

const Var kPrintExtraPrimaries([](const caf::SRSliceProxy* slc) -> float {
    int category = kCategory(slc);


    std::cout << std::endl << "Category: " << category << ", NoExtraMIP: " << kNoExtraMIP(slc) << ", NoExtraShower: " << kNoExtraShower(slc) << ", ExtraPrimaryLinFitLength: " << kExtraPrimaryLinFitLength(slc);
    if ( !kNoExtraMIP(slc) ) std::cout << ", ExtraMIPIdx: " << kExtraMIPIdx(slc);
    if ( !kNoExtraShower(slc) ) std::cout << ", ExtraShowerIdx: " << kExtraShowerIdx(slc);

    std::cout << std::endl << "~~~~~~~~~~~~~~~~~ G4 Primaries";
    for ( const auto &prim : slc->truth.prim ) {
      if ( prim.start_process != 0 ) continue;
      std::cout << std::endl << "G4ID: " << prim.G4ID << ", PDG: " << prim.pdg << ", momentum: " << sqrt( prim.startp.x*prim.startp.x + prim.startp.y*prim.startp.y + prim.startp.z*prim.startp.z ) << ", end_process: " << prim.end_process;
    }

    std::vector<unsigned> trackedPrimaries;
    if ( kRecoMuonIdx(slc) != -1 ) trackedPrimaries.push_back(kRecoMuonIdx(slc));
    if ( kRecoProtonIdx(slc) != -1 ) trackedPrimaries.push_back(kRecoProtonIdx(slc));
    if ( kScndProtonIdx(slc) != -1 ) trackedPrimaries.push_back(kScndProtonIdx(slc));

    std::cout << std::endl << "~~~~~~~~~~~~~~~~~ Pandora Primaries" << std::endl;
    for ( unsigned i = 1; i < slc->reco.pfp.size(); i++ ) {
      const auto &pfp = slc->reco.pfp.at(i);
      if ( pfp.parent_is_primary ) {
        if ( std::find(trackedPrimaries.begin(), trackedPrimaries.end(), i) != trackedPrimaries.end() ) std::cout << "Tracked Primary, ";
        else std::cout << "Extra Primary, ";
      }
      else std::cout << "Nonprimary, ";
      if ( !isnan(pfp.shw.truth.p.G4ID) ) std::cout << "G4ID: " << pfp.shw.truth.p.G4ID << ", ";
      if ( !isnan(pfp.shw.truth.p.pdg) ) std::cout << "pdg: " << pfp.shw.truth.p.pdg << ", ";
      if ( !isnan(pfp.shw.truth.p.parent) ) std::cout << "parent: " << pfp.shw.truth.p.parent << ", ";
      if ( !isnan(static_cast<int>(pfp.shw.truth.p.end_process)) ) std::cout << "end_process: " << pfp.shw.truth.p.end_process << ", ";
      if ( !isnan(pfp.shw.truth.p.length) ) std::cout << "true length: " << pfp.shw.truth.p.length << ", ";
      if ( !isnan(pfp.shw.truth.nmatches) ) std::cout << "nmatches: " << pfp.shw.truth.nmatches << ", ";
      if ( !isnan(pfp.pfochar.linfitlen) ) std::cout << "lin fit length: " << pfp.pfochar.linfitlen << ", ";
      if ( !isnan(pfp.shw.start.z) ) std::cout << "start pos: " << pfp.shw.start.x << ", " << pfp.shw.start.y << ", " << pfp.shw.start.z << ", " ;
      if ( !isnan(pfp.shw.bestplane_dEdx) ) std::cout << "bestplane_dEdx: " << pfp.shw.bestplane_dEdx << ", ";
      if ( !isnan(pfp.shw.bestplane_energy) ) std::cout << "bestplane_energy: " << pfp.shw.bestplane_energy << ", ";
      if ( !isnan(pfp.shw.conversion_gap) ) std::cout << "conversion_gap: " << pfp.shw.conversion_gap << ", ";
      if ( !isnan(pfp.shw.len) ) std::cout << "len: " << pfp.shw.len << ", ";
      if ( !isnan(pfp.shw.open_angle) ) std::cout << "open_angle: " << pfp.shw.open_angle << ", ";

      if ( !isnan(pfp.trk.len) ) std::cout << "len: " << pfp.trk.len << ", ";
      if ( !isnan(pfp.trk.chi2pid[2].chi2_proton) ) std::cout << "chi2_proton: " << pfp.trk.chi2pid[2].chi2_proton << ", ";
      if ( !isnan(pfp.trk.chi2pid[2].chi2_muon) ) std::cout << "chi2_muon: " << pfp.trk.chi2pid[2].chi2_muon << ", ";
      if ( !isnan(pfp.trk.chi2pid[2].pid_ndof) ) std::cout << "pid_ndof: " << pfp.trk.chi2pid[2].pid_ndof << ", ";
      std::cout << std::endl;
    }

    return 0.;
  });

const Var kContained([](const caf::SRSliceProxy* slc) -> int {
  if ( kRecoMuonContained(slc) ) return 1;
  else return 0;
  });

const Var kGENIEMode([](const caf::SRSliceProxy* slc) -> int {
    if ( slc->truth.index < 0 ) return -1;
    return slc->truth.genie_mode;
  });









std::vector<std::string> GetGENIEMultisigmaKnobNames(){

  return {
"RPA_CCQE",
"CoulombCCQE",
"NormCCMEC",
"NormNCMEC",
//"DecayAngMEC", --> MirrorSyst!
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

std::vector<std::string> GetDetectorKnobNames() {
  std::vector<std::string> knobs = {
"kNuMIXSecFrontIndPlaneGainSyst",
"kNuMIXSecFrontIndPlaneNoiseSyst",
"kNuMIXSecFrontIndPlaneNoiseSyst",
"kNuMIXSecMiddleIndPlaneTransparencySyst",
  };
  return knobs;
}

}//END NAMESPACE
