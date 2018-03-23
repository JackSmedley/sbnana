#include <TBranch.h>
#include <TFile.h>
#include <TTree.h>
#include "gallery/ValidHandle.h"
#include "gallery/Handle.h"
#include "canvas/Utilities/InputTag.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "nusimdata/SimulationBase/GTruth.h"
#include "json/json.h"
#include "Event.hh"
#include "Loader.hh"
#include "ProcessorBase.hh"

namespace core {

ProcessorBase::ProcessorBase() : fEventIndex(0), fOutputFilename("output.root") {}


ProcessorBase::~ProcessorBase() {}


void ProcessorBase::FillTree() {
  fTree->Fill();
  fEventIndex++;
}

void ProcessorBase::EventCleanup() {
  fEvent->interactions.clear();
}


void ProcessorBase::Initialize(char* config) {
  Json::Value* cfg = LoadConfig(config);
  Initialize(cfg);
}


void ProcessorBase::Setup(char* config) {
  Json::Value* cfg = LoadConfig(config);
  Setup(cfg);
}


void ProcessorBase::Setup(Json::Value* config) {
  // Load configuration parameters
  fTruthTag = { "generator" };

  if (config) {
    fTruthTag = { config->get("MCTruthTag", "generator").asString() };
    fOutputFilename = config->get("OutputFile", "output.root").asString();
  }

  // Open the output file and create the standard event tree
  fOutputFile = TFile::Open(fOutputFilename.c_str(), "recreate");
  fTree = new TTree("sbnana", "SBN Analysis Tree");
  fEvent = new Event();
  fTree->Branch("events", &fEvent);
}


void ProcessorBase::Teardown() {
  // Write the standard tree and close the output file
  fOutputFile->cd();
  fTree->Write();
  fOutputFile->Close();
}


void ProcessorBase::BuildEventTree(gallery::Event& ev) {
  // Get MCTruth information
  auto const& mctruths = *ev.getValidHandle<std::vector<simb::MCTruth> >(fTruthTag);
  
  gallery::Handle<std::vector<simb::GTruth> > gtruths_handle;
  ev.getByLabel(fTruthTag,gtruths_handle);
  bool genie_truth_is_valid = gtruths_handle.isValid();

  fTree->GetEntry(fEventIndex);

  // Populate event tree
  for (size_t i=0; i<mctruths.size(); i++) {
    Event::Interaction interaction;
    auto const& mctruth = mctruths.at(i);

    // Neutrino
    const simb::MCNeutrino& nu = mctruth.GetNeutrino();
    interaction.neutrino.isnc =   nu.CCNC()  && (nu.Mode()!=simb::kWeakMix);
    interaction.neutrino.iscc = (!nu.CCNC()) && (nu.Mode()!=simb::kWeakMix);
    interaction.neutrino.pdg = nu.Nu().PdgCode();
    interaction.neutrino.targetPDG = nu.Target();
    interaction.neutrino.genie_intcode = nu.Mode();
    interaction.neutrino.bjorkenX = nu.X();
    interaction.neutrino.inelasticityY = nu.Y();
    interaction.neutrino.q2 = nu.QSqr();
    interaction.neutrino.w = nu.W();
    interaction.neutrino.energy = nu.Nu().EndMomentum().Energy();
    interaction.neutrino.momentum = nu.Nu().EndMomentum().Vect();

    // Primary lepton
    const simb::MCParticle& lepton = nu.Lepton();
    interaction.lepton.pdg = lepton.PdgCode();
    interaction.lepton.energy = lepton.Momentum(0).Energy();
    interaction.lepton.momentum = lepton.Momentum(0).Vect();

    // Hadronic system
    for (int iparticle=0; iparticle<mctruth.NParticles(); iparticle++) {
      const simb::MCParticle& particle = mctruth.GetParticle(iparticle);

      if (particle.Process() != "primary") {
        continue;
      }

      Event::FinalStateParticle fsp;
      fsp.pdg = particle.PdgCode();
      fsp.energy = particle.Momentum(0).Energy();
      fsp.momentum = particle.Momentum(0).Vect();

      interaction.finalstate.push_back(fsp);
    }
    
    // GENIE specific
    if (genie_truth_is_valid) {
      auto const& gtruth = gtruths_handle->at(i);
      TLorentzVector q_nucframe = gtruth.fFShadSystP4-gtruth.fHitNucP4;
      interaction.neutrino.modq = q_nucframe.P();
      interaction.neutrino.q0   = q_nucframe.E();
    }
    
    fEvent->interactions.push_back(interaction);
  }
}

}  // namespace core

