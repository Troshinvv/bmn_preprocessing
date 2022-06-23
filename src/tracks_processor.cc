//
// Created by mikhail on 11/30/20.
//

#include "tracks_processor.h"

#include "TLorentzVector.h"
#include "TDatabasePDG.h"
#include "TMath.h"
#include "TaskManager.hpp"

#include <AnalysisTree/DataHeader.hpp>
#include <random>
#include <TH2F.h>

namespace AnalysisTree {

void TracksProcessor::Init() {
  auto man = TaskManager::GetInstance();
  auto chain = man->GetChain();

  AddInputBranch("TpcTracks");
  AddInputBranch("RecoEvent");
  in_sim_particles_ = chain->GetBranch("TpcTracks");
  in_sim_particles_.Freeze();

  auto in_sim_event_conf = chain->GetConfiguration()->GetBranchConfig("RecoEvent");
  auto out_sim_event_conf =in_sim_event_conf.Clone("RecoEventExt", DetType::kEventHeader);
  out_sim_event_conf.AddField<float>("centrality", "centrality, %");

  out_sim_event_ = Branch(out_sim_event_conf);
  out_sim_event_.SetMutable();

  man->AddBranch(&out_sim_event_);

}

void TracksProcessor::Exec() {
  using AnalysisTree::Particle;
  this->LoopRecTracks();
}

void TracksProcessor::LoopRecTracks() {
  using AnalysisTree::Particle;
  auto field_in_eta = in_sim_particles_.GetField("eta");
  auto field_in_pT = in_sim_particles_.GetField("pT");
  auto field_in_nhits = in_sim_particles_.GetField("nhits");
  auto field_in_dca_x = in_sim_particles_.GetField("dca_x");
  auto field_in_dca_y = in_sim_particles_.GetField("dca_y");
  auto field_in_dca_z = in_sim_particles_.GetField("dca_z");
  auto Mult=0;
  for (size_t i=0; i<in_sim_particles_.size(); ++i) { 
    auto in_particle = in_sim_particles_[i];
    auto eta = in_particle[field_in_eta];
    auto pT = in_particle[field_in_pT];
    auto nhits = in_particle[field_in_nhits];
    auto dca_x = in_particle[field_in_dca_x];
    auto dca_y = in_particle[field_in_dca_y];
    auto dca_z = in_particle[field_in_dca_z];
    if(nhits>16 && abs(eta)<0.5 && pT>0.15 && dca_x<1. && dca_y<1. && dca_z<1.)
    {
	    Mult=Mult+1;
    }
  }
  Float_t minCentPercent [14] = { 0, 5, 10, 15, 20, 25, 30, 35, 40, 50, 60, 70, 80, 90};
Float_t maxCentPercent [14] = { 5, 10, 15, 20, 25, 30, 35, 40, 50, 60, 70, 80, 90, 100};
Int_t minMult [14] = { 138, 116, 97, 82, 69, 58, 48, 40, 27, 17, 10, 6, 3, 1};
Int_t maxMult [14] = { 217, 138, 116, 97, 82, 69, 58, 48, 40, 27, 17, 10, 6, 3};
Int_t centBin = -1;
auto centrality = 0;
        for (Int_t i = 0; i < 14; i++)
        {
                if (Mult >= minMult[i] && Mult < maxMult[i])
                        centBin = i;
        }
        if (centBin == -1) centrality=-1;
	else 
		centrality=((maxCentPercent[centBin] - minCentPercent[centBin]) / 2. + minCentPercent[centBin]);
auto out_value=out_sim_event_.NewChannel();
auto field_centrality = out_sim_event_.GetField("centrality");
out_value.SetValue(field_centrality,(float) centrality);
}
} // namespace AnalysisTree
