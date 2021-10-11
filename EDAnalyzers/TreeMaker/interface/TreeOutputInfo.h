#ifndef TreeOutputInfo_H
#define TreeOutputInfo_H 1

#include <TH1F.h>
#include <TH2F.h>
#include <TMatrixD.h>
#include <TROOT.h>
#include <TTree.h>
#include <TVectorD.h>
#include <stdlib.h>

#include <iostream>
#include <map>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "EDAnalyzers/TreeMaker/interface/Constants.h"

namespace TreeOutputInfo {
  class TreeOutput {
  public:
    TTree *tree;

    // Run info //
    ULong64_t runNumber;
    ULong64_t eventNumber;
    ULong64_t luminosityNumber;
    ULong64_t bunchCrossingNumber;

    // Gen electron //
    int genEl_n;
    std::vector<float> v_genEl_E;
    std::vector<float> v_genEl_px;
    std::vector<float> v_genEl_py;
    std::vector<float> v_genEl_pz;
    std::vector<float> v_genEl_pT;
    std::vector<float> v_genEl_eta;
    std::vector<float> v_genEl_phi;

    // Gen photon //
    int genPh_n;
    std::vector<float> v_genPh_E;
    std::vector<float> v_genPh_px;
    std::vector<float> v_genPh_py;
    std::vector<float> v_genPh_pz;
    std::vector<float> v_genPh_pT;
    std::vector<float> v_genPh_eta;
    std::vector<float> v_genPh_phi;

    // Pileup //
    int pileup_n;

    // Rho //
    float rho;

    float simHit_n;
    std::vector<float> v_simHit_E;
    std::vector<float> v_simHit_x;
    std::vector<float> v_simHit_y;
    std::vector<float> v_simHit_z;
    std::vector<float> v_simHit_eta;
    std::vector<float> v_simHit_phi;
    std::vector<float> v_simHit_ET;
    std::vector<float> v_simHit_layer;
    std::vector<float> v_simHit_zside;
    std::vector<float> v_simHit_isCaloParticleMatched;
    std::vector<float> v_simHit_matchedSimClusIndex;
    std::vector<float> v_simHit_detector;

    float recHit_n;
    std::vector<float> v_recHit_E;
    std::vector<float> v_recHit_x;
    std::vector<float> v_recHit_y;
    std::vector<float> v_recHit_z;
    std::vector<float> v_recHit_eta;
    std::vector<float> v_recHit_phi;
    std::vector<float> v_recHit_ET;
    std::vector<float> v_recHit_layer;
    std::vector<float> v_recHit_zside;
    std::vector<float> v_recHit_detector;

    std::vector<float> v_recHit_matchedSimHitIndex;
    std::vector<float> v_recHit_matchedSimClusIndex;
    std::vector<float> v_recHit_isCaloParticleMatched;

    std::vector<float> v_recHit_iType;
    std::vector<float> v_recHit_iCell1;
    std::vector<float> v_recHit_iCell2;

    std::vector<float> v_recHit_SiThickness;

    float caloParticle_n;
    std::vector<float> v_caloParticle_E;
    std::vector<float> v_caloParticle_px;
    std::vector<float> v_caloParticle_py;
    std::vector<float> v_caloParticle_pz;
    std::vector<float> v_caloParticle_pT;
    std::vector<float> v_caloParticle_eta;
    std::vector<float> v_caloParticle_phi;
    std::vector<float> v_caloParticle_pdgid;

    char name[500];

    TreeOutput(std::string details, edm::Service<TFileService> fs) {
      //printf("Loading custom ROOT dictionaries. \n");
      //gROOT->ProcessLine(".L EDAnalyzers/TreeMaker/interface/CustomRootDict.cc+");
      //printf("Loaded custom ROOT dictionaries. \n");

      tree = fs->make<TTree>(details.c_str(), details.c_str());

      // Run info //
      tree->Branch("runNumber", &runNumber);
      tree->Branch("eventNumber", &eventNumber);
      tree->Branch("luminosityNumber", &luminosityNumber);
      tree->Branch("bunchCrossingNumber", &bunchCrossingNumber);

      // Gen electron //
      sprintf(name, "genEl_n");
      tree->Branch(name, &genEl_n);

      sprintf(name, "genEl_E");
      tree->Branch(name, &v_genEl_E);

      sprintf(name, "genEl_px");
      tree->Branch(name, &v_genEl_px);

      sprintf(name, "genEl_py");
      tree->Branch(name, &v_genEl_py);

      sprintf(name, "genEl_pz");
      tree->Branch(name, &v_genEl_pz);

      sprintf(name, "genEl_pT");
      tree->Branch(name, &v_genEl_pT);

      sprintf(name, "genEl_eta");
      tree->Branch(name, &v_genEl_eta);

      sprintf(name, "genEl_phi");
      tree->Branch(name, &v_genEl_phi);

      // Gen photon //
      sprintf(name, "genPh_n");
      tree->Branch(name, &genPh_n);

      sprintf(name, "genPh_E");
      tree->Branch(name, &v_genPh_E);

      sprintf(name, "genPh_px");
      tree->Branch(name, &v_genPh_px);

      sprintf(name, "genPh_py");
      tree->Branch(name, &v_genPh_py);

      sprintf(name, "genPh_pz");
      tree->Branch(name, &v_genPh_pz);

      sprintf(name, "genPh_pT");
      tree->Branch(name, &v_genPh_pT);

      sprintf(name, "genPh_eta");
      tree->Branch(name, &v_genPh_eta);

      sprintf(name, "genPh_phi");
      tree->Branch(name, &v_genPh_phi);

      // Pileup //
      sprintf(name, "pileup_n");
      tree->Branch(name, &pileup_n);

      // Rho //
      sprintf(name, "rho");
      tree->Branch(name, &rho);

      // Sim-hit //
      sprintf(name, "simHit_n");
      tree->Branch(name, &simHit_n);

      sprintf(name, "simHit_E");
      tree->Branch(name, &v_simHit_E);

      sprintf(name, "simHit_x");
      tree->Branch(name, &v_simHit_x);

      sprintf(name, "simHit_y");
      tree->Branch(name, &v_simHit_y);

      sprintf(name, "simHit_z");
      tree->Branch(name, &v_simHit_z);

      sprintf(name, "simHit_eta");
      tree->Branch(name, &v_simHit_eta);

      sprintf(name, "simHit_phi");
      tree->Branch(name, &v_simHit_phi);

      sprintf(name, "simHit_ET");
      tree->Branch(name, &v_simHit_ET);

      sprintf(name, "simHit_layer");
      tree->Branch(name, &v_simHit_layer);

      sprintf(name, "simHit_zside");
      tree->Branch(name, &v_simHit_zside);

      sprintf(name, "simHit_isCaloParticleMatched");
      tree->Branch(name, &v_simHit_isCaloParticleMatched);

      sprintf(name, "simHit_matchedSimClusIndex");
      tree->Branch(name, &v_simHit_matchedSimClusIndex);

      sprintf(name, "recHit_detector");
      tree->Branch(name, &v_recHit_detector);

      // Rec-hit //
      sprintf(name, "recHit_n");
      tree->Branch(name, &recHit_n);

      sprintf(name, "recHit_E");
      tree->Branch(name, &v_recHit_E);

      sprintf(name, "recHit_x");
      tree->Branch(name, &v_recHit_x);

      sprintf(name, "recHit_y");
      tree->Branch(name, &v_recHit_y);

      sprintf(name, "recHit_z");
      tree->Branch(name, &v_recHit_z);

      sprintf(name, "recHit_eta");
      tree->Branch(name, &v_recHit_eta);

      sprintf(name, "recHit_phi");
      tree->Branch(name, &v_recHit_phi);

      // sprintf(name, "recHit_ET");
      // tree->Branch(name, &v_recHit_ET);

      sprintf(name, "recHit_layer");
      tree->Branch(name, &v_recHit_layer);

      sprintf(name, "recHit_zside");
      tree->Branch(name, &v_recHit_zside);

      sprintf(name, "recHit_detector");
      tree->Branch(name, &v_recHit_detector);

      // sprintf(name, "recHit_matchedSimHitIndex");
      // tree->Branch(name, &v_recHit_matchedSimHitIndex);

      // sprintf(name, "recHit_matchedSimClusIndex");
      // tree->Branch(name, &v_recHit_matchedSimClusIndex);

      // sprintf(name, "recHit_isCaloParticleMatched");
      // tree->Branch(name, &v_recHit_isCaloParticleMatched);

      // sprintf(name, "recHit_iType");
      // tree->Branch(name, &v_recHit_iType);

      // sprintf(name, "recHit_iCell1");
      // tree->Branch(name, &v_recHit_iCell1);

      // sprintf(name, "recHit_iCell2");
      // tree->Branch(name, &v_recHit_iCell2);

      // sprintf(name, "recHit_SiThickness");
      // tree->Branch(name, &v_recHit_SiThickness);

      // // Calo-particle //
      // sprintf(name, "caloParticle_n");
      // tree->Branch(name, &caloParticle_n);

      // sprintf(name, "caloParticle_E");
      // tree->Branch(name, &v_caloParticle_E);

      // sprintf(name, "caloParticle_px");
      // tree->Branch(name, &v_caloParticle_px);

      // sprintf(name, "caloParticle_py");
      // tree->Branch(name, &v_caloParticle_py);

      // sprintf(name, "caloParticle_pz");
      // tree->Branch(name, &v_caloParticle_pz);

      // sprintf(name, "caloParticle_pT");
      // tree->Branch(name, &v_caloParticle_pT);

      // sprintf(name, "caloParticle_eta");
      // tree->Branch(name, &v_caloParticle_eta);

      // sprintf(name, "caloParticle_phi");
      // tree->Branch(name, &v_caloParticle_phi);

      // sprintf(name, "caloParticle_pdgid");
      // tree->Branch(name, &v_caloParticle_pdgid);
    }

    void fill() { tree->Fill(); }

    void clear() {
      // Gen electron //
      genEl_n = 0;
      v_genEl_E.clear();
      v_genEl_px.clear();
      v_genEl_py.clear();
      v_genEl_pz.clear();
      v_genEl_pT.clear();
      v_genEl_eta.clear();
      v_genEl_phi.clear();

      // Gen photon //
      genPh_n = 0;
      v_genPh_E.clear();
      v_genPh_px.clear();
      v_genPh_py.clear();
      v_genPh_pz.clear();
      v_genPh_pT.clear();
      v_genPh_eta.clear();
      v_genPh_phi.clear();

      // Pileup //
      pileup_n = 0;

      // Rho //
      rho = 0;

      // Sim-hit //
      simHit_n = 0;
      v_simHit_E.clear();
      v_simHit_x.clear();
      v_simHit_y.clear();
      v_simHit_z.clear();
      v_simHit_eta.clear();
      v_simHit_phi.clear();
      v_simHit_ET.clear();
      v_simHit_layer.clear();
      v_simHit_zside.clear();
      v_simHit_isCaloParticleMatched.clear();
      v_simHit_matchedSimClusIndex.clear();

      // Rec-hit //
      recHit_n = 0;
      v_recHit_E.clear();
      v_recHit_x.clear();
      v_recHit_y.clear();
      v_recHit_z.clear();
      v_recHit_eta.clear();
      v_recHit_phi.clear();
      v_recHit_ET.clear();
      v_recHit_layer.clear();
      v_recHit_zside.clear();
      v_recHit_matchedSimHitIndex.clear();
      v_recHit_matchedSimClusIndex.clear();
      v_recHit_isCaloParticleMatched.clear();
      v_recHit_iType.clear();
      v_recHit_iCell1.clear();
      v_recHit_iCell2.clear();
      v_recHit_detector.clear();

      v_recHit_SiThickness.clear();

      // Calo-particle
      caloParticle_n = 0;
      v_caloParticle_E.clear();
      v_caloParticle_px.clear();
      v_caloParticle_py.clear();
      v_caloParticle_pz.clear();
      v_caloParticle_pT.clear();
      v_caloParticle_eta.clear();
      v_caloParticle_phi.clear();
      v_caloParticle_pdgid.clear();
    }
  };
}  // namespace TreeOutputInfo

#endif
