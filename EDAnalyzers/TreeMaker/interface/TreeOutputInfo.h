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
    // ULong64_t bunchCrossingNumber;

    // Gen Particle //
    int genPart_n;
    std::vector<float> v_genPart_E;
    std::vector<float> v_genPart_px;
    std::vector<float> v_genPart_py;
    std::vector<float> v_genPart_pz;
    std::vector<float> v_genPart_pT;
    std::vector<float> v_genPart_eta;
    std::vector<float> v_genPart_phi;
    std::vector<int> v_genPart_pdgId;

    // // Pileup //
    // int pileup_n;

    // Rho //
    // float rho;

    float simHit_n;
    std::vector<float> v_simHit_E;
    std::vector<float> v_simHit_x;
    std::vector<float> v_simHit_y;
    std::vector<float> v_simHit_z;
    std::vector<float> v_simHit_eta;
    std::vector<float> v_simHit_phi;
    // std::vector<float> v_simHit_ET;
    std::vector<float> v_simHit_layer;
    std::vector<float> v_simHit_zside;
    // std::vector<float> v_simHit_isCaloParticleMatched;
    // std::vector<float> v_simHit_matchedSimClusIndex;
    std::vector<uint32_t> v_simHit_detId;
    std::vector<float> v_simHit_detector;

    float recHit_n;
    std::vector<float> v_recHit_E;
    std::vector<float> v_recHit_x;
    std::vector<float> v_recHit_y;
    std::vector<float> v_recHit_z;
    std::vector<float> v_recHit_eta;
    std::vector<float> v_recHit_phi;
    // std::vector<float> v_recHit_ET;
    std::vector<float> v_recHit_layer;
    std::vector<float> v_recHit_zside;
    std::vector<uint32_t> v_recHit_detId;
    std::vector<float> v_recHit_detector;

    // std::vector<float> v_recHit_matchedSimHitIndex;
    // std::vector<float> v_recHit_matchedSimClusIndex;
    // std::vector<float> v_recHit_isCaloParticleMatched;

    std::vector<float> v_recHit_iType;
    std::vector<float> v_recHit_iCell1;
    std::vector<float> v_recHit_iCell2;

    std::vector<float> v_recHit_SiThickness;

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
      // tree->Branch("bunchCrossingNumber", &bunchCrossingNumber);

      // Gen Particle //
      sprintf(name, "genPart_n");
      tree->Branch(name, &genPart_n);

      sprintf(name, "genPart_E");
      tree->Branch(name, &v_genPart_E);

      sprintf(name, "genPart_px");
      tree->Branch(name, &v_genPart_px);

      sprintf(name, "genPart_py");
      tree->Branch(name, &v_genPart_py);

      sprintf(name, "genPart_pz");
      tree->Branch(name, &v_genPart_pz);

      sprintf(name, "genPart_pT");
      tree->Branch(name, &v_genPart_pT);

      sprintf(name, "genPart_eta");
      tree->Branch(name, &v_genPart_eta);

      sprintf(name, "genPart_phi");
      tree->Branch(name, &v_genPart_phi);

      // // Pileup //
      // sprintf(name, "pileup_n");
      // tree->Branch(name, &pileup_n);

      // // Rho //
      // sprintf(name, "rho");
      // tree->Branch(name, &rho);

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

      // sprintf(name, "simHit_ET");
      // tree->Branch(name, &v_simHit_ET);

      sprintf(name, "simHit_layer");
      tree->Branch(name, &v_simHit_layer);

      sprintf(name, "simHit_zside");
      tree->Branch(name, &v_simHit_zside);

      // sprintf(name, "simHit_isCaloParticleMatched");
      // tree->Branch(name, &v_simHit_isCaloParticleMatched);

      // sprintf(name, "simHit_matchedSimClusIndex");
      // tree->Branch(name, &v_simHit_matchedSimClusIndex);

      sprintf(name, "simHit_detector");
      tree->Branch(name, &v_simHit_detector);

      sprintf(name, "simHit_detId");
      tree->Branch(name, &v_simHit_detId);

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

      sprintf(name, "recHit_detId");
      tree->Branch(name, &v_recHit_detId);
    }

    void fill() { tree->Fill(); }

    void clear() {
      // Gen Particle //
      genPart_n = 0;
      v_genPart_E.clear();
      v_genPart_px.clear();
      v_genPart_py.clear();
      v_genPart_pz.clear();
      v_genPart_pT.clear();
      v_genPart_eta.clear();
      v_genPart_phi.clear();
      v_genPart_pdgId.clear();

      // // Pileup //
      // pileup_n = 0;

      // // Rho //
      // rho = 0;

      // Sim-hit //
      simHit_n = 0;
      v_simHit_E.clear();
      v_simHit_x.clear();
      v_simHit_y.clear();
      v_simHit_z.clear();
      v_simHit_eta.clear();
      v_simHit_phi.clear();
      // v_simHit_ET.clear();
      v_simHit_layer.clear();
      v_simHit_zside.clear();
      // v_simHit_isCaloParticleMatched.clear();
      // v_simHit_matchedSimClusIndex.clear();
      v_simHit_detector.clear();
      v_simHit_detId.clear();

      // Rec-hit //
      recHit_n = 0;
      v_recHit_E.clear();
      v_recHit_x.clear();
      v_recHit_y.clear();
      v_recHit_z.clear();
      v_recHit_eta.clear();
      v_recHit_phi.clear();
      // v_recHit_ET.clear();
      v_recHit_layer.clear();
      v_recHit_zside.clear();
      // v_recHit_matchedSimHitIndex.clear();
      // v_recHit_matchedSimClusIndex.clear();
      // v_recHit_isCaloParticleMatched.clear();
      v_recHit_iType.clear();
      v_recHit_iCell1.clear();
      v_recHit_iCell2.clear();
      v_recHit_detector.clear();
      v_recHit_detId.clear();

      v_recHit_SiThickness.clear();
    }
  };
}  // namespace TreeOutputInfo

#endif
