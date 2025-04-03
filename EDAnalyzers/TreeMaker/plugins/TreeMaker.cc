// -*- C++ -*-
//
// Package:    EDAnalyzers/TreeMaker
// Class:      TreeMaker
//
/**\class TreeMaker TreeMaker.cc EDAnalyzers/TreeMaker/plugins/TreeMaker.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:
//         Created:  Sat, 11 May 2019 13:14:55 GMT
//
//

// system include files
#include <memory>

// user include files

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/CaloTowers/interface/CaloTowerDefs.h"
#include "DataFormats/Common/interface/MapOfVectors.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/ForwardDetId/interface/HGCEEDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCalDetId.h"
#include "DataFormats/FWLite/interface/ESHandle.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/HGCalReco/interface/Trackster.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHit.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/ParticleFlowReco/interface/HGCalMultiCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHit.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "Geometry/CaloTopology/interface/HGCalTopology.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"
#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"
#include "SimDataFormats/CaloAnalysis/interface/SimCluster.h"
#include "SimDataFormats/CaloHit/interface/PCaloHit.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "EDAnalyzers/TreeMaker/interface/Common.h"
#include "EDAnalyzers/TreeMaker/interface/Constants.h"
#include "EDAnalyzers/TreeMaker/interface/TreeOutputInfo.h"

#include <CLHEP/Matrix/Matrix.h>
#include <CLHEP/Vector/ThreeVector.h>
#include <CLHEP/Vector/ThreeVector.h>

#include <Compression.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TMatrixD.h>
#include <TTree.h>
#include <TVector2.h>
#include <TVectorD.h>
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

double HGCal_minEta = 1.479;
double HGCal_maxEta = 3.1;

double el_minPt = 0;     //15;
double el_maxPt = 99999;  //30;

double ph_minPt = 0;     //15;
double ph_maxPt = 99999;  //30;

double _largeVal = 999999999;

class TreeMaker : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit TreeMaker(const edm::ParameterSet &);
  ~TreeMaker();

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event &, const edm::EventSetup &) override;
  virtual void endJob() override;

  hgcal::RecHitTools recHitTools;

  int minLayer;
  int maxLayer;

  TreeOutputInfo::TreeOutput *treeOutput;
  // My stuff //
  void WriteParticleKinematicsToTree(reco::GenParticle &part,
                                     std::vector<CLHEP::HepLorentzVector> &);
  void populateRecHitMaps(edm::Handle<edm::SortedCollection<HGCRecHit, edm::StrictWeakOrdering<HGCRecHit>>> &v_hits,
                          std::vector<DetId> &v_caloHitId);
  void populateSimHitMaps(edm::Handle<std::vector<PCaloHit> > &v_hits,
                          std::vector<DetId> &v_caloHitId);
  bool isGunSample;

  bool storeSimHit;
  bool storeRecHit;

  // GenEventInfoProduct //
  edm::EDGetTokenT<GenEventInfoProduct> tok_generator;

  // Gen particles //
  edm::EDGetTokenT<std::vector<reco::GenParticle>> tok_genParticle;

  // RecHits //
  //Declare a mapping from the id to the hit
  std::map<DetId, const HGCRecHit *> m_recHit;
  edm::EDGetTokenT<edm::SortedCollection<HGCRecHit, edm::StrictWeakOrdering<HGCRecHit>>> tok_HGCEERecHit;
  edm::EDGetTokenT<edm::SortedCollection<HGCRecHit, edm::StrictWeakOrdering<HGCRecHit>>> tok_HGCHEFRecHit;
  edm::EDGetTokenT<edm::SortedCollection<HGCRecHit, edm::StrictWeakOrdering<HGCRecHit>>> tok_HGCHEBRecHit;

  // SimHits //
  //Declare a mapping from the id to the hit
  std::map<DetId, const PCaloHit *> m_simHit;
  edm::EDGetTokenT<std::vector<PCaloHit> > tok_HGCEESimHit;
  edm::EDGetTokenT<std::vector<PCaloHit> > tok_HGCHEFSimHit;
  edm::EDGetTokenT<std::vector<PCaloHit> > tok_HGCHEBSimHit;

  // Calo particles //
  //
  // edm::EDGetTokenT<std::vector<CaloParticle>> tok_caloParticle;
  edm::ESGetToken<CaloGeometry, CaloGeometryRecord> caloGeomToken_;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
TreeMaker::TreeMaker(const edm::ParameterSet &iConfig) :
  caloGeomToken_(esConsumes<CaloGeometry, CaloGeometryRecord>())
{
  usesResource("TFileService");
  edm::Service<TFileService> fs;

  // Compression
  //fs->file().SetCompressionAlgorithm(ROOT::kLZMA);
  //fs->file().SetCompressionLevel(8);

  treeOutput = new TreeOutputInfo::TreeOutput("tree", fs);

  minLayer = +9999;
  maxLayer = -9999;

  // My stuff //
  // bool debug = iConfig.getParameter<bool>("debug");
  isGunSample = iConfig.getParameter<bool>("isGunSample");

  // storeSimHit = iConfig.getParameter<bool>("storeSimHit");
  // storeRecHit = iConfig.getParameter<bool>("storeRecHit");

  // GenEventInfoProduct //
  tok_generator = consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("label_generator"));

  // Gen particles //
  tok_genParticle = consumes<std::vector<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("label_genParticle"));


  tok_HGCEERecHit = consumes<edm::SortedCollection<HGCRecHit, edm::StrictWeakOrdering<HGCRecHit>>>(iConfig.getParameter<edm::InputTag>("label_HGCEERecHit"));
  tok_HGCHEFRecHit = consumes<edm::SortedCollection<HGCRecHit, edm::StrictWeakOrdering<HGCRecHit>>>(iConfig.getParameter<edm::InputTag>("label_HGCHEFRecHit"));
  tok_HGCHEBRecHit = consumes<edm::SortedCollection<HGCRecHit, edm::StrictWeakOrdering<HGCRecHit>>>(iConfig.getParameter<edm::InputTag>("label_HGCHEBRecHit"));
  tok_HGCEESimHit = consumes<std::vector<PCaloHit> >(
      iConfig.getParameter<edm::InputTag>("label_HGCEESimHit"));
  tok_HGCHEFSimHit = consumes<std::vector<PCaloHit> >(
      iConfig.getParameter<edm::InputTag>("label_HGCHEFSimHit"));
  tok_HGCHEBSimHit = consumes<std::vector<PCaloHit> >(
      iConfig.getParameter<edm::InputTag>("label_HGCHEBSimHit"));


}

TreeMaker::~TreeMaker() {
  // do anything here that needs to be done at destruction time
  // (e.g. close files, deallocate resources etc.)
  delete treeOutput;
}

//
// member functions
//
void TreeMaker::WriteParticleKinematicsToTree(reco::GenParticle &part,
                                              std::vector<CLHEP::HepLorentzVector> &v_genPart_4mom) {
  int pdgId = part.pdgId();
  int status = part.status();

  // Check for eleectrons and (hard/promt) photons, otherwise skip
  bool validHardEl = abs(pdgId) == 11 && ((isGunSample && status == 1) || (!isGunSample && part.isHardProcess()));
  bool validHardPh = abs(pdgId) == 22 && ((isGunSample && status == 1) || (!isGunSample && part.isHardProcess()));
  bool validPromtPh = abs(pdgId) == 22 && Common::isPromptPhoton(part);
  if (!(validHardEl || validHardPh || validPromtPh)) {
    //return;
  }

  double &maxPt = (abs(pdgId) == 11) ? el_maxPt : ph_maxPt;
  double &minPt = (abs(pdgId) == 11) ? el_minPt : ph_minPt;

  // Check if the particle is in the HGCal, otherwise skip
  bool ptetaCut =
      fabs(part.eta()) > HGCal_minEta && fabs(part.eta()) < HGCal_maxEta && part.pt() > minPt && part.pt() < maxPt;
  if (!ptetaCut) {
    return;
  }

  printf("Gen part found: pdgId %d E %0.2f, pT %0.2f, eta %+0.2f, pz %+0.2f \n",
         pdgId,
         part.energy(),
         part.pt(),
         part.eta(),
         part.pz());

  std::vector<CLHEP::HepLorentzVector> &v_gen_4mom_ref = v_genPart_4mom;

  CLHEP::HepLorentzVector gen_4mom;
  gen_4mom.setT(part.energy());
  gen_4mom.setX(part.px());
  gen_4mom.setY(part.py());
  gen_4mom.setZ(part.pz());

  v_gen_4mom_ref.push_back(gen_4mom);
  treeOutput->v_genPart_pdgId.push_back(pdgId);

  //Attach the properties of the patricle to the relevant vector in the tree.
  treeOutput->v_genPart_E.push_back(gen_4mom.e());
  treeOutput->v_genPart_px.push_back(gen_4mom.px());
  treeOutput->v_genPart_py.push_back(gen_4mom.py());
  treeOutput->v_genPart_pz.push_back(gen_4mom.pz());
  treeOutput->v_genPart_pT.push_back(gen_4mom.perp());
  treeOutput->v_genPart_eta.push_back(gen_4mom.eta());
  treeOutput->v_genPart_phi.push_back(gen_4mom.phi());
  treeOutput->genPart_n++;

}

void TreeMaker::populateRecHitMaps(
    edm::Handle<edm::SortedCollection<HGCRecHit, edm::StrictWeakOrdering<HGCRecHit>>> &v_hits,
    std::vector<DetId> &v_caloHitId) {
  int nHits = v_hits->size();
  std::cout << "Populating RecHitMap with " << nHits << "hits.\n";
  for (int iRecHit = 0; iRecHit < nHits; iRecHit++) {
    const HGCRecHit *recHit = &(*v_hits)[iRecHit];

    m_recHit[recHit->id()] = recHit;
    v_caloHitId.push_back(recHit->id());
  }
  return;
}
void TreeMaker::populateSimHitMaps(
    edm::Handle<std::vector<PCaloHit>> &v_hits,
    std::vector<DetId> &v_caloHitId) {
  int nHits = v_hits->size();
  std::cout << "Populating SimHitMap with " << nHits << "hits.\n";
  for (int iSimHit = 0; iSimHit < nHits; iSimHit++) {
    const PCaloHit *simHit = &(*v_hits)[iSimHit];

    // Add the hitid to to the detector mapping
    m_simHit[simHit->id()] = simHit;
    // add the hit id to the vector witht the hits of the subdetector
    v_caloHitId.push_back(simHit->id());
  }
  return;
}
// ------------ method called for each event  ------------
void TreeMaker::analyze(const edm::Event &iEvent, const edm::EventSetup &es)
{
  using namespace edm;

  long long eventNumber = iEvent.id().event();
  printf("Event %llu \n", eventNumber);

  treeOutput->clear();

  // Load the geometry
  const CaloGeometry &geom = es.getData(caloGeomToken_);
  recHitTools.setGeometry(geom);

  //////////////////// Run info ////////////////////
  treeOutput->runNumber = iEvent.id().run();
  treeOutput->eventNumber = iEvent.id().event();
  treeOutput->luminosityNumber = iEvent.id().luminosityBlock();
  // treeOutput->bunchCrossingNumber = iEvent.bunchCrossing();

  // HGCal Topology
  /*
  edm::ESHandle<HGCalTopology> handle_topo_HGCalEE;
  //iSetup.get<IdealGeometryRecord>().get("HGCalEESensitive", handle_topo_HGCalEE);
  iSetup.get("HGCalEESensitive", handle_topo_HGCalEE);

  if (!handle_topo_HGCalEE.isValid()) {
    printf("Error: Invalid HGCalEE topology. \n");
    exit(EXIT_FAILURE);
  }
  */

  //////////////////// GenEventInfoProduct ////////////////////
  edm::Handle<GenEventInfoProduct> generatorHandle;
  iEvent.getByToken(tok_generator, generatorHandle);
  GenEventInfoProduct generator = *generatorHandle;

  printf("[%llu] Gen. evt. wt. %0.4g \n", eventNumber, generator.weight());

  //////////////////// Gen particle ////////////////////
  edm::Handle<std::vector<reco::GenParticle>> v_genParticle;
  iEvent.getByToken(tok_genParticle, v_genParticle);

  std::vector<CLHEP::HepLorentzVector> v_genPart_4mom;

  // Iterate over the generated Particles, filter for valid ones in the HGCal and write the kinematics to the tree
  for (int iPart = 0; iPart < (int)v_genParticle->size(); iPart++) {
    reco::GenParticle part = v_genParticle->at(iPart);
    WriteParticleKinematicsToTree(part, v_genPart_4mom);
  }

  // // Pileup
  // edm::Handle<std::vector<PileupSummaryInfo>> pileUps_reco;
  // iEvent.getByToken(tok_pileup, pileUps_reco);
  // treeOutput->pileup_n = Common::getPileup(pileUps_reco);

  // // Rho
  // edm::Handle<double> handle_rho;
  // iEvent.getByToken(tok_rho, handle_rho);
  // double rho = *handle_rho;

  // treeOutput->rho = rho;
  
  //Declare the vectors for EE, HEF, HEB

  typedef edm::Handle<std::vector<PCaloHit> > SimHitVector;
  SimHitVector v_HGCEESimHit;
  SimHitVector v_HGCHEFSimHit;
  SimHitVector v_HGCHEBSimHit;

  typedef edm::Handle<edm::SortedCollection<HGCRecHit, edm::StrictWeakOrdering<HGCRecHit>>> RecHitVector;
  RecHitVector v_HGCEERecHit;
  RecHitVector v_HGCHEFRecHit;
  RecHitVector v_HGCHEBRecHit;
  //Load them from the tokens
  iEvent.getByToken(tok_HGCEESimHit, v_HGCEESimHit);
  iEvent.getByToken(tok_HGCHEFSimHit, v_HGCHEFSimHit);
  iEvent.getByToken(tok_HGCHEBSimHit, v_HGCHEBSimHit);

  iEvent.getByToken(tok_HGCEERecHit, v_HGCEERecHit);
  iEvent.getByToken(tok_HGCHEFRecHit, v_HGCHEFRecHit);
  iEvent.getByToken(tok_HGCHEBRecHit, v_HGCHEBRecHit);
  //'Calnumber' -> List of hits by DetId
  std::map<int, std::vector<DetId>> m_SimHitsByCalo;
  std::map<int, std::vector<DetId>> m_RecoHitsByCalo;
  //EE = HGCalEE = 8, Electromagnetic Endcap (EE), a tungsten/copper-silicon sampling electromagnetic calorimeter with a depth of about 26 X0 and 1.5λ
  //HEF = HGCalHSi = 9,Forward Hadronic (HEF or FH), a stainless-steel-silicon hadron calorimeter 3.5λ deep
  //HEB HGCalHSc = 10, Backing Hadronic (HEB or BH), a 5λ stainless-steel-scintillator sampling backing calorimeter
  //HGCalTrigger = 11

  //Create instances of the vector
  m_SimHitsByCalo[8];
  m_SimHitsByCalo[9];
  m_SimHitsByCalo[10];

  m_RecoHitsByCalo[8];
  m_RecoHitsByCalo[9];
  m_RecoHitsByCalo[10];

  populateSimHitMaps(v_HGCEESimHit, m_SimHitsByCalo.at(8));
  populateSimHitMaps(v_HGCHEFSimHit, m_SimHitsByCalo.at(9));
  populateSimHitMaps(v_HGCHEBSimHit, m_SimHitsByCalo.at(10));

  populateRecHitMaps(v_HGCEERecHit, m_RecoHitsByCalo.at(8));
  populateRecHitMaps(v_HGCHEFRecHit, m_RecoHitsByCalo.at(9));
  populateRecHitMaps(v_HGCHEBRecHit, m_RecoHitsByCalo.at(10));

  std::vector<SimHitVector*> SimHitsLists {&v_HGCEESimHit, &v_HGCHEFSimHit,&v_HGCHEBSimHit};
  std::vector<RecHitVector*> RecHitsLists {&v_HGCEERecHit, &v_HGCHEFRecHit,&v_HGCHEBRecHit};

  for (SimHitVector* simhitslistpointer: SimHitsLists){
    for (PCaloHit simHit: **simhitslistpointer){
      int layer = recHitTools.getLayer(simHit.id()) - 1;  // Start from 0
      treeOutput->v_simHit_layer.push_back(layer);

      int zside = recHitTools.zside(simHit.id());
      treeOutput->v_simHit_zside.push_back(zside);

      uint32_t rawDetId = simHit.id();
      DetId detId(rawDetId);
      int detector = detId.det();
      treeOutput->v_simHit_detector.push_back(detector);
      treeOutput->v_simHit_detId.push_back(rawDetId);

      auto position = recHitTools.getPosition(simHit.id());

      treeOutput->v_simHit_E.push_back(simHit.energy());

      treeOutput->v_simHit_x.push_back(position.x());
      treeOutput->v_simHit_y.push_back(position.y());
      treeOutput->v_simHit_z.push_back(position.z());

      treeOutput->v_simHit_eta.push_back(position.eta());
      treeOutput->v_simHit_phi.push_back(position.phi());

      treeOutput->simHit_n++;
    }
  }
  // RecHits
  for (RecHitVector* rechitslistpointer: RecHitsLists){
    for (HGCRecHit recHit: (**rechitslistpointer)){
      int layer = recHitTools.getLayer(recHit.id()) - 1;  // Start from 0
      treeOutput->v_recHit_layer.push_back(layer);

      int zside = recHitTools.zside(recHit.id());
      treeOutput->v_recHit_zside.push_back(zside);

      uint32_t detector = recHit.id().det();
      int rawDetId = recHit.id().rawId();
      treeOutput->v_recHit_detector.push_back(detector);
      treeOutput->v_recHit_detId.push_back(rawDetId);

      auto position = recHitTools.getPosition(recHit.id());

      treeOutput->v_recHit_E.push_back(recHit.energy());

      treeOutput->v_recHit_x.push_back(position.x());
      treeOutput->v_recHit_y.push_back(position.y());
      treeOutput->v_recHit_z.push_back(position.z());

      treeOutput->v_recHit_eta.push_back(position.eta());
      treeOutput->v_recHit_phi.push_back(position.phi());

      treeOutput->recHit_n++;
    }
  }

  // Fill tree
  treeOutput->fill();

  printf("\n");
  fflush(stdout);
  fflush(stderr);
}

// ------------ method called once each job just before starting event loop  ------------
void TreeMaker::beginJob() {}

// ------------ method called once each job just after ending the event loop  ------------
void TreeMaker::endJob() {
  printf("minLayer = %d, maxLayer = %d \n", minLayer, maxLayer);

  fflush(stdout);
  fflush(stderr);
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void TreeMaker::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TreeMaker);
