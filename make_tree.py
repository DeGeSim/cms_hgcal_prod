import os

import FWCore.ParameterSet.Config as cms
from Configuration.AlCa.GlobalTag import GlobalTag
from Configuration.StandardSequences.Eras import eras

# from EDAnalyzers.TreeMaker.parseOptions_cff import options
from EDFilters.MyFilters.ApplyFilters import apply_filters

from Configuration.Eras.Era_Phase2C22I13M9_cff import Phase2C22I13M9

process = cms.Process('BuildTree',Phase2C22I13M9)

# import of standard configurations
process.load("Configuration.StandardSequences.Services_cff")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("Configuration.EventContent.EventContent_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.EndOfProcess_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")


process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T38', '')

process.load('Configuration.Geometry.GeometryExtended2026D113Reco_cff')


############################## Parse arguments ##############################
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = "INFO"
process.MessageLogger.cerr.FwkReport.reportEvery = 100
# process.MessageLogger.categories.append("BuildTree")
process.MessageLogger.cerr.INFO = cms.untracked.PSet(limit=cms.untracked.int32(-1))

process.options = cms.untracked.PSet(wantSummary=cms.untracked.bool(True))
process.maxEvents = cms.untracked.PSet(
    input=cms.untracked.int32(-1)
)

############################## File Paths ###################################
#### Check if inputfile is given, otherwise read filenames from the sourcefile

process.source = cms.Source(
    "PoolSource",
    fileNames=cms.untracked.vstring("/store/mc/Phase2Spring24DIGIRECOMiniAOD/SinglePhoton_E-1To1000_Eta-2_Phi-1p57_Z-321-CloseByParticleGun/GEN-SIM-DIGI-RAW-MINIAOD/noPU_AllTP_140X_mcRun4_realistic_v4-v1/2530000/b9a36d08-5642-4bdd-8dd0-161ac7e9fb06.root",
        "/store/mc/Phase2Spring24DIGIRECOMiniAOD/SinglePhoton_E-1To1000_Eta-2_Phi-1p57_Z-321-CloseByParticleGun/GEN-SIM-DIGI-RAW-MINIAOD/noPU_AllTP_140X_mcRun4_realistic_v4-v1/2530000/d96dcb5f-c916-48ce-abd8-0b31fb01b7ce.root"),
)



ed_analyzer_kwargs = {
    ############################## My stuff ##############################
    "debug": cms.bool(False),
    "isGunSample": cms.bool(True),
    "storeSimHit": cms.bool(True),
    "storeRecHit": cms.bool(True),
    ############################## GEN ##############################
    "label_generator": cms.InputTag("generator"),
    "label_genParticle": cms.InputTag("genParticles"),
    ############################## RECO ##############################
    "label_HGCEESimHit": cms.InputTag("g4SimHits", "HGCHitsEE"),
    "label_HGCHEFSimHit": cms.InputTag("g4SimHits", "HGCHitsHEfront"),
    "label_HGCHEBSimHit": cms.InputTag("g4SimHits", "HGCHitsHEback"),
    "label_HGCEERecHit": cms.InputTag("HGCalRecHit", "HGCEERecHits"),
    "label_HGCHEFRecHit": cms.InputTag("HGCalRecHit", "HGCHEFRecHits"),
    "label_HGCHEBRecHit": cms.InputTag("HGCalRecHit", "HGCHEBRecHits"),
}
process.treeMaker = cms.EDAnalyzer(
    "TreeMaker",
    **ed_analyzer_kwargs
)

settingsD = {
'isGunSample' : cms.bool(True),
'genEleFilter' : cms.bool(False),
'genPhoFilter' : cms.bool(True),
'genPartonFilter' : cms.bool(False),
}


## Apply filters to the generated particeles
process = apply_filters(process, settingsD)


outFile = "hgcal_tree.root"

# Remove old output file
#if os.path.isfile(outFile):
    #os.remove(outFile)

# Output
process.TFileService = cms.Service(
    "TFileService", fileName=cms.string("file:" + outFile)
)

process.schedule = cms.Schedule()


# Aging


process.reco_seq = cms.Sequence()



process.p = cms.Path(
    #process.filter_seq_genEle
    #* process.filter_seq_genPho
    #* process.filter_seq_genPart
    #* process.reco_seq
    process.treeMaker
)

process.schedule.insert(0, process.p)

process.MessageLogger = cms.Service(
    "MessageLogger",
    destinations=cms.untracked.vstring(
        "cerr",
    ),
    cerr=cms.untracked.PSet(
        # threshold  = cms.untracked.string("ERROR"),
        DEBUG=cms.untracked.PSet(limit=cms.untracked.int32(0)),
        WARNING=cms.untracked.PSet(limit=cms.untracked.int32(0)),
        ERROR=cms.untracked.PSet(limit=cms.untracked.int32(0)),
    ),
)
