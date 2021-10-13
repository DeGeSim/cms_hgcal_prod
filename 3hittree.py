import os

import FWCore.ParameterSet.Config as cms
from Configuration.AlCa.GlobalTag import GlobalTag
from Configuration.StandardSequences.Eras import eras

# from EDAnalyzers.TreeMaker.parseOptions_cff import options
from EDFilters.MyFilters.ApplyFilters import apply_filters
from settingsparse import cliargs, settingsD
from SLHCUpgradeSimulations.Configuration.aging import customise_aging_1000

process = cms.Process("BuildTree", eras.Phase2C9)

# import of standard configurations
process.load("Configuration.StandardSequences.Services_cff")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("Configuration.EventContent.EventContent_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.RawToDigi_cff")
process.load("Configuration.StandardSequences.L1Reco_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.RecoSim_cff")
process.load("Configuration.StandardSequences.EndOfProcess_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")


process.GlobalTag = GlobalTag(process.GlobalTag, "auto:phase2_realistic_T15", "")

process.load("Configuration.Geometry.GeometryExtended2026D49Reco_cff")
process.load("Configuration.Geometry.GeometryExtended2026D49_cff")


############################## Parse arguments ##############################
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = "INFO"
process.MessageLogger.cerr.FwkReport.reportEvery = 100
# process.MessageLogger.categories.append("BuildTree")
process.MessageLogger.cerr.INFO = cms.untracked.PSet(limit=cms.untracked.int32(-1))

process.options = cms.untracked.PSet(wantSummary=cms.untracked.bool(True))
process.maxEvents = cms.untracked.PSet(
    input=cms.untracked.int32(settingsD["maxEvents"].value())
)

############################## File Paths ###################################
#### Check if inputfile is given, otherwise read filenames from the sourcefile

process.source = cms.Source(
    "PoolSource",
    fileNames=cms.untracked.vstring(
        "file:{}/{}.root".format(
            settingsD["path"]["step2_output"].value(),
            cliargs.fileid,
        )
    ),
    # Run1:Event1 to Run2:Event2
    # eventsToProcess = cms.untracked.VEventRange("1:78722-1:78722"),
    # duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
)


# if len(settingsD["eventRange"]):
#     process.source.eventsToProcess = cms.untracked.VEventRange(settingsD["eventRange"])

# if settingsD["depGraph"]:
#     process.DependencyGraph = cms.Service("DependencyGraph")
#     process.source = cms.Source("EmptySource")
#     process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(0))


ed_analyzer_kwargs = {
    ############################## My stuff ##############################
    "debug": settingsD["debug"],
    "isGunSample": settingsD["isGunSample"],
    "storeSimHit": settingsD["storeSimHit"],
    "storeRecHit": settingsD["storeRecHit"],
    ############################## GEN ##############################
    "label_generator": cms.InputTag("generator"),
    "label_genParticle": cms.InputTag("genParticles"),
    ############################## RECO ##############################
    "label_HGCEESimHit": cms.InputTag("g4SimHits", "HGCHitsEE"),
    "label_HGCHEFSimHit": cms.InputTag("g4SimHits", "HGCHitsHEfront"),
    "label_HGCHEBSimHit": cms.InputTag("g4SimHits", "HGCHitsHEback"),
    # "label_HGCEERecHit": cms.InputTag("HGCalRecHit", "HGCEERecHits"),
    # "label_HGCHEFRecHit": cms.InputTag("HGCalRecHit", "HGCHEFRecHits"),
    # "label_HGCHEBRecHit": cms.InputTag("HGCalRecHit", "HGCHEBRecHits"),
}
process.treeMaker = cms.EDAnalyzer(
    "TreeMaker",
    **ed_analyzer_kwargs,
)

## Apply filters to the generated particeles
process = apply_filters(process, settingsD)


outFile = "{}/{}.root".format(
    settingsD["path"]["hittrees"].value(),
    cliargs.fileid,
)

# Remove old output file
if os.path.isfile(outFile):
    os.remove(outFile)

# Output
process.TFileService = cms.Service(
    "TFileService", fileName=cms.string("file:" + outFile)
)

process.schedule = cms.Schedule()


# Aging

customise_aging_1000(process)

process.reco_seq = cms.Sequence()

if settingsD["onRaw"]:
    process.reco_seq = cms.Sequence(
        process.RawToDigi * process.L1Reco * process.reconstruction_mod
    )


###### PixelCPE issue
process.TrackProducer.TTRHBuilder = "WithTrackAngle"
process.PixelCPEGenericESProducer.UseErrorsFromTemplates = False
process.PixelCPEGenericESProducer.LoadTemplatesFromDB = False
process.PixelCPEGenericESProducer.TruncatePixelCharge = False
process.PixelCPEGenericESProducer.IrradiationBiasCorrection = False
process.PixelCPEGenericESProducer.DoCosmics = False
process.PixelCPEGenericESProducer.Upgrade = cms.bool(True)
######


process.p = cms.Path(
    process.filter_seq_genEle
    * process.filter_seq_genPho
    * process.filter_seq_genPart
    * process.reco_seq
    * process.treeMaker
)

process.schedule.insert(0, process.p)

print("*" * 50)
print("process.schedule:", process.schedule)
print("*" * 50)

# Tracer
if settingsD["trace"]:
    process.Tracer = cms.Service("Tracer")


if settingsD["memoryCheck"]:
    process.SimpleMemoryCheck = cms.Service(
        "SimpleMemoryCheck",
        moduleMemorySummary=cms.untracked.bool(True),
    )


# Timing
if settingsD["printTime"]:
    process.Timing = cms.Service(
        "Timing",
        summaryOnly=cms.untracked.bool(False),
        useJobReport=cms.untracked.bool(True),
    )


# Debug
if settingsD["debugFile"]:

    process.out = cms.OutputModule(
        "PoolOutputModule", fileName=cms.untracked.string("debug.root")
    )

    process.output_step = cms.EndPath(process.out)
    process.schedule.extend([process.output_step])


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
