# Auto generated configuration file
# using:
# Revision: 1.19
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v
# with command line options: step2 --conditions auto:phase2_realistic_T15 -s DIGI:pdigi_valid,L1TrackTrigger,L1,DIGI2RAW,HLT:@fake2 --datatier GEN-SIM-DIGI-RAW -n 10 --geometry Extended2026D49 --era Phase2C9 --eventcontent FEVTDEBUGHLT --no_exec --filein file:step1.root --fileout file:step2.root
import os
import FWCore.ParameterSet.Config as cms
from Configuration.Eras.Era_Phase2C9_cff import Phase2C9
from settingsparse import settingsD, cliargs


process = cms.Process("HLT", Phase2C9)

# import of standard configurations
process.load("Configuration.StandardSequences.Services_cff")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.EventContent.EventContent_cff")
process.load("SimGeneral.MixingModule.mixNoPU_cfi")
process.load("Geometry.HGCalCommonData.testGeometryV14_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Digi_cff")
process.load("Configuration.StandardSequences.L1TrackTrigger_cff")
process.load("Configuration.StandardSequences.SimL1Emulator_cff")
process.load("Configuration.StandardSequences.DigiToRaw_cff")
process.load("HLTrigger.Configuration.HLT_Fake2_cff")
process.load("Configuration.StandardSequences.EndOfProcess_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

process.maxEvents = cms.untracked.PSet(
    output=cms.optional.untracked.allowed(cms.int32, cms.PSet)
)

# Input source
process.source = cms.Source(
    "PoolSource",
    dropDescendantsOfDroppedBranches=cms.untracked.bool(False),
    fileNames=cms.untracked.vstring(
        "file:{}/{}.root".format(
            settingsD["path"]["step1_output"].value(),
            cliargs.fileid,
        )
    ),
    inputCommands=cms.untracked.vstring(
        "keep *",
        # "drop *_genParticles_*_*",
        "drop *_genParticlesForJets_*_*",
        "drop *_kt4GenJets_*_*",
        "drop *_kt6GenJets_*_*",
        "drop *_iterativeCone5GenJets_*_*",
        "drop *_ak4GenJets_*_*",
        "drop *_ak7GenJets_*_*",
        "drop *_ak8GenJets_*_*",
        "drop *_ak4GenJetsNoNu_*_*",
        "drop *_ak8GenJetsNoNu_*_*",
        "drop *_genCandidatesForMET_*_*",
        "drop *_genParticlesForMETAllVisible_*_*",
        "drop *_genMetCalo_*_*",
        "drop *_genMetCaloAndNonPrompt_*_*",
        "drop *_genMetTrue_*_*",
        "drop *_genMetIC5GenJs_*_*",
    ),
    secondaryFileNames=cms.untracked.vstring(),
)

process.options = cms.untracked.PSet(
    FailPath=cms.untracked.vstring(),
    IgnoreCompletely=cms.untracked.vstring(),
    Rethrow=cms.untracked.vstring(),
    SkipEvent=cms.untracked.vstring(),
    allowUnscheduled=cms.obsolete.untracked.bool,
    canDeleteEarly=cms.untracked.vstring(),
    emptyRunLumiMode=cms.obsolete.untracked.string,
    eventSetup=cms.untracked.PSet(
        forceNumberOfConcurrentIOVs=cms.untracked.PSet(),
        numberOfConcurrentIOVs=cms.untracked.uint32(1),
    ),
    fileMode=cms.untracked.string("FULLMERGE"),
    forceEventSetupCacheClearOnNewRun=cms.untracked.bool(False),
    makeTriggerResults=cms.obsolete.untracked.bool,
    numberOfConcurrentLuminosityBlocks=cms.untracked.uint32(1),
    numberOfConcurrentRuns=cms.untracked.uint32(1),
    numberOfStreams=cms.untracked.uint32(0),
    numberOfThreads=cms.untracked.uint32(1),
    printDependencies=cms.untracked.bool(False),
    sizeOfStackForThreadsInKB=cms.optional.untracked.uint32,
    throwIfIllegalParameter=cms.untracked.bool(True),
    wantSummary=cms.untracked.bool(False),
)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation=cms.untracked.string("step2 nevts:10"),
    name=cms.untracked.string("Applications"),
    version=cms.untracked.string("$Revision: 1.19 $"),
)

# Output definition
# Output definition
os.system("mkdir -p %s" %(settingsD["path"]["step2_output"].value()))

process.FEVTDEBUGHLToutput = cms.OutputModule(
    "PoolOutputModule",
    dataset=cms.untracked.PSet(
        dataTier=cms.untracked.string("GEN-SIM-DIGI-RAW"),
        filterName=cms.untracked.string(""),
    ),
    fileName=cms.untracked.string(
        "file:{}/{}.root".format(
            settingsD["path"]["step2_output"].value(),
            cliargs.fileid,
        )
    ),
    outputCommands=process.FEVTDEBUGHLTEventContent.outputCommands,
    splitLevel=cms.untracked.int32(0),
)

# Additional output definition

# Other statements
from SimCalorimetry.HGCalSimProducers.hgcalDigitizer_cfi import hgceeDigitizer, hgchebackDigitizer, hgchefrontDigitizer, HGCAL_chargeCollectionEfficiencies, HGCAL_ileakParam_toUse, HGCAL_cceParams_toUse
from SimGeneral.MixingModule.caloTruthProducer_cfi import *

theDigitizers = cms.PSet(
    hgceeDigitizer = cms.PSet(hgceeDigitizer),
    hgchebackDigitizer = cms.PSet(hgchebackDigitizer),
    hgchefrontDigitizer = cms.PSet(hgchefrontDigitizer),
    calotruth = cms.PSet(caloParticles),
)

#process.mix.digitizers = cms.PSet(process.theDigitizersValid)
process.mix.digitizers = cms.PSet(cms.PSet(theDigitizers))
from Configuration.AlCa.GlobalTag import GlobalTag

process.GlobalTag = GlobalTag(process.GlobalTag, "auto:phase2_realistic_T15", "")

from SimCalorimetry.Configuration.SimCalorimetry_cff import *
from GeneratorInterface.Core.generatorSmeared_cfi import *

doAllDigiTask = cms.Task(generatorSmeared)#, calDigiTask)
pdigi_valid = cms.Sequence(doAllDigiTask)

#DigiToRawTask = cms.Task(L1TDigiToRawTask, siPixelRawData, SiStripDigiToRaw, ecalPacker, esDigiToRaw, hcalRawDataTask, cscpacker, dtpacker, rpcpacker, ctppsRawData, castorRawData, rawDataCollector)

process.load("Configuration.StandardSequences.RawToDigi_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")

from EventFilter.HGCalRawToDigi.HGCalRawToDigi_cfi import *
hgcalRawToDigiTask = cms.Task(hgcalDigis)

from RecoLocalCalo.Configuration.hgcalLocalReco_cff import *
hgcalLocalRecoTask = cms.Task(
    HGCalUncalibRecHit,
    HGCalRecHit,
    hgcalRecHitMapProducer,
    #hgcalLayerClusters,
    #hgcalMultiClusters,
    particleFlowRecHitHGC,
    #particleFlowClusterHGCal
)

# Path and EndPath definitions
process.digitisation_step = cms.Path(pdigi_valid)
#process.L1TrackTrigger_step = cms.Path(process.L1TrackTrigger)
#process.pL1TkPrimaryVertex = cms.Path(process.L1TkPrimaryVertex)
#process.pL1TkPhotonsCrystal = cms.Path(process.L1TkPhotonsCrystal)
#process.pL1TkIsoElectronsCrystal = cms.Path(process.L1TkIsoElectronsCrystal)
#process.pL1TkElectronsLooseCrystal = cms.Path(process.L1TkElectronsLooseCrystal)
#process.pL1TkElectronsHGC = cms.Path(process.L1TkElectronsHGC)
#process.pL1TkMuon = cms.Path(process.L1TkMuons + process.L1TkMuonsTP)
#process.pL1TkElectronsLooseHGC = cms.Path(process.L1TkElectronsLooseHGC)
#process.pL1TkElectronsEllipticMatchHGC = cms.Path(process.L1TkElectronsEllipticMatchHGC)
#process.pL1TkElectronsCrystal = cms.Path(process.L1TkElectronsCrystal)
#process.pL1TkPhotonsHGC = cms.Path(process.L1TkPhotonsHGC)
#process.pL1TkIsoElectronsHGC = cms.Path(process.L1TkIsoElectronsHGC)
#process.pL1TkElectronsEllipticMatchCrystal = cms.Path(
#    process.L1TkElectronsEllipticMatchCrystal
#)
#process.L1simulation_step = cms.Path(process.SimL1Emulator)
process.digi2raw_step = cms.Path(process.DigiToRaw)
#process.raw2digi_step = cms.Path(process.RawToDigi)
process.raw2digi_step = cms.Path(cms.Sequence(hgcalRawToDigiTask))
#process.reconstruction_step = cms.Path(process.reconstruction)
process.reconstruction_step = cms.Path(cms.Sequence(hgcalLocalRecoTask))
process.endjob_step = cms.EndPath(process.endOfProcess)
process.FEVTDEBUGHLToutput_step = cms.EndPath(process.FEVTDEBUGHLToutput)

# Schedule definition
process.schedule = cms.Schedule(
    process.digitisation_step,
    #process.L1TrackTrigger_step,
    #process.L1simulation_step,
    #process.digi2raw_step,
    process.raw2digi_step,
    process.reconstruction_step,
)
#process.schedule.extend(process.HLTSchedule)
process.schedule.extend([process.endjob_step, process.FEVTDEBUGHLToutput_step])
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask

associatePatAlgosToolsTask(process)

# customisation of the process.

# Automatic addition of the customisation function from HLTrigger.Configuration.customizeHLTforMC
#from HLTrigger.Configuration.customizeHLTforMC import customizeHLTforMC

# call to customisation function customizeHLTforMC imported from HLTrigger.Configuration.customizeHLTforMC
#process = customizeHLTforMC(process)

from SimCalorimetry.HGCalSimProducers.hgcalDigitizer_cfi import *

process = HGCal_disableNoise(process)

# End of customisation functions

# Customisation from command line

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete

process = customiseEarlyDelete(process)
# End adding early deletion
