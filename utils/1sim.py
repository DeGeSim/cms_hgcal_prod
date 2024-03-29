import os
from pprint import pprint
import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Phase2C9_cff import Phase2C9
from utils.settingsparse import settingsD, cliargs

from datetime import datetime
import time

process = cms.Process("SIM", Phase2C9)


# import of standard configurations
process.load("Configuration.StandardSequences.Services_cff")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.EventContent.EventContent_cff")
process.load("SimGeneral.MixingModule.mixNoPU_cfi")
#process.load("Geometry.HGCalCommonData.testGeometryV14_cff")
process.load('Configuration.Geometry.GeometryExtended2026D49Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2026D49_cff')
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Generator_cff")
process.load("IOMC.EventVertexGenerators.VtxSmearedHLLHC14TeV_cfi")
process.load("GeneratorInterface.Core.genFilterSummary_cff")
process.load("Configuration.StandardSequences.SimIdeal_cff")
process.load("Configuration.StandardSequences.EndOfProcess_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

process.maxEvents = cms.untracked.PSet(
    input=cms.untracked.int32(settingsD["maxEvents"].value()),
)

# Input source
process.source = cms.Source("EmptySource")

process.options = cms.untracked.PSet(
    FailPath=cms.untracked.vstring(),
    IgnoreCompletely=cms.untracked.vstring(),
    Rethrow=cms.untracked.vstring(),
    SkipEvent=cms.untracked.vstring(),
    allowUnscheduled=cms.obsolete.untracked.bool,
    canDeleteEarly=cms.untracked.vstring(),
    #    deleteNonConsumedUnscheduledModules = cms.untracked.bool(True),
    emptyRunLumiMode=cms.obsolete.untracked.string,
    eventSetup=cms.untracked.PSet(
        forceNumberOfConcurrentIOVs=cms.untracked.PSet(
            allowAnyLabel_=cms.required.untracked.uint32
        ),
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
    annotation=cms.untracked.string(str(settingsD)),
    name=cms.untracked.string("Applications"),
    version=cms.untracked.string("$Revision: 1.19 $"),
)

# Output definition
os.system("mkdir -p %s" % (settingsD["path"]["step1_output"].value()))

process.FEVTDEBUGoutput = cms.OutputModule(
    "PoolOutputModule",
    SelectEvents=cms.untracked.PSet(SelectEvents=cms.vstring("generation_step")),
    dataset=cms.untracked.PSet(
        dataTier=cms.untracked.string("GEN-SIM"), filterName=cms.untracked.string("")
    ),
    fileName=cms.untracked.string(
        "file:{}/{}.root".format(
            settingsD["path"]["step1_output"].value(),
            cliargs.fileid,
        )
    ),
    outputCommands=process.FEVTDEBUGEventContent.outputCommands,
    splitLevel=cms.untracked.int32(0),
)

# Additional output definition

# Other statements
process.genstepfilter.triggerConditions = cms.vstring("generation_step")
from Configuration.AlCa.GlobalTag import GlobalTag

process.GlobalTag = GlobalTag(process.GlobalTag, "auto:phase2_realistic_T15", "")

#dt = datetime.now()
#seed = int(time.mktime(dt.utctimetuple()) * 1000 + dt.microsecond / 1000)

#process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
#    # Include a PSet for each module label that needs a
#    # random engine.  The name is the module label.
#    # You must supply a seed or seeds.
#    # Optionally an engine type can be specified
#    generator = cms.PSet(
#        initialSeed = cms.untracked.uint32(cliargs.fileid)
#    ),
#)

process.RandomNumberGeneratorService.generator.initialSeed = cms.untracked.uint32(cliargs.fileid)


#process.generator = cms.EDProducer(
#    "FlatRandomEGunProducer",
#    PGunParameters=cms.PSet(
#        PartID=cms.vint32(settingsD["Gun"]["PartID"].value()),
#        MinEta=settingsD["Gun"]["MinEta"],
#        MaxEta=settingsD["Gun"]["MaxEta"],
#        MinE=settingsD["Gun"]["MinE"],
#        MaxE=settingsD["Gun"]["MaxE"],
#        MinPhi=settingsD["Gun"]["MinPhi"],
#        MaxPhi=settingsD["Gun"]["MaxPhi"],
#    ),
#    Verbosity=cms.untracked.int32(0),  ## set to 1 (or greater)  for printouts
#    psethack=cms.string("single gamma E 50"),  # Todo
#    AddAntiParticle=settingsD["Gun"]["AddAntiParticle"],
#    firstRun=cms.untracked.uint32(1),
#)
process.generator = cms.EDProducer("CloseByParticleGunProducer",
    AddAntiParticle = settingsD["Gun"]["AddAntiParticle"],
    PGunParameters = cms.PSet(
        Delta = cms.double(100),
        EnMax = settingsD["Gun"]["MaxE"],
        EnMin = settingsD["Gun"]["MinE"],
        MaxEta = settingsD["Gun"]["MaxEta"],
        MinEta = settingsD["Gun"]["MinEta"],
        MaxEnSpread = cms.bool(False),
        MaxPhi = settingsD["Gun"]["MaxPhi"],
        MinPhi = settingsD["Gun"]["MinPhi"],
        NParticles = cms.int32(1),
        Overlapping = cms.bool(False),
        PartID = cms.vint32(settingsD["Gun"]["PartID"].value()),
        Pointing = cms.bool(True),
        RandomShoot = cms.bool(False),
        RMax = cms.double(127),
        RMin = cms.double(41),
        ControlledByEta = cms.bool(True),
        ZMax = cms.double(321),
        ZMin = cms.double(320)
    ),
    Verbosity = cms.untracked.int32(0),
    firstRun = cms.untracked.uint32(1),
    psethack = cms.string('random particles in phi and r windows')
)



from SimCalorimetry.HGCalSimProducers.hgcalDigitizer_cfi import *
process = HGCal_disableNoise(process)

#Following Removes Mag Field
process.g4SimHits.UseMagneticField = False
process.g4SimHits.Physics.bField = cms.double(0.0)

process.ProductionFilterSequence = cms.Sequence(process.generator)


# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen)
process.simulation_step = cms.Path(process.psim)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.FEVTDEBUGoutput_step = cms.EndPath(process.FEVTDEBUGoutput)

# Schedule definition
process.schedule = cms.Schedule(
    process.generation_step,
    process.genfiltersummary_step,
    process.simulation_step,
    process.endjob_step,
    process.FEVTDEBUGoutput_step,
)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask

associatePatAlgosToolsTask(process)
# filter all path with the production filter sequence
for path in process.paths:
    getattr(process, path).insert(0, process.ProductionFilterSequence)


# Customisation from command line

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete

process = customiseEarlyDelete(process)
# End adding early deletion
