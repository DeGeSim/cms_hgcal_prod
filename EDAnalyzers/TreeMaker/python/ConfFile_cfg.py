import os
import FWCore.ParameterSet.Config as cms

processName = "Demo"

from Configuration.StandardSequences.Eras import eras

process = cms.Process(processName, eras.Phase2C9)

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

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, "auto:phase2_realistic_T15", "")

process.load("Configuration.Geometry.GeometryExtended2026D49Reco_cff")
process.load("Configuration.Geometry.GeometryExtended2026D49_cff")


############################## Parse arguments ##############################
from EDAnalyzers.TreeMaker.parseOptions_cff import options



process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.MessageLogger.categories.append('Demo')
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    limit = cms.untracked.int32(-1)
)

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(options.maxEvents))

############################## File Paths ###################################
#### Check if inputfile is given, otherwise read filenames from the sourcefile
fNames = []
if len(options.inputFiles):
    fNames = options.inputFiles
else:
    with open(options.sourceFile) as f:
        fNames = f.readlines()
#### Add file: prefix for local files as required by ROOT
for iFile, fName in enumerate(fNames):
    if "file:" not in fName and "root:" not in fName:
        fNames[iFile] = "file:%s" % (fName)

inputFileNames = cms.untracked.vstring(fNames)


def constructOutFilename():
    outFileSuffix = ""
    if options.onRaw:
        outFileSuffix = "%s_onRaw" % (outFileSuffix)

    if options.outFileNumber >= 0:
        outFileSuffix = "%s_%d" % (outFileSuffix, options.outFileNumber)

    return "ntupleTree%s.root" % (outFileSuffix)


outFile = constructOutFilename()

#### Create output folders, prefix path of the output file with the outdir
if len(options.outputDir):
    os.system("mkdir -p %s" % (options.outputDir))
    outFile = "%s/%s" % (options.outputDir, outFile)


process.source = cms.Source(
    "PoolSource",
    fileNames=inputFileNames,
    # Run1:Event1 to Run2:Event2
    # eventsToProcess = cms.untracked.VEventRange("1:78722-1:78722"),
    # duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
)


if len(options.eventRange):
    process.source.eventsToProcess = cms.untracked.VEventRange(options.eventRange)

if options.depGraph:
    process.DependencyGraph = cms.Service("DependencyGraph")
    process.source = cms.Source("EmptySource")
    process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(0))


ed_analyzer_kwargs = {
    ############################## My stuff ##############################
    "debug": cms.bool(False),
    "isGunSample": cms.bool(bool(options.isGunSample)),
    "storeSimHit": cms.bool(bool(options.storeSimHit)),
    "storeRecHit": cms.bool(bool(options.storeRecHit)),
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
process.treeMaker = cms.EDAnalyzer("TreeMaker", **ed_analyzer_kwargs)

## Apply filters to the generated particeles
from EDFilters.MyFilters.ApplyFilters import apply_filters

process = apply_filters(process, options)

# Remove old output file
if os.path.isfile(outFile):
    os.remove(outFile)

# Output file name modification
if outFile.startswith("/eos/cms"):
    outFile = outFile.replace("/eos/cms", "root://eoscms.cern.ch//eos/cms")


# Output
process.TFileService = cms.Service("TFileService", fileName=cms.string(outFile))

process.schedule = cms.Schedule()


# Aging
from SLHCUpgradeSimulations.Configuration.aging import customise_aging_1000

customise_aging_1000(process)


process.reco_seq = cms.Sequence()

if options.onRaw:
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

print "\n"
print "*" * 50
print "process.schedule:", process.schedule
print "*" * 50
print "\n"


# Tracer
if options.trace:
    process.Tracer = cms.Service("Tracer")


if options.memoryCheck:
    process.SimpleMemoryCheck = cms.Service(
        "SimpleMemoryCheck",
        moduleMemorySummary=cms.untracked.bool(True),
    )


# Timing
if options.printTime:
    process.Timing = cms.Service(
        "Timing",
        summaryOnly=cms.untracked.bool(False),
        useJobReport=cms.untracked.bool(True),
    )


# Debug
if options.debugFile:

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