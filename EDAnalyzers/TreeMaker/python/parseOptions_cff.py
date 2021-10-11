import FWCore.ParameterSet.VarParsing as VarParsing

options = VarParsing.VarParsing("analysis")


options.register(
    "sourceFile",
    "sourceFiles/SingleElectronFlatPtGun_fpantale_pT-0-200_eta-1p5-3p0_GEN-SIM-RECO/SingleElectronFlatPtGun_fpantale_pT-0-200_eta-1p5-3p0_GEN-SIM-RECO_mod.txt",  # Default value
    VarParsing.VarParsing.multiplicity.singleton,  # singleton or list
    VarParsing.VarParsing.varType.string,  # string, int, or float
    "File containing list of input files",  # Description
)

options.register(
    "outputDir",
    "",  # Default value
    VarParsing.VarParsing.multiplicity.singleton,  # singleton or list
    VarParsing.VarParsing.varType.string,  # string, int, or float
    "Output directory",  # Description
)

options.register(
    "outFileNumber",
    -1,  
    VarParsing.VarParsing.multiplicity.singleton,  # singleton or list
    VarParsing.VarParsing.varType.int,  # string, int, or float
    "File number (will be added to the filename if >= 0)",  # Description
)

options.register(
    "eventRange",
    "",  # Default value
    VarParsing.VarParsing.multiplicity.singleton,  # singleton or list
    VarParsing.VarParsing.varType.string,  # string, int, or float
    "Syntax: Run1:Event1-Run2:Event2 (includes both)",  # Description
)

# maxEvents = -1
# options.maxEvents = 15
# process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(options.maxEvents))
# options.register(
#     "maxEvents",
#     -1, # Default value
#     VarParsing.VarParsing.multiplicity.singleton,
#     VarParsing.VarParsing.varType.int,
#     "Maximum Number of included Events, -1 (default) removes the limit.",  # Description
# )

options.register(
    "debugFile",
    0,  # Default value
    VarParsing.VarParsing.multiplicity.singleton,  # singleton or list
    VarParsing.VarParsing.varType.int,  # string, int, or float
    "Create debug file",  # Description
)

options.register(
    "onRaw",
    0,  # Default value
    VarParsing.VarParsing.multiplicity.singleton,  # singleton or list
    VarParsing.VarParsing.varType.int,  # string, int, or float
    "Running on RAW",  # Description
)

options.register(
    "storeSimHit",
    0,  # Default value
    VarParsing.VarParsing.multiplicity.singleton,  # singleton or list
    VarParsing.VarParsing.varType.int,  # string, int, or float
    "Store sim-hits",  # Description
)

options.register(
    "storeRecHit",
    0,  # Default value
    VarParsing.VarParsing.multiplicity.singleton,  # singleton or list
    VarParsing.VarParsing.varType.int,  # string, int, or float
    "Store rec-hits",  # Description
)

options.register(
    "isGunSample",
    0,  # Default value
    VarParsing.VarParsing.multiplicity.singleton,  # singleton or list
    VarParsing.VarParsing.varType.int,  # string, int, or float
    "Is it a particle gun sample",  # Description
)

options.register(
    "genEleFilter",
    0,  # Default value
    VarParsing.VarParsing.multiplicity.singleton,  # singleton or list
    VarParsing.VarParsing.varType.int,  # string, int, or float
    "Apply gen-electron filter",  # Description
)

options.register(
    "genPhoFilter",
    0,  # Default value
    VarParsing.VarParsing.multiplicity.singleton,  # singleton or list
    VarParsing.VarParsing.varType.int,  # string, int, or float
    "Apply gen-photon filter",  # Description
)

options.register(
    "genPartonFilter",
    0,  # Default value
    VarParsing.VarParsing.multiplicity.singleton,  # singleton or list
    VarParsing.VarParsing.varType.int,  # string, int, or float
    "Apply gen-parton filter",  # Description
)

options.register(
    "trace",
    0,  # Default value
    VarParsing.VarParsing.multiplicity.singleton,  # singleton or list
    VarParsing.VarParsing.varType.int,  # string, int, or float
    "Trace modules",  # Description
)

options.register(
    "memoryCheck",
    0,  # Default value
    VarParsing.VarParsing.multiplicity.singleton,  # singleton or list
    VarParsing.VarParsing.varType.int,  # string, int, or float
    "Check memory usage",  # Description
)

options.register(
    "printTime",
    0,  # Default value
    VarParsing.VarParsing.multiplicity.singleton,  # singleton or list
    VarParsing.VarParsing.varType.int,  # string, int, or float
    "Print timing information",  # Description
)

options.register(
    "depGraph",
    0,  # Default value
    VarParsing.VarParsing.multiplicity.singleton,  # singleton or list
    VarParsing.VarParsing.varType.int,  # string, int, or float
    "Produce dependency graph only",  # Description
)

options.parseArguments()