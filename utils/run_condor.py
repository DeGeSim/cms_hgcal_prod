import argparse
import numpy
import os


cwd = os.getcwd()


# Argument parser
parser = argparse.ArgumentParser(formatter_class = argparse.RawTextHelpFormatter)

# TreeMaker_SingleElectron_PT2to100_PhaseIITDRSpring19DR-PU200_106X_upgrade2023_realistic_v3-v1_GEN-SIM-DIGI-RAW
parser.add_argument(
    "--processName",
    help = "Name of the process to be run",
    type = str,
    required = True,
)

# EDAnalyzers/TreeMaker/python/ConfFile_cfg.py
#parser.add_argument(
#    "--cmsRunFile",
#    help = "cmsRun file",
#    type = str,
#    required = True,
#)

#parser.add_argument(
#    "--cmsRunOptions",
#    help = "Options for the cmsRun file",
#    type = str,
#    required = False,
#    default = "",
#)

# /eos/cms/store/group/phys_egamma/sobhatta/HGCal_TreeMaker
#parser.add_argument(
#    "--outputDir",
#    help = "cmsRun output directory",
#    type = str,
#    required = True,
#)
#
#parser.add_argument(
#    "--outputFile",
#    help = "Output file name (WITHOUT .root)",
#    type = str,
#    required = False,
#    default = "output"
#)

parser.add_argument(
    "--suffix",
    help = "suffix",
    type = str,
    required = False,
    default = "",
)

parser.add_argument(
    "--nJob",
    help = "Numbers of jobs",
    type = int,
    required = True,
)

#parser.add_argument(
#    "--nUnitPerJob",
#    help = "Numbers of units to process per job",
#    type = int,
#    required = True,
#)

parser.add_argument(
    "--test",
    help = "Only create job files (do not submit)",
    default = False,
    action = "store_true",
)


# Parse arguments
args = parser.parse_args()


condorConfig = "utils/condor_config.sub"
condorScript = "utils/condor_script.sh"

condorConfig_name = condorConfig[condorConfig.rfind("/")+1: condorConfig.rfind(".")]
condorConfig_ext = condorConfig[condorConfig.rfind("."):]

condorScript_name = condorScript[condorScript.rfind("/")+1: condorScript.rfind(".")]
condorScript_ext = condorScript[condorScript.rfind("."):]


if (__name__ == "__main__") :
    
    processName = "%s%s" %(args.processName, args.suffix)
    
    #outputDir = "%s/%s" %(args.outputDir, processName)
    #command = "mkdir -p " + outputDir
    #print "Command:", command
    #os.system(command)
    #print ""
    
    #condorDir = "%s/condorJobs" %(outputDir)
    condorDir = "condorJobs/%s" %(processName)
    command = "mkdir -p " + condorDir
    print "Command:", command
    os.system(command)
    print ""
    
    
    print "****************************************************************************************************"
    print "****************************************************************************************************"
    print "Process:", args.processName
    #print "cmsRun file:", args.cmsRunFile
    #print "Output directory:", outputDir
    print "Condor directory:", condorDir
    print "# jobs:", args.nJob
    #print "# units per job:", args.nUnitPerJob
    print "****************************************************************************************************"
    print "****************************************************************************************************"
    print ""
    
    
    condorConfig_content = ""
    
    with open(condorConfig, "r") as f :
        
        condorConfig_content = f.read()
    
    
    condorScript_content = ""
    
    with open(condorScript, "r") as f :
        
        condorScript_content = f.read()
    
    
    nDigit = len(str(args.nJob))
    
    for iJob in range(0, args.nJob) :
        
        jobNumberStr = str(iJob+1)
        
        
        condorConfig_mod = condorConfig_name + "_" + jobNumberStr + condorConfig_ext
        condorConfig_mod = "%s/%s" %(condorDir, condorConfig_mod)
        
        condorScript_mod = condorScript_name + "_" + jobNumberStr + condorScript_ext
        condorScript_mod = "%s/%s" %(condorDir, condorScript_mod)
        
        # Condor config
        condorConfig_content_mod = condorConfig_content
        condorConfig_content_mod = condorConfig_content_mod.replace("@exe@", condorScript_mod)
        condorConfig_content_mod = condorConfig_content_mod.replace("@log@", condorDir + "/" + "job_%s" %(jobNumberStr) + ".log")
        condorConfig_content_mod = condorConfig_content_mod.replace("@out@", condorDir + "/" + "job_%s" %(jobNumberStr) + ".out")
        condorConfig_content_mod = condorConfig_content_mod.replace("@err@", condorDir + "/" + "job_%s" %(jobNumberStr) + ".err")
        
        print "Writing: %s" %(condorConfig_mod)
        
        with open(condorConfig_mod, "w") as f :
            
            f.write(condorConfig_content_mod)
        
        
        #outputDir = "%s/%s" %(args.outputDir, processName)
        #outputFile = "%s_%d.root" %(args.outputFile, iJob+1)
        #
        # Condor script
        #cmsRun_cmd = ("cmsRun "
        #    "%s print "
        #    "outputDir=%s "
        #    "outputFile=%s "
        #    "seed=%d "
        #    "nEvent=%d "
        #    "%s"
        #%(
        #    args.cmsRunFile,
        #    outputDir,
        #    outputFile,
        #    iJob+1,
        #    args.nUnitPerJob,
        #    args.cmsRunOptions,
        #))
        
        cmd = "./utils/production_chain.sh %d" %(iJob+1)
        
        condorScript_content_mod = condorScript_content
        condorScript_content_mod = condorScript_content_mod.replace("@dir@", cwd)
        condorScript_content_mod = condorScript_content_mod.replace("@cmd@", cmd)
        
        print "Writing: %s" %(condorScript_mod)
        
        with open(condorScript_mod, "w") as f :
            
            f.write(condorScript_content_mod)
        
        command = "chmod +x %s" %(condorScript_mod)
        print "Command:", command
        os.system(command)
        
        
        # Submit job
        command = "condor_submit %s" %(condorConfig_mod)
        print "Command:", command
        
        commandReturn = 1
        
        if (not args.test) :
            
            # Repeat until job is submission is successful (returns 0)
            while (commandReturn) :
                
                commandReturn = os.system(command)
        
        
        print "\n"
        
    
    
    print "\n"
    print "****************************************************************************************************"
    print "****************************************************************************************************"
    print "Total # unit:", args.nJob
    print "Total # job:", args.nJob
    print "****************************************************************************************************"
    print "****************************************************************************************************"
    print "\n"
