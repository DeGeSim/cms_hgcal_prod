import argparse
import glob
import os


# Argument parser
parser = argparse.ArgumentParser()

# List of directories
parser.add_argument(
    "-d",
    "--directories",
    help = "List of directories",
    nargs = "*",
    type = str,
    required = True
)

parser.add_argument(
    "--submit",
    help = "Submit the failed jobs",
    action = "store_true"
)

parser.add_argument(
    "--success",
    help = "Will only show the jobs that have succeeded",
    action = "store_true"
)

# Parse arguments
args = parser.parse_args()


l_directory = args.directories

print "You have entered:"
print "\n".join(l_directory)
print ""

if (args.submit) :
    
    print "Will SUBMIT the failed jobs."

else :
    
    print "Will NOT submit the failed jobs."

print ""
print "Enter YES to proceed, or anything else to cancel."
choice = raw_input("Input: ")

if (choice != "YES") :
    
    print "Exiting."
    exit()


extension_log = ".log"
extension_err = ".err"
extension_out = ".out"

# Case insentitive
#l_errorStr = ["error", "segmentation", "killed", "SCRAM fatal"]
l_errorStr = ["Fatal Exception", "FatalRootError", "fatal system signal has occurred"]
l_errorStrSkip = ["unknown branch"]

# The format is: (return value <value>)
retValStr_start = "(return value "
retValStr_end = ")"

jobSubStr = "Job submitted"
jobExeStr = "Job executing"
jobRmvStr = "Job removed"
jobAbortStr = "Job was aborted"
jobTermStr = "Job terminated."
jobRequeueStr = "Job terminated and was requeued"

condorConfig_nameTemplate = "condor_config_%s.sub"
command_template = "condor_submit %s"


nDir = len(l_directory)

l_nJob_total = [0] * nDir
l_nJob_succ = [0] * nDir
l_nJob_fail = [0] * nDir
l_nJob_rmvd = [0] * nDir
l_nJob_abrt = [0] * nDir


# Check if all the directories exist
for iDir in range(0, nDir) :
    
    directory = l_directory[iDir]
    
    if (not os.path.isdir(directory)) :
        
        print "Directory does not exist: %s" %(directory)
        exit(1)


for iDir in range(0, nDir) :
    
    directory = l_directory[iDir]
    
    print "\n"
    print "*"*50
    print "Directory %d/%d: %s" %(iDir+1, nDir, directory)
    print "*"*50
    
    nJob_total = 0
    nJob_succ = 0
    nJob_fail = 0
    nJob_rmvd = 0
    nJob_abrt = 0
    
    l_fileName = glob.glob(directory + "/**")
    l_fileName = [f for f in l_fileName if extension_log in f]
    l_fileName.sort(key = lambda l: int(l[l.rfind("_")+1: l.rfind(extension_log)]))
    
    nJob_total += len(l_fileName)
    
    count = 0
    
    for iFile in range(0, len(l_fileName)) :
        
        fileName = l_fileName[iFile]
        
        fileContent = []
        
        with open(fileName) as f :
            
            fileContent = f.readlines()
        
        
        hasCompleted = False
        hasError = False
        retVal = 0
        
        isRemoved = False
        rmvLine = ""
        
        isAborted = False
        isRequeued = False
        
        
        # N.B. Log files are appended to (not overwritten) on resubmitting the job
        # Hence use the last occurence of the job summary
        for line in fileContent :
            
            if (jobTermStr in line) :
                
                # Final termination
                isRequeued = False
            
            if (jobRequeueStr in line) :
                
                isRequeued = True
            
            if (retValStr_start in line) :
                
                if (not isRequeued) :
                    
                    hasCompleted = True
                    retValStr = line
            
            # Check for removed jobs
            if (jobRmvStr in line) :
                
                hasCompleted = True
                isRemoved = True
                rmvLine = line.strip()
            
            # Check for aborted jobs
            if (jobAbortStr in line) :
                
                hasCompleted = True
                isRemoved = False
                isAborted = True
            
            # The job must not have been resubmitted
            if (jobSubStr in line or jobExeStr in line) :
                
                hasCompleted = False
                isRemoved = False
                isAborted = False
            
        
        
        if (not hasCompleted) :
            
            #print "#"*10, fileName
            #exit(0)
            
            continue
        
        
        if (not isRemoved and not isAborted) :
            
            retValStr = retValStr[retValStr.find(retValStr_start) + len(retValStr_start):]
            retValStr = retValStr[: retValStr.find(retValStr_end)]
            
            retVal = int(retValStr)
        
        hasError = bool(retVal) or isRemoved or isAborted
        
        
        nJob_succ += int(not hasError)
        nJob_fail += int(hasError and not isRemoved and not isAborted)
        nJob_rmvd += int(isRemoved)
        nJob_abrt += int(isAborted)
        
        hasPrinted = False
        
        if ((hasError and not args.success) or (not hasError and args.success)) :
            
            count += 1
            
            number = fileName[fileName.rfind("_")+1: fileName.rfind(extension_log)]
            
            condorConfig_name = condorConfig_nameTemplate %(number)
            condorConfig_name = "%s/%s" %(directory, condorConfig_name)
            command = command_template %(condorConfig_name)
            
            print "Count:", count
            print "Total count:", nJob_fail+nJob_rmvd+nJob_abrt
            print "File:", fileName
            
            if (isRemoved) :
                
                print rmvLine
            
            elif (isAborted) :
                
                print "Job aborted."
            
            else :
                
                print "Return value: %d" %(retVal)
            
            print "Condor config:", condorConfig_name
            
            hasPrinted = True
        
        if (hasError and args.submit) :
            
            # Delete the err and out files
            fileName_err = fileName.replace(extension_log, extension_err)
            print "Deleting error file: %s" %(fileName_err)
            os.system("rm %s" %(fileName_err))
            
            fileName_out = fileName.replace(extension_log, extension_out)
            print "Deleting output file: %s" %(fileName_out)
            os.system("rm %s" %(fileName_out))
            
            print "Command:", command
            os.system(command)
            
            hasPrinted = True
        
        if (hasPrinted) :
            
            print "\n"
    
    
    l_nJob_total[iDir] = nJob_total
    l_nJob_succ[iDir]  = nJob_succ
    l_nJob_fail[iDir]  = nJob_fail
    l_nJob_rmvd[iDir]  = nJob_rmvd
    l_nJob_abrt[iDir]  = nJob_abrt


print "\n"
print "_"*100
print ""

for iDir in range(0, nDir) :
    
    directory = l_directory[iDir]
    
    nDigit = len(str(l_nJob_total[iDir]))
    
    print "Directory %d/%d: %s" %(iDir+1, nDir, directory)
    print "Number of succeeded jobs: %0*d/%d" %(nDigit, l_nJob_succ[iDir], l_nJob_total[iDir])
    print "Number of failed jobs   : %0*d/%d" %(nDigit, l_nJob_fail[iDir], l_nJob_total[iDir])
    print "Number of removed jobs  : %0*d/%d" %(nDigit, l_nJob_rmvd[iDir], l_nJob_total[iDir])
    print "Number of aborted jobs  : %0*d/%d" %(nDigit, l_nJob_abrt[iDir], l_nJob_total[iDir])
    print ""


nJob_total = sum(l_nJob_total)
nJob_succ  = sum(l_nJob_succ)
nJob_fail  = sum(l_nJob_fail)
nJob_rmvd  = sum(l_nJob_rmvd)
nJob_abrt  = sum(l_nJob_abrt)

nDigit = len(str(nJob_total))

print "\n"
print "_"*100
print ""
print "Total number of succeeded jobs: %0*d/%d" %(nDigit, nJob_succ, nJob_total)
print "Total number of failed jobs   : %0*d/%d" %(nDigit, nJob_fail, nJob_total)
print "Total number of removed jobs  : %0*d/%d" %(nDigit, nJob_rmvd, nJob_total)
print "Total number of aborted jobs  : %0*d/%d" %(nDigit, nJob_abrt, nJob_total)
print "\n"
