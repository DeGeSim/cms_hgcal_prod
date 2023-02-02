#!/bin/sh

DIR=@dir@

echo "Working directory: "$DIR
cd $DIR
echo "Moved to working directory."

#export SCRAM_ARCH=slc7_amd64_gcc820
#export CPATH=$CPATH:$DIR
#export ROOT_INCLUDE_PATH=$ROOT_INCLUDE_PATH:$DIR

# Proxy path must be common to all nodes
# Default path is /tmp/, but thatis not common for all nodes
export X509_USER_PROXY=`voms-proxy-info`

source /cvmfs/cms.cern.ch/cmsset_default.sh

eval `scramv1 runtime -sh`

#echo "Creating voms-proxy..."
#eval vpxy
echo "Proxy info:"
voms-proxy-info -all
echo ""

@cmd@
