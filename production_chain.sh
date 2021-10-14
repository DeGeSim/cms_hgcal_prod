#!/bin/bash
source ~/.profile
source ~/utils/bashFunctionCollection.sh 

export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch
source $VO_CMS_SW_DIR/cmsset_default.sh
export CMSVER=CMSSW_11_3_0_pre5
export SCRAM_ARCH=slc7_amd64_gcc900
eval `scramv1 runtime -sh`

NCORES=1

logandrun cmsRun --numThreads $NCORES 1sim.py fileid=$1

logandrun cmsRun --numThreads $NCORES 2reco.py fileid=$1

logandrun cmsRun --numThreads $NCORES 3hittree.py fileid=$1
