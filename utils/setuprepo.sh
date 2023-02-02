#!/bin/bash
#export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch
#source $VO_CMS_SW_DIR/cmsset_default.sh
export CMSVER=CMSSW_11_3_0_pre5
#export SCRAM_ARCH=slc7_amd64_gcc900
cmsrel $CMSVER
cd $CMSVER/src
# alias cmsenv='eval `scramv1 runtime -sh`'
cmsenv
#git cms-init
#git cms-merge-topic WilliamKorcari:HGCal_GAN_production
ln -s ../../IOMC ./
ln -s ../../EDFilters ./
ln -s ../../EDAnalyzers ./
ln -s ../../utils ./

scram b -j8
