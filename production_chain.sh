#!/bin/bash
set -xe
Nevents=10000000
NCORES=6
cd hgcal_gen/
cmsRun step1_SingleElectronPt15Eta1p7_2p7_cfi_GEN_SIM.py \
    maxEvents=$Nevents \
    --numThreads $NCORES

cmsRun step2_DIGI_L1TrackTrigger_L1_DIGI2RAW_HLT.py \
    maxEvents=$Nevents \
    --numThreads $NCORES

cmsRun step3_RAW2DIGI_L1Reco_RECO_RECOSIM_PAT_VALIDATION_DQM.py \
    maxEvents=$Nevents \
    --numThreads $NCORES

cd ..
cmsRun EDAnalyzers/TreeMaker/python/ConfFile_cfg.py \
    inputFiles=step1.root \
    genEleFilter=0 \
    genPartonFilter=0 \
    isGunSample=1 \
    onRaw=0 \
    storeSimHit=1 \
    storeRecHit=1 \
    debugFile=0 \
    maxEvents=$Nevents