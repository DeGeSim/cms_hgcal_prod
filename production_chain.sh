#!/bin/bash
set -xe
Nevents=2
NCORES=6

cmsRun --numThreads $NCORES hgcal_gen/1_gen.py

cmsRun --numThreads $NCORES hgcal_gen/2_digi.py

cmsRun --numThreads $NCORES hgcal_gen/3_reco.py

cmsRun --numThreads $NCORES\
    EDAnalyzers/TreeMaker/python/ConfFile_cfg.py \
    inputFiles=step1.root \
    genEleFilter=0 \
    genPartonFilter=0 \
    isGunSample=1 \
    onRaw=0 \
    storeSimHit=1 \
    storeRecHit=1 \
    debugFile=0 \
    maxEvents=$Nevents