#!/bin/bash
set -xe
Nevents=20
NCORES=6

cmsRun --numThreads $NCORES hgcal_gen/1_gen.py fileid=1

cmsRun --numThreads $NCORES hgcal_gen/2_digi.py fileid=1

#cmsRun --numThreads $NCORES hgcal_gen/3_reco.py fileid=1

cmsRun --numThreads $NCORES \
   EDAnalyzers/TreeMaker/python/ConfFile_cfg.py \
   fileid=1