#!/bin/bash
set -xe
Nevents=20
NCORES=6

cmsRun --numThreads $NCORES 1sim.py fileid=1

cmsRun --numThreads $NCORES 2reco.py fileid=1

cmsRun --numThreads $NCORES 3hittree.py fileid=1