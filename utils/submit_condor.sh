#!/usr/bin/env python


import os


cmd = (
    "python utils/run_condor.py "
    #"--cmsRunFile=gammaGun_noTracker_cfg.py "
    "--processName=photon-gun "
    #"--outputDir=output "
    "--nJob=2000 "
    #"--nUnitPerJob=10 "
    #"--cmsRunOptions=\"minE=49.99 maxE=50.01 minEta=1.99 maxEta=2.01\" "
    #"--test "
)

print cmd


os.system(cmd)

