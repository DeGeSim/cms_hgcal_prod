########## Filters ##########

# from EDFilters.MyFilters.GenParticleFilter_cfi import *
import FWCore.ParameterSet.Config as cms
from .GenParticleFilter_cfi import *

def apply_filters(process, options):
    # Gen-ele filter
    process.GenParticleFilter_ele = GenParticleFilter.clone()
    process.GenParticleFilter_ele.atLeastN = cms.int32(1)
    process.GenParticleFilter_ele.pdgIds = cms.vint32(11)
    process.GenParticleFilter_ele.minPt = cms.double(10)
    process.GenParticleFilter_ele.minEta = cms.double(1.479)
    process.GenParticleFilter_ele.maxEta = cms.double(3.1)
    process.GenParticleFilter_ele.isGunSample = cms.bool(bool(options.isGunSample))
    # process.GenParticleFilter_ele.debug = cms.bool(True)

    process.filter_seq_genEle = cms.Sequence()

    if options.genEleFilter:
        process.filter_seq_genEle = cms.Sequence(process.GenParticleFilter_ele)

    # Gen-pho filter
    process.GenParticleFilter_pho = GenParticleFilter.clone()
    process.GenParticleFilter_pho.atLeastN = cms.int32(1)
    process.GenParticleFilter_pho.pdgIds = cms.vint32(22)
    process.GenParticleFilter_pho.minPt = cms.double(10)
    process.GenParticleFilter_pho.minEta = cms.double(1.479)
    process.GenParticleFilter_pho.maxEta = cms.double(3.1)
    process.GenParticleFilter_pho.isGunSample = cms.bool(bool(options.isGunSample))
    # process.GenParticleFilter_pho.debug = cms.bool(True)

    process.filter_seq_genPho = cms.Sequence()

    if options.genPhoFilter:
        process.filter_seq_genPho = cms.Sequence(process.GenParticleFilter_pho)

    # Gen-parton filter
    process.GenParticleFilter_part = GenParticleFilter.clone()
    process.GenParticleFilter_part.atLeastN = cms.int32(1)
    process.GenParticleFilter_part.pdgIds = cms.vint32(1, 2, 3, 4, 5, 21)
    process.GenParticleFilter_part.minPt = cms.double(10)
    process.GenParticleFilter_part.minEta = cms.double(1.479)
    process.GenParticleFilter_part.maxEta = cms.double(3.1)
    process.GenParticleFilter_part.isGunSample = cms.bool(bool(options.isGunSample))
    # process.GenParticleFilter_part.debug = cms.bool(True)

    process.filter_seq_genPart = cms.Sequence()

    if options.genPartonFilter:
        process.filter_seq_genPart = cms.Sequence(process.GenParticleFilter_part)
    return process
