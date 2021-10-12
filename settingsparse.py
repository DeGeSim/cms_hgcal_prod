from pprint import pprint

import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
import yaml

cliargs = VarParsing.VarParsing("analysis")

cliargs.register(
    "fileid",
    "1",  # Default value
    VarParsing.VarParsing.multiplicity.singleton,  # singleton or list
    VarParsing.VarParsing.varType.string,  # string, int, or float
    "file identification string",  # Description
)

cliargs.parseArguments()

with open("conf.yaml", "r") as f:
    settingsD = yaml.load(f, Loader=yaml.SafeLoader)


def maptocmstype(element):
    if type(element) is int:
        return cms.int32(element)
    if type(element) is float:
        return cms.double(element)
    if type(element) is bool:
        return cms.bool(element)
    if type(element) is str:
        return cms.string(element)
    if type(element) is dict:
        return {key: maptocmstype(val) for key, val in element.items()}
    if type(element) in (list, tuple):
        return type(element)(maptocmstype(sube) for sube in element)


settingsD = maptocmstype(settingsD)

pprint(settingsD)
