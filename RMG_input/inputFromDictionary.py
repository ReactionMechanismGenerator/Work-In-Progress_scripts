#!/usr/bin/env python
# encoding: utf-8

"""
This script imports an RMG input file and a species dictionary file and merges them by
adding all of the species in the dictionary to the species list in the input file.
"""

import argparse
from rmgpy.rmg.input import readInputFile
from rmgpy.rmg.main import RMG
from rmgpy.chemkin import loadSpeciesDictionary


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()
    parser.add_argument('inputPath', metavar='INPUT', type=str, nargs=1,
        help='The path of the RMG input file')
    parser.add_argument('dictionaryPath', metavar='DICTIONARY', type=str, nargs=1,
        help='The path of the RMG dictionary file')
    parser.add_argument('outputPath', metavar='OUTPUT', type=str, nargs=1,
        help='The path where the output will be saved')
    
    args = parser.parse_args()
    inputPath = args.inputPath[0]
    dictionaryPath = args.dictionaryPath[0]
    outputPath = args.outputPath[0]
    
    # Read RMG input file
    rmg = RMG()
    readInputFile(inputPath, rmg)
    # Read species dictionary
    speciesDict = loadSpeciesDictionary(dictionaryPath)
    
    addedMoleFractions = {}
    for key, spec in speciesDict.items():
        rmg.initialSpecies.append(spec)
        addedMoleFractions[spec] = 0
    
    for system in rmg.reactionSystems:
        system.initialMoleFractions.update(addedMoleFractions)
    
    # Save new input file
    rmg.saveInput(outputPath)
    
    
    