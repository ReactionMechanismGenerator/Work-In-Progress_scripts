#!/usr/bin/env python
# encoding: utf-8

"""
This script takes a species dictionary and converts it to a text file with 2 columns:
The first column is the SMILES for each species, and the second is the species name.
"""

import argparse
from rmgpy.species import Species
from rmgpy.chemkin import loadSpeciesDictionary, saveSpeciesDictionary, ChemkinError


def main(inputPath, outputPath):

    # Load species dictionary
    speciesDict = loadSpeciesDictionary(inputPath)

    """
    Save the given list of `species` as adjacency lists in a text file `path` 
    on disk.

    If `oldStyle==True` then it saves it in the old RMG-Java syntax.
    """
    with open(outputPath, 'w') as f:
        for label, spec in speciesDict.iteritems():
            try:
                f.write(spec.molecule[0].toSMILES())
                f.write(' ' + label)
            except:
                raise ChemkinError('Ran into error saving dictionary for species {0}. Please check your files.'.format(
                    label))
            f.write('\n')
if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('inputPath', metavar='INPUT', type=str, nargs=1,
        help='The path of the input species dictionary')
    parser.add_argument('outputPath', metavar='OUTPUT', type=str, nargs=1,
        help='The path to save the output file')

    args = parser.parse_args()
    inputPath = args.inputPath[0]
    outputPath = args.outputPath[0]

    main(inputPath, outputPath)
