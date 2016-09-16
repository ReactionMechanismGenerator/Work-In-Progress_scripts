#!/usr/bin/env python
# encoding: utf-8

"""
This script reads a text file with a list of species names and
SMILES strings and converts it to an RMG dictionary.

The text file should have one species per line, with the species
name/identifier first, followed by the SMILES string.

Comments marked with '!' are ignored.
"""

import argparse
import re
from rmgpy.species import Species
from rmgpy.chemkin import saveSpeciesDictionary


def main(inputPath, outputPath):

    species = []

    with open(inputPath, 'r+b') as f:
        for line0 in f:
            line1 = line0.strip()
            if line1 and line1[0] != '!':
                # Parse line using regex
                match = re.match(r"([^\s!]+)\s+([^\s!]+)", line1.strip())

                # Get label and SMILES string
                label = match.group(1)
                smiles = match.group(2)

                # Save in species dictionary
                spec = Species().fromSMILES(smiles)
                spec.label = label
                species.append(spec)

    # Write to dictionary file
    saveSpeciesDictionary(outputPath, species)


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('inputPath', metavar='INPUT', type=str, nargs=1,
        help='The path of the SMILES input file')
    parser.add_argument('outputPath', metavar='OUTPUT', type=str, nargs=1,
        help='The path to save the output file')

    args = parser.parse_args()
    inputPath = args.inputPath[0]
    outputPath = args.outputPath[0]

    main(inputPath, outputPath)
