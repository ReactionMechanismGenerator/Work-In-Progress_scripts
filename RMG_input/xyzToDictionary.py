#!/usr/bin/env python
# encoding: utf-8

"""
This script reads a text file with a list of OpenBabel format xyz coordinates
and converts it to an RMG dictionary.

Separate entries should be separated by at least one empty line

OpenBabel format:
    First line is the number of atoms
    Second line is for species name
    Remaining lines are for Cartesian coordinates (Element  X  Y  Z)
"""

import argparse
import pybel

from rmgpy.molecule.parser import fromOBMol
from rmgpy.molecule.molecule import Molecule
from rmgpy.species import Species
from rmgpy.chemkin import saveSpeciesDictionary


def main(inputPath, outputPath):

    species = []

    with open(inputPath, 'r+b') as f:
        xyz = ''
        label = ''
        lineNum = 0
        specNum = 0
        for line0 in f:
            line1 = line0.strip()
            if line1:
                # The line is not blank, so add it to the xyz string
                xyz += line1 + '\n'
                lineNum += 1
                # In openbabel xyz format, the second line is a comment/name
                if lineNum == 2:
                    label = line1
            else:
                # Check if we've found a valid entry
                if not xyz:
                    continue
                lines = xyz.strip().split('\n')
                # The first line is the number of atoms, which should match the number of lines minus the first two
                if int(lines[0]) != len(lines) - 2:
                    print 'Invalid entry found'
                    xyz = ''
                    lineNum = 0
                    continue

                # Assume that we finished reading this entry
                specNum += 1

                # Make sure the label is appropriate
                label = label.split()[0]
                label = label.strip()
                if not label:
                    label = 'SPC({0})'.format(specNum)

                # Read xyz string using openbabel
                obmol = pybel.readstring('xyz', xyz)

                # Convert to RMG molecule object
                mol = fromOBMol(Molecule(), obmol.OBMol)

                # Save in species dictionary
                spec = Species(molecule=[mol])
                spec.label = label
                species.append(spec)

                # Reset our counter
                xyz = ''
                lineNum = 0

    # Write to dictionary file
    saveSpeciesDictionary(outputPath, species)


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('inputPath', metavar='INPUT', type=str, nargs=1,
        help='The path of the xyz input file')
    parser.add_argument('outputPath', metavar='OUTPUT', type=str, nargs=1,
        help='The path to save the output file')

    args = parser.parse_args()
    inputPath = args.inputPath[0]
    outputPath = args.outputPath[0]

    main(inputPath, outputPath)
