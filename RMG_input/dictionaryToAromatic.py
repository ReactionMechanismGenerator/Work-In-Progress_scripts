#!/usr/bin/env python
# encoding: utf-8

"""
This script takes a species dictionary and converts adjacency lists
of aromatic species to the aromatic (benzene bond) form.
"""

import argparse
from rmgpy.species import Species
from rmgpy.chemkin import loadSpeciesDictionary, saveSpeciesDictionary
from rmgpy.molecule.resonance import generateAromaticResonanceStructures


def main(inputPath, outputPath):

    # Replicate rmgpy.chemkin.loadSpeciesDictionary here with some modifcations
    # We don't need to compare against inerts, and we want to maintain the same species order
    speciesList = []
    with open(inputPath, 'r+b') as f:
        adjlist = ''
        for line in f:
            if line.strip() == '' and adjlist.strip() != '':
                # Finish this adjacency list
                species = Species().fromAdjacencyList(adjlist)

                # Try to generate the aromatic form of the molecule
                aromatic = generateAromaticResonanceStructures(species.molecule[0])
                if aromatic:
                    species.molecule = aromatic

                speciesList.append(species)
                adjlist = ''
            else:
                if "InChI" in line:
                    line = line.split()[0] + '\n'
                if '//' in line:
                    index = line.index('//')
                    line = line[0:index]
                adjlist += line

    # Write to dictionary file
    saveSpeciesDictionary(outputPath, speciesList)


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
