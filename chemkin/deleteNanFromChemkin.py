#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This script removes Chebyshev reactions that contain at least one NaN coefficient. 

To use, pass the paths of the Chemkin file on the command-line, e.g.

    $ python removeDuplicates.py /path/to/chem_annotated.inp

"""

import os
import argparse
import numpy

import rmgpy.chemkin
from rmgpy.chemkin import loadChemkinFile, writeKineticsEntry, saveSpeciesDictionary, getSpeciesIdentifier, writeThermoEntry
from rmgpy.rmg.model import CoreEdgeReactionModel
from rmgpy.data.kinetics.library import LibraryReaction
from rmgpy.kinetics import Chebyshev

################################################################################

def main(chemkin, dictionary):
    # Load Chemkin file
    reactionModel = CoreEdgeReactionModel()
    reactionModel.core.species, reactionModel.core.reactions = loadChemkinFile(chemkin, dictionary)

    # Identify reactions to be removed
    reactionList = reactionModel.core.reactions
    nanChebyshevReactionsToRemove = []
    for index1 in range(len(reactionList)):
        reaction1 = reactionList[index1]
        if isinstance(reaction1.kinetics, Chebyshev):
            for t in range(reaction1.kinetics.degreeT):
                for p in range(reaction1.kinetics.degreeP):
                    if numpy.isnan(reaction1.kinetics.coeffs.value_si[t, p]):
                        nanChebyshevReactionsToRemove.append(reaction1)
                        break
                else:
                    continue
                break


    # Remove the identified reactions
    for reaction in nanChebyshevReactionsToRemove:
        reactionList.remove(reaction)

    # Write new Chemkin file
    path = 'overall_chem_annotated_new.inp'
    speciesList = reactionModel.core.species + reactionModel.outputSpeciesList
    rxnList = reactionModel.core.reactions + reactionModel.outputReactionList
    with open(chemkin, 'rb') as old, open(path, 'wb') as new:
        # Copy species and reactions blocks to new Chemkin file
        line = old.readline()
        while 'SPECIES' not in line.upper():
            new.write(line)
            line = old.readline()
        # Species section
        new.write('SPECIES\n')
        for spec in speciesList:
            label = getSpeciesIdentifier(spec)
            new.write('    {0!s:<16}    ! {1}\n'.format(label, str(spec)))
        new.write('END\n\n\n\n')

        while 'THERM' not in line.upper():
            line = old.readline()
        # Thermodynamics section
        new.write('THERM ALL\n')
        new.write('    300.000  1000.000  5000.000\n\n')
        for spec in speciesList:
            new.write(writeThermoEntry(spec))
            new.write('\n')
        new.write('END\n\n\n\n')

        while 'REACTIONS' not in line.upper():
            line = old.readline()
        new.write('REACTIONS    KCAL/MOLE   MOLES\n\n')
        rmgpy.chemkin.__chemkin_reaction_count = 0
        for rxn in rxnList:
            new.write(writeKineticsEntry(rxn, speciesList=speciesList, verbose=True))
            new.write('\n')
        new.write('END\n\n')
    print "New Chemkin file contains {0} reactions.".format(rmgpy.chemkin.__chemkin_reaction_count)

    # Write new species dictionary
    saveSpeciesDictionary('overall_species_dictionary_new.txt', speciesList)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('chemkin', metavar='CHEMKIN', type=str, nargs=1,
                        help='the Chemkin input file to visualize')
    parser.add_argument('dictionary', metavar='DICTIONARY', type=str, nargs=1,
                        help='the species dictionary')
    args = parser.parse_args()

    chemkin = os.path.abspath(args.chemkin[0])
    dictionary = os.path.abspath(args.dictionary[0])
    main(chemkin, dictionary)

