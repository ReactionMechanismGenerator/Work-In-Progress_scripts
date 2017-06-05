#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This script generates transport data given a species dictionary.
All libraries are used by default (GRI-Mech and PrimaryTransportLibrary).

To use, pass the path of the species dictionary on the command-line, e.g.

    $ python generateTransportData.py /path/to/species_dictionary.txt

The resulting tran.dat file is placed in the execution directory.
"""

import os
import argparse

from rmgpy import settings
from rmgpy.data.rmg import RMGDatabase
from rmgpy.chemkin import loadSpeciesDictionary, saveTransportFile

################################################################################

def main(dictionary):

    database = RMGDatabase()
    database.loadTransport(path=os.path.join(settings['database.directory'], 'transport'), transportLibraries=None)

    speciesDict = loadSpeciesDictionary(dictionary)

    saveTransportFile('tran.dat', speciesDict.values())


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('dictionary', metavar='DICTIONARY', type=str, nargs=1,
                        help='the RMG species dictionary')

    args = parser.parse_args()

    dictionary = os.path.abspath(args.dictionary[0])

    main(dictionary)

