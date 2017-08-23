#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This script generates transport data given a species dictionary.
Results from all possible sources are retrieved.

To use, pass the path of the species dictionary on the command-line, e.g.

    $ python generateTransportData.py /path/to/species_dictionary.txt

The resulting tran.dat file is placed in the execution directory.
"""

import os
import argparse

from rmgpy import settings
from rmgpy.data.transport import TransportDatabase
from rmgpy.chemkin import loadSpeciesDictionary, saveTransportFile
import rmgpy.constants as constants

################################################################################

def main(dictionary):

    database = TransportDatabase()
    database.load(path=os.path.join(settings['database.directory'], 'transport'), libraries=None)

    speciesDict = loadSpeciesDictionary(dictionary)

    path = 'tran.dat'

    with open(path, 'w') as f:
        f.write("! {0:15} {1:8} {2:9} {3:9} {4:9} {5:9} {6:9} {7:9}\n".format('Species', 'Shape', 'LJ-depth', 'LJ-diam',
                                                                              'DiplMom', 'Polzblty', 'RotRelaxNum',
                                                                              'Data'))
        f.write(
            "! {0:15} {1:8} {2:9} {3:9} {4:9} {5:9} {6:9} {7:9}\n".format('Name', 'Index', 'epsilon/k_B', 'sigma', 'mu',
                                                                          'alpha', 'Zrot', 'Source'))
        for label, spec in speciesDict.iteritems():
            transportDataList = database.getAllTransportProperties(spec)
            for transportData, _, _ in transportDataList:
                if (not transportData):
                    missingData = True
                else:
                    missingData = False

                if missingData:
                    f.write('! {0:19s} {1!r}\n'.format(label, transportData))
                else:
                    f.write('{0:19} {1:d}   {2:9.3f} {3:9.3f} {4:9.3f} {5:9.3f} {6:9.3f}    ! {7:s}\n'.format(
                        label,
                        transportData.shapeIndex,
                        transportData.epsilon.value_si / constants.R,
                        transportData.sigma.value_si * 1e10,
                        (transportData.dipoleMoment.value_si * constants.c * 1e21 if transportData.dipoleMoment else 0),
                        (transportData.polarizability.value_si * 1e30 if transportData.polarizability else 0),
                        (transportData.rotrelaxcollnum if transportData.rotrelaxcollnum else 0),
                        transportData.comment,
                    ))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('dictionary', metavar='DICTIONARY', type=str, nargs=1,
                        help='the RMG species dictionary')

    args = parser.parse_args()

    dictionary = os.path.abspath(args.dictionary[0])

    main(dictionary)

