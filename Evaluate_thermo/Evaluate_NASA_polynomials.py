#!/usr/bin/env python
# encoding: utf-8

"""
This script takes a chemkin input file that includes both SPECIES and THERMO blocks, and outputs the thermo properties
(heat capacity, enthalpy, entropy and gibb's free energy) of each species by evaluating their NASA polynomials at specified
temperatures. The evaluated thermo properties are output in both tabular format (in "evaluated_NASA_polynomial.py"),
and as plots (in individual .pdf files for each species).

Run using the  following commands:

python Evaluate_NASA_polynomials.py "path_to_chemkin_file"
"""

import argparse
import os
from rmgpy.chemkin import loadChemkinFile, getSpeciesIdentifier
from rmgpy import settings
from rmgpy.cantherm.output import prettify
import numpy.linalg

#Specify the temperatures at which to evaluate the NASA polynomials
T_list = [300, 400, 500, 600, 800, 1000, 1500, 2000, 2400, 3000]

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()
    parser.add_argument('chemkinPath', metavar='CHEMKIN', type=str, nargs=1,
        help='The path of the chemkin file')
    
    args = parser.parse_args()
    chemkinPath = args.chemkinPath[0]

    #Save output in same directory as input
    outputDirectory = os.path.dirname(os.path.abspath(chemkinPath))

    # Initialize (and clear!) the output files for the job
    outputFile = os.path.join(outputDirectory, 'evaluated_NASA_polynomial.py')
    with open(outputFile, 'w') as f:
        pass

    #Load the checmkin file
    speciesList, reactionList = loadChemkinFile(chemkinPath)

    # Make full species identifier the species labels
    for species in speciesList:
        species.label = getSpeciesIdentifier(species)
        species.index = -1

    # Evaluate the NASA polynomials for thermo properties of each species at the specified T's and save in tabular format
    # to output file
    for i in range(len(speciesList)): 
        species = speciesList[i]
        if species.thermo:

            f = open(outputFile, 'a')

            f.write('# Thermodynamics for {0}:\n'.format(species.label))
            H298 = species.getThermoData().getEnthalpy(298) / 4184.
            S298 = species.getThermoData().getEntropy(298) / 4.184
            f.write('#   Enthalpy of formation (298 K)   = {0:9.3f} kcal/mol\n'.format(H298))
            f.write('#   Entropy of formation (298 K)    = {0:9.3f} cal/(mol*K)\n'.format(S298))
            f.write('#    =========== =========== =========== =========== ===========\n')
            f.write('#    Temperature Heat cap.   Enthalpy    Entropy     Free energy\n')
            f.write('#    (K)         (cal/mol*K) (kcal/mol)  (cal/mol*K) (kcal/mol)\n')
            f.write('#    =========== =========== =========== =========== ===========\n')
            for T in T_list:
                Cp = species.getThermoData().getHeatCapacity(T) / 4.184
                H = species.getThermoData().getEnthalpy(T) / 4184.
                S = species.getThermoData().getEntropy(T) / 4.184
                G = species.getThermoData().getFreeEnergy(T) / 4184.
                f.write('#    {0:11g} {1:11.3f} {2:11.3f} {3:11.3f} {4:11.3f}\n'.format(T, Cp, H, S, G))
            f.write('#    =========== =========== =========== =========== ===========\n')

            string = 'thermo(label={0!r}, thermo={1!r})'.format(species.label, species.getThermoData())
            f.write('{0}\n\n'.format(prettify(string)))

            f.close()
        else:
            print ('Species {0} did not contain any thermo data and was omitted from the output file.'.format(str(species)))

    #Convert T_list into an array
    T_list = numpy.asarray(T_list)

    # Plot the thermo properties of each species at the specified T's and save as .pdf files
    for i in range(len(speciesList)):
        species = speciesList[i]
        if species.thermo:

            # Skip this step if matplotlib is not installed
            try:
                import pylab
            except ImportError:
                print('Pylab not available so no thermo plots generated.')
                break

            Cplist1 = numpy.zeros_like(T_list)
            Hlist1 = numpy.zeros_like(T_list)
            Slist1 = numpy.zeros_like(T_list)
            Glist1 = numpy.zeros_like(T_list)

            thermo = species.getThermoData()
            for i in range(T_list.shape[0]):
                Cplist1[i] = thermo.getHeatCapacity(T_list[i])
                Slist1[i] = thermo.getEntropy(T_list[i])
                Hlist1[i] = thermo.getEnthalpy(T_list[i]) * 0.001
                Glist1[i] = thermo.getFreeEnergy(T_list[i]) * 0.001

            fig = pylab.figure(figsize=(10, 8))

            pylab.subplot(2, 2, 1)
            pylab.plot(T_list, Cplist1 / 4.184, '-b')
            pylab.xlabel('Temperature (K)')
            pylab.ylabel('Heat capacity (cal/mol*K)')
            pylab.legend(['statmech', 'thermo'], loc=4)

            pylab.subplot(2, 2, 2)
            pylab.plot(T_list, Slist1 / 4.184, '-b')
            pylab.xlabel('Temperature (K)')
            pylab.ylabel('Entropy (cal/mol*K)')

            pylab.subplot(2, 2, 3)
            pylab.plot(T_list, Hlist1 / 4.184, '-b')
            pylab.xlabel('Temperature (K)')
            pylab.ylabel('Enthalpy (kcal/mol)')

            pylab.subplot(2, 2, 4)
            pylab.plot(T_list, Glist1 / 4.184, '-b')
            pylab.xlabel('Temperature (K)')
            pylab.ylabel('Gibbs free energy (kcal/mol)')

            fig.subplots_adjust(left=0.10, bottom=0.08, right=0.95, top=0.95, wspace=0.35, hspace=0.20)
            pylab.savefig(os.path.join(outputDirectory, '{0}_thermo.pdf'.format(str(species))))
            pylab.close()

        else:
            print (
            'Species {0} did not contain any thermo data and was omitted from the thermo plots.'.format(str(species)))