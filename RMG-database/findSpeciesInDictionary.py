'''
This script is used to get species names from a library for a molecule of interest.
If using a dictionary from training reaction, you could have multiple names associated with a molecule
because of different atom labelling.

Last modified by Nathan Yee on Oct 24 2016
'''
import os.path
from rmgpy import settings
from rmgpy.molecule import Molecule
from rmgpy.chemkin import loadSpeciesDictionary
from rmgpy.species import Species


databaseDirectory = settings['database.directory']
dictionaryPath = os.path.join(databaseDirectory,"kinetics/families/H_Abstraction/training/dictionary.txt")
speciesDict = loadSpeciesDictionary(dictionaryPath)

outputPath = os.path.join(databaseDirectory, "output.txt") #Where to print out file of results


targetSmiles = ["c1ccccc1C", "c1ccccc1[CH2]", "[OH]", 'O', ]
nameDict = {}
for target in targetSmiles:
    nameDict[target] = []
    targetSpecies = Species().fromSMILES(target)
    for name, species in speciesDict.iteritems():
        if species.isIsomorphic(targetSpecies): nameDict[target].append(name)

#write out the output
with open (outputPath, 'wb') as outFile:
    for target in targetSmiles:
        outFile.write(target+":\n")
        for name in nameDict[target]:
            outFile.write(name+"\n")
        outFile.write("\n")

