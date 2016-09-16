"""
purpose: 
helps with inputing large data sets into training reactions.  the reactions must be in the form A*T^n*exp(-Ea/RT).

details: 
You give it a reaction family, path to the directory holding RMG-database, reactant SMILES (as an array), 
and rate information (as an array).  Optional arguments include 
references, short description, long description, and rank.  The script will append 
all the reactions to dictionary.txt, and reactions.py in the corresponding family's training reaction folder 

to use: 
import the method addTraining in python, and call the method using the arguments in the function.  
An example code exists with the name addTraining-run.py (watch out: the example will change your dictionary)

limitations:
1. current units supported are : A is in cm3/mols (or s-1 for unimolecular rxns) and E is in kcal/mol.
2. written to work with bimolecular reactions but never tested. let me know if it doesn't work

last edited by Mark on July 29th
"""
import math
import re
import csv
import numpy as np
from rmgpy.data.kinetics.database import *
from rmgpy.data.rmg import RMGDatabase
from rmgpy.molecule.molecule import *
from rmgpy.data.rmg import RMGDatabase
from rmgpy.data.base import ForbiddenStructures
from rmgpy.data.reference import *
from rmgpy.data.kinetics.family import KineticsFamily
import os.path

	
def addTraining(reactionFamily, A, n, Ea, smiles1, 
			smiles2=None, 
			rank=3, 
			shortDesc='', 
			longDesc='',
			reference=None,
			directory='/home/mark/',
			startingIndex=1,):

	"""
this method inputs a reaction into training reactions using just the reactants and reaction family.  
A is in cm3/mols and E is in kcal/mol.  A, n, Ea, smiles1, and smiles2 are lists containting multiple reactions 
from the family.  

This currently does not if multiple reaction paths found, or with bimolecular reactions. 
	"""
	
	path=os.path.join(directory,'RMG-database/input/kinetics/families')
	#throw errors
	if len(A)!=len(n) or len(n)!=len(Ea) or len(n)!=len(smiles1):
# 		error('lengths of matrixes not equal')
		print('lengths of matrixes not equal')
	#load databases
# 	kdatabase=KineticsDatabase()
	rdatabase=RMGDatabase()
	rdatabase.load(os.path.join(directory,'RMG-database/input'),kineticsFamilies=[reactionFamily])
# 	kdatabase.loadFamilies(path,families=[reactionFamily])
	#convert the smiles into molecule objects
	reactant1=[]
	reactant2=[]
	
	for index, smile in enumerate(smiles1):
		reactant1.append(Molecule(SMILES=smiles1[index]))
		if smiles2:
			reactant2.append(SMILES=Molecule(smiles2[index]))
#find possible reaction paths for reactants
	reactionList=[]
	for index, smile in enumerate(smiles1):
		if smiles2:
			reactantList=[reactant1[index], reactant2[index]]
		else:
			reactantList=[reactant1[index]]
		reactions=rdatabase.kinetics.families[reactionFamily].generateReactions(reactantList)
		
		#throw error if multiple paths exist
		if len(reactions)>1:
			print "multiple reactions found for "+smile
		reaction=reactions[0]
		
		
		#####now we have to do some workaround to obtain the graph values for the reactants
		####note for BIMOLECULAR reactions: this part may mix up molecule and species lists. fix
		#convert to list of molecules from list of species to use 'KineticsFamily.applyRecipe'
		products_with_molecules=[]
		for product in reaction.products:
			products_with_molecules.append(product.molecule[0])
		#find reactants with the atom labels using 'applyRecipe'
		newreactants=rdatabase.kinetics.families[reactionFamily].applyRecipe(products_with_molecules,False)
		print newreactants
		#convert back to species list from molecule list
		react_with_species=[]
		for reactant in newreactants:
			react_with_species.append(Species(molecule=[reactant]))
		reaction.reactants=react_with_species


		#make labels for species...the = in species labels leads to erros parsing the reaction, so it is replaced with !
		reaction.reactants[0].label=reaction.reactants[0].molecule[0].toSMILES().replace('=','!')+''+str(index)+'r1'
		reaction.products[0].label=reaction.products[0].molecule[0].toSMILES().replace('=','!')+''+str(index)+'p1'
		if len(reaction.reactants)>1:
			reaction.reactants[1].label=reaction.reactants[1].molecule[0].toSMILES().replace('=','!')+''+str(index)+'r2'
		if len(reactions[0].products)>1:
			reaction.products[1].label=reaction.products[1].molecule[0].toSMILES().replace('=','!')+''+str(index)+'p2'
			
		print reaction.__str__()
		reactionList.append(reaction)
		#print reactionList	
	
	#write to dictionary file
	with open(os.path.join(path,reactionFamily,'training','dictionary.txt'), 'ab') as f:
		f.write('\n')
		for index, reaction in enumerate(reactionList):
			for species in reaction.reactants:
				f.write(species.toAdjacencyList())
				f.write('\n')
			for species in reaction.products:
				f.write(species.toAdjacencyList())
				f.write('\n')
	print 'finished writing dictionary to file'
	#create strings for putting together a reactions.py entry
	startstr="""
entry(
    index = """
	afterindex=""",
    label = \""""
	afterlabel="""\",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = ("""
	if smiles2:		#if using bimolecular reactions, need to change the units for prefactor
		afterprefactor=""", 'cm^3/(mol*s), '*|/', 2.51189),
        n = """
	else:
		afterprefactor=""", 's^-1', '*|/', 2.51189),
        n = """
	afterexponential=""",
        Ea = ("""
	afterEa=""", 'kcal/mol', '+|-', 1.5),
        T0 = (1, 'K'),
    ),
    reference = """
	afterreference=""",
    rank = """
	afterrank=""",
    referenceType = \"theory\",
    shortDesc =u\' """
	aftershortdesc="""\',
    longDesc = u\"\"\""""
	endstr="""\"\"\",
)"""


	#write to reactions.py file - 'ab' appends to bottom of the file
	with open(os.path.join(path,reactionFamily,'training','reactions.py'), 'ab') as f:
		for index, reaction in enumerate(reactionList):
			f.write('\n\n')
			f.write(startstr)
			f.write(str(startingIndex+index))
			f.write(afterindex)
			f.write(reaction.__str__())
			f.write(afterlabel)
			f.write(str(A[index]))
			f.write(afterprefactor)
			f.write(str(n[index]))
			f.write(afterexponential)
			f.write(str(Ea[index]))
			f.write(afterEa)
			if reference:
				f.write(reference.toPrettyRepr())
			else:
				f.write('None')
			f.write(afterreference)
			f.write(str(rank))
			f.write(afterrank)
			f.write(shortDesc)
			f.write(aftershortdesc)
			f.write(longDesc)
			f.write(endstr)
	
	print 'finished writing reactions.py to file'
