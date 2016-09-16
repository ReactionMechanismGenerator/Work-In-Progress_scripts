#runs addTraining method to the specifications in this file- updated by mark July 29th
import math
from rmgpy.data.reference import *
from addTraining import addTraining

#input optional parameters to put in RMG
shortDesc0=u"""CBS-QB3 calculation with 1-d rotor treatment at B3LYP/631G(d)"""
longDesc0=u"""
Quantum chemistry calculations CBS-QB3 calculation with 1-d rotor treatment at 
B3LYP/631G(d)" using Gaussian 03 and Gaussian 09. High-pressure-limit rate 
coefficient computed TST with Eckart Tunnelling"
"""
reference0=Article(
        authors = ["K. Wang", "S. Villano", "A. Dean"],
        title = u'Reactions of allylic radicals that impact molecular weight growth kinetics',
        journal = "Phys. Chem. Chem. Phys.",
        volume = "17",
        pages = """6255-6273""",
        year = "2015",
)
#required conditions from paper
family='H_Abstraction'
A=[18.06,26.3,18.4,64.3,29.7,19.17,16.9,23.0,120.0,80.50,84.7,275.00,333.0,1110.0,429.5,427.5,945.0,890.0,21.0,299.7,10.5,8.32,667.5,570.0,1500,590.0,1.52e3,6.07]
n=[3.27, 3.24,3.27,3.13,3.22,3.28,3.30,3.22,3.09,3.12,3.12,3.08,2.92,2.73,2.95,2.93,2.84,2.86,3.21,2.84,3.31,3.29,2.90,2.91,2.80,2.90,2.93,3.42]
Ea=[6.85,7.03,7.15,7.32,7.04,6.73,6.48,5.48,5.73,5.82,5.78,4.37,4.91,5.32,5.24,4.01,3.72,3.77,5.77,5.63,5.46,5.72,3.45,4.67,3.38,5.06,2.93,8.66]
reactant1smiles=['CC=CC', 'CC=CCC', 'C=CC', 'C=C(C)CC', 'C=C(C)C', 'CC=C(C)C', 'CC=C(C)C', 
                 'C=CCCC=C', 'CC=CCC', 'C=CCC', 'C=CCCC', 'C=C(C)CC', 
                 'C=CC(C)C','C=CC(C)CC', 
                 'C1C=CCC1', 'C1C=CCCC1', 
                 'CC1=CCCC1', 'CC1=CCCCC1', 
                 'C=CCC=C', 'C=CC=CC', 'C=CC(C)=CC', 'C=CC=C(C)C', 
                 'C1=CCC=CC1', 'C1=CC=CCC1', 'CC1=CC=CC1', 'C1=CC=C1', 'c1ccccc1', 'c1cccc1']
reactant2smiles=['[CH3]']*28
degeneracy=[6,3,3,3,6,3,3,4,2,2,2,2,1,1,4,4,1,1,2,3,3,6,4,4,1,2,6,3]

# print len(A)
# print len(n)
# print len(Ea)
# print len(reactant1smiles)
# print len(reactant2smiles)
# print len(degeneracy)

addTraining(family,A,n,Ea,reactant1smiles,
        degeneracy,
		smiles2=reactant2smiles,			#for bimolecular reactions
        rank=3,
		reference=reference0,
		shortDesc=shortDesc0,
		longDesc=longDesc0,
		directory=r'c:/')
