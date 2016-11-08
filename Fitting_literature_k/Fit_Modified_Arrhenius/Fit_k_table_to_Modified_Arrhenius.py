#!/usr/bin/env python
# encoding: utf-8

"""
This script takes a table of high-P k(T) values in a csv file with specific units (k[=]cm^3/mol*s, T[=]K), 
and in a specific format (see example csv file '2012_Kislov_Phenyl_Propene_k_table_example.csv') 
and fits a Modified Arrhenius expression to each column using the built-in features of Cantherm.

The fits are output in CHEMKIN format, and if pylab is available a plot of each fit is also generated.
"""
import math
import numpy
import logging
import csv
import re
import os
import os.path

from rmgpy.kinetics.arrhenius import Arrhenius

###########Provide inputs#################################################################################
#Specify location of csv file with table of T and k data
#Note that the units must be correctly defined such that:
#T is in K
#k is in cm, mol, s units

k_data_path = '2012_Kislov_Phenyl_Propene_k_table_example.csv'

#Overall name of system being fit
system_name = 'Phenyl + Propene'

###########Read table of k data (must be in mol, cm, s, K units)########################################
f = csv.reader(open(k_data_path, 'r'))

columns = zip(*f)

kdata = numpy.zeros((len(columns[0][:])-3,len(columns)-1,), numpy.float64)

Tdata = numpy.array(columns[0][3:],dtype=numpy.float64)

reaction_string = []

order = numpy.zeros(len(columns)-1,dtype=numpy.float64)

for col_num, col in enumerate(columns[1:]):
    reaction_string.append(col[0])
    order[col_num] = (col[1])
    for row_num, row in enumerate(col[3:]):
	kdata[row_num,col_num] = numpy.array(row,dtype=numpy.float64)

########################Fit Modified Arrhenius expression to each column of kdata###################################################

kinetics = []
kunits = []

for rxn_number in range(kdata.shape[1]):
    kunits.append({1: 's^-1', 2: 'cm^3/(mol*s)', 3: 'cm^6/(mol^2*s)'}[order[rxn_number]])
    kinetics.append(Arrhenius().fitToData(Tdata, kdata[:,rxn_number], kunits=kunits[rxn_number]))

########################Output Modified Arrhenius fit to CHEMKIN input file format#####################################
# Initialize (and clear!) the output files for the job
outputDirectory = os.path.dirname(os.path.abspath(k_data_path))

chemkinFile = os.path.join(outputDirectory, system_name.replace(" ", "_") + '.inp')
with open(chemkinFile, 'w') as f:
    pass


Tcount = Tdata.shape[0]
reactioncount = kdata.shape[1]

f = open(chemkinFile, 'a')

first_line = 'REACTIONS MOL KCAL/MOLE'

f.write('{0}\n'.format(first_line))

for rxn_number in range(reactioncount):

    factor = 1e6 ** (order[rxn_number]-1)
 
    string = '{0!s:51} {1:9.3e} {2:9.3f} {3:9.3f}\n'.format(
	reaction_string[rxn_number], 
	kinetics[rxn_number].A.value_si * factor, 
	kinetics[rxn_number].n.value_si, 
	kinetics[rxn_number].Ea.value_si / 4184.)
                
    f.write('{0}\n'.format(string))
      
last_line = 'END'

f.write('{0}\n'.format(last_line))
      
f.close()

########################Evaluate Modified Arrhenius Fit at T of kdata##############################################
kfit = numpy.zeros((Tcount, reactioncount))
for rxn in range(reactioncount):
    for t in range(Tcount):
        kfit[t,rxn] = kinetics[rxn].getRateCoefficient(Tdata[t])
	kfit[t,rxn] *= 1e6 ** (order[rxn]-1)

########################Plot comparison of kdata and kfit if pylab is available##############################
# Skip this step if matplotlib is not installed
try:
    import matplotlib.pyplot as pylab
    import matplotlib.ticker
    import matplotlib.cm
    cm = matplotlib.cm.jet


    for rxn_number in range(reactioncount):

    	pylab.semilogy(1000.0 / Tdata, kdata[:,rxn_number], 'ok')
    	pylab.semilogy(1000.0 / Tdata, kfit[:,rxn_number], '-k')
    	pylab.xlabel('1000 / Temperature (1000/K)')
    	pylab.ylabel('Rate coefficient ({0})'.format(kunits[rxn_number]))
    	pylab.title(reaction_string[rxn_number])
    	pylab.savefig(os.path.join(outputDirectory, reaction_string[rxn_number].replace(" ", "_") + '_Modified_Arrhenius_fit.pdf'))
    	pylab.close()
    
except ImportError:
    print('Pylab cannot be imported, so fit plots were not generated')