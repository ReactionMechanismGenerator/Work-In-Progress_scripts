#!/usr/bin/env python
# encoding: utf-8

"""
This script takes a table of k(T,P) values in a csv file with specific units (k[=]cm^3/mol*s, T[=]K, P[=]atm), 
and in a specific format (see example csv file '2006_Joshi_OH_CO_k1_table_example.csv') 
and fits a Chebyshev polynomial of specified order to it using the built-in features of Cantherm.

The fit is output in CHEMKIN format, and if pylab is available a plot of the fit is also generated.
"""
import math
import numpy
import logging
import csv
import re
import os
import os.path

from rmgpy.kinetics import Chebyshev

###########Provide inputs#################################################################################
#Reaction name using Chemkin format (be sure to include (+M) on both sides of reaction)
reaction_string = 'OH + CO (+M) <=> Products (+M)'

#Order of reaction (determines units)
order = 2

#Specify location of csv file with table of T, P and k data
#Note that the units must be correctly defined such that:
#T is in K
#P is in atm
#k is in cm^3/mol*s

k_data_path = '2006_Joshi_OH_CO_k1_table_example.csv'

#Details of Chebyshev fit
num_Cheb_coeff_1 = 6
num_Cheb_coeff_2= 4

###########Read table of k data (must be in mol, cm, s, atm, K units)########################################
f = csv.reader(open(k_data_path, 'r'))

columns = zip(*f)

kdata = numpy.zeros((len(columns[0][:])-1,len(columns)-1,), numpy.float64)

Tdata = numpy.array(columns[0][1:],dtype=numpy.float64)

Pdata = numpy.zeros(len(columns)-1, numpy.float64)

for col_num, col in enumerate(columns[1:]):
    Pdata[col_num] = col[0]
    for row_num, row in enumerate(col[1:]):
	kdata[row_num,col_num] = numpy.array(row,dtype=numpy.float64)

########################Fit Chebyshev polynomial to kdata###################################################
kunits = {1: 's^-1', 2: 'cm^3/(mol*s)', 3: 'cm^6/(mol^2*s)'}[order]

kinetics = Chebyshev().fitToData(Tdata, Pdata*101325., kdata, kunits, num_Cheb_coeff_1, num_Cheb_coeff_2, Tdata[0], Tdata[-1], Pdata[0]*101325., Pdata[-1]*101325.)

########################Output Chebyshev fit to CHEMKIN input file format#####################################
# Initialize (and clear!) the output files for the job
outputDirectory = os.path.dirname(os.path.abspath(k_data_path))

chemkinFile = os.path.join(outputDirectory, reaction_string.replace(" ", "_") + '.inp')
with open(chemkinFile, 'w') as f:
    pass
        
Tcount = Tdata.shape[0]
Pcount = Pdata.shape[0]

f = open(chemkinFile, 'a')
 
string = '{0!s:51} 1.0 0.0 0.0\n'.format(reaction_string)
       
coeffs = kinetics.coeffs.value_si.copy()
coeffs[0,0] += 6 * (order - 1)
string += 'TCHEB/ {0:<9.3f} {1:<9.3f}/\n'.format(Tdata[0], Tdata[-1])
string += 'PCHEB/ {0:<9.3f} {1:<9.3f}/\n'.format(Pdata[0], Pdata[-1])
string += 'CHEB/ {0:d} {1:d}/\n'.format(num_Cheb_coeff_1, num_Cheb_coeff_2)
if num_Cheb_coeff_2 < 6:
    for i in range(num_Cheb_coeff_1):
        string += 'CHEB/'
        for j in range(num_Cheb_coeff_2):
            string += ' {0:<12.3e}'.format(coeffs[i,j])
        string += '/\n'
else:
    coeffs_list = []
    for i in range(num_Cheb_coeff_1):
        for j in range(num_Cheb_coeff_2):
            coeffs_list.append(coeffs[i,j])
	coeffs_list[0] += 6 * (order - 1)
        for i in range(len(coeffs_list)):
            if i % 5 == 0: string += '    CHEB/'
            string += ' {0:<12.3e}'.format(coeffs_list[i])
            if i % 5 == 4: string += '/\n'  

f.write('{0}\n'.format(string))
            
f.close()

########################Evaluate Chebyshev Fit at T,P of kdata##############################################
kfit = numpy.zeros((Tcount, Pcount))
for t in range(Tcount):
    for p in range(Pcount):
        kfit[t,p] = kinetics.getRateCoefficient(Tdata[t], Pdata[p]*101325.)

########################Plot comparison of kdata and kfit if pylab is available##############################
try:
    import matplotlib.pyplot as pylab
    import matplotlib.ticker
    import matplotlib.cm
    cm = matplotlib.cm.jet

    fig = pylab.figure(figsize=(10,6))

    kfit *= 1e6 ** (order-1)

    pylab.subplot(1,2,1)
    for p in range(Pcount):
        pylab.semilogy(1000.0 / Tdata, kdata[:,p], color=cm(1.*p/(Pcount-1)), marker='o', linestyle='')
        pylab.semilogy(1000.0 / Tdata, kfit[:,p], color=cm(1.*p/(Pcount-1)), marker='', linestyle='-')
    pylab.xlabel('1000 / Temperature (1000/K)')
    pylab.ylabel('Rate coefficient ({0})'.format(kunits))
    pylab.title(reaction_string)
                
    pylab.subplot(1,2,2)
    for t in range(Tcount):
        pylab.loglog(Pdata, kdata[t,:], color=cm(1.*t/(Tcount-1)), marker='o', linestyle='')
        pylab.loglog(Pdata, kfit[t,:], color=cm(1.*t/(Tcount-1)), marker='', linestyle='-')
    pylab.xlabel('Pressure (atm)')
    pylab.ylabel('Rate coefficient ({0})'.format(kunits))
    pylab.title(reaction_string)
                
    fig.subplots_adjust(left=0.10, bottom=0.13, right=0.95, top=0.92, wspace=0.3, hspace=0.3)

    pylab.savefig(os.path.join(outputDirectory, reaction_string.replace(" ", "_") + '_Chebyshev_fit.pdf'))
    pylab.close()

except ImportError:
    print('Pylab cannot be imported, so fit plots were not generated')