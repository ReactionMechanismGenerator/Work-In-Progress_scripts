#!/usr/bin/env python
# encoding: utf-8

import math
import numpy
import matplotlib
import matplotlib.ticker
import logging

import rmgpy.constants as constants
#from rmgpy.cantherm.main import *

from rmgpy.statmech.conformer import Conformer
from rmgpy.cantherm.statmech import projectRotors

"""
Return the optimum geometry of the molecular configuration from the
Gaussian log file. If multiple such geometries are identified, only the
last is returned.
"""


# Geometry Input: [atomic number] [x] [y] [z]
geometry = """
6 1.182654 1.317879 0.000000
6 -0.152787 -0.748354 0.000000
6 0.000000 0.705075 0.000000
6 -1.295173 -1.390475 0.000000
1 2.106737 0.755782 0.000000
1 1.258597 2.394457 0.000000
1 0.772822 -1.328363 0.000000
1 -0.913461 1.285149 0.000000
1 -1.632862 -2.411773 0.000000
"""

number = []; coord = []
Natoms = 0

lines = geometry.strip().split('\n')
for line in lines:
    data = line.split()
    number.append(int(data[0]))
    coord.append([float(data[1]), float(data[2]), float(data[3])])
    # Automatically determine the number of atoms
    Natoms += 1


coord = numpy.array(coord, numpy.float64)
number = numpy.array(number, numpy.int)
mass = numpy.zeros(len(number), numpy.float64)
# Use the atomic mass of the most common isotope rather than the
# average atomic mass
# These values were taken from "Atomic Weights and Isotopic Compositions" v3.0 (July 2010) from NIST
for i in range(len(number)):
    if number[i] == 1:
        mass[i] = 1.00782503207
    elif number[i] == 6:
        mass[i] = 12.0
    elif number[i] == 7:
        mass[i] = 14.0030740048
    elif number[i] == 8:
        mass[i] = 15.99491461956
    elif number[i] == 15:
        mass[i] = 30.97376163
    elif number[i] == 16:
        mass[i] = 31.97207100
    elif number[i] == 17:
        mass[i] = 35.4527
    else:
        print 'Atomic number {0:d} not yet supported in loadGeometry().'.format(number[i])

conformer = Conformer(number=number, mass = (mass,"amu"), coordinates = (coord,"angstroms") )


"""
Define rotors
"""

pivots=[2,3]

top=[2,4,7,9]

symmetry=1

linear=False

TS=False

rotors = [
        ('fake.log', pivots, top, symmetry, 'best'),
    ]

"""
Calculate reduced moment of inertia for hindered rotor
"""

reduced_moment_of_inertia = conformer.getInternalReducedMomentOfInertia(pivots,top) * constants.Na * 1e23 #Units are amu*angstrom^2

print reduced_moment_of_inertia
