from rmgpy.thermo.wilhoit import Wilhoit
from rmgpy.thermo.thermodata import ThermoData
import rmgpy.constants as constants
import numpy
from rmgpy.chemkin import writeThermoEntry
from rmgpy.molecule.element import getElement

SpeciesIdentifier = 'C6H5I'

elementCounts = {'C': 6, 'H': 5, 'I': 1}

Tlist = numpy.array([300.0,400.0,500.0,600.0,800.0,1000.0])
Cplist = numpy.array([24.2, 31.1, 36.9, 41.4, 48, 52.55])*4.184 #J/mol*K
H298 = 40.5*4184 #J/mol
S298 = 81.35*4.184 #J/mol*K

# Polyatomic species
linear = False
Nfreq = 30
Nrotors = 0
Cp0 = (3.5 if linear else 4.0) * constants.R
CpInf = Cp0 + (Nfreq + 0.5 * Nrotors) * constants.R

wilhoit = Wilhoit()

wilhoit.fitToData(Tlist, Cplist, Cp0, CpInf, H298, S298, B0=500.0)\

NASA_fit = wilhoit.toNASA(Tlist[0], Tlist[-1], Tint=500.0)

string = ''

f = open('chem.inp', 'a')

# Line 1
string += '{0:<16}        '.format(SpeciesIdentifier)
if len(elementCounts) <= 4:
    # Use the original Chemkin syntax for the element counts
    for key, count in elementCounts.iteritems():
        if isinstance(key, tuple):
            symbol, isotope = key
            chemkinName = getElement(symbol, isotope=isotope).chemkinName
        else:
            chemkinName = key
        string += '{0!s:<2}{1:>3d}'.format(chemkinName, count)
    string += '     ' * (4 - len(elementCounts))
else:
    string += '     ' * 4
string += 'G{0:>10.3f}{1:>10.3f}{2:>8.2f}      1'.format(NASA_fit.polynomials[0].Tmin.value_si,
                                                         NASA_fit.polynomials[1].Tmax.value_si,
                                                         NASA_fit.polynomials[0].Tmax.value_si)
if len(elementCounts) > 4:
    string += '&\n'
    # Use the new-style Chemkin syntax for the element counts
    # This will only be recognized by Chemkin 4 or later
    for key, count in elementCounts.iteritems():
        if isinstance(key, tuple):
            symbol, isotope = key
            chemkinName = getElement(symbol, isotope=isotope).chemkinName
        else:
            chemkinName = key
        string += '{0!s:<2}{1:>3d}'.format(chemkinName, count)
string += '\n'

# Line 2
string += '{0:< 15.8E}{1:< 15.8E}{2:< 15.8E}{3:< 15.8E}{4:< 15.8E}    2\n'.format(NASA_fit.polynomials[1].c0,
                                                                                  NASA_fit.polynomials[1].c1,
                                                                                  NASA_fit.polynomials[1].c2,
                                                                                  NASA_fit.polynomials[1].c3,
                                                                                  NASA_fit.polynomials[1].c4)

# Line 3
string += '{0:< 15.8E}{1:< 15.8E}{2:< 15.8E}{3:< 15.8E}{4:< 15.8E}    3\n'.format(NASA_fit.polynomials[1].c5,
                                                                                  NASA_fit.polynomials[1].c6,
                                                                                  NASA_fit.polynomials[0].c0,
                                                                                  NASA_fit.polynomials[0].c1,
                                                                                  NASA_fit.polynomials[0].c2)

# Line 4
string += '{0:< 15.8E}{1:< 15.8E}{2:< 15.8E}{3:< 15.8E}                   4\n'.format(NASA_fit.polynomials[0].c3,
                                                                                      NASA_fit.polynomials[0].c4,
                                                                                      NASA_fit.polynomials[0].c5,
                                                                                      NASA_fit.polynomials[0].c6)

f.write('{0}\n'.format(string))
f.close()