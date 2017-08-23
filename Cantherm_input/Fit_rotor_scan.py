#!/usr/bin/env python
# encoding: utf-8

"""
Fit a hindered rotor scan (in 2-column text file format) to both a Fourier series and a single cosine function. 
Outputs the best fit.

Run using the  following commands:

python Fit_rotor_scan.py "path_to_scan_log.txt" symmetry_of_rotor

"""

import argparse
import os
import rmgpy.cantherm.statmech
from rmgpy.cantherm.statmech import HinderedRotor, ScanLog
from rmgpy import settings
from rmgpy.cantherm.output import prettify
import numpy
import rmgpy.constants as constants

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()

    parser.add_argument('scanLogPath', metavar='SCAN', type=str, nargs=1,
        help='The path of scanlog text file')

    parser.add_argument('symmetry',type=int, nargs=1,
        help='Symmetry of rotor scan')

    args = parser.parse_args()
    scanLogPath = args.scanLogPath[0]
    symmetry = args.symmetry[0]

    directory = os.path.abspath(os.path.dirname(scanLogPath))

    #Load the scanlog file
    scanlog=ScanLog(scanLogPath)
    angle, Vlist = scanlog.load()

    #Initialize HinderedRotor class
    fourierRotor = HinderedRotor(inertia=(1, "amu*angstrom^2"), symmetry=symmetry)
    fourierRotor.fitFourierPotentialToData(angle, Vlist)

    #Fit scan to both cosine and fourier series
    cosineRotor = HinderedRotor(inertia=(1, "amu*angstrom^2"), symmetry=symmetry)
    cosineRotor.fitCosinePotentialToData(angle, Vlist)
    fourierRotor = HinderedRotor(inertia=(1, "amu*angstrom^2"), symmetry=symmetry)
    fourierRotor.fitFourierPotentialToData(angle, Vlist)

    #Evaluate fits
    Vlist_cosine = numpy.zeros_like(angle)
    Vlist_fourier = numpy.zeros_like(angle)
    for i in range(angle.shape[0]):
        Vlist_cosine[i] = cosineRotor.getPotential(angle[i])
        Vlist_fourier[i] = fourierRotor.getPotential(angle[i])

    #Determine which fit is better
    rms_cosine = numpy.sqrt(numpy.sum((Vlist_cosine - Vlist) * (Vlist_cosine - Vlist)) / (len(Vlist) - 1)) / 4184.
    rms_fourier = numpy.sqrt(numpy.sum((Vlist_fourier - Vlist) * (Vlist_fourier - Vlist)) / (len(Vlist) - 1)) / 4184.

    # Keep the rotor with the most accurate potential
    rotor = cosineRotor if rms_cosine < rms_fourier else fourierRotor
    # However, keep the cosine rotor if it is accurate enough, the
    # fourier rotor is not significantly more accurate, and the cosine
    # rotor has the correct symmetry
    if rms_cosine < 0.05 and rms_cosine / rms_fourier < 2.0 and rms_cosine / rms_fourier < 4.0 and symmetry == cosineRotor.symmetry:
        rotor = cosineRotor

    """
    Plot the potential for the rotor, along with its cosine and Fourier
    series potential fits. The plot is saved to a set of files of the form
    ``hindered_rotor_1.pdf``.
    """

    try:
        import pylab
    except ImportError:
        error

    fig = pylab.figure(figsize=(6, 5))
    pylab.plot(angle, Vlist / 4184., 'ok')
    linespec = '-r' if rotor is cosineRotor else '--r'
    pylab.plot(angle, Vlist_cosine / 4184., linespec)
    linespec = '-b' if rotor is fourierRotor else '--b'
    pylab.plot(angle, Vlist_fourier / 4184., linespec)
    pylab.legend(['scan', 'cosine', 'fourier'], loc=1)
    pylab.xlim(0, 2 * constants.pi)
    pylab.xlabel('Angle')
    pylab.ylabel('Potential (kcal/mol)')

    axes = fig.get_axes()[0]
    axes.set_xticks([float(j * constants.pi / 4) for j in range(0, 9)])
    axes.set_xticks([float(j * constants.pi / 8) for j in range(0, 17)], minor=True)
    axes.set_xticklabels(
        ['$0$', '$\pi/4$', '$\pi/2$', '$3\pi/4$', '$\pi$', '$5\pi/4$', '$3\pi/2$', '$7\pi/4$', '$2\pi$'])

    pylab.savefig(os.path.join(directory, 'rotor_fit.pdf'))
    pylab.close()

    print(rotor)