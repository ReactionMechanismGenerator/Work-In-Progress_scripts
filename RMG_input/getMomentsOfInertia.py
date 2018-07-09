#!/usr/bin/env python
# encoding: utf-8

import numpy

masses = numpy.array([1,1,1])
coordinates = numpy.array([[1,0,0],[0,1,0],[0,0,1]])

def getMomentsOfInertia(masses,coordinates):
    """
    moment of inertia about the center of mass in the coordinate system defined in 
    coordinates
    coordinates and masses are numpy arrays
    
    C: 12.0107 amu, H: 1.0079 amu, O: 15.9994 amu
    """
    assert coordinates.shape[1] == 3, 'orientation of coordinates matrix is wrong try with transpose: arr.T'
    Ixx = 0.0
    Iyy = 0.0
    Izz = 0.0
    cm = 0.0
    
    tmass = sum(masses)
    for i,c in enumerate(coordinates):
        cm += masses[i]*c
    cm = cm/tmass

    for i,m in enumerate(masses):
        coord = coordinates[i,:]
        Ixx += m*(coord[0]-cm[0])**2
    for i,m in enumerate(masses):
        coord = coordinates[i,:]
        Iyy += m*(coord[1]-cm[1])**2
    for i,m in enumerate(masses):
        coord = coordinates[i,:]
        Izz += m*(coord[2]-cm[2])**2
    return [Ixx,Iyy,Izz]

print getMomentsOfInertia(masses,coordinates)


