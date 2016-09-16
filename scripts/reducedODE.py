import math
import re
import csv
import numpy as np
import matplotlib.pyplot as plt
from getRateCoefficients import getRateCoefficients
from getHeatsReaction import getHeatsReaction
from calcconc import calcconc
from scipy import integrate as int

def f(t,y,Tinit,Pinit,stoic):
    '''
    y is vector containing a T, HO2,H2O2,formal
    Tinit is initial temperature (K), pinit is initial pressure (atm)
    stoic is the stoiciometric ratio of the system >1 is fuel rich
    '''
    T=y[0]
    ho2=y[1]
    h2o2=y[2]
    formal=y[3]
    Pinit=Pinit*101325#pascal conversion
    R=8.314#m3Pa/molK
    molardensity=Pinit/R/Tinit#initial density in mol/m3
    
    o2init=1/(5+1.5*stoic)*molardensity
    fuelinit=1.5*stoic/(5+1.5*stoic)*molardensity
    n2init=4/(5+1.5*stoic)*molardensity
#     fuelinit=20.4820165855#mol/m3 at 700k 10atm
#     o2init=30.7230161734#mol/m3 at 700k 10atm
    Cv=1.25#poor estimate in (J/gK) of heat capacity of stoiciometric conditions at 700K
    M=(28.*n2init+32.*o2init+32.*fuelinit)/molardensity#g/mol molecular mass of stoiciometric mixture (may write explicitly later)
    massdensity=M*molardensity #g/m3
    ratecoefficient=getRateCoefficients(
           '/home/mark/Research/Alcohols/Methanol/ChemKin/RMGupdated/meth_core_v3.inp',
           '/home/mark/Research/Alcohols/Methanol/ChemKin/RMGupdated/species_dictionary.txt',
           [1,27,26,22,55],
           [True,True,True,True,True],
           [T],
           10*101325,
           {'N2': 12./17, 'O2':3./17,'methanol':2./17})
    heats=getHeatsReaction(
            '/home/mark/Research/Alcohols/Methanol/ChemKin/RMGupdated/meth_core_v3.inp',
            '/home/mark/Research/Alcohols/Methanol/ChemKin/RMGupdated/species_dictionary.txt',
            [4,70],#4 is OH+fuel and 70 is O2 +HCO
            [True,True],
            [T])
    r=[ratecoefficient[0][0]*o2init*fuelinit,
       ratecoefficient[1][0]*ho2*fuelinit,
       ratecoefficient[2][0]*h2o2,
       ratecoefficient[3][0]*ho2**2,
       ratecoefficient[4][0]*ho2*formal]
    return [(-2*r[2]*heats[0][0]-r[4]*heats[1][0])/massdensity/Cv,
            2*r[0]+2*r[2]-2*r[3],
            r[1]-r[2]+r[3]+r[4],
            r[0]+r[1]+2*r[2]-r[4]]
# initial conditions
initial_temp=700
initial_pres=10#atm
stoic_ratio=1
y0=[initial_temp,0.,0,0]
stepsize=2
tfinal=6
# solve ode
solver=int.ode(f).set_integrator('vode',method='bdf',order=3)
solver.set_initial_value(y0, 0).set_f_params(initial_temp,initial_pres,stoic_ratio)
tval=np.array((0), dtype=float)
yval=np.array(y0,dtype=float)
print yval.shape
while solver.successful() and solver.t < tfinal:
    solver.integrate(solver.t+stepsize)
    tval=np.append(tval, solver.t)
    print yval.shape
    yval=np.append(yval,solver.y,axis=1)
    print solver.t
print yval

print tval

# T=[:, 0]
# ho2=[:, 1]
# h2o2=[:, 2]
# formal=[:, 3]
# 
# plt.figure()
# plt.plot(time_grid,ho2,)
