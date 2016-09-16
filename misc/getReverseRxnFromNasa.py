'''
Created on Dec 2, 2014

This gets reverse Rxn kinetics from NASA polynomials.

Be careful! right now it only uses the the low temperature junction from the NASA polynomial

@author: Nathan Yee
'''
import re
import math
from decimal import Decimal

nasaPath="C:/Users/User1/Dropbox/Research/RMG/thermo/aromatics_ATCT.txt"
libraryPath="C:/Users/User1/Dropbox/Research/RMG/thermo/Burcat_library.txt"

nasaList=[]
libraryList=[]

R=1.9872041e-3
Trange=[300,400,500,600,800,1000,1500]
TWOPLACES = Decimal(10) ** -2

r1="""
CH2O(17)                C 1  H 2  O 1       G298.000   3000.000  1000.00       1
 1.21253087E+00 1.03100703E-02-5.32426860E-06 1.32521523E-09-1.29097877E-13    2
-1.39186182E+04 1.66174479E+01 4.02550517E+00-4.06770880E-03 2.17177398E-05    3
-2.08782820E-08 6.59719571E-12-1.43510901E+04 3.76247699E+00                   4"""

r2="""
C2H5(39)                C 2  H 5            G298.000   3000.000  1000.00       1
 1.39424662E+00 1.85742180E-02-9.00112446E-06 2.12877239E-09-1.99100582E-13    2
 1.34011406E+04 1.66581828E+01 2.01543665E+00 1.36856127E-02 3.09981682E-06    3
-9.89324355E-09 3.98938932E-12 1.33584121E+04 1.41655049E+01                   4"""

p1="""
C3H7O(76)               C 3  H 7  O 1       G298.000   3000.000  1000.00       1
 3.71165592E+00 2.98888018E-02-1.50650214E-05 3.68054175E-09-3.53437601E-13    2
-6.55169654E+03 8.27016784E+00 1.28341389E+00 3.57266194E-02-1.76969950E-05    3
 1.14202428E-09 1.40747781E-12-5.88259252E+03 2.09279823E+01                   4"""

O2="""
O2(2)                   O 2                 G100.000   5000.000  1074.55       1
 3.15382081E+00 1.67804370E-03-7.69974231E-07 1.51275461E-10-1.08782413E-14    2
-1.04081728E+03 6.16755829E+00 3.53732243E+00-1.21571646E-03 5.31620250E-06    3
-4.89446430E-09 1.45846256E-12-1.03858849E+03 4.68368183E+00                   4"""

meth="""
CH4O(1)                 C 1  H 4  O 1       G100.000   5000.000  995.37        1
 3.01638671E+00 1.12129856E-02-4.28041856E-06 7.85995770E-10-5.52590086E-14    2
-2.58401877E+04 7.97393606E+00 3.91002253E+00-1.31507446E-03 2.80668730E-05    3
-2.98995317E-08 9.91740210E-12-2.55753744E+04 5.89067671E+00                   4"""

H2O="""
H2O(7)                  H 2  O 1            G100.000   5000.000  1130.24       1
 2.84325037E+00 2.75108559E-03-7.81031566E-07 1.07243658E-10-5.79392397E-15    2
-2.99586127E+04 5.91042036E+00 4.05763585E+00-7.87936058E-04 2.90877559E-06    3
-1.47518940E-09 2.12843251E-13-3.02815866E+04-3.11364125E-01                   4"""

CO="""
CO(27)                  C 1  O 1            G100.000   5000.000  1669.93       1
 2.92794119E+00 1.81932724E-03-8.35318755E-07 1.51271329E-10-9.88888522E-15    2
-1.42826359E+04 6.50662559E+00 3.59709701E+00-1.02424113E-03 2.83337188E-06    3
-1.75825933E-09 3.42589052E-13-1.43331255E+04 3.45318822E+00                   4"""

HO2="""
HO2(8)                  H 1  O 2            G100.000   5000.000  932.08        1
 4.05113698E+00 2.15394221E-03-5.68214775E-07 9.60359750E-11-7.17412654E-15    2
 1.09371051E+02 3.58811323E+00 4.03106589E+00-2.55637542E-03 1.47311108E-05    3
-1.63677597E-08 5.88954024E-12 3.21462891E+02 4.80119907E+00                   4"""

OH="""
OH(5)                   H 1  O 1            G100.000   5000.000  1145.75       1
 3.07194127E+00 6.04013329E-04-1.39769778E-08-2.13449117E-11 2.48068109E-15    2
 3.57938611E+03 4.57799176E+00 3.51456784E+00 2.92796151E-05-5.32170670E-07    3
 1.01949890E-09-3.85948479E-13 3.41425421E+03 2.10434956E+00                   4"""



#Gets low temperature G in kcal/mol
def getH(nasa,T):
    nasaList=re.split('\n', nasa)
    cleanNasaList=[]
    for line in nasaList:
        if line=="": pass
        else: cleanNasaList.append(line)
    
    #get temperature range
    TempRanges=re.split('\s+',cleanNasaList[0])
    TempRanges=TempRanges[0:-1]
    while len(TempRanges)>2:
        TempRanges=TempRanges[1:]
    TempRanges=[float(x) for x in TempRanges]
    
    #Get nasa polynomials
    nasapoly=[]
    for x in range(1,4):
        reLine=re.findall('\-?[0-9].*?E[+-][0-9][0-9]', cleanNasaList[x])
        nasapoly.extend(reLine)
    if not len(nasapoly) == 14:
        print "did not parse as 14 nasa polynomials"
        return
    
    nasapolyNew=[]
    for num in nasapoly:
        nasapolyNew.append(float(num))
    nasapoly=nasapolyNew
    
    if T>TempRanges[0]: c=7
    else: c=0
    H=nasapoly[7-c]+nasapoly[8-c]/2*T+(nasapoly[9-c]*T**2)/3+(nasapoly[10-c]*T**3)/4+(nasapoly[11-c]*T**4)/5+nasapoly[12-c]/T
    H=H*T*R
    return H

def getG(nasa, T):
    nasaList=re.split('\n', nasa)
    cleanNasaList=[]
    for line in nasaList:
        if line=="": pass
        else: cleanNasaList.append(line)
    
    #get temperature range
    TempRanges=re.split('\s+',cleanNasaList[0])
    TempRanges=TempRanges[0:-1]
    while len(TempRanges)>2:
        TempRanges=TempRanges[1:]
    TempRanges=[float(x) for x in TempRanges]
    
    #Get nasa polynomials
    nasapoly=[]
    for x in range(1,4):
        reLine=re.findall('\-?[0-9].*?E[+-][0-9][0-9]', cleanNasaList[x])
        nasapoly.extend(reLine)
    if not len(nasapoly) == 14:
        print "did not parse as 14 nasa polynomials"
        return
    
    nasapolyNew=[]
    for num in nasapoly:
        nasapolyNew.append(float(num))
    nasapoly=nasapolyNew
    
    if T>TempRanges[0]: c=7
    else: c=0
    H=nasapoly[7-c]+nasapoly[8-c]/2*T+(nasapoly[9-c]*T**2)/3+(nasapoly[10-c]*T**3)/4+(nasapoly[11-c]*T**4)/5+nasapoly[12-c]/T
    H=H*T*R
    S=nasapoly[7-c]*math.log(T)+nasapoly[8-c]*T+(nasapoly[9-c]*T**2)/2+(nasapoly[10-c]*T**3)/3+(nasapoly[11-c]*T**4)/4+nasapoly[13-c]
    S=S*R
    
    G=H-T*S
                
        
    return G

#Put in the NASA's for the forward reaction. If not that many products, put in an empty string
def getReverseRxnRate(kf, T, rNasa1, pNasa1, rNasa2='', pNasa2=''):
    order=0
    rG1=getG(rNasa1,T)
    pG1=getG(pNasa1,T)
    if rNasa2=='': rG2=0
    else: 
        rG2=getG(rNasa2,T)
        order=order-1
    if pNasa2=='': pG2=0
    else: 
        pG2=getG(pNasa2,T)
        order+=1
    
    deltaG=pG1+pG2-rG1-rG2
    
    K=math.exp(-deltaG/R/T)*1/(83.144598*T)**order
    
    return kf/K

def print_listofNASA(nasa):
    nasaList=re.split('\n', nasa)
    cleanNasaList=[]
    for line in nasaList:
        if line=="": pass
        else: cleanNasaList.append(line)
    
    #get temperature range
    TempRanges=re.split('\s+',cleanNasaList[0])
    TempRanges=TempRanges[0:-1]
    while len(TempRanges)>2:
        TempRanges=TempRanges[1:]
    TempRanges=[float(x) for x in TempRanges]
    
    #Get nasa polynomials
    nasapoly=[]
    for x in range(1,4):
        reLine=re.findall('\-?[0-9].*?E[+-][0-9][0-9]', cleanNasaList[x])
        nasapoly.extend(reLine)
    if not len(nasapoly) == 14:
        print "did not parse as 14 nasa polynomials"
        return
    
    nasapolyNew=[]
    for num in nasapoly:
        nasapolyNew.append(float(num))
    nasapoly=nasapolyNew
    
    print nasapoly

if __name__ == '__main__':
    Tlist=[x for x in xrange(600,1000,10)]
    # kforwardList=[1.0e11*T**0*math.exp(-3.5/T/R) for T in Tlist]
    #
    # for index, T in enumerate(Tlist):
    #     kback=getReverseRxnRate(kforwardList[index], T, r1, p1, r2)
    #     print T, kback, math.log10(kback)

    for index, T in enumerate(Tlist):
        meth_H=getH(meth,T)
        O2_H=getH(O2,T)
        H2O_H=getH(H2O, T)
        CO_H=getH(CO, T)
        HO2_H=getH(HO2, T)
        OH_H=getH(OH, T)

        deltaH=H2O_H+CO_H+OH_H+HO2_H-meth_H-2*O2_H
        print T, deltaH