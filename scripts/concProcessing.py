import math
import re
import csv
import matplotlib.pyplot as plt

#change information here

#creating filenames to analyze

filename='/home/mark/Research/Alcohols/Methanol/ChemKin/RMGupdated/750K_stoic.csv'
fuelname='CH4O'

T=750

colors=((0,0,0),(1,0,0), (0,1,0),(0,0,1),(1,0,1),(0,1,1),
        (.5,0,0),(0,.5,0),(0,0,.5),(0,.5,.5),(.5,.5,0),(.5,0,.5),
        (1,.5,0),(1,0,.5),(.5,1,0),(0,1,.5),(.5,0,1),(0,.5,1),
        (0,0,0),(1,0,0), (0,1,0),(0,0,1),(1,0,1),(0,1,1),
        (.5,0,0),(0,.5,0),(0,0,.5),(0,.5,.5),(.5,.5,0),(.5,0,.5),
        (1,.5,0),(1,0,.5),(.5,1,0),(0,1,.5),(.5,0,1),(0,.5,1),
        (0,0,0),(1,0,0), (0,1,0),(0,0,1),(1,0,1),(0,1,1),
        (.5,0,0),(0,.5,0),(0,0,.5),(0,.5,.5),(.5,.5,0),(.5,0,.5),
        (1,.5,0),(1,0,.5),(.5,1,0),(0,1,.5),(.5,0,1),(0,.5,1))
#filepath='/home/mark/Research/ChemKin/Methanol/updated model/'
#filenames=['675kstoic.csv', '700kstoic.csv','725kstoic.csv',
#           '750kstoic.csv','800kstoic.csv','850kstoic.csv',
#          '950kstoic.csv']
#T=[675, 700,725,750,800,850,950]#kelvin
#plotting info
#plt.axis([0,15,-20,1]) #range of plot. first two are time, second two 
                        #are log concentration. comment out to plot everything





#constants used in function
#from chemkin file [prefactor(moles/s) n kcal/mol
rxn1rates=[2.823e4, 2.483, 9.514] #ho2 plus fuel
rxn2rates=[4.2e14, 0, 11.982,1.3e11,0,-1.629] #ho2 squared - chemkin#22&23
Rkcal=1.98720 #kcal/molK
Rgas=0.082056 #L atm/molK
conckey='Mole fraction' #finds mole fraction concentrations





#read csv file
reactionIndex={} #column number of each component
with open(filename, 'rb') as csvfile:
    reader=csv.reader(csvfile)
    #find column number of each component (reactionIndex)
    for line in reader:
        for index, cell in enumerate(line):
            if re.search('Pressure', cell):
                pressureColumn=index
            if re.search('Time', cell):
                timeColumn=index
            if re.search(conckey, cell)and (re.search('H2O2',cell) or re.search('HO2',cell) or re.search('CH2O',cell)):
                print cell
                reactionIndex[cell]=index #cell is text
        break
    
    timeList=[]     #stores time values
    pressureList=[] #stores pressure values
    ConcDict={}   #stores concentration values 
    for key in reactionIndex:   #key is the string of text
        ConcDict[key]=[]
    
    
    #stores data from the file and converts mole fractions to Mol/L concentration
    for line in reader:
        timeList.append(float(line[timeColumn]))
        pressureList.append(float(line[pressureColumn]))
        for key in reactionIndex:
            ConcDict[key].append(float(line[reactionIndex[key]])*pressureList[-1]/Rgas/T)
#end stuff with open reader
    
    


#find fuel concentration.
for key in ConcDict:
    if re.search(fuelname,key):
        fuelconc=ConcDict[key][0]
        
#plot data
print len(ConcDict)
index=-1
for key in ConcDict:
    index+=1
    plt.plot(timeList[1:], ConcDict[key][1:], color=colors[index],label=re.sub('Mole fraction ', '', re.sub('Soln.*', '', key)))

#find characteristic HO2 concentration
#HO2knee=
plt.yscale('log')
# plt.ylim([1e-16,4e-14])
# plt.xlim([0,.3])
plt.legend()
plt.xlabel(r'time $(s)$')
plt.ylabel(r'concentration $(\frac{mol}{L})$')
plt.show()
        