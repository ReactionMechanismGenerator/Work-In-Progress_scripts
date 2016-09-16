import re

#Fill out these inputs
inpath = '/Users/Nate/Desktop/chem_annotated.inp' #path to chemkin file
outpath = '/Users/Nate/Desktop/chem_annotated_new.inp' #path to new chemkin file
speciesName = 'CH3(17)' #Name of species with parenthesis

#Import as list into python
textList=[]
with open(inpath, 'rb') as inputFile:
    for line in inputFile:
        textList.append(line)

newTextList= []

pastThermo = False
deleteMode = False
for line in textList:
    #This block deletes reactions or thermo based on the fact that every reaction/thermo block is separated by blank lines
    if pastThermo:
        if deleteMode:
            #This means we reached the end of the thermo block for this species
            if line.strip() == '':
                #stop watching to delete
                deleteMode = False
                #count how many lines to delete
                for index, reverseLine in enumerate(reversed(newTextList)):
                    if reverseLine.strip() == '':
                        break
                #delete the lines
                newTextList=newTextList[0:len(newTextList)-index]
            else: newTextList.append(line)
        #ignore lines that are comments
        elif len(line.strip()) == 0 or line.strip()[0]=='!': newTextList.append(line)
        elif re.search(re.escape(speciesName), line):
            #start watching to delete
            deleteMode = True
        else: newTextList.append(line)
    #This block is for skipping the species declaration, each species should only occur in one line
    else:
        if re.search('THERM', line):
            newTextList.append(line)
            pastThermo = True
            continue
        #ignore lines that are comments
        elif len(line.strip()) == 0 or line.strip()[0]=='!': newTextList.append(line)
        elif re.search(re.escape(speciesName), line): pass
        else: newTextList.append(line)

#rewrite out chemkin list
with open(outpath, 'wb') as outFile:
    for line in newTextList:
        outFile.write(line)