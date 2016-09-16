"""
This script creates a combinatorial list of simple reactors to copy and paste into an inpute file
"""


outPath = "/Users/Nate/Dropbox (MIT)/Research/butanolSurfaces/RMG models/reactors.txt"
temperature = range(700,2300,200) #K
pressure = [1,5,10] #bar
#key is label and value is list of concentrations.
# Every list must have equal length, as each list element corresponds to a specific condition
concentrations={
    "iBuOH": [0.0654, 0.0338, 0.0172, 0.00867],
    "O2": [0.196, 0.203, 0.206, 0.208],
    "N2": [0.738, 0.763, 0.776, 0.783]
}

#These are not combinatorially added, they apply to every simple reactor
terminationConversion = {"iBuOH": 0.99}
terminationTime = 3600 #s

textList=[]

concKeys = concentrations.keys() #ensures that the order is preserved, although the one ordering is still random
for tx in temperature:
    for px in pressure:
        for mixIndex in range(len(concentrations.values()[0])):
            textList.append("simpleReactor(")
            textList.append("    temperature=({0}),'K'),".format(tx))
            textList.append("    pressure=(10.0,'bar'),".format(px))
            textList.append("    initialMoleFractions={")
            for cx in concKeys:
                textList.append("        '{0}': {1},".format(cx, concentrations[cx][mixIndex]))
            textList.append("    },")
            textList.append("    terminationConversion={")
            for termx in terminationConversion:
                textList.append("        '{0}': {1},".format(termx, terminationConversion[termx]))
            textList.append("    },")
            textList.append("    terminationTime=({0},'s'),".format(terminationTime))
            textList.append(")")
            textList.append("")

#add new line to every
for index in range(len(textList)):
    textList[index] = textList[index]+"\n"

with open(outPath, 'wb') as output:
    for line in textList:
        output.write(line)
