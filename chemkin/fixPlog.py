"""
Previously, we had a bug in RMG where we printed out chemkin.inp with PLOG reactions marked as duplicates.
See RMG-PY issue #147 for more information.

The correct way to represent duplicate PLOGs is with one reaction where all Plog lines (with the same pressure
more than once) all under this one reaction.

This script takes a chemkin file with the incorrect format and gives a new chemkin file with the correct
format.

WIP notes:
If there are comments surrounding duplicate PLOGS, they will not be handled correctly.
This script uses '=' signs to parse reactions. If there are '=' in comments, this script will break

Last updated by nyee on 9/18/2016
"""

import re

def getReactionBlock(textList, index):
    """

    :param textList: full list of lines from the chemkin file
    :param index: line that contains part of a reaction
    :return: list of all lines associated with the reaction

    The output will include any black lines after the reaction block
    """
    reactionBlock = [textList[index]]
    #get all lines above index until we find a reaction string
    newIndex = index
    while not re.search('\=', reactionBlock[0]) and newIndex >0:
        newIndex += -1
        reactionBlock.insert(0, textList[newIndex])

    #get all lines below index until we find a reaction string. At the end
    #of this while loop, we should have 2 lines with rection strings in them
    #of hit the end of the file
    numberOfReactions=0
    newIndex = index
    for line in reactionBlock:
        if re.search('\=', line): numberOfReactions+=1
    while numberOfReactions < 2 and newIndex < len(textList)-1:
        newIndex += 1
        reactionBlock.append(textList[newIndex])

        numberOfReactions=0
        for line in reactionBlock:
            if re.search('\=', line): numberOfReactions+=1

    #remove the last reaction
    if re.search('\=', reactionBlock[-1]):
        reactionBlock = reactionBlock[:-1]

    return reactionBlock

def readKineticsReaction(line):

    """
    Parse the first line of of a Chemkin reaction entry.
    Copied and modified method from rmgpy/chemkin.py
    """
    tokens = line.split()

    rmg = True
    try:
        float(tokens[-6])
    except (ValueError, IndexError):
        rmg = False
    AuncertaintyType = '+|-'
    if rmg:
        A = float(tokens[-6])
        n = float(tokens[-5])
        Ea = float(tokens[-4])
        try:
            dA = float(tokens[-3])
        except ValueError:
            AuncertaintyType = '*|/'
            dA = float(tokens[-3][1:])
        dn = float(tokens[-2])
        dEa = float(tokens[-1])
        reaction = ''.join(tokens[:-6])
    else:
        # A = Ffloat(tokens[-3])
        # n = Ffloat(tokens[-2])
        # Ea = Ffloat(tokens[-1])
        dA = 0.0
        dn = 0.0
        dEa = 0.0
        reaction = ''.join(tokens[:-3])
    thirdBody = False

    # Split the reaction equation into reactants and products
    reversible = True
    reactants, products = reaction.split('=')
    if '<=>' in reaction:
        reactants = reactants[:-1]
        products = products[1:]
    elif '=>' in reaction:
        products = products[1:]
        reversible = False
    if '(+M)' in reactants: reactants = reactants.replace('(+M)','')
    if '(+m)' in reactants: reactants = reactants.replace('(+m)','')
    if '(+M)' in products:  products = products.replace('(+M)','')
    if '(+m)' in products:  products = products.replace('(+m)','')

    reactantSet = []
    productSet = []
    for reactant in reactants.split('+'):
        reactant = reactant.strip()
        reactantSet.append(reactant)
    reactantSet = set(reactantSet)

    for product in products.split('+'):
        product = product.strip()
        productSet.append(product)
    productSet = set(productSet)

    return (reactantSet, productSet)

if __name__ == '__main__':
    inpath = '/Users/Nate/Dropbox (MIT)/Research/butanolSurfaces/RMG models/publishedModel/chem.inp' #path to chemkin file
    outpath = '/Users/Nate/Dropbox (MIT)/Research/butanolSurfaces/RMG models/publishedModel/published_fixedPlog.inp' #path to new chemkin file
    speciesName = 'CH3(17)' #Name of species with parenthesis

    #Import as list into python
    textList=[]
    with open(inpath, 'rb') as inputFile:
        for line in inputFile:
            textList.append(line)

    #collect all duplicate PLOG reactions
    reactionSectionIndex = 0
    duplicatePlogs = [] #list of duplicate plogs.
    for index, line in enumerate(textList):
        if reactionSectionIndex == 0:
            if re.match("REACTIONS", line):
                reactionSectionIndex =index
        else:
            if re.match("DUPLICATE", line):
                reactionBlock = getReactionBlock(textList, index)
                plog = False
                for reactionLine in reactionBlock:
                    if re.match("PLOG", reactionLine):
                        duplicatePlogs.append(reactionBlock)
                        break
    #merge and process the new duplicates
    mergedPlogs = []
    for index, reaction1 in enumerate(duplicatePlogs):
        for reaction2 in duplicatePlogs[index+1:]:
            (r1, p1) = readKineticsReaction(reaction1[0])
            (r2, p2) = readKineticsReaction(reaction2[0])
            if r1 == r2 and p1 == p2:
                #merge the two reactions
                print "found duplicate reactions starting with {0}".format(reaction1[0].strip())
                duplicateIndex = 0
                for reactionIndex, line in enumerate(reaction1):
                    if re.match("DUPLICATE", line):  duplicateIndex = reactionIndex
                merge = reaction1[:duplicateIndex]
                merge.extend(reaction2[1:])
                #remove the duplicate line
                merge = [mergeLine for mergeLine in merge if not re.match("DUPLICATE", mergeLine)]

                #plogs must be in ascending order for chemkin to read, so reorder the lines
                sortingDictionary = {}
                for mergeLine in merge:
                    #next line stupidly unreadable...might want to fix
                    if re.match("PLOG", mergeLine.strip()):
                        tokens = mergeLine.split()
                        for token in tokens:
                            if re.search('[0-9]', token): break
                        plogValue = float(re.sub("[^0-9.]", "", token))
                        sortingDictionary[mergeLine] = plogValue
                # merge = sorted(merge, key=float(re.sub("[^0-9.]", "", re.search("PLOG.*?/s").group(0).strip())))
                sortedPlogLines = sorted(sortingDictionary.keys(), key=lambda mergeLine: sortingDictionary[mergeLine])
                sortedMerge = []
                pastPlog = False
                for mergeLine in merge:
                    if re.match("PLOG", mergeLine.strip()):
                        if not pastPlog:
                            pastPlog = True
                            sortedMerge.extend(sortedPlogLines)
                    else: sortedMerge.append(mergeLine)
                mergedPlogs.append(sortedMerge)

    #now create newTextList, replacing any old duplicate plogs
    newTextList = []
    alreadyAdded = [] #These merged Plogs have already been added
    ignore = -1 # ignore up to end of reaction block
    for index, line in enumerate(textList):
        #don't do anything outside reaction section
        if index > ignore:
            if index > reactionSectionIndex and not re.match('\!', line.strip()):
                if re.search('\=', line):
                    # print index
                    reactionBlock = getReactionBlock(textList, index)
                    duplicate = False
                    plog = False
                    for reactionLine in reactionBlock:
                        if re.match("PLOG", reactionLine): plog = True
                        elif re.match("DUPLICATE", reactionLine): duplicate = True
                    #This is both a duplicate and a plog so probably will match something we have in mergedPlogs
                    if duplicate and plog:
                        #check and then extend
                        (r1, p1) = readKineticsReaction(line) #line we are on
                        for reaction2 in mergedPlogs:
                            (r2, p2) = readKineticsReaction(reaction2[0])
                            # we found an old duplicate plog matching a plog in merged
                            if r1 == r2 and p1 == p2:
                                if reaction2 not in alreadyAdded:
                                    newTextList.extend(reaction2)
                                    alreadyAdded.append(reaction2)
                                ignore = index + len(reactionBlock) - 1 #don't write any of the rest of the lines of the original
                                break
                        else: newTextList.append(line)
                    else: newTextList.append(line)
                else: newTextList.append(line)
            else: newTextList.append(line)

    # print "old text length", len(textList)
    # print "new text length", len(newTextList)
    with open(outpath, 'wb') as outFile:
        for line in newTextList:
            outFile.write(line)

    test ="""iBuOH=CH3+C3H7O-N2     1.000e+00 0.000     0.000
PLOG/0.0013 3.86E+123 -32.0 138678/
PLOG/0.0132 2.93E+112 -28.32 136273/
PLOG/0.1316 8.13E+104 -25.73 136287/
PLOG/1.0 7.15E+98 -23.67 136787.4/
PLOG/10.0 9.95E+24 -2.35 90619.7/
PLOG/100.0 5.40E+63 -13.26 120314.53/
DUPLICATE

iBuOH=CH3+C3H7O-N2     1.000e+00 0.000     0.000
PLOG/0.0013 1.33E+76 -18.54 107646.6/
PLOG/0.0132 2.47E+70 -16.65 106177/
PLOG/0.1316 8.83E+65 -15.18 105114.6/
PLOG/1.0 9.67E+56 -12.32 102256.9/
PLOG/10.0 -1.01E+61 -11.88 132857.36/
PLOG/100.0 2.71E+38 -6.52 95614.2/
DUPLICATE

!iBuOH(+M)=iC4H8+H2O(+M)                             1.000e+00 0.000     0.000
!TCHEB / 290.0 3000.0 /	PCHEB / 0.009869232667160128 98.69232667160128 /
!CHEB / 6	4 /
!CHEB / -1.2898000e+01  3.6717000e-01 -2.1813000e-02 -1.9835000e-03 /
!CHEB /  1.9638000e+01  6.7713000e-01 -3.5512000e-02 -4.1935000e-03 /
!CHEB /  5.7880000e-01  5.2604000e-01 -1.5684000e-02 -4.3274000e-03 /
!CHEB / -1.1510000e-01  3.3238000e-01  5.4661000e-03 -3.3240000e-03 /
!CHEB / -1.8781000e-01  1.5367000e-01  1.8605000e-02 -8.4913000e-04 /
!CHEB / -1.1955000e-01  2.9990000e-02  2.0763000e-02  2.0734000e-03 /
"""
    # testList = re.split("\n", test)
    # reactionBlock = getReactionBlock(testList, 12)
    #
    # # for line in reactionBlock: print line
    #
    # (reactantSet, productSet) = readKineticsReaction(reactionBlock[0])
    # print reactantSet
    # print productSet

