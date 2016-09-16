'''
This script is used to export the RMG-database to a human readable csv file. 
To use: input as the first argument the path to a library or depository folder, and a second argument path to an output csv file:
eg: python generateDatabaseCsv ~/RMG-database/input/kinetics/families/H_Abstraction/training ~/RMG-database/test.csv

Last edited by Nate before July 29th 2015.
'''
import os.path
import re
import csv
import argparse
from rmgpy.data.rmg import RMGDatabase
from rmgpy import settings
from rmgpy.kinetics import Arrhenius, Chebyshev, ThirdBody,Troe, KineticsData


def getDatabaseType(directory):
    """
    This function takes a directory name and returns a string giving the type of database entry.
    """
    if os.path.isdir(directory):
        if os.path.exists(os.path.join(directory, 'rules.py')):
            return 'rules'
        #This includes depositories and libraries
        if re.search(re.escape('kinetics\libraries'), directory):
            return 'kineticsLibrary'
        if re.search('families.*training', directory):
            return 'trainingSet'

if __name__ == '__main__':
    # load database
    parser = argparse.ArgumentParser()
    parser.add_argument('file', metavar='FILE', type=str, nargs='?',
        help='the directory for the library or rules to export')
    parser.add_argument('output', metavar='OUTPUT', type=str, nargs='?', help='the csv where we give the ouput')
    
    args = parser.parse_args()
    
    path = settings['database.directory']
    database = RMGDatabase()
    databaseType=getDatabaseType(args.file)
    
    if databaseType=='kineticsLibrary':
        
        database.load(path,
             thermoLibraries=None,
             transportLibraries=None,
             reactionLibraries=[args.file],
             seedMechanisms=None,
             kineticsFamilies='none',
             kineticsDepositories=[],
             statmechLibraries=None,
             depository=True,
             solvation=True,
             )
        with open(args.output, 'wb') as csvfile:
            csvwriter=csv.writer(csvfile)
            entries=database.kinetics.libraries.values()[0].entries
            
            headingWritten=False
            for index in entries:
                x=entries[index].data
                if type(x) is Arrhenius:
                    if not headingWritten:
                        headingWritten=True
                        csvwriter.writerow(['Label', 'index', 'A (mol-sec-m^3)', 'n', 'Ea (kJ/mol)', 'T0 (K)', 'degeneracy'])           
                    csvwriter.writerow([entries[index].label, 
                                        index, 
                                        x._A.value_si,
                                        x._n.value_si, 
                                        x._Ea.value_si,
                                        x._T0.value_si,
                                        entries[index].item.degeneracy])
    
            #Need another for loop here
            #headingWritten=False
    if databaseType=='trainingSet':
        familyName=re.sub('.*families', '', args.file)
        familyName=re.sub('training', '', familyName)
        familyName=re.sub(re.escape('/'), '', familyName)
        familyName=re.sub('\\\\', '', familyName)
        database.load(path,
             thermoLibraries=None,
             transportLibraries=None,
             reactionLibraries=[],
             seedMechanisms=None,
             kineticsFamilies=[familyName],
             kineticsDepositories=[],
             statmechLibraries=None,
             depository=True,
             solvation=True,
             )
        
        templateDict={}
        with open(args.output, 'wb') as csvfile:
            csvwriter=csv.writer(csvfile)
            entries=database.kinetics.families[familyName].depositories[0].entries
            
            headingWritten=False
            for index in entries:
                x=entries[index].data
                template=database.kinetics.families[familyName].getReactionTemplate(entries[index].item)
                templateStr=''
                for group in template:
                    if templateStr!='': templateStr+=';'
                    templateStr+=group.label
                
                
                for key, value in templateDict.iteritems():
                    if templateStr==value:
                        print entries[index].label, 'has matching template to', key
                templateDict[entries[index].label]=templateStr
                
                
                if type(x) is Arrhenius:
                    if not headingWritten:
                        headingWritten=True
                        csvwriter.writerow(['Label', 'index', 'template', 'A (mol-sec-m^3)', 'n', 'Ea (kJ/mol)', 'T0 (K)', 'degeneracy'])
                    csvwriter.writerow([entries[index].label, 
                                        index, 
                                        templateStr,
                                        x._A.value_si,
                                        x._n.value_si, 
                                        x._Ea.value_si,
                                        x._T0.value_si,
                                        entries[index].item.degeneracy])
    

