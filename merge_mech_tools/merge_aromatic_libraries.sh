#!/bin/bash
source deactivate
source activate rmg_env

#Create chemkin mechanisms (with thermo) for each aromatic library and merge with main mechanism

#As opposed to simply "appending" reaction libraries by including ('aromatic_library_name',True) in reactionLibraries list of RMG input file,
#which will only append reactions involving all species existing in the mechanism already,
#this script will append every reaction in the library to the main mechanism, including reactions involving species not appearing in the main mechanism.
#The thermo for these species  will be estimated by RMG according to the input to exportKineticsLibraryToChemkin.py

declare -a list_of_aromatic_libraries=('phenyl_diacetylene_benzofulvenyl_effective' 'C10H11' 'C3' 
'Fulvene_H' 'naphthalene_H' 'vinylCPD_H' 'biCPD_H_shift' '2016_Mebel_C6H4C2H_C2H2_High_P'
'2006_Park_Phenyl_Propene' '2012_Kislov_Phenyl_Propene' '2016_Mebel_C9H9' '2005_Ismail_C6H5_C4H6' 
'2016_Mebel_Indene_CH3' '2015_Wang_K_C6H9' '2015_Buras_Vinyl_1_3_Butadiene' 'New_Phenyl_Propene_Pathway' 
'Methyl_CPD_H_abs' 'C3H3_C7H7_Matsugi' 'C6H5_C4H4_all_TST_rates' 'c-C5H5_CH3_Sharma_TST' 
'C10H9_Mebel_TST' 'kislovB' 'C10H8_HACA')

number_of_aromatic_libraries=${#list_of_aromatic_libraries[@]}

append_aromatic_libraries='yes'

if [ “$append_aromatic_libraries” == “yes” ]; then 
    mkdir aromatic_libraries
    cd aromatic_libraries
    for (( i=0; i<${number_of_aromatic_libraries}; i++ ));
    do
        mkdir ${list_of_aromatic_libraries[$i]}
	cd ${list_of_aromatic_libraries[$i]}
	python /home/zjburas/RMG-database/scripts/exportKineticsLibraryToChemkin.py ${list_of_aromatic_libraries[$i]}
	cd ../../
	python /home/zjburas/RMG-Py/scripts/mergeModels.py --model1 overall_chem_annotated.inp overall_species_dictionary.txt --model2 aromatic_libraries/${list_of_aromatic_libraries[$i]}/chem.inp aromatic_libraries/${list_of_aromatic_libraries[$i]}/species_dictionary.txt
	mv chem.inp overall_chem_annotated.inp
    	mv species_dictionary.txt overall_species_dictionary.txt	
	cd aromatic_libraries
    done
    cd ../
fi