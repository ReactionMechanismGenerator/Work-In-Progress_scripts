from rmgpy.tools.plot import SimulationPlot
import os.path
import numpy
import copy
from rmgpy.tools.data import GenericData

########################################################################################################################
#Inputs

#Experimental Conditions
T = [605.0, 605.0] #K
P = [10.0, 10.0] #Torr
Total_Flow = [320.0, 320.0] #sccm
CalMix_Flow = [10.0,  10.0] #sccm
Photolysis_Power = [25.0, 25.0]#mJ
C6H5I_Flow = [100, 100] #sccm of He through bubbler
C3H6_Flow = [0, 15]

#Specify location of simulated and experimental species profiles

simulation_csv_file=[r'C:\Users\User1\Documents\Research Data\Phenyl + Propene\RMG_runs\No_Pdep\range_of_T_P\forbid_infeasible_Intra_R_Add_or_H_Mig\Phenyl_Propene_Training_Rxns\GAV_thermo\0_1_tol_18_maxC\solver\simulation_1_332.csv',
                     r'C:\Users\User1\Documents\Research Data\Phenyl + Propene\RMG_runs\No_Pdep\range_of_T_P\forbid_infeasible_Intra_R_Add_or_H_Mig\Phenyl_Propene_Training_Rxns\GAV_thermo\0_1_tol_18_maxC\solver\simulation_2_332.csv',
                     ]

experiment_csv_file=[r'C:\Users\User1\Documents\Research Data\Phenyl + Propene\MS_experiments\CSV_files\AutoIntegratedPeakAreas_20170916_600K_10Torr_no_propene.csv',
                     r'C:\Users\User1\Documents\Research Data\Phenyl + Propene\MS_experiments\CSV_files\AutoIntegratedPeakAreas_20170916_600K_10Torr_extra_low_propene.csv',
                     ]

#Specify which simulated and experimental species to plot

simulated_species_to_plot={'C6H5(1)_obs': 'C6H5(1)_obs',
                           'C6H6(49)_obs': 'C6H6(49)_obs',
                            'C7H7(86)_obs': 'C7H7(86)_obs',
                            'C8H8(125)_obs': 'C8H8(125)_obs',
                           'C9H10(62)_obs': 'C9H10(62)_obs',
                           'C9H10(63)_obs': 'C9H10(63)_obs',
                            'C9H11(51)_obs': 'C9H11(51)_obs',
                           'I_obs': 'I_obs',
                           'HI_obs': 'HI_obs',
                           'CH3I_obs': 'CH3I_obs',
                           'C12H10(48)_obs': 'C12H10(48)_obs',
                           'C3H5I_obs': 'C3H5I_obs',
                           'C7H7I_obs': 'C7H7I_obs',
                           'C9H11I_obs': 'C9H11I_obs',
                           }

experiment_species_to_plot={'77': '77',
                           '78': '78',
                           '91': '91',
                           '100': '100',
                           '104': '104',
                           '118': '118',
                           '119': '119',
                            '127': '127',
                            '128': '128',
                            '142': '142',
                            '154': '154',
                            '168': '168',
                           '246': '246',
                           }

#Specify total cross sections (MB, at 10.5 eV) to use
PICS_list ={'C6H5(1)_obs': 17.0,
       'C6H6(49)_obs': 31.8,
       'C9H11(51)_obs': 40.0,
       'C7H7(86)_obs': 25.5,
       'C9H10(62)_obs': 38.8,
       'C9H10(63)_obs': 38.8,
       'C8H8(125)_obs': 42.9,
       'C9H11I_obs': 50.0,
       'CH3I_obs': 48.2,
       'C3H5I_obs': 50.0,
       'C7H7I_obs': 50.0,
       'I_obs': 74.0,
       'HI_obs': 44.0,
       '100': 9.9,
        'C12H10(48)_obs': 64,
       }

#Specify fragmentation patterns of iodide species
C9H11I_frag = {'118': 0.45, '119': 0.45, '246': 0.1}
C7H7I_frag = {'91': 1.0, '218': 0.0}
C3H5I_frag = {'41': 0.5, '168': 0.5}

#Specify internal standards for either "absolute" or "relative" normalization
absolute_internal_standard = '100'

internal_standard_mass_dictionary = {'I_obs': '127',
                                      'C8H8(125)_obs': '104',
                                      }

relative_internal_standard  = ['I_obs', 'C8H8(125)_obs']

########################################################################################################################

#Iterate over all experimental conditions
for i, run in enumerate(experiment_csv_file):

    ########################################################################################################################
    #Load simulated and experimental data

    total_conc = (P[i]/760.0)/(0.08206*T[i])*6.022e23/1000.0 #molecules/cm3
    absolute_internal_standard_conc = CalMix_Flow[i]/Total_Flow[i]*.0001*total_conc #molecules/cm3

    output_directory = os.path.join(os.path.dirname(simulation_csv_file[i]),
                                                '{0}_K_{1}_Torr_{2}_C3H6_{3}_C6H5I_{4}_photolysis'.format(int(T[i]), int(P[i]),int(C3H6_Flow[i]),
                                                                                                 int(C6H5I_Flow[i]),
                                                                                                 int(Photolysis_Power[i])))

    if not os.path.exists(output_directory):
        os.makedirs(output_directory)


    simulated_data = SimulationPlot(ylabel='Mole Fraction',csvFile=simulation_csv_file[i], species=simulated_species_to_plot)

    experiment_data = SimulationPlot(ylabel='Mole Fraction',csvFile=run, species=experiment_species_to_plot)

    simulated_data.load()

    experiment_data.load()

    #Weight the simulated data by cross section

    for dataobject in simulated_data.yVar:
        dataobject.data *= PICS_list[dataobject.label]

    #Normalize on an "absolute" scale by using calibration gas as internal standard

    for dataobject in experiment_data.yVar:
        if dataobject.label == absolute_internal_standard:
            absolute_internal_standard_signal = numpy.mean(dataobject.data)

    for label, PICS in PICS_list.iteritems():
        if label == absolute_internal_standard:
            absolute_internal_standard_PICS = PICS
            break

    yVar_absolute_scaled = []

    for dataobject in simulated_data.yVar:
        dataobject_absolute_scaled = copy.deepcopy(dataobject)
        dataobject_absolute_scaled.data *= absolute_internal_standard_signal/(absolute_internal_standard_PICS*absolute_internal_standard_conc)*total_conc
        yVar_absolute_scaled.append(dataobject_absolute_scaled)

    simulated_data_absolute_scaled = SimulationPlot(xVar = simulated_data.xVar, yVar = yVar_absolute_scaled,ylabel='Absolute Signal', species=simulated_species_to_plot)

    #Subtract an average of experimental signal before t=0 from all experimental signal
    time = experiment_data.xVar.data[0]
    average_experiment_data_background = numpy.zeros(len(experiment_species_to_plot))
    counter = 0

    while time <= 0.0:
        for j, dataobject in enumerate(experiment_data.yVar):
            average_experiment_data_background[j] += dataobject.data[counter]
        counter += 1
        time = experiment_data.xVar.data[counter]

    average_experiment_data_background /= counter

    for j, dataobject in enumerate(experiment_data.yVar):
        dataobject.data -= average_experiment_data_background[j]

    #Normalize on an "relative" scale by using one of the products as internal standard

    for dataobject in simulated_data.yVar:
        if dataobject.label == relative_internal_standard[i]:
            relative_simulated_internal_standard_signal = numpy.max(dataobject.data)
            break

    for dataobject in experiment_data.yVar:
        if dataobject.label == internal_standard_mass_dictionary[relative_internal_standard[i]]:
            relative_experiment_internal_standard_signal = numpy.max(dataobject.data)
            break

    yVar_relative_scaled = []

    for dataobject in simulated_data.yVar:
        dataobject_relative_scaled = copy.deepcopy(dataobject)
        dataobject_relative_scaled.data /= relative_simulated_internal_standard_signal
        yVar_relative_scaled.append(dataobject_relative_scaled)

    simulated_data_relative_scaled = SimulationPlot(xVar = simulated_data.xVar, yVar = yVar_relative_scaled,ylabel='Relative Signal', species=simulated_species_to_plot)

    yVar_relative_scaled = []

    for dataobject in experiment_data.yVar:
        dataobject_relative_scaled = copy.deepcopy(dataobject)
        dataobject_relative_scaled.data /= relative_experiment_internal_standard_signal
        yVar_relative_scaled.append(dataobject_relative_scaled)

    experiment_data_relative_scaled = SimulationPlot(xVar = experiment_data.xVar, yVar = yVar_relative_scaled,ylabel='Relative Signal', species=experiment_species_to_plot)
    ########################################################################################################################

    ########################################################################################################################
    #Plot all simulated profiles together
    simulated_data.plot(filename=os.path.join(output_directory, 'all_species.png'))

    ########################################################################################################################

    ########################################################################################################################
    #Plot I and HI together
    simulated_I_species = {'I_obs': 'I_obs',
                           'HI_obs': 'HI_obs',
                           }

    experiment_I_species = {'127': '127',
                            '128': '128',
                            }

    simulated_data_absolute_scaled_I_species = copy.deepcopy(simulated_data_absolute_scaled)
    simulated_data_relative_scaled_I_species = copy.deepcopy(simulated_data_relative_scaled)
    experiment_data_I_species = copy.deepcopy(experiment_data)
    experiment_data_relative_scaled_I_species = copy.deepcopy(experiment_data_relative_scaled)

    simulated_data_absolute_scaled_I_species.species = simulated_I_species
    simulated_data_relative_scaled_I_species.species = simulated_I_species
    experiment_data_I_species.species = experiment_I_species
    experiment_data_relative_scaled_I_species.species = experiment_I_species

    simulated_data_absolute_scaled_I_species.comparePlot(filename=os.path.join(output_directory,'absolute_I.png'), otherSimulationPlot = experiment_data_I_species)

    simulated_data_relative_scaled_I_species.comparePlot(filename=os.path.join(output_directory,'relative_I.png'), otherSimulationPlot = experiment_data_relative_scaled_I_species)
    ########################################################################################################################

    ########################################################################################################################
    #Plot 104 (styrene)
    simulated_104_species = {'C8H8(125)_obs': 'C8H8(125)_obs',
                           }

    experiment_104_species = {'104': '104',
                            }

    simulated_data_absolute_scaled_104_species = copy.deepcopy(simulated_data_absolute_scaled)
    simulated_data_relative_scaled_104_species = copy.deepcopy(simulated_data_relative_scaled)
    experiment_data_104_species = copy.deepcopy(experiment_data)
    experiment_data_relative_scaled_104_species = copy.deepcopy(experiment_data_relative_scaled)

    simulated_data_absolute_scaled_104_species.species = simulated_104_species
    simulated_data_relative_scaled_104_species.species = simulated_104_species
    experiment_data_104_species.species = experiment_104_species
    experiment_data_relative_scaled_104_species.species = experiment_104_species

    simulated_data_absolute_scaled_104_species.comparePlot(filename=os.path.join(output_directory,'absolute_104.png'), otherSimulationPlot = experiment_data_104_species)

    simulated_data_relative_scaled_104_species.comparePlot(filename=os.path.join(output_directory,'relative_104.png'), otherSimulationPlot = experiment_data_relative_scaled_104_species)
    ########################################################################################################################

    ########################################################################################################################
    #Plot 119
    simulated_119_species = {'C9H11(51)_obs': 'C9H11(51)_obs',
                             'C9H11I_obs': 'C9H11I_obs',
                             'Summed Simulated 119 amu': 'Summed Simulated 119 amu',
                           }

    experiment_119_species = {'119': '119',
                            }

    simulated_data_absolute_scaled_119_species = copy.deepcopy(simulated_data_absolute_scaled)
    simulated_data_relative_scaled_119_species = copy.deepcopy(simulated_data_relative_scaled)
    experiment_data_119_species = copy.deepcopy(experiment_data)
    experiment_data_relative_scaled_119_species = copy.deepcopy(experiment_data_relative_scaled)

    simulated_data_absolute_scaled_119_species.species = simulated_119_species
    simulated_data_relative_scaled_119_species.species = simulated_119_species
    experiment_data_119_species.species = experiment_119_species
    experiment_data_relative_scaled_119_species.species = experiment_119_species

    summed_simulated_data_absolute_scaled_119_species = GenericData(label='Summed Simulated 119 amu',data=numpy.zeros_like(simulated_data_absolute_scaled_119_species.xVar.data))
    for dataobject in simulated_data_absolute_scaled_119_species.yVar:
        if dataobject.label in 'C9H11I_obs':
            dataobject.data *= C9H11I_frag['119']
        if dataobject.label in simulated_119_species:
            summed_simulated_data_absolute_scaled_119_species.data += dataobject.data

    simulated_data_absolute_scaled_119_species.yVar.append(summed_simulated_data_absolute_scaled_119_species)

    summed_simulated_data_relative_scaled_119_species = GenericData(label='Summed Simulated 119 amu',data=numpy.zeros_like(simulated_data_relative_scaled_119_species.xVar.data))
    for dataobject in simulated_data_relative_scaled_119_species.yVar:
        if dataobject.label in 'C9H11I_obs':
            dataobject.data *= C9H11I_frag['119']
        if dataobject.label in simulated_119_species:
            summed_simulated_data_relative_scaled_119_species.data += dataobject.data

    simulated_data_relative_scaled_119_species.yVar.append(summed_simulated_data_relative_scaled_119_species)

    simulated_data_absolute_scaled_119_species.comparePlot(filename=os.path.join(output_directory,'absolute_119.png'), otherSimulationPlot = experiment_data_119_species)

    simulated_data_relative_scaled_119_species.comparePlot(filename=os.path.join(output_directory,'relative_119.png'), otherSimulationPlot = experiment_data_relative_scaled_119_species)
    ########################################################################################################################

    ########################################################################################################################
    #Plot 118
    simulated_118_species = {'C9H11I_obs': 'C9H11I_obs',
                             'C9H10(62)_obs': 'C9H10(62)_obs',
                             'C9H10(63)_obs': 'C9H10(63)_obs',
                             'Summed Simulated 118 amu': 'Summed Simulated 118 amu',
                           }

    experiment_118_species = {'118': '118',
                            }

    simulated_data_absolute_scaled_118_species = copy.deepcopy(simulated_data_absolute_scaled)
    simulated_data_relative_scaled_118_species = copy.deepcopy(simulated_data_relative_scaled)
    experiment_data_118_species = copy.deepcopy(experiment_data)
    experiment_data_relative_scaled_118_species = copy.deepcopy(experiment_data_relative_scaled)

    simulated_data_absolute_scaled_118_species.species = simulated_118_species
    simulated_data_relative_scaled_118_species.species = simulated_118_species
    experiment_data_118_species.species = experiment_118_species
    experiment_data_relative_scaled_118_species.species = experiment_118_species

    summed_simulated_data_absolute_scaled_118_species = GenericData(label='Summed Simulated 118 amu',data=numpy.zeros_like(simulated_data_absolute_scaled_118_species.xVar.data))
    for dataobject in simulated_data_absolute_scaled_118_species.yVar:
        if dataobject.label in 'C9H11I_obs':
            dataobject.data *= C9H11I_frag['118']
        if dataobject.label in simulated_118_species:
            summed_simulated_data_absolute_scaled_118_species.data += dataobject.data

    simulated_data_absolute_scaled_118_species.yVar.append(summed_simulated_data_absolute_scaled_118_species)

    summed_simulated_data_relative_scaled_118_species = GenericData(label='Summed Simulated 118 amu',data=numpy.zeros_like(simulated_data_relative_scaled_118_species.xVar.data))
    for dataobject in simulated_data_relative_scaled_118_species.yVar:
        if dataobject.label in 'C9H11I_obs':
            dataobject.data *= C9H11I_frag['118']
        if dataobject.label in simulated_118_species:
            summed_simulated_data_relative_scaled_118_species.data += dataobject.data

    simulated_data_relative_scaled_118_species.yVar.append(summed_simulated_data_relative_scaled_118_species)

    simulated_data_absolute_scaled_118_species.comparePlot(filename=os.path.join(output_directory,'absolute_118.png'), otherSimulationPlot = experiment_data_118_species)

    simulated_data_relative_scaled_118_species.comparePlot(filename=os.path.join(output_directory,'relative_118.png'), otherSimulationPlot = experiment_data_relative_scaled_118_species)
    ########################################################################################################################

    ########################################################################################################################
    #Plot 91
    simulated_91_species = {'C7H7I_obs': 'C7H7I_obs',
                             'C7H7(86)_obs': 'C7H7(86)_obs',
                             'Summed Simulated 91 amu': 'Summed Simulated 91 amu',
                           }

    experiment_91_species = {'91': '91',
                            }

    simulated_data_absolute_scaled_91_species = copy.deepcopy(simulated_data_absolute_scaled)
    simulated_data_relative_scaled_91_species = copy.deepcopy(simulated_data_relative_scaled)
    experiment_data_91_species = copy.deepcopy(experiment_data)
    experiment_data_relative_scaled_91_species = copy.deepcopy(experiment_data_relative_scaled)

    simulated_data_absolute_scaled_91_species.species = simulated_91_species
    simulated_data_relative_scaled_91_species.species = simulated_91_species
    experiment_data_91_species.species = experiment_91_species
    experiment_data_relative_scaled_91_species.species = experiment_91_species

    summed_simulated_data_absolute_scaled_91_species = GenericData(label='Summed Simulated 91 amu',data=numpy.zeros_like(simulated_data_absolute_scaled_91_species.xVar.data))
    for dataobject in simulated_data_absolute_scaled_91_species.yVar:
        if dataobject.label in 'C7H7I_obs':
            dataobject.data *= C7H7I_frag['91']
        if dataobject.label in simulated_91_species:
            summed_simulated_data_absolute_scaled_91_species.data += dataobject.data

    simulated_data_absolute_scaled_91_species.yVar.append(summed_simulated_data_absolute_scaled_91_species)

    summed_simulated_data_relative_scaled_91_species = GenericData(label='Summed Simulated 91 amu',data=numpy.zeros_like(simulated_data_relative_scaled_91_species.xVar.data))
    for dataobject in simulated_data_relative_scaled_91_species.yVar:
        if dataobject.label in 'C7H7I_obs':
            dataobject.data *= C7H7I_frag['91']
        if dataobject.label in simulated_91_species:
            summed_simulated_data_relative_scaled_91_species.data += dataobject.data

    simulated_data_relative_scaled_91_species.yVar.append(summed_simulated_data_relative_scaled_91_species)

    simulated_data_absolute_scaled_91_species.comparePlot(filename=os.path.join(output_directory,'absolute_91.png'), otherSimulationPlot = experiment_data_91_species)

    simulated_data_relative_scaled_91_species.comparePlot(filename=os.path.join(output_directory,'relative_91.png'), otherSimulationPlot = experiment_data_relative_scaled_91_species)
    ########################################################################################################################

    ########################################################################################################################
    #Plot 77 and 78 together
    simulated_77_species = {'C6H5(1)_obs': 'C6H5(1)_obs',
                           'C6H6(49)_obs': 'C6H6(49)_obs',
                            'C12H10(48)_obs': 'C12H10(48)_obs'
                           }

    experiment_77_species = {'77': '77',
                            '78': '78',
                             '154': '154',
                            }

    simulated_data_absolute_scaled_77_species = copy.deepcopy(simulated_data_absolute_scaled)
    simulated_data_relative_scaled_77_species = copy.deepcopy(simulated_data_relative_scaled)
    experiment_data_77_species = copy.deepcopy(experiment_data)
    experiment_data_relative_scaled_77_species = copy.deepcopy(experiment_data_relative_scaled)

    simulated_data_absolute_scaled_77_species.species = simulated_77_species
    simulated_data_relative_scaled_77_species.species = simulated_77_species
    experiment_data_77_species.species = experiment_77_species
    experiment_data_relative_scaled_77_species.species = experiment_77_species

    simulated_data_absolute_scaled_77_species.comparePlot(filename=os.path.join(output_directory,'absolute_77.png'), otherSimulationPlot = experiment_data_77_species)

    simulated_data_relative_scaled_77_species.comparePlot(filename=os.path.join(output_directory,'relative_77.png'), otherSimulationPlot = experiment_data_relative_scaled_77_species)
    ########################################################################################################################

    ########################################################################################################################
    #Plot Alkyl Iodide species
    simulated_alkyl_I_species = {'C9H11I_obs': 'C9H11I_obs',
                             'CH3I_obs': 'CH3I_obs',
                             'C3H5I_obs': 'C3H5I_obs',
                           }

    experiment_alkyl_I_species = {'246': '246',
                                  '142': '142',
                                  '168': '168',
                            }

    simulated_data_absolute_scaled_alkyl_I_species = copy.deepcopy(simulated_data_absolute_scaled)
    simulated_data_relative_scaled_alkyl_I_species = copy.deepcopy(simulated_data_relative_scaled)
    experiment_data_alkyl_I_species = copy.deepcopy(experiment_data)
    experiment_data_relative_scaled_alkyl_I_species = copy.deepcopy(experiment_data_relative_scaled)

    simulated_data_absolute_scaled_alkyl_I_species.species = simulated_alkyl_I_species
    simulated_data_relative_scaled_alkyl_I_species.species = simulated_alkyl_I_species
    experiment_data_alkyl_I_species.species = experiment_alkyl_I_species
    experiment_data_relative_scaled_alkyl_I_species.species = experiment_alkyl_I_species

    for dataobject in simulated_data_absolute_scaled_alkyl_I_species.yVar:
        if dataobject.label in 'C9H11I_obs':
            dataobject.data *= C9H11I_frag['246']
        if dataobject.label in 'C3H5I_obs':
            dataobject.data *= C3H5I_frag['168']

    for dataobject in simulated_data_relative_scaled_alkyl_I_species.yVar:
        if dataobject.label in 'C9H11I_obs':
            dataobject.data *= C9H11I_frag['246']
        if dataobject.label in 'C3H5I_obs':
            dataobject.data *= C3H5I_frag['168']

    simulated_data_absolute_scaled_alkyl_I_species.comparePlot(filename=os.path.join(output_directory,'absolute_alkyl_I.png'), otherSimulationPlot = experiment_data_alkyl_I_species)

    simulated_data_relative_scaled_alkyl_I_species.comparePlot(filename=os.path.join(output_directory,'relative_alkyl_I.png'), otherSimulationPlot = experiment_data_relative_scaled_alkyl_I_species)
    ########################################################################################################################