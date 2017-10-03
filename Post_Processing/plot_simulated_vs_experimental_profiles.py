from rmgpy.tools.plot import SimulationPlot
import os.path
import numpy as np
import copy
from rmgpy.tools.data import GenericData
import matplotlib.pyplot as plt


def sort_profiles(species_list, unsorted_profiles):

    #Load data to plot and then sort before making plots
    unsorted_profiles.load()

    sorted_yVar = []

    for species_name in species_list:
        for dataobject in unsorted_profiles.yVar:
            if dataobject.label == species_name:
                sorted_yVar.append(dataobject)
                break
        else:
            sorted_yVar.append(GenericData(label=species_name, data=np.zeros_like(unsorted_profiles.xVar.data)))


    unsorted_profiles.yVar = sorted_yVar

def fragment(species_data, fragmentation_pattern, parent_species_name):
    fragment_species = []

    for dataobject in species_data.yVar:
        if parent_species_name in dataobject.label:
            for fragment, fragment_fraction in fragmentation_pattern.iteritems():
                fragment_label = dataobject.label + '_frag_' + str(fragment)
                fragment_data = dataobject.data*fragment_fraction
                fragment_species.append(GenericData(label=fragment_label,data=fragment_data))
            dataobject = None

    species_data.yVar += fragment_species

def sum_species(species_to_sum, species_data, summed_species_name):
    summed_species = GenericData(label=summed_species_name,data=np.zeros_like(species_data.xVar.data))
    for dataobject in species_data.yVar:
        if dataobject.label in species_to_sum:
            summed_species.data += dataobject.data

    species_data.yVar.append(summed_species)

########################################################################################################################
#Inputs

#Arrays of Experimental Conditions
T = [605.0,
     605.0,
     707.0,
     707.0,
     700.0,
     700.0,
     687.0,
     687.0,
     707.0,
     707.0,
     707.0,
     605.0,
     605.0,
     605.0,
     707.0,
     707.0,
     ] #K

P = [10.0,
     10.0,
     10.0,
     10.0,
     25.0,
     25.0,
     50.0,
     50.0,
     10.0,
     10.0,
     10.0,
     10.0,
     10.0,
     10.0,
     10.0,
     10.0,
     ] #Torr

Total_Flow = [320.0,
              320.0,
              280.0,
              280.0,
              700.0,
              700.0,
              1400.0,
              1400.0,
              280.0,
              280.0,
              280.0,
              320.0,
              320.0,
              320.0,
              280.0,
              280.0,
              ] #sccm

CalMix_Flow = [5.0,
               10.0,
               5.0,
               10.0,
               5.0,
               10.0,
               10.0,
               10.0,
               10.0,
               10.0,
               10.0,
               10.0,
               10.0,
               10.0,
               10.0,
               10.0,
               ] #sccm

Photolysis_Power = [25.0,
                    25.0,
                    25.0,
                    25.0,
                    25.0,
                    25.0,
                    25.0,
                    25.0,
                    25.0,
                    25.0,
                    25.0,
                    25.0,
                    25.0,
                    25.0,
                    36.0,
                    12.5,
                    ]#mJ

C6H5I_Flow = [100,
              100,
              100,
              100,
              100,
              100,
              100,
              100,
              100,
              100,
              200,
              100,
              100,
              200,
              70,
              200,
              ] #sccm of He through bubbler

C3H6_Flow = [0,
             15,
             0,
             15,
             0,
             15,
             0,
             15,
             30,
             60,
             30,
             30,
             60,
             30,
             60,
             30,
             ]#sccm

#Specify location of simulated and experimental species profiles

simulation_csv_file=[r'C:\Users\User1\Documents\Research Data\Phenyl + Propene\RMG_runs\W_Pdep_13_maxatoms\solver\simulation_1_240.csv',
                     r'C:\Users\User1\Documents\Research Data\Phenyl + Propene\RMG_runs\W_Pdep_13_maxatoms\solver\simulation_2_240.csv',
                     r'C:\Users\User1\Documents\Research Data\Phenyl + Propene\RMG_runs\W_Pdep_13_maxatoms\solver\simulation_3_240.csv',
                     r'C:\Users\User1\Documents\Research Data\Phenyl + Propene\RMG_runs\W_Pdep_13_maxatoms\solver\simulation_4_240.csv',
                     r'C:\Users\User1\Documents\Research Data\Phenyl + Propene\RMG_runs\W_Pdep_13_maxatoms\solver\simulation_5_240.csv',
                     r'C:\Users\User1\Documents\Research Data\Phenyl + Propene\RMG_runs\W_Pdep_13_maxatoms\solver\simulation_6_240.csv',
                     r'C:\Users\User1\Documents\Research Data\Phenyl + Propene\RMG_runs\W_Pdep_13_maxatoms\solver\simulation_7_240.csv',
                     r'C:\Users\User1\Documents\Research Data\Phenyl + Propene\RMG_runs\W_Pdep_13_maxatoms\solver\simulation_8_240.csv',
                     r'C:\Users\User1\Documents\Research Data\Phenyl + Propene\RMG_runs\W_Pdep_13_maxatoms\solver\simulation_9_240.csv',
                     r'C:\Users\User1\Documents\Research Data\Phenyl + Propene\RMG_runs\W_Pdep_13_maxatoms\solver\simulation_10_240.csv',
                     r'C:\Users\User1\Documents\Research Data\Phenyl + Propene\RMG_runs\W_Pdep_13_maxatoms\solver\simulation_11_240.csv',
                     r'C:\Users\User1\Documents\Research Data\Phenyl + Propene\RMG_runs\W_Pdep_13_maxatoms\solver\simulation_12_240.csv',
                     r'C:\Users\User1\Documents\Research Data\Phenyl + Propene\RMG_runs\W_Pdep_13_maxatoms\solver\simulation_13_240.csv',
                     r'C:\Users\User1\Documents\Research Data\Phenyl + Propene\RMG_runs\W_Pdep_13_maxatoms\solver\simulation_14_240.csv',
                     r'C:\Users\User1\Documents\Research Data\Phenyl + Propene\RMG_runs\W_Pdep_13_maxatoms\solver\simulation_15_240.csv',
                     r'C:\Users\User1\Documents\Research Data\Phenyl + Propene\RMG_runs\W_Pdep_13_maxatoms\solver\simulation_16_240.csv',
                     ]

experiment_csv_file=[r'C:\Users\User1\Documents\Research Data\Phenyl + Propene\MS_experiments\CSV_files\ManuallyIntegratedPeakAreas_20170916_600K_10Torr_no_propene.csv',
                     r'C:\Users\User1\Documents\Research Data\Phenyl + Propene\MS_experiments\CSV_files\ManuallyIntegratedPeakAreas_20170916_600K_10Torr_extra_low_propene.csv',
                     r'C:\Users\User1\Documents\Research Data\Phenyl + Propene\MS_experiments\CSV_files\ManuallyIntegratedPeakAreas_20170912_700K_10Torr_no_propene.csv',
                     r'C:\Users\User1\Documents\Research Data\Phenyl + Propene\MS_experiments\CSV_files\ManuallyIntegratedPeakAreas_20170912_700K_10Torr_extra_low_propene.csv',
                     r'C:\Users\User1\Documents\Research Data\Phenyl + Propene\MS_experiments\CSV_files\ManuallyIntegratedPeakAreas_20170916_700K_25Torr_no_propene.csv',
                     r'C:\Users\User1\Documents\Research Data\Phenyl + Propene\MS_experiments\CSV_files\ManuallyIntegratedPeakAreas_20170916_700K_25Torr_extra_low_propene.csv',
                     r'C:\Users\User1\Documents\Research Data\Phenyl + Propene\MS_experiments\CSV_files\ManuallyIntegratedPeakAreas_20170916_700K_50Torr_no_propene.csv',
                     r'C:\Users\User1\Documents\Research Data\Phenyl + Propene\MS_experiments\CSV_files\ManuallyIntegratedPeakAreas_20170916_700K_50Torr_extra_low_propene.csv',
                     r'C:\Users\User1\Documents\Research Data\Phenyl + Propene\MS_experiments\CSV_files\ManuallyIntegratedPeakAreas_20170912_700K_10Torr_low_propene.csv',
                     r'C:\Users\User1\Documents\Research Data\Phenyl + Propene\MS_experiments\CSV_files\ManuallyIntegratedPeakAreas_20170912_700K_10Torr_medium_propene.csv',
                     r'C:\Users\User1\Documents\Research Data\Phenyl + Propene\MS_experiments\CSV_files\ManuallyIntegratedPeakAreas_20170912_700K_10Torr_low_propene_double_radical.csv',
                     r'C:\Users\User1\Documents\Research Data\Phenyl + Propene\MS_experiments\CSV_files\ManuallyIntegratedPeakAreas_20170916_600K_10Torr_low_propene.csv',
                     r'C:\Users\User1\Documents\Research Data\Phenyl + Propene\MS_experiments\CSV_files\ManuallyIntegratedPeakAreas_20170916_600K_10Torr_medium_propene.csv',
                     r'C:\Users\User1\Documents\Research Data\Phenyl + Propene\MS_experiments\CSV_files\ManuallyIntegratedPeakAreas_20170916_600K_10Torr_low_propene_double_radical.csv',
                     r'C:\Users\User1\Documents\Research Data\Phenyl + Propene\MS_experiments\CSV_files\ManuallyIntegratedPeakAreas_20170912_700K_10Torr_medium_propene_max_laser.csv',
                     r'C:\Users\User1\Documents\Research Data\Phenyl + Propene\MS_experiments\CSV_files\ManuallyIntegratedPeakAreas_20170912_700K_10Torr_low_propene_double_precursor.csv',
                     ]

#Specify which simulated and experimental species to plot

simulated_species_to_plot={'C6H5(1)_obs': 'C6H5(1)_obs',
                           'C6H6(49)_obs': 'C6H6(49)_obs',
                            'C7H7(190)_obs': 'C7H7(190)_obs',
                            'C8H8(114)_obs': 'C8H8(114)_obs',
                           'C9H10(59)_obs': 'C9H10(59)_obs',
                           'C9H10(60)_obs': 'C9H10(60)_obs',
                           'C9H10(113)_obs': 'C9H10(113)_obs',
                            'C9H11(51)_obs': 'C9H11(51)_obs',
                            'C9H11(50)_obs': 'C9H11(50)_obs',
                           'C9H12(74)_obs': 'C9H12(74)_obs',
                           'I_obs': 'I_obs',
                           'HI_obs': 'HI_obs',
                           'S(298)': 'S(298)',
                           'CH3I_obs': 'CH3I_obs',
                           'C12H10(48)_obs': 'C12H10(48)_obs',
                           'mz_160_obs': 'mz_160_obs',
                           'C3H5I_obs': 'C3H5I_obs',
                           'C7H7I_obs': 'C7H7I_obs',
                           'C18H22(77)': 'C18H22(77)',
                           'C9H11I_obs': 'C9H11I_obs',
                           'I2_obs': 'I2_obs',
                           'CH3(18)_obs': 'CH3(18)_obs',
                           'C3H5(45)_obs': 'C3H5(45)_obs',
                           }

experiment_species_to_plot={'42': '42',
                            '54': '54',
                            '68': '68',
                            '77': '77',
                            '78': '78',
                            '84': '84',
                            '91': '91',
                            '100': '100',
                            '104': '104',
                            '118': '118',
                            '119': '119',
                            '120': '120',
                            '127': '127',
                            '128': '128',
                            '134': '134',
                            '142': '142',
                            '154': '154',
                            '160': '160',
                            '168': '168',
                            '218': '218',
                            '238': '238',
                            '246': '246',
                            '254': '254',
                           }

#Specify total cross sections (MB, at 10.5 eV) to use
PICS_list ={'C6H5(1)_obs': 17.0,
       'C6H6(49)_obs': 31.8,
       'C9H11(51)_obs': 10,  #40.0,
       'C9H11(50)_obs': 10,  # 40.0,
       'C7H7(190)_obs': 25.5,
       'C9H10(59)_obs': 38.8,
       'C9H10(60)_obs': 38.8,
       'C9H10(113)_obs': 38.8,
       'C8H8(114)_obs': 42.9,
       'C9H11I_obs': 80,  #50.0,
       'CH3I_obs': 48.2,
       'C3H5I_obs': 50.0,
       'C7H7I_obs': 50.0,
       'I_obs': 74.0,
       'HI_obs': 44.0,
       '100': 9.9,
        'C12H10(48)_obs': 64.0,
        'C9H12(74)_obs': 30.0,
        'C18H22(77)': 60.0,
        'I2_obs': 50,
        'S(298)': 30,
        'mz_160_obs': 40,
        '54': 16.3,
        '68': 14.4,
        '84': 21.3,
        '42': 10.7,
        'CH3(18)_obs': 6.5,
        'C3H5(45)_obs': 6.09,
            }

#Specify dictionary of masses
mass_dictionary = {'C6H5(1)_obs': 77.0,
       'C6H6(49)_obs': 78.0,
       'C9H11(51)_obs': 119.0,
       'C9H11(50)_obs': 119.0,
       'C7H7(190)_obs': 91.0,
       'C9H10(59)_obs': 118.0,
       'C9H10(60)_obs': 118.0,
       'C9H10(113)_obs': 118.0,
       'C8H8(114)_obs': 104.0,
       'C9H11I_obs': 246.0,
       'CH3I_obs': 142.0,
       'C3H5I_obs': 168.0,
       'C7H7I_obs': 218.0,
       'I_obs': 127.0,
       'HI_obs': 128.0,
       '100': 100.0,
        'C12H10(48)_obs': 154.0,
        'C9H12(74)_obs': 120.0,
        'C18H22(77)': 238.0,
        'I2_obs': 254.0,
        'S(298)': 134.0,
        'mz_160_obs': 160.0,
        '54': 54.0,
        '68': 68.0,
        '84': 84.0,
        '42': 42.0,
        'CH3(18)_obs': 15,
        'C3H5(45)_obs': 41,
        }

#Specify fragmentation patterns
#Key is the name of the parent species
#values are dictionaries where key is the fragment mass and value is the fraction of fragmentation to that mass
fragmentation_dict = {'C9H11I_obs': {'118': 0.3, '119': 0.6, '246': 0.1},
                      'C7H7I_obs': {'91': 0.9, '218': 0.1},
                      'C3H5I_obs': {'41': 0.5, '168': 0.5},
                      'C9H11(51)_obs': {'91': 0.5, '119': 0.5},
                      }

#Specify internal standards for either "absolute" or "relative" normalization
absolute_internal_standard = ['54', '68', '84', '100']

relative_internal_standard  = ['I_obs', 'C8H8(114)_obs', 'I_obs', 'C8H8(114)_obs', 'I_obs', 'C8H8(114)_obs', 'I_obs', 'C8H8(114)_obs', 'C8H8(114)_obs', 'C8H8(114)_obs', 'C8H8(114)_obs', 'C8H8(114)_obs', 'C8H8(114)_obs', 'C8H8(114)_obs', 'C8H8(114)_obs', 'C8H8(114)_obs']

relative_internal_standard_type = ['max', 'final', 'max', 'final', 'max', 'final', 'max', 'final', 'final', 'final', 'final', 'final', 'final', 'final', 'final', 'final']

#Specify dictionary of isotopes to include in simulation
#First entry in value list is the label that will be given to the simulated isotopologue
#Second entry is the fraction of the isotopologue
#Third entry is the mass of the isotopologue
simulated_isotopes_dict = {'C6H5(1)_obs': ['C6H5(1)_obs_C13_isotope', 6*0.011, 78.0],
                           'C9H11(51)_obs': ['C9H11(51)_obs_C13_isotope', 9*0.011, 120.0],
                           'C9H11(50)_obs': ['C9H11(50)_obs_C13_isotope', 9*0.011, 120.0],
                           'C9H11I_obs': ['C9H11I_obs_C13_isotope', 9*0.011, 247.0],
                           'C9H10(59)_obs': ['C9H10(59)_obs_C13_isotope', 9*0.011, 119.0],
                           'C9H10(60)_obs': ['C9H10(60)_obs_C13_isotope', 9*0.011, 119.0],
                           'C9H10(113)_obs': ['C9H10(113)_obs_C13_isotope', 9*0.011, 119.0],
                           }

#Specify species to sum
#Key will be the name of the summed species
#Values are lists of species to sum (they should all have the same mass)
species_to_sum_dictionary = {'Summed Simulated 78 amu': ['C6H6(49)_obs','C6H5(1)_obs_C13_isotope'],
                             'Summed Simulated 91 amu': ['C7H7(190)_obs', 'C7H7I_obs_frag_91', 'C9H11(51)_obs_frag_91'],
                             'Summed Simulated 118 amu': ['C9H10(59)_obs', 'C9H10(60)_obs', 'C9H10(113)_obs', 'C9H11I_obs_frag_118'],
                             'Summed Simulated 119 amu': ['C9H11(51)_obs_frag_119', 'C9H11(50)_obs', 'C9H11I_obs_frag_119', 'C9H10(59)_obs_C13_isotope', 'C9H10(60)_obs_C13_isotope', 'C9H10(113)_obs_C13_isotope', 'C9H11I_obs_C13_isotope_frag_118'],
                             'Summed Simulated 120 amu': ['C9H12(74)_obs', 'C9H11(51)_obs_C13_isotope_frag_119', 'C9H11(50)_obs_C13_isotope', 'C9H11I_obs_C13_isotope_frag_119'],
                             }

########################################################################################################################

#Iterate over all experimental conditions
for i, run in enumerate(experiment_csv_file):

    ########################################################################################################################
    #Load simulated and experimental data

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

    ########################################################################################################################

    ########################################################################################################################
    #Create C13 isotope versions of specified simulated species

    simulated_isotopes_list = []

    for dataobject in simulated_data.yVar:
        if dataobject.label in simulated_isotopes_dict.keys():
            simulated_isotope = GenericData(label=simulated_isotopes_dict[dataobject.label][0],
                                            data=dataobject.data*simulated_isotopes_dict[dataobject.label][1])
            simulated_isotopes_list.append(simulated_isotope)
            mass_dictionary[simulated_isotopes_dict[dataobject.label][0]] = simulated_isotopes_dict[dataobject.label][2]
            PICS_list[simulated_isotopes_dict[dataobject.label][0]] = PICS_list[dataobject.label]
            dataobject.data *= (1 - simulated_isotopes_dict[dataobject.label][1])

    simulated_data.yVar += simulated_isotopes_list

    ########################################################################################################################

    ########################################################################################################################
    # Fit a Mass-dependent conversion factor from simulated concentration to signal (Mass Discrimination Factor)
    total_conc = (P[i]/760.0)/(0.08206*T[i])*6.022e23/1000.0 #molecules/cm3
    absolute_internal_standard_conc = CalMix_Flow[i]/Total_Flow[i]*.0001*total_conc #molecules/cm3

    absolute_internal_standard_signal = {}

    for dataobject in experiment_data.yVar:
        if dataobject.label in absolute_internal_standard:
            absolute_internal_standard_signal[dataobject.label] = np.mean(dataobject.data)

    absolute_internal_standard_PICS = {}

    for label, PICS in PICS_list.iteritems():
        if label in absolute_internal_standard:
            absolute_internal_standard_PICS[label] = PICS

    mass_discrimination_factors = []
    absolute_internal_standard_masses = []

    for label, signal in absolute_internal_standard_signal.iteritems():
        mass_discrimination_factor = signal/(absolute_internal_standard_PICS[label]*absolute_internal_standard_conc)
        mass_discrimination_factors.append(mass_discrimination_factor)
        absolute_internal_standard_masses.append(mass_dictionary[label])

    mass_discrimination_fit = np.polyfit(np.log10(absolute_internal_standard_masses), np.log10(mass_discrimination_factors), 1)

    if mass_discrimination_fit[0] > 1 or mass_discrimination_fit[0] < 0:
        mass_discrimination_fit[0] = 0.0
        mass_discrimination_fit[1] = np.log10(np.mean(mass_discrimination_factors))

    mass_discrimination_fit_x = np.arange(15,250,1)
    mass_discrimination_fit_y = 10**(np.polyval(mass_discrimination_fit, np.log10(mass_discrimination_fit_x)))

    fig = plt.figure()

    ax = fig.add_subplot(111)

    ax.plot(absolute_internal_standard_masses, mass_discrimination_factors, '.')

    ax.plot(mass_discrimination_fit_x, mass_discrimination_fit_y, '.')

    plt.xlabel('m/z (amu)')

    plt.ylabel('Mass Discrimination Factor')

    ax.grid('on')
    handles, labels = ax.get_legend_handles_labels()

    filename = os.path.join(output_directory, 'mz_factors.png')

    fig.savefig(filename, bbox_inches='tight')

    ########################################################################################################################

    ########################################################################################################################
    #Scale Data

    #Weight the simulated data by cross section and mass discrimination factor

    for dataobject in simulated_data.yVar:
        dataobject.data *= PICS_list[dataobject.label]*10**(np.polyval(mass_discrimination_fit, np.log10(mass_dictionary[dataobject.label])))*total_conc

    # Normalize on an "absolute" scale by using calibration gas as internal standard

    yVar_absolute_scaled = copy.deepcopy(simulated_data.yVar)

    simulated_data_absolute_scaled = SimulationPlot(xVar = simulated_data.xVar, yVar = yVar_absolute_scaled,ylabel='Absolute Signal', species=simulated_species_to_plot)

    #Subtract an average of experimental signal before t=0 from all experimental signal
    time = experiment_data.xVar.data[0]
    average_experiment_data_background = np.zeros(len(experiment_species_to_plot))
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
            if relative_internal_standard_type[i] == 'max':
                relative_simulated_internal_standard_signal = np.max(dataobject.data)
            elif relative_internal_standard_type[i] == 'final':
                relative_simulated_internal_standard_signal = dataobject.data[-1]
            break

    for dataobject in experiment_data.yVar:
        if dataobject.label == str(int(mass_dictionary[relative_internal_standard[i]])):
            if relative_internal_standard_type[i] == 'max':
                relative_experiment_internal_standard_signal = np.max(dataobject.data)
            elif relative_internal_standard_type[i] == 'final':
                relative_experiment_internal_standard_signal = np.mean(dataobject.data[-3:])
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
    #Apply fragmentation patterns

    for parent_species_name, fragmentation_pattern in fragmentation_dict.iteritems():
        fragment(species_data = simulated_data_absolute_scaled, fragmentation_pattern = fragmentation_pattern, parent_species_name = parent_species_name)
        fragment(species_data=simulated_data_relative_scaled, fragmentation_pattern=fragmentation_pattern,
                 parent_species_name=parent_species_name)

    ########################################################################################################################

    ########################################################################################################################
    #Sum species of the same mass

    for summed_species_name, species_to_sum in species_to_sum_dictionary.iteritems():
        sum_species(species_to_sum = species_to_sum, species_data = simulated_data_absolute_scaled, summed_species_name = summed_species_name)
        sum_species(species_to_sum=species_to_sum, species_data=simulated_data_relative_scaled,
                    summed_species_name=summed_species_name)

    ########################################################################################################################

    ########################################################################################################################
    #Plot all simulated profiles together
    simulated_data.plot(filename=os.path.join(output_directory, 'all_species.png'))

    ########################################################################################################################


    ########################################################################################################################
    # Plot summed major species together
    simulated_major_species = ['Summed Simulated 78 amu',
                            'Summed Simulated 91 amu',
                            'C8H8(114)_obs',
                            'Summed Simulated 118 amu',
                            'Summed Simulated 119 amu',
                            ]

    experiment_major_species = ['78',
                             '91',
                             '104',
                             '118',
                             '119',
                             ]

    simulated_data_absolute_scaled_major_species = copy.deepcopy(simulated_data_absolute_scaled)
    simulated_data_relative_scaled_major_species = copy.deepcopy(simulated_data_relative_scaled)
    experiment_data_major_species = copy.deepcopy(experiment_data)
    experiment_data_relative_scaled_major_species = copy.deepcopy(experiment_data_relative_scaled)

    simulated_data_absolute_scaled_major_species.species = dict(zip(simulated_major_species, simulated_major_species))
    simulated_data_relative_scaled_major_species.species = dict(zip(simulated_major_species, simulated_major_species))
    experiment_data_major_species.species = dict(zip(experiment_major_species, experiment_major_species))
    experiment_data_relative_scaled_major_species.species = dict(zip(experiment_major_species, experiment_major_species))

    # Sort the profiles so that they will appear in a predictable order in plot
    sort_profiles(simulated_major_species, simulated_data_absolute_scaled_major_species)
    sort_profiles(simulated_major_species, simulated_data_relative_scaled_major_species)
    sort_profiles(experiment_major_species, experiment_data_major_species)
    sort_profiles(experiment_major_species, experiment_data_relative_scaled_major_species)

    simulated_data_absolute_scaled_major_species.comparePlot(filename=os.path.join(output_directory, 'absolute_major.png'),
                                                          otherSimulationPlot=experiment_data_major_species,
                                                          dataAlreadyLoaded=True)

    simulated_data_relative_scaled_major_species.comparePlot(filename=os.path.join(output_directory, 'relative_major.png'),
                                                          otherSimulationPlot=experiment_data_relative_scaled_major_species,
                                                          dataAlreadyLoaded=True)
    ########################################################################################################################

    # Plot 15 amu
    simulated_15_species = ['CH3(18)_obs',
                            ]

    experiment_15_species = ['15',
                             ]

    simulated_data_absolute_scaled_15_species = copy.deepcopy(simulated_data_absolute_scaled)
    simulated_data_relative_scaled_15_species = copy.deepcopy(simulated_data_relative_scaled)
    experiment_data_15_species = copy.deepcopy(experiment_data)
    experiment_data_relative_scaled_15_species = copy.deepcopy(experiment_data_relative_scaled)

    simulated_data_absolute_scaled_15_species.species = dict(zip(simulated_15_species, simulated_15_species))
    simulated_data_relative_scaled_15_species.species = dict(zip(simulated_15_species, simulated_15_species))
    experiment_data_15_species.species = dict(zip(experiment_15_species, experiment_15_species))
    experiment_data_relative_scaled_15_species.species = dict(zip(experiment_15_species, experiment_15_species))

    # Sort the profiles so that they will appear in a predictable order in plot
    sort_profiles(simulated_15_species, simulated_data_absolute_scaled_15_species)
    sort_profiles(simulated_15_species, simulated_data_relative_scaled_15_species)
    sort_profiles(experiment_15_species, experiment_data_15_species)
    sort_profiles(experiment_15_species, experiment_data_relative_scaled_15_species)

    simulated_data_absolute_scaled_15_species.comparePlot(filename=os.path.join(output_directory, 'absolute_15.png'),
                                                          otherSimulationPlot=experiment_data_15_species,
                                                          dataAlreadyLoaded=True)

    simulated_data_relative_scaled_15_species.comparePlot(filename=os.path.join(output_directory, 'relative_15.png'),
                                                          otherSimulationPlot=experiment_data_relative_scaled_15_species,
                                                          dataAlreadyLoaded=True)
    ########################################################################################################################

    ########################################################################################################################
    # Plot allyl and phenyl together
    simulated_41_species = ['C3H5(45)_obs',
                            'C6H5(1)_obs',
                            ]

    experiment_41_species = ['41',
                             '77'
                             ]

    simulated_data_absolute_scaled_41_species = copy.deepcopy(simulated_data_absolute_scaled)
    simulated_data_relative_scaled_41_species = copy.deepcopy(simulated_data_relative_scaled)
    experiment_data_41_species = copy.deepcopy(experiment_data)
    experiment_data_relative_scaled_41_species = copy.deepcopy(experiment_data_relative_scaled)

    simulated_data_absolute_scaled_41_species.species = dict(zip(simulated_41_species, simulated_41_species))
    simulated_data_relative_scaled_41_species.species = dict(zip(simulated_41_species, simulated_41_species))
    experiment_data_41_species.species = dict(zip(experiment_41_species, experiment_41_species))
    experiment_data_relative_scaled_41_species.species = dict(zip(experiment_41_species, experiment_41_species))

    # Sort the profiles so that they will appear in a predictable order in plot
    sort_profiles(simulated_41_species, simulated_data_absolute_scaled_41_species)
    sort_profiles(simulated_41_species, simulated_data_relative_scaled_41_species)
    sort_profiles(experiment_41_species, experiment_data_41_species)
    sort_profiles(experiment_41_species, experiment_data_relative_scaled_41_species)

    simulated_data_absolute_scaled_41_species.comparePlot(filename=os.path.join(output_directory, 'absolute_41.png'),
                                                          otherSimulationPlot=experiment_data_41_species,
                                                          dataAlreadyLoaded=True)

    simulated_data_relative_scaled_41_species.comparePlot(filename=os.path.join(output_directory, 'relative_41.png'),
                                                          otherSimulationPlot=experiment_data_relative_scaled_41_species,
                                                          dataAlreadyLoaded=True)
    ########################################################################################################################

    ########################################################################################################################
    #Plot 77 and 154 together
    simulated_77_species = ['C6H5(1)_obs',
                            'C12H10(48)_obs',
                           ]

    experiment_77_species = ['77',
                             '154',
                            ]

    simulated_data_absolute_scaled_77_species = copy.deepcopy(simulated_data_absolute_scaled)
    simulated_data_relative_scaled_77_species = copy.deepcopy(simulated_data_relative_scaled)
    experiment_data_77_species = copy.deepcopy(experiment_data)
    experiment_data_relative_scaled_77_species = copy.deepcopy(experiment_data_relative_scaled)

    simulated_data_absolute_scaled_77_species.species = dict(zip(simulated_77_species, simulated_77_species))
    simulated_data_relative_scaled_77_species.species = dict(zip(simulated_77_species, simulated_77_species))
    experiment_data_77_species.species = dict(zip(experiment_77_species, experiment_77_species))
    experiment_data_relative_scaled_77_species.species = dict(zip(experiment_77_species, experiment_77_species))

    #Sort the profiles so that they will appear in a predictable order in plot
    sort_profiles(simulated_77_species, simulated_data_absolute_scaled_77_species)
    sort_profiles(simulated_77_species, simulated_data_relative_scaled_77_species)
    sort_profiles(experiment_77_species, experiment_data_77_species)
    sort_profiles(experiment_77_species, experiment_data_relative_scaled_77_species)

    simulated_data_absolute_scaled_77_species.comparePlot(filename=os.path.join(output_directory,'absolute_77.png'), otherSimulationPlot = experiment_data_77_species, dataAlreadyLoaded=True)

    simulated_data_relative_scaled_77_species.comparePlot(filename=os.path.join(output_directory,'relative_77.png'), otherSimulationPlot = experiment_data_relative_scaled_77_species, dataAlreadyLoaded=True)
    ########################################################################################################################

    ########################################################################################################################
    # Plot 77, 78 and 154 together
    simulated_no_propene_species = ['C6H5(1)_obs',
                            'Summed Simulated 78 amu',
                            'C12H10(48)_obs',
                            ]

    experiment_no_propene_species = ['77',
                             '78',
                             '154',
                             ]

    simulated_data_absolute_scaled_no_propene_species = copy.deepcopy(simulated_data_absolute_scaled)
    simulated_data_relative_scaled_no_propene_species = copy.deepcopy(simulated_data_relative_scaled)
    experiment_data_no_propene_species = copy.deepcopy(experiment_data)
    experiment_data_relative_scaled_no_propene_species = copy.deepcopy(experiment_data_relative_scaled)

    simulated_data_absolute_scaled_no_propene_species.species = dict(zip(simulated_no_propene_species, simulated_no_propene_species))
    simulated_data_relative_scaled_no_propene_species.species = dict(zip(simulated_no_propene_species, simulated_no_propene_species))
    experiment_data_no_propene_species.species = dict(zip(experiment_no_propene_species, experiment_no_propene_species))
    experiment_data_relative_scaled_no_propene_species.species = dict(zip(experiment_no_propene_species, experiment_no_propene_species))

    # Sort the profiles so that they will appear in a predictable order in plot
    sort_profiles(simulated_no_propene_species, simulated_data_absolute_scaled_no_propene_species)
    sort_profiles(simulated_no_propene_species, simulated_data_relative_scaled_no_propene_species)
    sort_profiles(experiment_no_propene_species, experiment_data_no_propene_species)
    sort_profiles(experiment_no_propene_species, experiment_data_relative_scaled_no_propene_species)

    simulated_data_absolute_scaled_no_propene_species.comparePlot(filename=os.path.join(output_directory, 'absolute_noC3H6.png'),
                                                          otherSimulationPlot=experiment_data_no_propene_species,
                                                          dataAlreadyLoaded=True)

    simulated_data_relative_scaled_no_propene_species.comparePlot(filename=os.path.join(output_directory, 'relative_noC3H6.png'),
                                                          otherSimulationPlot=experiment_data_relative_scaled_no_propene_species,
                                                          dataAlreadyLoaded=True)
    ########################################################################################################################

    ########################################################################################################################
    # Plot 78
    simulated_78_species = ['Summed Simulated 78 amu'] + species_to_sum_dictionary['Summed Simulated 78 amu']

    experiment_78_species = ['78']

    simulated_data_absolute_scaled_78_species = copy.deepcopy(simulated_data_absolute_scaled)
    simulated_data_relative_scaled_78_species = copy.deepcopy(simulated_data_relative_scaled)
    experiment_data_78_species = copy.deepcopy(experiment_data)
    experiment_data_relative_scaled_78_species = copy.deepcopy(experiment_data_relative_scaled)

    simulated_data_absolute_scaled_78_species.species = dict(zip(simulated_78_species, simulated_78_species))
    simulated_data_relative_scaled_78_species.species = dict(zip(simulated_78_species, simulated_78_species))
    experiment_data_78_species.species = dict(zip(experiment_78_species, experiment_78_species))
    experiment_data_relative_scaled_78_species.species = dict(zip(experiment_78_species, experiment_78_species))

    #Sort the profiles so that they will appear in a predictable order in plot
    sort_profiles(simulated_78_species, simulated_data_absolute_scaled_78_species)
    sort_profiles(simulated_78_species, simulated_data_relative_scaled_78_species)
    sort_profiles(experiment_78_species, experiment_data_78_species)
    sort_profiles(experiment_78_species, experiment_data_relative_scaled_78_species)

    simulated_data_absolute_scaled_78_species.comparePlot(filename=os.path.join(output_directory, 'absolute_78.png'),
                                                          otherSimulationPlot=experiment_data_78_species, dataAlreadyLoaded=True)

    simulated_data_relative_scaled_78_species.comparePlot(filename=os.path.join(output_directory, 'relative_78.png'),
                                                          otherSimulationPlot=experiment_data_relative_scaled_78_species, dataAlreadyLoaded=True)
    ########################################################################################################################

    ########################################################################################################################
    #Plot 91

    simulated_91_species = ['Summed Simulated 91 amu'] + species_to_sum_dictionary['Summed Simulated 91 amu']

    experiment_91_species = ['91']


    simulated_data_absolute_scaled_91_species = copy.deepcopy(simulated_data_absolute_scaled)
    simulated_data_relative_scaled_91_species = copy.deepcopy(simulated_data_relative_scaled)
    experiment_data_91_species = copy.deepcopy(experiment_data)
    experiment_data_relative_scaled_91_species = copy.deepcopy(experiment_data_relative_scaled)

    simulated_data_absolute_scaled_91_species.species = dict(zip(simulated_91_species, simulated_91_species))
    simulated_data_relative_scaled_91_species.species = dict(zip(simulated_91_species, simulated_91_species))
    experiment_data_91_species.species = dict(zip(experiment_91_species, experiment_91_species))
    experiment_data_relative_scaled_91_species.species = dict(zip(experiment_91_species, experiment_91_species))

    #Sort the profiles so that they will appear in a predictable order in plot
    sort_profiles(simulated_91_species, simulated_data_absolute_scaled_91_species)
    sort_profiles(simulated_91_species, simulated_data_relative_scaled_91_species)
    sort_profiles(experiment_91_species, experiment_data_91_species)
    sort_profiles(experiment_91_species, experiment_data_relative_scaled_91_species)

    simulated_data_absolute_scaled_91_species.comparePlot(filename=os.path.join(output_directory,'absolute_91.png'), otherSimulationPlot = experiment_data_91_species, dataAlreadyLoaded=True)

    simulated_data_relative_scaled_91_species.comparePlot(filename=os.path.join(output_directory,'relative_91.png'), otherSimulationPlot = experiment_data_relative_scaled_91_species, dataAlreadyLoaded=True)
    ########################################################################################################################

    ########################################################################################################################
    #Plot 104 (styrene)
    simulated_104_species = ['C8H8(114)_obs',
                           ]

    experiment_104_species = ['104',
                            ]

    simulated_data_absolute_scaled_104_species = copy.deepcopy(simulated_data_absolute_scaled)
    simulated_data_relative_scaled_104_species = copy.deepcopy(simulated_data_relative_scaled)
    experiment_data_104_species = copy.deepcopy(experiment_data)
    experiment_data_relative_scaled_104_species = copy.deepcopy(experiment_data_relative_scaled)

    simulated_data_absolute_scaled_104_species.species = dict(zip(simulated_104_species, simulated_104_species))
    simulated_data_relative_scaled_104_species.species = dict(zip(simulated_104_species, simulated_104_species))
    experiment_data_104_species.species = dict(zip(experiment_104_species, experiment_104_species))
    experiment_data_relative_scaled_104_species.species = dict(zip(experiment_104_species, experiment_104_species))

    #Sort the profiles so that they will appear in a predictable order in plot
    sort_profiles(simulated_104_species, simulated_data_absolute_scaled_104_species)
    sort_profiles(simulated_104_species, simulated_data_relative_scaled_104_species)
    sort_profiles(experiment_104_species, experiment_data_104_species)
    sort_profiles(experiment_104_species, experiment_data_relative_scaled_104_species)

    simulated_data_absolute_scaled_104_species.comparePlot(filename=os.path.join(output_directory,'absolute_104.png'), otherSimulationPlot = experiment_data_104_species, dataAlreadyLoaded=True)

    simulated_data_relative_scaled_104_species.comparePlot(filename=os.path.join(output_directory,'relative_104.png'), otherSimulationPlot = experiment_data_relative_scaled_104_species, dataAlreadyLoaded=True)
    ########################################################################################################################

    ########################################################################################################################
    #Plot 118
    simulated_118_species = ['Summed Simulated 118 amu'] + species_to_sum_dictionary['Summed Simulated 118 amu']

    experiment_118_species = ['118']

    simulated_data_absolute_scaled_118_species = copy.deepcopy(simulated_data_absolute_scaled)
    simulated_data_relative_scaled_118_species = copy.deepcopy(simulated_data_relative_scaled)
    experiment_data_118_species = copy.deepcopy(experiment_data)
    experiment_data_relative_scaled_118_species = copy.deepcopy(experiment_data_relative_scaled)

    simulated_data_absolute_scaled_118_species.species = dict(zip(simulated_118_species, simulated_118_species))
    simulated_data_relative_scaled_118_species.species = dict(zip(simulated_118_species, simulated_118_species))
    experiment_data_118_species.species = dict(zip(experiment_118_species, experiment_118_species))
    experiment_data_relative_scaled_118_species.species = dict(zip(experiment_118_species, experiment_118_species))

    #Sort the profiles so that they will appear in a predictable order in plot
    sort_profiles(simulated_118_species, simulated_data_absolute_scaled_118_species)
    sort_profiles(simulated_118_species, simulated_data_relative_scaled_118_species)
    sort_profiles(experiment_118_species, experiment_data_118_species)
    sort_profiles(experiment_118_species, experiment_data_relative_scaled_118_species)

    simulated_data_absolute_scaled_118_species.comparePlot(filename=os.path.join(output_directory,'absolute_118.png'), otherSimulationPlot = experiment_data_118_species, dataAlreadyLoaded=True)

    simulated_data_relative_scaled_118_species.comparePlot(filename=os.path.join(output_directory,'relative_118.png'), otherSimulationPlot = experiment_data_relative_scaled_118_species, dataAlreadyLoaded=True)

    ########################################################################################################################

    ########################################################################################################################
    #Plot 119
    simulated_119_species = ['Summed Simulated 119 amu'] + species_to_sum_dictionary['Summed Simulated 119 amu']

    experiment_119_species = ['119']

    simulated_data_absolute_scaled_119_species = copy.deepcopy(simulated_data_absolute_scaled)
    simulated_data_relative_scaled_119_species = copy.deepcopy(simulated_data_relative_scaled)
    experiment_data_119_species = copy.deepcopy(experiment_data)
    experiment_data_relative_scaled_119_species = copy.deepcopy(experiment_data_relative_scaled)

    simulated_data_absolute_scaled_119_species.species = dict(zip(simulated_119_species, simulated_119_species))
    simulated_data_relative_scaled_119_species.species = dict(zip(simulated_119_species, simulated_119_species))
    experiment_data_119_species.species = dict(zip(experiment_119_species, experiment_119_species))
    experiment_data_relative_scaled_119_species.species = dict(zip(experiment_119_species, experiment_119_species))

    #Sort the profiles so that they will appear in a predictable order in plot
    sort_profiles(simulated_119_species, simulated_data_absolute_scaled_119_species)
    sort_profiles(simulated_119_species, simulated_data_relative_scaled_119_species)
    sort_profiles(experiment_119_species, experiment_data_119_species)
    sort_profiles(experiment_119_species, experiment_data_relative_scaled_119_species)

    simulated_data_absolute_scaled_119_species.comparePlot(filename=os.path.join(output_directory,'absolute_119.png'), otherSimulationPlot = experiment_data_119_species, dataAlreadyLoaded=True)

    simulated_data_relative_scaled_119_species.comparePlot(filename=os.path.join(output_directory,'relative_119.png'), otherSimulationPlot = experiment_data_relative_scaled_119_species, dataAlreadyLoaded=True)
    ########################################################################################################################

    ########################################################################################################################
    #Plot 120 amu
    simulated_120_species = ['Summed Simulated 120 amu'] + species_to_sum_dictionary['Summed Simulated 120 amu']

    experiment_120_species = ['120']

    simulated_data_absolute_scaled_120_species = copy.deepcopy(simulated_data_absolute_scaled)
    simulated_data_relative_scaled_120_species = copy.deepcopy(simulated_data_relative_scaled)
    experiment_data_120_species = copy.deepcopy(experiment_data)
    experiment_data_relative_scaled_120_species = copy.deepcopy(experiment_data_relative_scaled)

    simulated_data_absolute_scaled_120_species.species = dict(zip(simulated_120_species, simulated_120_species))
    simulated_data_relative_scaled_120_species.species = dict(zip(simulated_120_species, simulated_120_species))
    experiment_data_120_species.species = dict(zip(experiment_120_species, experiment_120_species))
    experiment_data_relative_scaled_120_species.species = dict(zip(experiment_120_species, experiment_120_species))

    #Sort the profiles so that they will appear in a predictable order in plot
    sort_profiles(simulated_120_species, simulated_data_absolute_scaled_120_species)
    sort_profiles(simulated_120_species, simulated_data_relative_scaled_120_species)
    sort_profiles(experiment_120_species, experiment_data_120_species)
    sort_profiles(experiment_120_species, experiment_data_relative_scaled_120_species)

    simulated_data_absolute_scaled_120_species.comparePlot(filename=os.path.join(output_directory,'absolute_120.png'), otherSimulationPlot = experiment_data_120_species, dataAlreadyLoaded=True)

    simulated_data_relative_scaled_120_species.comparePlot(filename=os.path.join(output_directory,'relative_120.png'), otherSimulationPlot = experiment_data_relative_scaled_120_species, dataAlreadyLoaded=True)
    ########################################################################################################################

    ########################################################################################################################
    #Plot remaining 119 bimolecular products together (134, 160 and 238 amu)
    simulated_119_bi_products_species = ['S(298)',
                                         'mz_160_obs',
                                         'C18H22(77)',
                                         ]

    experiment_119_bi_products_species = ['134',
                                          '160',
                                          '238',
                                         ]

    simulated_data_absolute_scaled_119_bi_products_species = copy.deepcopy(simulated_data_absolute_scaled)
    simulated_data_relative_scaled_119_bi_products_species = copy.deepcopy(simulated_data_relative_scaled)
    experiment_data_119_bi_products_species = copy.deepcopy(experiment_data)
    experiment_data_relative_scaled_119_bi_products_species = copy.deepcopy(experiment_data_relative_scaled)

    simulated_data_absolute_scaled_119_bi_products_species.species = dict(zip(simulated_119_bi_products_species, simulated_119_bi_products_species))
    simulated_data_relative_scaled_119_bi_products_species.species = dict(zip(simulated_119_bi_products_species, simulated_119_bi_products_species))
    experiment_data_119_bi_products_species.species = dict(zip(experiment_119_bi_products_species, experiment_119_bi_products_species))
    experiment_data_relative_scaled_119_bi_products_species.species = dict(zip(experiment_119_bi_products_species, experiment_119_bi_products_species))

    #Sort the profiles so that they will appear in a predictable order in plot
    sort_profiles(simulated_119_bi_products_species, simulated_data_absolute_scaled_119_bi_products_species)
    sort_profiles(simulated_119_bi_products_species, simulated_data_relative_scaled_119_bi_products_species)
    sort_profiles(experiment_119_bi_products_species, experiment_data_119_bi_products_species)
    sort_profiles(experiment_119_bi_products_species, experiment_data_relative_scaled_119_bi_products_species)

    simulated_data_absolute_scaled_119_bi_products_species.comparePlot(filename=os.path.join(output_directory,'absolute_119_bi.png'), otherSimulationPlot = experiment_data_119_bi_products_species, dataAlreadyLoaded=True)

    simulated_data_relative_scaled_119_bi_products_species.comparePlot(filename=os.path.join(output_directory,'relative_119_bi.png'), otherSimulationPlot = experiment_data_relative_scaled_119_bi_products_species, dataAlreadyLoaded=True)
    ########################################################################################################################

    ########################################################################################################################
    #Plot I, HI together
    simulated_IHI_species = ['I_obs',
                           'HI_obs',
                           ]

    experiment_IHI_species = ['127',
                            '128',
                            ]

    simulated_data_absolute_scaled_IHI_species = copy.deepcopy(simulated_data_absolute_scaled)
    simulated_data_relative_scaled_IHI_species = copy.deepcopy(simulated_data_relative_scaled)
    experiment_data_IHI_species = copy.deepcopy(experiment_data)
    experiment_data_relative_scaled_IHI_species = copy.deepcopy(experiment_data_relative_scaled)

    simulated_data_absolute_scaled_IHI_species.species = dict(zip(simulated_IHI_species, simulated_IHI_species))
    simulated_data_relative_scaled_IHI_species.species = dict(zip(simulated_IHI_species, simulated_IHI_species))
    experiment_data_IHI_species.species = dict(zip(experiment_IHI_species, experiment_IHI_species))
    experiment_data_relative_scaled_IHI_species.species = dict(zip(experiment_IHI_species, experiment_IHI_species))

    #Sort the profiles so that they will appear in a predictable order in plot
    sort_profiles(simulated_IHI_species, simulated_data_absolute_scaled_IHI_species)
    sort_profiles(simulated_IHI_species, simulated_data_relative_scaled_IHI_species)
    sort_profiles(experiment_IHI_species, experiment_data_IHI_species)
    sort_profiles(experiment_IHI_species, experiment_data_relative_scaled_IHI_species)

    simulated_data_absolute_scaled_IHI_species.comparePlot(filename=os.path.join(output_directory,'absolute_IHI.png'), otherSimulationPlot = experiment_data_IHI_species, dataAlreadyLoaded=True)

    simulated_data_relative_scaled_IHI_species.comparePlot(filename=os.path.join(output_directory,'relative_IHI.png'), otherSimulationPlot = experiment_data_relative_scaled_IHI_species, dataAlreadyLoaded=True)
    ########################################################################################################################

    ########################################################################################################################
    #Plot I by itself
    simulated_I_species = ['I_obs',
                           ]

    experiment_I_species = ['127',
                            ]

    simulated_data_absolute_scaled_I_species = copy.deepcopy(simulated_data_absolute_scaled)
    simulated_data_relative_scaled_I_species = copy.deepcopy(simulated_data_relative_scaled)
    experiment_data_I_species = copy.deepcopy(experiment_data)
    experiment_data_relative_scaled_I_species = copy.deepcopy(experiment_data_relative_scaled)

    simulated_data_absolute_scaled_I_species.species = dict(zip(simulated_I_species, simulated_I_species))
    simulated_data_relative_scaled_I_species.species = dict(zip(simulated_I_species, simulated_I_species))
    experiment_data_I_species.species = dict(zip(experiment_I_species, experiment_I_species))
    experiment_data_relative_scaled_I_species.species = dict(zip(experiment_I_species, experiment_I_species))

    #Sort the profiles so that they will appear in a predictable order in plot
    sort_profiles(simulated_I_species, simulated_data_absolute_scaled_I_species)
    sort_profiles(simulated_I_species, simulated_data_relative_scaled_I_species)
    sort_profiles(experiment_I_species, experiment_data_I_species)
    sort_profiles(experiment_I_species, experiment_data_relative_scaled_I_species)

    # #Pad t<0 with zeros
    # simulated_data_absolute_scaled_I_species.xVar.data = np.insert(simulated_data_absolute_scaled_I_species.xVar.data, 0, [-0.002, 0.0])
    # for dataobject in simulated_data_absolute_scaled_I_species.yVar:
    #     dataobject.data = np.insert(dataobject.data, 0, [0.0, 0.0])
    #
    # simulated_data_relative_scaled_I_species.xVar.data = np.insert(simulated_data_relative_scaled_I_species.xVar.data, 0, [-0.002, 0.0])
    # for dataobject in simulated_data_relative_scaled_I_species.yVar:
    #     dataobject.data = np.insert(dataobject.data, 0, [0.0, 0.0])

    simulated_data_absolute_scaled_I_species.comparePlot(filename=os.path.join(output_directory,'absolute_I.png'), otherSimulationPlot = experiment_data_I_species, dataAlreadyLoaded=True)

    simulated_data_relative_scaled_I_species.comparePlot(filename=os.path.join(output_directory,'relative_I.png'), otherSimulationPlot = experiment_data_relative_scaled_I_species, dataAlreadyLoaded=True)
    ########################################################################################################################

    ########################################################################################################################
    # Plot Alkyl Iodide species
    simulated_alkyl_I_species = ['CH3I_obs',
                                 'C3H5I_obs_frag_168',
                                 'C7H7I_obs_frag_218',
                                 'C9H11I_obs_frag_246',
                                 ]

    experiment_alkyl_I_species = ['142',
                                  '168',
                                  '218',
                                  '246',
                                  ]

    simulated_data_absolute_scaled_alkyl_I_species = copy.deepcopy(simulated_data_absolute_scaled)
    simulated_data_relative_scaled_alkyl_I_species = copy.deepcopy(simulated_data_relative_scaled)
    experiment_data_alkyl_I_species = copy.deepcopy(experiment_data)
    experiment_data_relative_scaled_alkyl_I_species = copy.deepcopy(experiment_data_relative_scaled)

    simulated_data_absolute_scaled_alkyl_I_species.species = dict(zip(simulated_alkyl_I_species, simulated_alkyl_I_species))
    simulated_data_relative_scaled_alkyl_I_species.species = dict(zip(simulated_alkyl_I_species, simulated_alkyl_I_species))
    experiment_data_alkyl_I_species.species = dict(zip(experiment_alkyl_I_species, experiment_alkyl_I_species))
    experiment_data_relative_scaled_alkyl_I_species.species = dict(zip(experiment_alkyl_I_species, experiment_alkyl_I_species))

    #Sort the profiles so that they will appear in a predictable order in plot
    sort_profiles(simulated_alkyl_I_species, simulated_data_absolute_scaled_alkyl_I_species)
    sort_profiles(simulated_alkyl_I_species, simulated_data_relative_scaled_alkyl_I_species)
    sort_profiles(experiment_alkyl_I_species, experiment_data_alkyl_I_species)
    sort_profiles(experiment_alkyl_I_species, experiment_data_relative_scaled_alkyl_I_species)

    simulated_data_absolute_scaled_alkyl_I_species.comparePlot(
        filename=os.path.join(output_directory, 'absolute_alkyl_I.png'),
        otherSimulationPlot=experiment_data_alkyl_I_species, dataAlreadyLoaded=True)

    simulated_data_relative_scaled_alkyl_I_species.comparePlot(
        filename=os.path.join(output_directory, 'relative_alkyl_I.png'),
        otherSimulationPlot=experiment_data_relative_scaled_alkyl_I_species, dataAlreadyLoaded=True)
    ########################################################################################################################

    ########################################################################################################################
    #Plot I2
    simulated_I2_species = ['I2_obs']

    experiment_I2_species = ['254']

    simulated_data_absolute_scaled_I2_species = copy.deepcopy(simulated_data_absolute_scaled)
    simulated_data_relative_scaled_I2_species = copy.deepcopy(simulated_data_relative_scaled)
    experiment_data_I2_species = copy.deepcopy(experiment_data)
    experiment_data_relative_scaled_I2_species = copy.deepcopy(experiment_data_relative_scaled)

    simulated_data_absolute_scaled_I2_species.species = dict(zip(simulated_I2_species, simulated_I2_species))
    simulated_data_relative_scaled_I2_species.species = dict(zip(simulated_I2_species, simulated_I2_species))
    experiment_data_I2_species.species = dict(zip(experiment_I2_species, experiment_I2_species))
    experiment_data_relative_scaled_I2_species.species = dict(zip(experiment_I2_species, experiment_I2_species))

    #Sort the profiles so that they will appear in a predictable order in plot
    sort_profiles(simulated_I2_species, simulated_data_absolute_scaled_I2_species)
    sort_profiles(simulated_I2_species, simulated_data_relative_scaled_I2_species)
    sort_profiles(experiment_I2_species, experiment_data_I2_species)
    sort_profiles(experiment_I2_species, experiment_data_relative_scaled_I2_species)

    simulated_data_absolute_scaled_I2_species.comparePlot(filename=os.path.join(output_directory,'absolute_I2.png'), otherSimulationPlot = experiment_data_I2_species, dataAlreadyLoaded=True)

    simulated_data_relative_scaled_I2_species.comparePlot(filename=os.path.join(output_directory,'relative_I2.png'), otherSimulationPlot = experiment_data_relative_scaled_I2_species, dataAlreadyLoaded=True)
    ########################################################################################################################