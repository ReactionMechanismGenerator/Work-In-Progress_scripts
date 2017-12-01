from rmgpy.tools.plot import SimulationPlot
import os.path
import numpy as np
import copy
from rmgpy.tools.data import GenericData
import matplotlib as mpl
# Force matplotlib to not use any Xwindows backend.
# This must be called before pylab, matplotlib.pyplot, or matplotlib.backends is imported
mpl.use('Agg')
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

def sum_species(species_to_sum, species_data, summed_species_name, weights = None):
    summed_species = GenericData(label=summed_species_name,data=np.zeros_like(species_data.xVar.data))
    for dataobject in species_data.yVar:
        if dataobject.label in species_to_sum:
            if weights:
                summed_species.data += weights[dataobject.label]*dataobject.data
            else:
                summed_species.data += dataobject.data

    species_data.yVar.append(summed_species)

def comparePlot(data, otherData, linecolors, linewidths, markers, filename='', title='', xlabel='', ylabel='',
                legendloc='upper left', xlimits=None, ylimits=None, legend_title='m/z (amu) = ', legend_labels=None,
                logy=False, offset=None):

    mpl.rc('font', family='sans-serif')

    fig = plt.figure()

    ax = fig.add_subplot(111)

    # Plot the sets of data
    for i, plot in enumerate([data, otherData]):
        # Reset the color cycle per plot to get matching colors in each set
        plt.gca().set_prop_cycle(None)

        xVar = plot.xVar
        yVar = plot.yVar
        # Convert yVar to a list if it wasn't one already
        if isinstance(yVar, GenericData):
            yVar = [yVar]

        for j, y in enumerate(yVar):
            if i == 1:
                if legend_labels:
                    label = legend_labels[j]
                else:
                    label = y.label
            else:
                # label = y.label
                label = ''

            if np.any(np.nonzero(y.data)):
                style = markers[i][j] # + '-'
                if logy is True:
                    if offset is not None:
                        ax.semilogy(xVar.data, y.data + offset[j], style, label=label, mfc='none', mec=linecolors[i][j],
                            linewidth=linewidths[i][j], mew=1, c=linecolors[i][j])
                    else:
                        ax.semilogy(xVar.data, y.data, style, label=label, mfc='none', mec=linecolors[i][j],
                                    linewidth=linewidths[i][j], mew=1, c=linecolors[i][j])
                else:
                    if offset is not None:
                        ax.plot(xVar.data, y.data + offset[j], style, label=label, mfc='none', mec = linecolors[i][j], linewidth = linewidths[i][j], mew = 1, c = linecolors[i][j])
                    else:
                        ax.plot(xVar.data, y.data, style, label=label, mfc='none', mec=linecolors[i][j],
                                linewidth=linewidths[i][j], mew=1, c=linecolors[i][j])

                # Plot the second set of data

    # Prioritize using the function's x and y labels, otherwise the labels from this data object
    if xlabel:
        plt.xlabel(xlabel)
    elif data.xlabel:
        plt.xlabel(data.xlabel)
    elif data.xVar.label:
        xlabel = data.xVar.label
        if data.xVar.units: xlabel += ' ({0})'.format(data.xVar.units)
        plt.xlabel(xlabel)

    if ylabel:
        plt.ylabel(ylabel)
    elif data.ylabel:
        plt.ylabel(data.ylabel)

    # Use user inputted title
    if title:
        plt.title(title)

    ax.grid('on')

    # # set the x-spine (see below for more info on `set_position`)
    # ax.spines['left'].set_position('zero')
    #
    # # set the y-spine
    # ax.spines['bottom'].set_position('zero')

    if offset is not None:
        for j, hline in enumerate(offset):
            ax.axhline(y=offset[j], color='k')
    else:
        ax.axhline(y=0, color='k')

    ax.axvline(x=0, color='k')

    ax.set_xlim(xlimits)
    ax.set_ylim(ylimits)

    handles, labels = ax.get_legend_handles_labels()
    if labels:
        # Create a legend outside the plot and adjust width based off of longest legend label
        maxStringLength = max([len(label) for label in labels])
        width = 1.2 + .011 * maxStringLength * 2
        legend = ax.legend(handles, labels, loc=legendloc, numpoints=1, #bbox_to_anchor=(width, 1),
                           ncol=1, title=legend_title)  # bbox_to_anchor=(1.01,.9)
        fig.savefig(filename, bbox_inches='tight') #, bbox_extra_artists=(legend,))
    else:
        fig.savefig(filename, bbox_inches='tight')

def compareSubPlot(data, otherData, linecolors, linewidths, markers, numsubplotrows, rowspans, filename='', title='', xlabel='', ylabel='',
                   legendloc='upper left', xlimits=None, ylimits=None, legend_title='m/z (amu) = ', legend_labels=None,
                   logy=False):

    mpl.rc('font', family='sans-serif')

    fig = plt.figure()

    rowindex = 0

    # Plot the sets of data
    for j, rowspan in enumerate(rowspans):

        # Reset the color cycle per plot to get matching colors in each set
        plt.gca().set_prop_cycle(None)

        for i, plot in enumerate([data, otherData]):

            xVar = plot.xVar
            yVar = plot.yVar

            # Convert yVar to a list if it wasn't one already
            if isinstance(yVar, GenericData):
                yVar = [yVar]

            y = yVar[j]

            if i == 1:
                if legend_labels:
                    label = legend_labels[j]
                else:
                    label = y.label

                ax.text(0.01, 0.95,label+' amu',ha='left',va='top', transform=ax.transAxes)

            else:
                # label = y.label
                label = ''

                ax = plt.subplot2grid((numsubplotrows, 1), (rowindex, 0), rowspan=rowspan)

                rowindex += rowspan

                if j == 0:
                    # Use user inputted title
                    if title:
                        plt.title(title)

            if np.any(np.nonzero(y.data)):
                style = markers[i][j]  # + '-'
                if logy is True:
                    ax.semilogy(xVar.data, y.data, style, label=label, mfc='none', mec=linecolors[i][j],
                                    linewidth=linewidths[i][j], mew=1, c=linecolors[i][j])
                else:
                    ax.plot(xVar.data, y.data, style, label=label, mfc='none', mec=linecolors[i][j],
                                linewidth=linewidths[i][j], mew=1, c=linecolors[i][j])

                ax.set_xlim(xlimits)
                ax.set_ylim(ylimits[j])

        plt.yticks(np.arange(0.0, ylimits[j][1], 0.05))
        plt.tick_params(axis='x', which='both', top='off')

    # Prioritize using the function's x and y labels, otherwise the labels from this data object
    if xlabel:
        plt.xlabel(xlabel)
    elif data.xlabel:
        plt.xlabel(data.xlabel)
    elif data.xVar.label:
        xlabel = data.xVar.label
        if data.xVar.units: xlabel += ' ({0})'.format(data.xVar.units)
        plt.xlabel(xlabel)

    if ylabel:
        fig.text(0.1, 0.5, ylabel, ha='center', va='center', rotation='vertical')
    elif data.ylabel:
        fig.text(0.1, 0.5, data.ylabel, ha='center', va='center', rotation='vertical')

    fig.subplots_adjust(hspace=0)
    plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)
    plt.setp([a.get_yticklabels() for a in fig.axes[:]], visible=False)

    fig.savefig(filename, bbox_inches='tight')

########################################################################################################################
#Inputs

#Arrays of Experimental Conditions
T = [
     605.0,
     605.0,
     605.0,
     605.0,
     605.0,
     707.0,
     707.0,
     707.0,
     707.0,
     707.0,
     707.0,
     700.0,
     700.0,
     687.0,
     687.0,
     605.0,
     605.0,
     605.0,
     707.0,
     707.0,
     707.0,
     605.0,
     605.0,
     605.0,
     707.0,
     707.0,
     707.0,
     ] #K

P = [
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
     10.0,
     10.0,
     10.0,
     10.0,
     ] #Torr

Total_Flow = [
              320.0,
              320.0,
              320.0,
              320.0,
              320.0,
              280.0,
              280.0,
              280.0,
              280.0,
              280.0,
              280.0,
              700.0,
              700.0,
              1400.0,
              1400.0,
              320.0,
              320.0,
              320.0,
              280.0,
              280.0,
              280.0,
              320.0,
              320.0,
              320.0,
              280.0,
              280.0,
              280.0,
              ] #sccm

CalMix_Flow = [
               5.0,
               10.0,
               10.0,
               10.0,
               10.0,
               5.0,
               10.0,
               10.0,
               10.0,
               10.0,
               10.0,
               5.0,
               10.0,
               10.0,
               10.0,
               0.0,
               0.0,
               0.0,
               0.0,
               0.0,
               0.0,
               0.0,
               0.0,
               0.0,
               0.0,
               0.0,
               0.0,
               ] #sccm

Photolysis_Power = [
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
                    25.0,
                    25.0,
                    25.0,
                    25.0,
                    ]#mJ

C6H5I_Flow = [
              100,
              100,
              100,
              100,
              200,
              100,
              100,
              100,
              70,
              200,
              200,
              100,
              100,
              100,
              100,
              100,
              100,
              100,
              100,
              100,
              100,
              100,
              100,
              100,
              100,
              100,
              100,
              ] #sccm of He through bubbler

C3H6_Flow = [
             0,
             15,
             30,
             60,
             30,
             0,
             15,
             30,
             60,
             30,
             30,
             0,
             15,
             0,
             15,
             15,
             30,
             60,
             15,
             30,
             60,
             15,
             30,
             60,
             15,
             30,
             60,
             ]#sccm

I0 = [
      1.93,
      1.99,
      2.02,
      2.02,
      1.94,
      1.97,
      2.03,
      2.11,
      2.11,
      2.01,
      2.02,
      2.0,
      1.94,
      2.05,
      2.07,
      2.0,
      2.0,
      1.99,
      1.88,
      1.93,
      1.92,
      2.2,
      2.2,
      2.14,
      1.91,
      1.97,
      2.12,
     ]#Volts

Wavelength = [
              505.3,
              505.3,
              505.3,
              505.3,
              505.3,
              505.3,
              505.3,
              505.3,
              505.3,
              505.3,
              505.3,
              505.3,
              505.3,
              505.3,
              505.3,
              408.4,
              408.4,
              408.4,
              408.4,
              408.4,
              408.4,
              447.7,
              447.7,
              447.7,
              447.7,
              447.7,
              447.7,
     ] #nm

#Specify location of simulated and experimental species profiles

simulation_csv_file=[
                     r'C:\Users\User1\Documents\Research Data\Phenyl + Propene\RMG_runs\W_Pdep_13_maxatoms\solver\simulation_1_240.csv',
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
                     r'C:\Users\User1\Documents\Research Data\Phenyl + Propene\RMG_runs\W_Pdep_13_maxatoms\solver\simulation_17_240.csv',
                     r'C:\Users\User1\Documents\Research Data\Phenyl + Propene\RMG_runs\W_Pdep_13_maxatoms\solver\simulation_18_240.csv',
                     r'C:\Users\User1\Documents\Research Data\Phenyl + Propene\RMG_runs\W_Pdep_13_maxatoms\solver\simulation_19_240.csv',
                     r'C:\Users\User1\Documents\Research Data\Phenyl + Propene\RMG_runs\W_Pdep_13_maxatoms\solver\simulation_20_240.csv',
                     r'C:\Users\User1\Documents\Research Data\Phenyl + Propene\RMG_runs\W_Pdep_13_maxatoms\solver\simulation_21_240.csv',
                     r'C:\Users\User1\Documents\Research Data\Phenyl + Propene\RMG_runs\W_Pdep_13_maxatoms\solver\simulation_16_240.csv',
                     r'C:\Users\User1\Documents\Research Data\Phenyl + Propene\RMG_runs\W_Pdep_13_maxatoms\solver\simulation_17_240.csv',
                     r'C:\Users\User1\Documents\Research Data\Phenyl + Propene\RMG_runs\W_Pdep_13_maxatoms\solver\simulation_18_240.csv',
                     r'C:\Users\User1\Documents\Research Data\Phenyl + Propene\RMG_runs\W_Pdep_13_maxatoms\solver\simulation_19_240.csv',
                     r'C:\Users\User1\Documents\Research Data\Phenyl + Propene\RMG_runs\W_Pdep_13_maxatoms\solver\simulation_20_240.csv',
                     r'C:\Users\User1\Documents\Research Data\Phenyl + Propene\RMG_runs\W_Pdep_13_maxatoms\solver\simulation_21_240.csv',
                     ]

experiment_csv_file=[
                     r'C:\Users\User1\Documents\Research Data\Phenyl + Propene\MS_experiments\CSV_files\ManuallyIntegratedPeakAreas_20170916_600K_10Torr_no_propene.csv',
                     r'C:\Users\User1\Documents\Research Data\Phenyl + Propene\MS_experiments\CSV_files\ManuallyIntegratedPeakAreas_20170916_600K_10Torr_extra_low_propene.csv',
                     r'C:\Users\User1\Documents\Research Data\Phenyl + Propene\MS_experiments\CSV_files\ManuallyIntegratedPeakAreas_20170916_600K_10Torr_low_propene.csv',
                     r'C:\Users\User1\Documents\Research Data\Phenyl + Propene\MS_experiments\CSV_files\ManuallyIntegratedPeakAreas_20170916_600K_10Torr_medium_propene.csv',
                     r'C:\Users\User1\Documents\Research Data\Phenyl + Propene\MS_experiments\CSV_files\ManuallyIntegratedPeakAreas_20170916_600K_10Torr_low_propene_double_radical.csv',
                     r'C:\Users\User1\Documents\Research Data\Phenyl + Propene\MS_experiments\CSV_files\ManuallyIntegratedPeakAreas_20170912_700K_10Torr_no_propene.csv',
                     r'C:\Users\User1\Documents\Research Data\Phenyl + Propene\MS_experiments\CSV_files\ManuallyIntegratedPeakAreas_20170912_700K_10Torr_extra_low_propene.csv',
                     r'C:\Users\User1\Documents\Research Data\Phenyl + Propene\MS_experiments\CSV_files\ManuallyIntegratedPeakAreas_20170912_700K_10Torr_low_propene.csv',
                     r'C:\Users\User1\Documents\Research Data\Phenyl + Propene\MS_experiments\CSV_files\ManuallyIntegratedPeakAreas_20170912_700K_10Torr_medium_propene_max_laser.csv',
                     r'C:\Users\User1\Documents\Research Data\Phenyl + Propene\MS_experiments\CSV_files\ManuallyIntegratedPeakAreas_20170912_700K_10Torr_low_propene_double_precursor.csv',
                     r'C:\Users\User1\Documents\Research Data\Phenyl + Propene\MS_experiments\CSV_files\ManuallyIntegratedPeakAreas_20170912_700K_10Torr_low_propene_double_radical.csv',
                     r'C:\Users\User1\Documents\Research Data\Phenyl + Propene\MS_experiments\CSV_files\ManuallyIntegratedPeakAreas_20170916_700K_25Torr_no_propene.csv',
                     r'C:\Users\User1\Documents\Research Data\Phenyl + Propene\MS_experiments\CSV_files\ManuallyIntegratedPeakAreas_20170916_700K_25Torr_extra_low_propene.csv',
                     r'C:\Users\User1\Documents\Research Data\Phenyl + Propene\MS_experiments\CSV_files\ManuallyIntegratedPeakAreas_20170916_700K_50Torr_no_propene.csv',
                     r'C:\Users\User1\Documents\Research Data\Phenyl + Propene\MS_experiments\CSV_files\ManuallyIntegratedPeakAreas_20170916_700K_50Torr_extra_low_propene.csv',
                     None,
                     None,
                     None,
                     None,
                     None,
                     None,
                     None,
                     None,
                     None,
                     None,
                     None,
                     None,
                     ]

abs_experiment_csv_file=[
                     r'C:\Users\User1\Documents\Research Data\Phenyl + Propene\Abs_experiments\CSV_files\505_Absorbance_20170916_600K_10Torr_no_propene.csv',
                     r'C:\Users\User1\Documents\Research Data\Phenyl + Propene\Abs_experiments\CSV_files\505_Absorbance_20170916_600K_10Torr_extra_low_propene.csv',
                     r'C:\Users\User1\Documents\Research Data\Phenyl + Propene\Abs_experiments\CSV_files\505_Absorbance_20170916_600K_10Torr_low_propene.csv',
                     r'C:\Users\User1\Documents\Research Data\Phenyl + Propene\Abs_experiments\CSV_files\505_Absorbance_20170916_600K_10Torr_medium_propene.csv',
                     r'C:\Users\User1\Documents\Research Data\Phenyl + Propene\Abs_experiments\CSV_files\505_Absorbance_20170916_600K_10Torr_low_propene_double_radical.csv',
                     r'C:\Users\User1\Documents\Research Data\Phenyl + Propene\Abs_experiments\CSV_files\505_Absorbance_20170912_700K_10Torr_no_propene.csv',
                     r'C:\Users\User1\Documents\Research Data\Phenyl + Propene\Abs_experiments\CSV_files\505_Absorbance_20170912_700K_10Torr_extra_low_propene.csv',
                     r'C:\Users\User1\Documents\Research Data\Phenyl + Propene\Abs_experiments\CSV_files\505_Absorbance_20170912_700K_10Torr_low_propene.csv',
                     r'C:\Users\User1\Documents\Research Data\Phenyl + Propene\Abs_experiments\CSV_files\505_Absorbance_20170912_700K_10Torr_medium_propene_max_laser.csv',
                     r'C:\Users\User1\Documents\Research Data\Phenyl + Propene\Abs_experiments\CSV_files\505_Absorbance_20170912_700K_10Torr_low_propene_double_precursor.csv',
                     r'C:\Users\User1\Documents\Research Data\Phenyl + Propene\Abs_experiments\CSV_files\505_Absorbance_20170912_700K_10Torr_low_propene_double_radical.csv',
                     r'C:\Users\User1\Documents\Research Data\Phenyl + Propene\Abs_experiments\CSV_files\505_Absorbance_20170916_700K_25Torr_no_propene.csv',
                     r'C:\Users\User1\Documents\Research Data\Phenyl + Propene\Abs_experiments\CSV_files\505_Absorbance_20170916_700K_25Torr_extra_low_propene.csv',
                     r'C:\Users\User1\Documents\Research Data\Phenyl + Propene\Abs_experiments\CSV_files\505_Absorbance_20170916_700K_50Torr_no_propene.csv',
                     r'C:\Users\User1\Documents\Research Data\Phenyl + Propene\Abs_experiments\CSV_files\505_Absorbance_20170916_700K_50Torr_extra_low_propene.csv',
                     r'C:\Users\User1\Documents\Research Data\Phenyl + Propene\Abs_experiments\CSV_files\408_Absorbance_20171014_600K_10Torr_extra_low_propene.csv',
                     r'C:\Users\User1\Documents\Research Data\Phenyl + Propene\Abs_experiments\CSV_files\408_Absorbance_20171014_600K_10Torr_low_propene.csv',
                     r'C:\Users\User1\Documents\Research Data\Phenyl + Propene\Abs_experiments\CSV_files\408_Absorbance_20171014_600K_10Torr_medium_propene.csv',
                     r'C:\Users\User1\Documents\Research Data\Phenyl + Propene\Abs_experiments\CSV_files\408_Absorbance_20171014_700K_10Torr_extra_low_propene.csv',
                     r'C:\Users\User1\Documents\Research Data\Phenyl + Propene\Abs_experiments\CSV_files\408_Absorbance_20171014_700K_10Torr_low_propene.csv',
                     r'C:\Users\User1\Documents\Research Data\Phenyl + Propene\Abs_experiments\CSV_files\408_Absorbance_20171014_700K_10Torr_medium_propene.csv',
                     r'C:\Users\User1\Documents\Research Data\Phenyl + Propene\Abs_experiments\CSV_files\447_Absorbance_20171015_600K_10Torr_extra_low_propene.csv',
                     r'C:\Users\User1\Documents\Research Data\Phenyl + Propene\Abs_experiments\CSV_files\447_Absorbance_20171015_600K_10Torr_low_propene.csv',
                     r'C:\Users\User1\Documents\Research Data\Phenyl + Propene\Abs_experiments\CSV_files\447_Absorbance_20171015_600K_10Torr_medium_propene.csv',
                     r'C:\Users\User1\Documents\Research Data\Phenyl + Propene\Abs_experiments\CSV_files\447_Absorbance_20171015_700K_10Torr_extra_low_propene.csv',
                     r'C:\Users\User1\Documents\Research Data\Phenyl + Propene\Abs_experiments\CSV_files\447_Absorbance_20171015_700K_10Torr_low_propene.csv',
                     r'C:\Users\User1\Documents\Research Data\Phenyl + Propene\Abs_experiments\CSV_files\447_Absorbance_20171015_700K_10Torr_medium_propene.csv',

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
                           'I': 'I',
                           'C6H5(1)': 'C6H5(1)',
                           'C3H5(45)': 'C3H5(45)',
                           'C7H7(190)': 'C7H7(190)',
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
        'CH3(18)_obs': 6.7,
        'C3H5(45)_obs': 6.09,
        'I': 74.0,
        'C6H5(1)': 17.0,
        'C3H5(45)': 6.09,
        'C7H7(190)': 25.5,
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
        'I': 127.0,
        'C6H5(1)': 77,
        'C3H5(45)': 41,
        'C7H7(190)': 91.0,
        }

#Specify fragmentation patterns
#Key is the name of the parent species
#values are dictionaries where key is the fragment mass and value is the fraction of fragmentation to that mass
fragmentation_dict = {'C9H11I_obs': {'118': 0.3, '119': 0.6, '246': 0.1},
                      'C7H7I_obs': {'91': 0.91, '218': 0.09},
                      'C3H5I_obs': {'41': 0.55, '168': 0.45},
                      'C9H11(51)_obs': {'91': 0.5, '119': 0.5},
                      }

#Specify internal standards for either "absolute" or "relative" normalization
absolute_internal_standard = ['54', '68', '84', '100']
absolute_internal_standards_plot_symbols = ['s', 'o', '^', 'D', '*', 'x']
absolute_internal_standard_names = ['1,3-Butadiene', 'Furan', 'Cyclohexane', 'Heptane']

relative_internal_standard  = ['I_obs', 'C8H8(114)_obs', 'C8H8(114)_obs', 'C8H8(114)_obs', 'C8H8(114)_obs',
                               'I_obs', 'C8H8(114)_obs', 'C8H8(114)_obs', 'C8H8(114)_obs', 'C8H8(114)_obs', 'C8H8(114)_obs',
                               'I_obs', 'C8H8(114)_obs',
                               'I_obs', 'C8H8(114)_obs']

relative_internal_standard_type = ['max', 'final', 'final', 'final', 'final',
                                   'max', 'final', 'final', 'final', 'final', 'final',
                                   'max', 'final',
                                   'max', 'final']

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
species_to_sum_dictionary = {'Summed Simulated 78': ['C6H6(49)_obs','C6H5(1)_obs_C13_isotope'],
                             'Summed Simulated 91': ['C7H7(190)_obs', 'C7H7I_obs_frag_91', 'C9H11(51)_obs_frag_91'],
                             'Summed Simulated 118': ['C9H10(59)_obs', 'C9H10(60)_obs', 'C9H10(113)_obs', 'C9H11I_obs_frag_118'],
                             'Summed Simulated 119': ['C9H11(51)_obs_frag_119', 'C9H11(50)_obs', 'C9H11I_obs_frag_119', 'C9H10(59)_obs_C13_isotope', 'C9H10(60)_obs_C13_isotope', 'C9H10(113)_obs_C13_isotope', 'C9H11I_obs_C13_isotope_frag_118'],
                             'Summed Simulated 120': ['C9H12(74)_obs', 'C9H11(51)_obs_C13_isotope_frag_119', 'C9H11(50)_obs_C13_isotope', 'C9H11I_obs_C13_isotope_frag_119'],
                             }

#Absorbance weights for allyl and phenyl radical at 408 nm and different temperatures
Absorbance_weights_408 = {600: {'C3H5(45)': 0.93, 'C6H5(1)': 1.0},
                          700: {'C3H5(45)': 1.2, 'C6H5(1)': 1.0},
                         }

#Absorbance weights for benzyl and phenyl radical at 447.7 nm and different temperatures
Absorbance_weights_447 = {600: {'C7H7(190)': 1.0, 'C6H5(1)': 1.0},
                          700: {'C7H7(190)': 1.0, 'C6H5(1)': 1.0},
                         }

#Specify which experiment #'s to plot
expts_to_plot = [7]
########################################################################################################################

#Iterate over all experimental conditions
for i, run in enumerate(experiment_csv_file):

    #Plot both MS and Absorbance results if both are provided
    if run is not None and i + 1 in expts_to_plot:

        ########################################################################################################################
        #Load simulated and experimental data

        output_directory = os.path.join(os.path.dirname(simulation_csv_file[i]),
                                                    '{0}_K_{1}_Torr_{2}_C3H6_{3}_C6H5I_{4}_photolysis_{5}_nm'.format(int(T[i]), int(P[i]),int(C3H6_Flow[i]),
                                                                                                     int(C6H5I_Flow[i]),
                                                                                                     int(Photolysis_Power[i]), int(Wavelength[i])))

        if not os.path.exists(output_directory):
            os.makedirs(output_directory)


        simulated_data = SimulationPlot(ylabel='Mole Fraction',csvFile=simulation_csv_file[i], species=simulated_species_to_plot)

        experiment_data = SimulationPlot(ylabel='Mole Fraction',csvFile=run, species=experiment_species_to_plot)

        abs_experiment_data = SimulationPlot(ylabel='Absorbance', csvFile=abs_experiment_csv_file[i], species={'Intensity': 'Intensity'})

        simulated_data.load()

        experiment_data.load()

        abs_experiment_data.load()

        ########################################################################################################################

        ########################################################################################################################
        #Switch x axis from seconds to milliseconds
        simulated_data.xVar.data *= 1000
        simulated_data.xVar.label = 'Time'
        simulated_data.xVar.units = 'ms'

        experiment_data.xVar.data *= 1000
        experiment_data.xVar.label = 'Time'
        experiment_data.xVar.units = 'ms'

        abs_experiment_data.xVar.data *= 1000
        abs_experiment_data.xVar.label = 'Time'
        abs_experiment_data.xVar.units = 'ms'
        ########################################################################################################################

        ########################################################################################################################
        #Normalize simulation and absorbance experiments to initial Phenyl concentration for comparison

        #Normalized simulated data
        yVar_normalized = []

        for dataobject in simulated_data.yVar:
            if dataobject.label == 'C6H5(1)':
                C6H5_0 = dataobject.data[0]
                break

        for dataobject in simulated_data.yVar:
            dataobject_normalized = copy.deepcopy(dataobject)
            dataobject_normalized.data /= C6H5_0
            yVar_normalized.append(dataobject_normalized)

        simulated_data_normalized = SimulationPlot(xVar=simulated_data.xVar, yVar=yVar_normalized,
                                                        ylabel='Normalized Concentration', species=simulated_species_to_plot)

        #Add a new data series that is the normalized phenyl + allyl simulated data
        sum_species(['C3H5(45)','C6H5(1)'], simulated_data_normalized, 'combined_phenyl_allyl')

        # Add a new data series that is the normalized phenyl + benzyl simulated data
        sum_species(['C7H7(190)', 'C6H5(1)'], simulated_data_normalized, 'combined_phenyl_benzyl',
                    weights=Absorbance_weights_447[round(T[i], -2)])

        #Normalize Absorbance experiments

        yVar_normalized = []

        for dataobject in abs_experiment_data.yVar:
            dataobject_normalized = copy.deepcopy(dataobject)
            dataobject_normalized.data = np.log(I0[i]*1000/(I0[i]*1000+dataobject_normalized.data))
            linear_fit = np.polyfit(abs_experiment_data.xVar.data[np.where((abs_experiment_data.xVar.data >= 0.01) & (abs_experiment_data.xVar.data <=0.1))],
                                    dataobject_normalized.data[np.where((abs_experiment_data.xVar.data >= 0.01) & (abs_experiment_data.xVar.data <=0.1))], 1)
            dataobject_normalized.data /= np.polyval(linear_fit, 0)
            #Do a moving average to smooth data
            dataobject_normalized.data = np.convolve(dataobject_normalized.data, np.ones((5,)) / 5, mode = 'same')
            yVar_normalized.append(dataobject_normalized)

        abs_experiment_data_normalized = SimulationPlot(xVar=abs_experiment_data.xVar, yVar=yVar_normalized,
                                                        ylabel='Normalized Absorbance')

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

        print i+1
        print mass_discrimination_fit[0]
        print 10**mass_discrimination_fit[1]

        mass_discrimination_fit_x = np.arange(15,250,1)
        mass_discrimination_fit_y = 10**(np.polyval(mass_discrimination_fit, np.log10(mass_discrimination_fit_x)))

        fig = plt.figure()

        ax = fig.add_subplot(111)

        #Sort internal standards by mass
        absolute_internal_standard_masses, mass_discrimination_factors = zip(*sorted(zip(absolute_internal_standard_masses, mass_discrimination_factors)))

        for j, absolute_internal_standard_mass in enumerate(absolute_internal_standard_masses):
            ax.errorbar(absolute_internal_standard_mass, mass_discrimination_factors[j], yerr = 0.15*mass_discrimination_factors[j],
                        label=absolute_internal_standard_names[j], marker = absolute_internal_standards_plot_symbols[j], mfc='none', mec = 'k', mew = 1, ecolor = 'k',
                        c='w')

        ax.plot(mass_discrimination_fit_x, mass_discrimination_fit_y, 'g-', lw = 2)

        plt.xlabel('m/z (amu)')

        plt.ylabel('Mass Discrimination Factor')

        ax.grid('off')
        ax.legend(loc='upper left', title='Internal Standards:', numpoints=1)
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

        simulated_data_absolute_scaled = SimulationPlot(xVar = simulated_data.xVar, yVar = yVar_absolute_scaled,ylabel='Integrated MS Signal', species=simulated_species_to_plot)

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
        simulated_major_species = [
                                'Summed Simulated 78',
                                'Summed Simulated 91',
                                'C8H8(114)_obs',
                                'Summed Simulated 118',
                                'Summed Simulated 119',
                                ]

        experiment_major_species = [
                                 '78',
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

        #Define appearance of plot
        linecolors = [['b','g','r','c','m'],['b','g','r','c','m']]
        linewidths = [[2,2,2,2,2], [0.5,0.5,0.5,0.5,0.5]]
        markers = [['-','-','-','-','-'], ['o','s','^','*','D']]
        title = '{0} K, {1} Torr, MS Experiment {2}'.format(int(T[i]), int(P[i]), int(i+1))
        # offset = [0, 0.05, 0.1, 0.3, 0.35]
        rowspans = [8,8,21,8,7]

        comparePlot(data=simulated_data_absolute_scaled_major_species, filename=os.path.join(output_directory, 'absolute_major.png'),
                                                              otherData=experiment_data_major_species,
                    linecolors = linecolors, linewidths = linewidths, markers = markers, title=title, xlimits=(-1, 50), #, ylimits=(-0.01, 0.2),
                    legendloc='best') #, offset = offset)

        comparePlot(data=simulated_data_relative_scaled_major_species, filename=os.path.join(output_directory, 'relative_major.png'),
                                                              otherData=experiment_data_relative_scaled_major_species,
                    linecolors=linecolors, linewidths=linewidths, markers=markers, title=title)

        compareSubPlot(data=simulated_data_absolute_scaled_major_species, otherData=experiment_data_major_species,
                       linecolors = linecolors, linewidths = linewidths, markers = markers, numsubplotrows = 52, rowspans = rowspans,
                       filename=os.path.join(output_directory, 'absolute_major_sublplots.png'), title=title,
                       legendloc='lower right', xlimits=(-1,46), ylimits=[(0,0.07),(0,0.07),(0,0.2),(0,0.07),(0,0.06)])
        ########################################################################################################################

        ########################################################################################################################
        # Plot summed major species together on a log scale
        simulated_major_species = [
                                'Summed Simulated 78',
                                'Summed Simulated 91',
                                'C8H8(114)_obs',
                                'Summed Simulated 118',
                                'Summed Simulated 119',
                                ]

        experiment_major_species = [
                                 '78',
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

        #Define appearance of plot
        linecolors = [['b','g','r','c','m'],['b','g','r','c','m']]
        linewidths = [[2,2,2,2,2], [0.5,0.5,0.5,0.5,0.5]]
        markers = [['','','','',''], ['o','s','^','*','D']]
        title = '{0} K, {1} Torr, MS Experiment {2}'.format(int(T[i]), int(P[i]), int(i+1))

        comparePlot(data=simulated_data_absolute_scaled_major_species, filename=os.path.join(output_directory, 'absolute_major_log.png'),
                                                              otherData=experiment_data_major_species,
                    linecolors = linecolors, linewidths = linewidths, markers = markers, title=title, logy=True, ylimits=(0.001, 1.0))

        comparePlot(data=simulated_data_relative_scaled_major_species, filename=os.path.join(output_directory, 'relative_major_log.png'),
                                                              otherData=experiment_data_relative_scaled_major_species,
                    linecolors=linecolors, linewidths=linewidths, markers=markers, title=title, logy=True)

        ########################################################################################################################
        # Plot summed major species relative to time-dependent styrene
        simulated_major_species = [
                                'Summed Simulated 78',
                                'Summed Simulated 91',
                                'C8H8(114)_obs',
                                'Summed Simulated 118',
                                'Summed Simulated 119',
                                ]

        experiment_major_species = [
                                 '78',
                                 '91',
                                 '104',
                                 '118',
                                 '119',
                                 ]

        simulated_data_absolute_scaled_major_species = copy.deepcopy(simulated_data_absolute_scaled)
        experiment_data_major_species = copy.deepcopy(experiment_data)

        simulated_data_absolute_scaled_major_species.species = dict(zip(simulated_major_species, simulated_major_species))
        experiment_data_major_species.species = dict(zip(experiment_major_species, experiment_major_species))

        # Sort the profiles so that they will appear in a predictable order in plot
        sort_profiles(simulated_major_species, simulated_data_absolute_scaled_major_species)
        sort_profiles(experiment_major_species, experiment_data_major_species)

        #Normalize all major species signals to t-dependent styrene signal
        yVar_relative_scaled = []

        for dataobject in simulated_data_absolute_scaled_major_species.yVar:
            dataobject_relative_scaled = copy.deepcopy(dataobject)
            dataobject_relative_scaled.data /= simulated_data_absolute_scaled_major_species.yVar[2].data
            yVar_relative_scaled.append(dataobject_relative_scaled)

        simulated_data_relative_tdep = SimulationPlot(xVar=simulated_data_absolute_scaled_major_species.xVar, yVar=yVar_relative_scaled,
                                                        ylabel='Signal Relative to Styrene', species=simulated_major_species)

        yVar_relative_scaled = []

        for dataobject in experiment_data_major_species.yVar:
            dataobject_relative_scaled = copy.deepcopy(dataobject)
            dataobject_relative_scaled.data /= experiment_data_major_species.yVar[2].data
            yVar_relative_scaled.append(dataobject_relative_scaled)

        experiment_data_relative_tdep = SimulationPlot(xVar=experiment_data_major_species.xVar, yVar=yVar_relative_scaled,
                                                         ylabel='Signal Relative to Styrene', species=experiment_major_species)

        # Sort the profiles so that they will appear in a predictable order in plot
        # sort_profiles(simulated_major_species, simulated_data_relative_tdep)
        # sort_profiles(experiment_major_species, experiment_data_relative_tdep)

        #Define appearance of plot
        linecolors = [['b','g','r','c','m'],['b','g','r','c','m']]
        linewidths = [[2,2,2,2,2], [0.5,0.5,0.5,0.5,0.5]]
        markers = [['','','','',''], ['o','s','^','*','D']]
        title = '{0} K, {1} Torr, MS Experiment {2}'.format(int(T[i]), int(P[i]), int(i+1))

        comparePlot(data=simulated_data_relative_tdep, filename=os.path.join(output_directory, 'relative_tdep_major.png'),
                                                              otherData=experiment_data_relative_tdep,
                    linecolors = linecolors, linewidths = linewidths, markers = markers, title=title, ylimits=(-0.1, 1.1), legendloc = 'best')
        ########################################################################################################################

        ########################################################################################################################
        # Plot summed 78 amu relative to time-dependent styrene
        simulated_major_species = [
                                'Summed Simulated 78',
                                'C8H8(114)_obs',
                                ]

        experiment_major_species = [
                                 '78',
                                 '104',
                                 ]

        simulated_data_absolute_scaled_major_species = copy.deepcopy(simulated_data_absolute_scaled)
        experiment_data_major_species = copy.deepcopy(experiment_data)

        simulated_data_absolute_scaled_major_species.species = dict(zip(simulated_major_species, simulated_major_species))
        experiment_data_major_species.species = dict(zip(experiment_major_species, experiment_major_species))

        # Sort the profiles so that they will appear in a predictable order in plot
        sort_profiles(simulated_major_species, simulated_data_absolute_scaled_major_species)
        sort_profiles(experiment_major_species, experiment_data_major_species)

        #Normalize all major species signals to t-dependent styrene signal
        yVar_relative_scaled = []

        for dataobject in simulated_data_absolute_scaled_major_species.yVar:
            dataobject_relative_scaled = copy.deepcopy(dataobject)
            dataobject_relative_scaled.data /= simulated_data_absolute_scaled_major_species.yVar[1].data
            yVar_relative_scaled.append(dataobject_relative_scaled)

        simulated_data_relative_tdep = SimulationPlot(xVar=simulated_data_absolute_scaled_major_species.xVar, yVar=yVar_relative_scaled,
                                                        ylabel='Signal Relative to Styrene (104 amu)', species=simulated_major_species)

        yVar_relative_scaled = []

        for dataobject in experiment_data_major_species.yVar:
            dataobject_relative_scaled = copy.deepcopy(dataobject)
            dataobject_relative_scaled.data /= experiment_data_major_species.yVar[1].data
            yVar_relative_scaled.append(dataobject_relative_scaled)

        experiment_data_relative_tdep = SimulationPlot(xVar=experiment_data_major_species.xVar, yVar=yVar_relative_scaled,
                                                         ylabel='Signal Relative to Styrene (104 amu)', species=experiment_major_species)

        # Sort the profiles so that they will appear in a predictable order in plot
        # sort_profiles(simulated_major_species, simulated_data_relative_tdep)
        # sort_profiles(experiment_major_species, experiment_data_relative_tdep)

        #Define appearance of plot
        linecolors = [['b','r'],['b','r']]
        linewidths = [[2,2,2,2,2], [0.5,0.5,0.5,0.5,0.5]]
        markers = [['','','','',''], ['o','^']]
        title = '{0} K, {1} Torr, MS Experiment {2}'.format(int(T[i]), int(P[i]), int(i+1))

        comparePlot(data=simulated_data_relative_tdep, filename=os.path.join(output_directory, 'relative_tdep_78.png'),
                                                              otherData=experiment_data_relative_tdep,
                    linecolors = linecolors, linewidths = linewidths, markers = markers, title=title, ylimits=(-0.1, 1.1), legendloc = 'best')
        ########################################################################################################################

        ########################################################################################################################
        # Plot summed 91 amu relative to time-dependent styrene
        simulated_major_species = [
                                'Summed Simulated 91',
                                'C8H8(114)_obs',
                                ]

        experiment_major_species = [
                                 '91',
                                 '104',
                                 ]

        simulated_data_absolute_scaled_major_species = copy.deepcopy(simulated_data_absolute_scaled)
        experiment_data_major_species = copy.deepcopy(experiment_data)

        simulated_data_absolute_scaled_major_species.species = dict(zip(simulated_major_species, simulated_major_species))
        experiment_data_major_species.species = dict(zip(experiment_major_species, experiment_major_species))

        # Sort the profiles so that they will appear in a predictable order in plot
        sort_profiles(simulated_major_species, simulated_data_absolute_scaled_major_species)
        sort_profiles(experiment_major_species, experiment_data_major_species)

        #Normalize all major species signals to t-dependent styrene signal
        yVar_relative_scaled = []

        for dataobject in simulated_data_absolute_scaled_major_species.yVar:
            dataobject_relative_scaled = copy.deepcopy(dataobject)
            dataobject_relative_scaled.data /= simulated_data_absolute_scaled_major_species.yVar[1].data
            yVar_relative_scaled.append(dataobject_relative_scaled)

        simulated_data_relative_tdep = SimulationPlot(xVar=simulated_data_absolute_scaled_major_species.xVar, yVar=yVar_relative_scaled,
                                                        ylabel='Signal Relative to Styrene (104 amu)', species=simulated_major_species)

        yVar_relative_scaled = []

        for dataobject in experiment_data_major_species.yVar:
            dataobject_relative_scaled = copy.deepcopy(dataobject)
            dataobject_relative_scaled.data /= experiment_data_major_species.yVar[1].data
            yVar_relative_scaled.append(dataobject_relative_scaled)

        experiment_data_relative_tdep = SimulationPlot(xVar=experiment_data_major_species.xVar, yVar=yVar_relative_scaled,
                                                         ylabel='Signal Relative to Styrene (104 amu)', species=experiment_major_species)

        # Sort the profiles so that they will appear in a predictable order in plot
        # sort_profiles(simulated_major_species, simulated_data_relative_tdep)
        # sort_profiles(experiment_major_species, experiment_data_relative_tdep)

        #Define appearance of plot
        linecolors = [['g','r'],['g','r']]
        linewidths = [[2,2,2,2,2], [0.5,0.5,0.5,0.5,0.5]]
        markers = [['','','','',''], ['s','^']]
        title = '{0} K, {1} Torr, MS Experiment {2}'.format(int(T[i]), int(P[i]), int(i+1))

        comparePlot(data=simulated_data_relative_tdep, filename=os.path.join(output_directory, 'relative_tdep_91.png'),
                                                              otherData=experiment_data_relative_tdep,
                    linecolors = linecolors, linewidths = linewidths, markers = markers, title=title, ylimits=(-0.1, 1.1), legendloc = 'best')
        ########################################################################################################################

        ########################################################################################################################
        # Plot summed 118 amu relative to time-dependent styrene
        simulated_major_species = [
                                'C8H8(114)_obs',
                                'Summed Simulated 118',
                                ]

        experiment_major_species = [
                                 '104',
                                 '118',
                                 ]

        simulated_data_absolute_scaled_major_species = copy.deepcopy(simulated_data_absolute_scaled)
        experiment_data_major_species = copy.deepcopy(experiment_data)

        simulated_data_absolute_scaled_major_species.species = dict(zip(simulated_major_species, simulated_major_species))
        experiment_data_major_species.species = dict(zip(experiment_major_species, experiment_major_species))

        # Sort the profiles so that they will appear in a predictable order in plot
        sort_profiles(simulated_major_species, simulated_data_absolute_scaled_major_species)
        sort_profiles(experiment_major_species, experiment_data_major_species)

        #Normalize all major species signals to t-dependent styrene signal
        yVar_relative_scaled = []

        for dataobject in simulated_data_absolute_scaled_major_species.yVar:
            dataobject_relative_scaled = copy.deepcopy(dataobject)
            dataobject_relative_scaled.data /= simulated_data_absolute_scaled_major_species.yVar[0].data
            yVar_relative_scaled.append(dataobject_relative_scaled)

        simulated_data_relative_tdep = SimulationPlot(xVar=simulated_data_absolute_scaled_major_species.xVar, yVar=yVar_relative_scaled,
                                                        ylabel='Signal Relative to Styrene (104 amu)', species=simulated_major_species)

        yVar_relative_scaled = []

        for dataobject in experiment_data_major_species.yVar:
            dataobject_relative_scaled = copy.deepcopy(dataobject)
            dataobject_relative_scaled.data /= experiment_data_major_species.yVar[0].data
            yVar_relative_scaled.append(dataobject_relative_scaled)

        experiment_data_relative_tdep = SimulationPlot(xVar=experiment_data_major_species.xVar, yVar=yVar_relative_scaled,
                                                         ylabel='Signal Relative to Styrene (104 amu)', species=experiment_major_species)

        # Sort the profiles so that they will appear in a predictable order in plot
        # sort_profiles(simulated_major_species, simulated_data_relative_tdep)
        # sort_profiles(experiment_major_species, experiment_data_relative_tdep)

        #Define appearance of plot
        linecolors = [['r','c'],['r','c']]
        linewidths = [[2,2,2,2,2], [0.5,0.5,0.5,0.5,0.5]]
        markers = [['','','','',''], ['^','*']]
        title = '{0} K, {1} Torr, MS Experiment {2}'.format(int(T[i]), int(P[i]), int(i+1))

        comparePlot(data=simulated_data_relative_tdep, filename=os.path.join(output_directory, 'relative_tdep_118.png'),
                                                              otherData=experiment_data_relative_tdep,
                    linecolors = linecolors, linewidths = linewidths, markers = markers, title=title, ylimits=(-0.1, 1.1), legendloc = 'best')
        ########################################################################################################################

        ########################################################################################################################
        # Plot summed 119 amu relative to time-dependent styrene
        simulated_major_species = [
                                'C8H8(114)_obs',
                                'Summed Simulated 119',
                                ]

        experiment_major_species = [
                                 '104',
                                 '119',
                                 ]

        simulated_data_absolute_scaled_major_species = copy.deepcopy(simulated_data_absolute_scaled)
        experiment_data_major_species = copy.deepcopy(experiment_data)

        simulated_data_absolute_scaled_major_species.species = dict(zip(simulated_major_species, simulated_major_species))
        experiment_data_major_species.species = dict(zip(experiment_major_species, experiment_major_species))

        # Sort the profiles so that they will appear in a predictable order in plot
        sort_profiles(simulated_major_species, simulated_data_absolute_scaled_major_species)
        sort_profiles(experiment_major_species, experiment_data_major_species)

        #Normalize all major species signals to t-dependent styrene signal
        yVar_relative_scaled = []

        for dataobject in simulated_data_absolute_scaled_major_species.yVar:
            dataobject_relative_scaled = copy.deepcopy(dataobject)
            dataobject_relative_scaled.data /= simulated_data_absolute_scaled_major_species.yVar[0].data
            yVar_relative_scaled.append(dataobject_relative_scaled)

        simulated_data_relative_tdep = SimulationPlot(xVar=simulated_data_absolute_scaled_major_species.xVar, yVar=yVar_relative_scaled,
                                                        ylabel='Signal Relative to Styrene (104 amu)', species=simulated_major_species)

        yVar_relative_scaled = []

        for dataobject in experiment_data_major_species.yVar:
            dataobject_relative_scaled = copy.deepcopy(dataobject)
            dataobject_relative_scaled.data /= experiment_data_major_species.yVar[0].data
            yVar_relative_scaled.append(dataobject_relative_scaled)

        experiment_data_relative_tdep = SimulationPlot(xVar=experiment_data_major_species.xVar, yVar=yVar_relative_scaled,
                                                         ylabel='Signal Relative to Styrene (104 amu)', species=experiment_major_species)

        # Sort the profiles so that they will appear in a predictable order in plot
        # sort_profiles(simulated_major_species, simulated_data_relative_tdep)
        # sort_profiles(experiment_major_species, experiment_data_relative_tdep)

        #Define appearance of plot
        linecolors = [['r','m'],['r','m']]
        linewidths = [[2,2,2,2,2], [0.5,0.5,0.5,0.5,0.5]]
        markers = [['','','','',''], ['^','D']]
        title = '{0} K, {1} Torr, MS Experiment {2}'.format(int(T[i]), int(P[i]), int(i+1))

        comparePlot(data=simulated_data_relative_tdep, filename=os.path.join(output_directory, 'relative_tdep_119.png'),
                                                              otherData=experiment_data_relative_tdep,
                    linecolors = linecolors, linewidths = linewidths, markers = markers, title=title, ylimits=(-0.1, 1.1), legendloc = 'best')
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

        #Define appearance of plot
        linecolors = [['k'],['k']]
        linewidths = [[2], [0.5]]
        markers = [[''], ['o']]
        title = '{0} K, {1} Torr, MS Experiment {2}'.format(int(T[i]), int(P[i]), int(i+1))

        comparePlot(data=simulated_data_absolute_scaled_15_species, filename=os.path.join(output_directory, 'absolute_15.png'),
                                                              otherData=experiment_data_15_species,
                    linecolors = linecolors, linewidths = linewidths, markers = markers, title=title)

        comparePlot(data=simulated_data_relative_scaled_15_species, filename=os.path.join(output_directory, 'relative_15.png'),
                                                              otherData=experiment_data_relative_scaled_15_species,
                    linecolors=linecolors, linewidths=linewidths, markers=markers, title=title)

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

        #Define appearance of plot
        linecolors = [['b','g'],['b','g']]
        linewidths = [[2,2], [0.5,0.5]]
        markers = [['',''], ['o','s']]
        title = '{0} K, {1} Torr, MS Experiment {2}'.format(int(T[i]), int(P[i]), int(i+1))


        comparePlot(data=simulated_data_absolute_scaled_41_species, filename=os.path.join(output_directory, 'absolute_41.png'),
                                                              otherData=experiment_data_41_species,
                    linecolors = linecolors, linewidths = linewidths, markers = markers, title=title, legendloc = 'upper right')

        comparePlot(data=simulated_data_relative_scaled_41_species, filename=os.path.join(output_directory, 'relative_41.png'),
                                                              otherData=experiment_data_relative_scaled_41_species,
                    linecolors=linecolors, linewidths=linewidths, markers=markers, title=title, legendloc = 'upper right')

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

        #Define appearance of plot
        linecolors = [['b','g','r','c','m'],['b','g','r','c','m']]
        linewidths = [[2,2,2,2,2], [0.5,0.5,0.5,0.5,0.5]]
        markers = [['','','','',''], ['o','s','^','*','D']]
        title = '{0} K, {1} Torr, MS Experiment {2}'.format(int(T[i]), int(P[i]), int(i+1))

        comparePlot(data=simulated_data_absolute_scaled_77_species, filename=os.path.join(output_directory, 'absolute_77.png'),
                                                              otherData=experiment_data_77_species,
                    linecolors = linecolors, linewidths = linewidths, markers = markers, title=title, legendloc = 'upper right')

        comparePlot(data=simulated_data_relative_scaled_77_species, filename=os.path.join(output_directory, 'relative_77.png'),
                                                              otherData=experiment_data_relative_scaled_77_species,
                    linecolors=linecolors, linewidths=linewidths, markers=markers, title=title, legendloc = 'upper right')

        ########################################################################################################################

        ########################################################################################################################
        # Plot 77, 78 and 154 together
        simulated_no_propene_species = [
                                'C6H5(1)_obs',
                                'Summed Simulated 78',
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

        #Define appearance of plot
        linecolors = [['b','g','r','c','m'],['b','g','r','c','m']]
        linewidths = [[2,2,2,2,2], [0.5,0.5,0.5,0.5,0.5]]
        markers = [['','','','',''], ['o','s','^','*','D']]
        title = '{0} K, {1} Torr, MS Experiment {2}'.format(int(T[i]), int(P[i]), int(i+1))

        comparePlot(data=simulated_data_absolute_scaled_no_propene_species, filename=os.path.join(output_directory, 'absolute_noC3H6.png'),
                                                              otherData=experiment_data_no_propene_species,
                    linecolors = linecolors, linewidths = linewidths, markers = markers, title=title, legendloc = 'upper right')

        comparePlot(data=simulated_data_relative_scaled_no_propene_species, filename=os.path.join(output_directory, 'relative_noC3H6.png'),
                                                              otherData=experiment_data_relative_scaled_no_propene_species,
                    linecolors=linecolors, linewidths=linewidths, markers=markers, title=title, legendloc = 'upper right')
        ########################################################################################################################

        ########################################################################################################################
        # Plot 78
        simulated_78_species = ['Summed Simulated 78'] + species_to_sum_dictionary['Summed Simulated 78']

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

        #Define appearance of plot
        linecolors = [['b','k','g','r','c','m'],['b','k','g','r','c','m']]
        linewidths = [[2,1,1,1,1,1], [0.5,0.5,0.5,0.5,0.5,0.5]]
        markers = [['','','','',''], ['o','o','s','^','*','D']]
        title = '{0} K, {1} Torr, MS Experiment {2}'.format(int(T[i]), int(P[i]), int(i+1))

        comparePlot(data=simulated_data_absolute_scaled_78_species, filename=os.path.join(output_directory, 'absolute_78.png'),
                                                              otherData=experiment_data_78_species,
                    linecolors = linecolors, linewidths = linewidths, markers = markers, title=title)

        comparePlot(data=simulated_data_relative_scaled_78_species, filename=os.path.join(output_directory, 'relative_78.png'),
                                                              otherData=experiment_data_relative_scaled_78_species,
                    linecolors=linecolors, linewidths=linewidths, markers=markers, title=title)

        ########################################################################################################################

        ########################################################################################################################
        # Plot 78 and 104 amu together (only major species that don't have fragmentation/C13 contributions)
        simulated_nofrag_species = ['Summed Simulated 78',
                                'C8H8(114)_obs',
                                 ]

        experiment_nofrag_species = ['78',
                                '104',
                                ]

        simulated_data_absolute_scaled_nofrag_species = copy.deepcopy(simulated_data_absolute_scaled)
        simulated_data_relative_scaled_nofrag_species = copy.deepcopy(simulated_data_relative_scaled)
        experiment_data_nofrag_species = copy.deepcopy(experiment_data)
        experiment_data_relative_scaled_nofrag_species = copy.deepcopy(experiment_data_relative_scaled)

        simulated_data_absolute_scaled_nofrag_species.species = dict(zip(simulated_nofrag_species, simulated_nofrag_species))
        simulated_data_relative_scaled_nofrag_species.species = dict(zip(simulated_nofrag_species, simulated_nofrag_species))
        experiment_data_nofrag_species.species = dict(zip(experiment_nofrag_species, experiment_nofrag_species))
        experiment_data_relative_scaled_nofrag_species.species = dict(zip(experiment_nofrag_species, experiment_nofrag_species))

        # Sort the profiles so that they will appear in a predictable order in plot
        sort_profiles(simulated_nofrag_species, simulated_data_absolute_scaled_nofrag_species)
        sort_profiles(simulated_nofrag_species, simulated_data_relative_scaled_nofrag_species)
        sort_profiles(experiment_nofrag_species, experiment_data_nofrag_species)
        sort_profiles(experiment_nofrag_species, experiment_data_relative_scaled_nofrag_species)

        #Define appearance of plot
        linecolors = [['b','r'],['b','r']]
        linewidths = [[2,2], [0.5,0.5]]
        markers = [['',''], ['o','^']]
        title = '{0} K, {1} Torr, MS Experiment {2}'.format(int(T[i]), int(P[i]), int(i+1))

        comparePlot(data=simulated_data_absolute_scaled_nofrag_species, filename=os.path.join(output_directory, 'absolute_nofrag.png'),
                                                              otherData=experiment_data_nofrag_species,
                    linecolors = linecolors, linewidths = linewidths, markers = markers, title=title)

        comparePlot(data=simulated_data_relative_scaled_nofrag_species, filename=os.path.join(output_directory, 'relative_nofrag.png'),
                                                              otherData=experiment_data_relative_scaled_nofrag_species,
                    linecolors=linecolors, linewidths=linewidths, markers=markers, title=title)
        ########################################################################################################################

        ########################################################################################################################
        #Plot 91

        simulated_91_species = ['Summed Simulated 91'] + species_to_sum_dictionary['Summed Simulated 91']

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

        #Define appearance of plot
        linecolors = [['g','b','k','r','c','m'],['g','b','k','r','c','m']]
        linewidths = [[2,1,1,1,1,1], [0.5,0.5,0.5,0.5,0.5,0.5]]
        markers = [['','','','',''], ['s','o','s','^','*','D']]
        title = '{0} K, {1} Torr, MS Experiment {2}'.format(int(T[i]), int(P[i]), int(i+1))

        comparePlot(data=simulated_data_absolute_scaled_91_species, filename=os.path.join(output_directory, 'absolute_91.png'),
                                                              otherData=experiment_data_91_species,
                    linecolors = linecolors, linewidths = linewidths, markers = markers, title=title)

        comparePlot(data=simulated_data_relative_scaled_91_species, filename=os.path.join(output_directory, 'relative_91.png'),
                                                              otherData=experiment_data_relative_scaled_91_species,
                    linecolors=linecolors, linewidths=linewidths, markers=markers, title=title)

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

        #Define appearance of plot
        linecolors = [['r'],['r']]
        linewidths = [[2], [0.5]]
        markers = [[''], ['^']]
        title = '{0} K, {1} Torr, MS Experiment {2}'.format(int(T[i]), int(P[i]), int(i+1))

        comparePlot(data=simulated_data_absolute_scaled_104_species, filename=os.path.join(output_directory, 'absolute_104.png'),
                                                              otherData=experiment_data_104_species,
                    linecolors = linecolors, linewidths = linewidths, markers = markers, title=title)

        comparePlot(data=simulated_data_relative_scaled_104_species, filename=os.path.join(output_directory, 'relative_104.png'),
                                                              otherData=experiment_data_relative_scaled_104_species,
                    linecolors=linecolors, linewidths=linewidths, markers=markers, title=title)

        ########################################################################################################################

        ########################################################################################################################
        #Plot 118
        simulated_118_species = ['Summed Simulated 118'] + species_to_sum_dictionary['Summed Simulated 118']

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

        #Define appearance of plot
        linecolors = [['c','b','g','r','k','m'],['c','b','g','r','k','m']]
        linewidths = [[2,1,1,1,1,1], [0.5,0.5,0.5,0.5,0.5,0.5]]
        markers = [['','','','',''], ['*','o','s','^','*','D']]
        title = '{0} K, {1} Torr, MS Experiment {2}'.format(int(T[i]), int(P[i]), int(i+1))

        comparePlot(data=simulated_data_absolute_scaled_118_species, filename=os.path.join(output_directory, 'absolute_118.png'),
                                                              otherData=experiment_data_118_species,
                    linecolors = linecolors, linewidths = linewidths, markers = markers, title=title)

        comparePlot(data=simulated_data_relative_scaled_118_species, filename=os.path.join(output_directory, 'relative_118.png'),
                                                              otherData=experiment_data_relative_scaled_118_species,
                    linecolors=linecolors, linewidths=linewidths, markers=markers, title=title)

        ########################################################################################################################

        ########################################################################################################################
        #Plot 119
        simulated_119_species = ['Summed Simulated 119'] + species_to_sum_dictionary['Summed Simulated 119']

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

        #Define appearance of plot
        linecolors = [['m','b','g','r','c','k','gold','orange','brown'],['m','b','g','r','c','k']]
        linewidths = [[2,1,1,1,1,1,1,1,1], [0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5]]
        markers = [['','','','','','','',''], ['D','o','s','^','*','D','^','*','D']]
        title = '{0} K, {1} Torr, MS Experiment {2}'.format(int(T[i]), int(P[i]), int(i+1))

        comparePlot(data=simulated_data_absolute_scaled_119_species, filename=os.path.join(output_directory, 'absolute_119.png'),
                                                              otherData=experiment_data_119_species,
                    linecolors = linecolors, linewidths = linewidths, markers = markers, title=title)

        comparePlot(data=simulated_data_relative_scaled_119_species, filename=os.path.join(output_directory, 'relative_119.png'),
                                                              otherData=experiment_data_relative_scaled_119_species,
                    linecolors=linecolors, linewidths=linewidths, markers=markers, title=title)

        ########################################################################################################################

        ########################################################################################################################
        #Plot 120
        simulated_120_species = ['Summed Simulated 120'] + species_to_sum_dictionary['Summed Simulated 120']

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

        #Define appearance of plot
        linecolors = [['k','b','g','r','c','m'],['k','b','g','r','c','m']]
        linewidths = [[2,1,1,1,1,1], [0.5,0.5,0.5,0.5,0.5,0.5]]
        markers = [['','','','',''], ['o','o','s','^','*','D']]
        title = '{0} K, {1} Torr, MS Experiment {2}'.format(int(T[i]), int(P[i]), int(i+1))

        comparePlot(data=simulated_data_absolute_scaled_120_species, filename=os.path.join(output_directory, 'absolute_120.png'),
                                                              otherData=experiment_data_120_species,
                    linecolors = linecolors, linewidths = linewidths, markers = markers, title=title) #, ylimits=(-0.002, 0.012))

        comparePlot(data=simulated_data_relative_scaled_120_species, filename=os.path.join(output_directory, 'relative_120.png'),
                                                              otherData=experiment_data_relative_scaled_120_species,
                    linecolors=linecolors, linewidths=linewidths, markers=markers, title=title)

        ########################################################################################################################

        ########################################################################################################################
        #Plot remaining 119 bimolecular products together (134, 160 and 238 amu)
        simulated_119_bi_products_species = ['S(298)',
                                             'mz_160_obs',
                                             # 'C18H22(77)',
                                             ]

        experiment_119_bi_products_species = ['134',
                                              '160',
                                              # '238',
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

        #Define appearance of plot
        linecolors = [['b','g','r','c','m'],['b','g','r','c','m']]
        linewidths = [[2,2,2,2,2], [0.5,0.5,0.5,0.5,0.5]]
        markers = [['','','','',''], ['o','s','^','*','D']]
        title = '{0} K, {1} Torr, MS Experiment {2}'.format(int(T[i]), int(P[i]), int(i+1))

        comparePlot(data=simulated_data_absolute_scaled_119_bi_products_species, filename=os.path.join(output_directory, 'absolute_119_bi.png'),
                                                              otherData=experiment_data_119_bi_products_species,
                    linecolors = linecolors, linewidths = linewidths, markers = markers, title=title) #, ylimits=(-0.001, 0.0035))

        comparePlot(data=simulated_data_relative_scaled_119_bi_products_species, filename=os.path.join(output_directory, 'relative_119_bi.png'),
                                                              otherData=experiment_data_relative_scaled_119_bi_products_species,
                    linecolors=linecolors, linewidths=linewidths, markers=markers, title=title)

        ########################################################################################################################

        ########################################################################################################################
        #Plot I, HI together
        simulated_IHI_species = [
                                'I_obs',
                               'HI_obs',
                               ]

        experiment_IHI_species = [
                                '127',
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

        #Define appearance of plot
        linecolors = [['k','b','g','r','c','m'],['k','b','g','r','c','m']]
        linewidths = [[2,2,2,2,2,2], [0.5,0.5,0.5,0.5,0.5,0.5]]
        markers = [['','','','','',''], ['o','s','^','*','D']]
        title = '{0} K, {1} Torr, MS Experiment {2}'.format(int(T[i]), int(P[i]), int(i+1))

        comparePlot(data=simulated_data_absolute_scaled_IHI_species, filename=os.path.join(output_directory, 'absolute_IHI.png'),
                                                              otherData=experiment_data_IHI_species,
                    linecolors = linecolors, linewidths = linewidths, markers = markers, title=title, legendloc = 'upper right')

        comparePlot(data=simulated_data_relative_scaled_IHI_species, filename=os.path.join(output_directory, 'relative_IHI.png'),
                                                              otherData=experiment_data_relative_scaled_IHI_species,
                    linecolors=linecolors, linewidths=linewidths, markers=markers, title=title, legendloc = 'upper right')

        ########################################################################################################################

        ########################################################################################################################
        #Plot I by itself
        simulated_I_species = ['I_obs',
                               'I'
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

        #Pad t<0 with zeros
        simulated_data_absolute_scaled_I_species.xVar.data = np.insert(simulated_data_absolute_scaled_I_species.xVar.data, 0, [-5, 0.0])
        for dataobject in simulated_data_absolute_scaled_I_species.yVar:
            dataobject.data = np.insert(dataobject.data, 0, [0.0, 0.0])

        simulated_data_relative_scaled_I_species.xVar.data = np.insert(simulated_data_relative_scaled_I_species.xVar.data, 0, [-5, 0.0])
        for dataobject in simulated_data_relative_scaled_I_species.yVar:
            dataobject.data = np.insert(dataobject.data, 0, [0.0, 0.0])

        #Define appearance of plot
        linecolors = [['k','b','g','r','c','m'],['k','b','g','r','c','m']]
        linewidths = [[2,2,2,2,2,2], [0.5,0.5,0.5,0.5,0.5,0.5]]
        markers = [['','','','','',''], ['o','s','^','*','D']]
        title = '{0} K, {1} Torr, MS Experiment {2}'.format(int(T[i]), int(P[i]), int(i+1))

        comparePlot(data=simulated_data_absolute_scaled_I_species, filename=os.path.join(output_directory, 'absolute_I.png'),
                                                              otherData=experiment_data_I_species,
                    linecolors = linecolors, linewidths = linewidths, markers = markers, legendloc = 'upper right',) # title=title, )

        comparePlot(data=simulated_data_relative_scaled_I_species, filename=os.path.join(output_directory, 'relative_I.png'),
                                                              otherData=experiment_data_relative_scaled_I_species,
                    linecolors=linecolors, linewidths=linewidths, markers=markers, title=title, legendloc = 'upper right')

        ########################################################################################################################

        ########################################################################################################################
        # Plot Alkyl Iodide species
        simulated_alkyl_I_species = ['CH3I_obs',
                                     'C3H5I_obs_frag_168',
                                     # 'C7H7I_obs_frag_218',
                                     'C9H11I_obs_frag_246',
                                     ]

        experiment_alkyl_I_species = ['142',
                                      '168',
                                      # '218',
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

        #Define appearance of plot
        linecolors = [['b','g','r','c','m'],['b','g','r','c','m']]
        linewidths = [[2,2,2,2,2], [0.5,0.5,0.5,0.5,0.5]]
        markers = [['','','','',''], ['o','s','^','*','D']]
        title = '{0} K, {1} Torr, MS Experiment {2}'.format(int(T[i]), int(P[i]), int(i+1))

        comparePlot(data=simulated_data_absolute_scaled_alkyl_I_species, filename=os.path.join(output_directory, 'absolute_alkyl_I.png'),
                                                              otherData=experiment_data_alkyl_I_species,
                    linecolors = linecolors, linewidths = linewidths, markers = markers, title=title)

        comparePlot(data=simulated_data_relative_scaled_alkyl_I_species, filename=os.path.join(output_directory, 'relative_alkyl_I.png'),
                                                              otherData=experiment_data_relative_scaled_alkyl_I_species,
                    linecolors=linecolors, linewidths=linewidths, markers=markers, title=title)

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

        #Define appearance of plot
        linecolors = [['k','b','g','r','c','m'],['k','b','g','r','c','m']]
        linewidths = [[2,2,2,2,2,2], [0.5,0.5,0.5,0.5,0.5,0.5]]
        markers = [['','','','','',''], ['o','s','^','*','D']]
        title = '{0} K, {1} Torr, MS Experiment {2}'.format(int(T[i]), int(P[i]), int(i+1))

        comparePlot(data=simulated_data_absolute_scaled_I2_species, filename=os.path.join(output_directory, 'absolute_I2.png'),
                                                              otherData=experiment_data_I2_species,
                    linecolors = linecolors, linewidths = linewidths, markers = markers, title=title)

        comparePlot(data=simulated_data_relative_scaled_I2_species, filename=os.path.join(output_directory, 'relative_I2.png'),
                                                              otherData=experiment_data_relative_scaled_I2_species,
                    linecolors=linecolors, linewidths=linewidths, markers=markers, title=title)

        ########################################################################################################################

        ########################################################################################################################
        # Plot Normalized phenyl and Absorbance together
        simulated_77_species = ['C6H5(1)',
                                ]

        simulated_data_normalized_77_species = copy.deepcopy(simulated_data_normalized)
        simulated_data_normalized_77_species = copy.deepcopy(simulated_data_normalized)

        simulated_data_normalized_77_species.species = dict(zip(simulated_77_species, simulated_77_species))

        # Sort the profiles so that they will appear in a predictable order in plot
        sort_profiles(simulated_77_species, simulated_data_normalized_77_species)

        # Define appearance of plot
        linecolors = [['k', 'g', 'r', 'c', 'm'], ['r', 'c', 'm']]
        linewidths = [[0.5, 0.5, 0.5, 0.5, 0.5], [2, 2, 2]]
        markers = [['o', 's', '^', '*', 'D'], ['', '', '', '', '']]
        title = '{0} K, {1} Torr, {2} nm Absorbance Experiment {3}'.format(int(T[i]), int(P[i]), Wavelength[i], int(i + 1))

        comparePlot(data=abs_experiment_data_normalized, filename=os.path.join(output_directory, '505_abs.png'),
                    otherData=simulated_data_normalized_77_species,
                    linecolors=linecolors, linewidths=linewidths, markers=markers, title=title, legendloc='upper right',
                    xlimits=(-0.1, 2), ylimits=(-.1, 1.1), legend_title='Simulated Species:', ylabel='A/A$_0$', legend_labels=['Phenyl Radical'])

        ########################################################################################################################

        ########################################################################################################################
        # Plot Normalized phenyl+allyl and Absorbance together
        simulated_phenyl_allyl_species = ['C6H5(1)',
                                          'C3H5(45)',
                                          'combined_phenyl_allyl',
                                          ]

        simulated_data_combined_phenyl_allyl_species = copy.deepcopy(simulated_data_normalized)
        simulated_data_combined_phenyl_allyl_species = copy.deepcopy(simulated_data_normalized)

        simulated_data_combined_phenyl_allyl_species.species = dict(
            zip(simulated_phenyl_allyl_species, simulated_phenyl_allyl_species))

        # Sort the profiles so that they will appear in a predictable order in plot
        sort_profiles(simulated_phenyl_allyl_species, simulated_data_combined_phenyl_allyl_species)

        # Define appearance of plot
        linecolors = [['k', 'g', 'r', 'c', 'm'], ['g', 'b', 'r']]
        linewidths = [[0.5, 0.5, 0.5, 0.5, 0.5], [2, 2, 2]]
        markers = [['o', 's', '^', '*', 'D'], ['-', '-', '', '', '']]
        title = '{0} K, {1} Torr, {2} nm Absorbance Experiment {3}'.format(int(T[i]), int(P[i]), Wavelength[i], int(i + 1))

        comparePlot(data=abs_experiment_data_normalized, filename=os.path.join(output_directory, '408_abs.png'),
                    otherData=simulated_data_combined_phenyl_allyl_species,
                    linecolors=linecolors, linewidths=linewidths, markers=markers, title=title, legendloc='upper right',
                    xlimits=(-0.1, 2), ylimits=(-.1, 1.1), legend_title='Simulated Species:', ylabel='$\mathregular{A/A_0}$', legend_labels=['Phenyl Radical', 'Allyl Radical', 'Total'])

        ########################################################################################################################

        ########################################################################################################################
        # Plot Normalized phenyl+benzyl and Absorbance together

        simulated_phenyl_benzyl_species = ['C6H5(1)',
                                           'C7H7(190)',
                                           'combined_phenyl_benzyl',
                                           ]

        simulated_data_combined_phenyl_benzyl_species = copy.deepcopy(simulated_data_normalized)
        simulated_data_combined_phenyl_benzyl_species = copy.deepcopy(simulated_data_normalized)

        simulated_data_combined_phenyl_benzyl_species.species = dict(
            zip(simulated_phenyl_benzyl_species, simulated_phenyl_benzyl_species))

        # Sort the profiles so that they will appear in a predictable order in plot
        sort_profiles(simulated_phenyl_benzyl_species, simulated_data_combined_phenyl_benzyl_species)

        # Define appearance of plot
        linecolors = [['k', 'g', 'r', 'c', 'm'], ['g', 'b', 'r']]
        linewidths = [[0.5, 0.5, 0.5, 0.5, 0.5], [2, 2, 2]]
        markers = [['o', 's', '^', '*', 'D'], ['-', '-', '', '', '']]
        title = '{0} K, {1} Torr, {2} nm Absorbance Experiment {3}'.format(int(T[i]), int(P[i]), Wavelength[i], int(i + 1))

        comparePlot(data=abs_experiment_data_normalized, filename=os.path.join(output_directory, '447_abs.png'),
                    otherData=simulated_data_combined_phenyl_benzyl_species,
                    linecolors=linecolors, linewidths=linewidths, markers=markers, title=title, legendloc='upper right',
                    xlimits=(-0.1, 2), ylimits=(-.1, 1.1), legend_title='Simulated Species:', ylabel='$\mathregular{A/A_0}$', legend_labels=['Phenyl Radical', 'Benzyl Radical', 'Total'])

        ########################################################################################################################

    else:
        ########################################################################################################################
        # Load simulated and experimental data

        output_directory = os.path.join(os.path.dirname(simulation_csv_file[i]),
                                                    '{0}_K_{1}_Torr_{2}_C3H6_{3}_C6H5I_{4}_photolysis_{5}_nm'.format(int(T[i]), int(P[i]),int(C3H6_Flow[i]),
                                                                                                     int(C6H5I_Flow[i]),
                                                                                                     int(Photolysis_Power[i]), int(Wavelength[i])))

        if not os.path.exists(output_directory):
            os.makedirs(output_directory)

        simulated_data = SimulationPlot(ylabel='Mole Fraction', csvFile=simulation_csv_file[i],
                                        species=simulated_species_to_plot)

        abs_experiment_data = SimulationPlot(ylabel='Absorbance', csvFile=abs_experiment_csv_file[i],
                                             species={'Intensity': 'Intensity'})

        simulated_data.load()

        abs_experiment_data.load()

        ########################################################################################################################

        ########################################################################################################################
        # Switch x axis from seconds to milliseconds
        simulated_data.xVar.data *= 1000
        simulated_data.xVar.label = 'Time'
        simulated_data.xVar.units = 'ms'

        abs_experiment_data.xVar.data *= 1000
        abs_experiment_data.xVar.label = 'Time'
        abs_experiment_data.xVar.units = 'ms'
        ########################################################################################################################

        ########################################################################################################################
        # Normalize simulation and absorbance experiments to initial Phenyl concentration for comparison

        # Normalized simulated data
        yVar_normalized = []

        for dataobject in simulated_data.yVar:
            if dataobject.label == 'C6H5(1)':
                C6H5_0 = dataobject.data[0]
                break

        for dataobject in simulated_data.yVar:
            dataobject_normalized = copy.deepcopy(dataobject)
            dataobject_normalized.data /= C6H5_0
            yVar_normalized.append(dataobject_normalized)

        simulated_data_normalized = SimulationPlot(xVar=simulated_data.xVar, yVar=yVar_normalized,
                                                   ylabel='Normalized Concentration', species=simulated_species_to_plot)

        # Add a new data series that is the normalized phenyl + allyl simulated data
        sum_species(['C3H5(45)', 'C6H5(1)'], simulated_data_normalized, 'combined_phenyl_allyl', weights = Absorbance_weights_408[round(T[i],-2)])

        # Add a new data series that is the normalized phenyl + benzyl simulated data
        sum_species(['C7H7(190)', 'C6H5(1)'], simulated_data_normalized, 'combined_phenyl_benzyl',
                    weights=Absorbance_weights_447[round(T[i], -2)])

        # Normalize Absorbance experiments

        yVar_normalized = []

        for dataobject in abs_experiment_data.yVar:
            dataobject_normalized = copy.deepcopy(dataobject)
            dataobject_normalized.data = np.log(I0[i] * 1000 / (I0[i] * 1000 + dataobject_normalized.data))
            linear_fit = np.polyfit(abs_experiment_data.xVar.data[np.where(
                (abs_experiment_data.xVar.data >= 0.01) & (abs_experiment_data.xVar.data <= 0.1))],
                                    dataobject_normalized.data[np.where((abs_experiment_data.xVar.data >= 0.01) & (
                                    abs_experiment_data.xVar.data <= 0.1))], 1)
            dataobject_normalized.data /= np.polyval(linear_fit, 0)
            # Do a moving average to smooth data
            dataobject_normalized.data = np.convolve(dataobject_normalized.data, np.ones((5,)) / 5, mode='same')
            yVar_normalized.append(dataobject_normalized)

        abs_experiment_data_normalized = SimulationPlot(xVar=abs_experiment_data.xVar, yVar=yVar_normalized,
                                                        ylabel='Normalized Absorbance')

        ########################################################################################################################

        ########################################################################################################################
        # Plot Normalized phenyl and Absorbance together
        simulated_77_species = ['C6H5(1)',
                                ]

        simulated_data_normalized_77_species = copy.deepcopy(simulated_data_normalized)
        simulated_data_normalized_77_species = copy.deepcopy(simulated_data_normalized)

        simulated_data_normalized_77_species.species = dict(zip(simulated_77_species, simulated_77_species))

        # Sort the profiles so that they will appear in a predictable order in plot
        sort_profiles(simulated_77_species, simulated_data_normalized_77_species)

        # Define appearance of plot
        linecolors = [['k', 'g', 'r', 'c', 'm'], ['r', 'c', 'm']]
        linewidths = [[0.5, 0.5, 0.5, 0.5, 0.5], [2, 2, 2]]
        markers = [['o', 's', '^', '*', 'D'], ['', '', '', '', '']]
        title = '{0} K, {1} Torr, {2} nm Absorbance Experiment {3}'.format(int(T[i]), int(P[i]), Wavelength[i], int(i + 1))

        comparePlot(data=abs_experiment_data_normalized, filename=os.path.join(output_directory, '505_abs.png'),
                    otherData=simulated_data_normalized_77_species,
                    linecolors=linecolors, linewidths=linewidths, markers=markers, title=title, legendloc='upper right',
                    xlimits=(-0.1, 2), ylimits=(-.1, 1.1), legend_title='Simulated Species:', ylabel='$\mathregular{A/A_0}$', legend_labels=['Phenyl Radical'])

        ########################################################################################################################

        ########################################################################################################################
        # Plot Normalized phenyl+allyl and Absorbance together
        simulated_phenyl_allyl_species = ['C6H5(1)',
                                          'C3H5(45)',
                                          'combined_phenyl_allyl',
                                          ]

        simulated_data_combined_phenyl_allyl_species = copy.deepcopy(simulated_data_normalized)
        simulated_data_combined_phenyl_allyl_species = copy.deepcopy(simulated_data_normalized)

        simulated_data_combined_phenyl_allyl_species.species = dict(
            zip(simulated_phenyl_allyl_species, simulated_phenyl_allyl_species))

        # Sort the profiles so that they will appear in a predictable order in plot
        sort_profiles(simulated_phenyl_allyl_species, simulated_data_combined_phenyl_allyl_species)

        # Define appearance of plot
        linecolors = [['k', 'g', 'r', 'c', 'm'], ['g', 'b', 'r']]
        linewidths = [[0.5, 0.5, 0.5, 0.5, 0.5], [2, 2, 2]]
        markers = [['o', 's', '^', '*', 'D'], ['-', '-', '', '', '']]
        title = '{0} K, {1} Torr, {2} nm Absorbance Experiment {3}'.format(int(T[i]), int(P[i]), Wavelength[i], int(i + 1))

        comparePlot(data=abs_experiment_data_normalized, filename=os.path.join(output_directory, '408_abs.png'),
                    otherData=simulated_data_combined_phenyl_allyl_species,
                    linecolors=linecolors, linewidths=linewidths, markers=markers, title=title, legendloc='upper right',
                    xlimits=(-0.1, 2), ylimits=(-.1, 1.1), legend_title='Simulated Species:', ylabel='$\mathregular{A/A_0}$', legend_labels=['Phenyl Radical', 'Allyl Radical', 'Total'])

        ########################################################################################################################

        ########################################################################################################################
        # Plot Normalized phenyl+benzyl and Absorbance together

        simulated_phenyl_benzyl_species = ['C6H5(1)',
                                          'C7H7(190)',
                                          'combined_phenyl_benzyl',
                                          ]

        simulated_data_combined_phenyl_benzyl_species = copy.deepcopy(simulated_data_normalized)
        simulated_data_combined_phenyl_benzyl_species = copy.deepcopy(simulated_data_normalized)

        simulated_data_combined_phenyl_benzyl_species.species = dict(
            zip(simulated_phenyl_benzyl_species, simulated_phenyl_benzyl_species))

        # Sort the profiles so that they will appear in a predictable order in plot
        sort_profiles(simulated_phenyl_benzyl_species, simulated_data_combined_phenyl_benzyl_species)

        # Define appearance of plot
        linecolors = [['k', 'g', 'r', 'c', 'm'], ['g', 'b', 'r']]
        linewidths = [[0.5, 0.5, 0.5, 0.5, 0.5], [2, 2, 2]]
        markers = [['o', 's', '^', '*', 'D'], ['-', '-', '', '', '']]
        title = '{0} K, {1} Torr, {2} nm Absorbance Experiment {3}'.format(int(T[i]), int(P[i]), Wavelength[i], int(i + 1))

        comparePlot(data=abs_experiment_data_normalized, filename=os.path.join(output_directory, '447_abs.png'),
                    otherData=simulated_data_combined_phenyl_benzyl_species,
                    linecolors=linecolors, linewidths=linewidths, markers=markers, title=title, legendloc='upper right',
                    xlimits=(-0.1, 2), ylimits=(-.1, 1.1), legend_title='Simulated Species:', ylabel='$\mathregular{A/A_0}$', legend_labels=['Phenyl Radical', 'Benzyl Radical', 'Total'])

        ########################################################################################################################