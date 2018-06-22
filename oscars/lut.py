from __future__ import print_function

import numpy as np
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt

class lut1d:
    """Class for 1D lookup tables for insertion devices as a function of gap"""

    def __init__(self, ifile=None, name=None):

        self.lut1d_filename = ifile
        self.name = name
        self.splines_gap_vs_energy = dict()
        self.splines_flux_vs_energy = dict()
        self.energy_range = dict()
        self.splines_energy_vs_gap = dict()
        self.splines_flux_vs_gap = dict()
        self.gap_range = dict()
        
        if self.lut1d_filename is not None:
            self.read_file_lut1d(self.lut1d_filename)
            
        return


    def clear (self):
        """Clear some internal data"""

        self.splines_gap_vs_energy = dict()
        self.splines_flux_vs_energy = dict()
        self.energy_range = dict()
        self.splines_energy_vs_gap = dict()
        self.splines_flux_vs_gap = dict()
        self.gap_range = dict()

        return

        
    def read_file_lut1d(self, ifile):
        """read a file and setup data accordingly for the lut1d.  The file format is a simple text file.  Comments may be included on any
           line by beginning the line with '#'.  Entries for the table begin with the 'harmonic' keyword followed by whitespace and the harmonic number.
           This is then followed by multiple lines each containing energy gap flux (whitespace separated) for as many points as desired.
           The units for these are arbitrary, but suggested are: [eV, mm, a.u.].  The lines do not need to be any any particular order, but
           must be grouped by harmonic number.
           For example this is a valid file:

           # A comment line - U42 planar undulator

           harmonic 1
            313.39 11.5 175950363894016.0
            536.12 16.0 266545674797824.0
            874.82 21.0 347655414526208.0
           1752.20 35.0 180941839137408.0
           2008.33 55.0  10694510247968.0

           harmonic 3
            941.47 11.5 212151858914304.0
           1360.22 14.5 252854000205824.0
           1890.43 17.5 262521650105344.0
           3089.79 23.0 184171317321728.0

        """
        
        self.clear()
        self.lut1d_filename = ifile

        with open(ifile) as fi:
            for line in fi:
                ll = line.rstrip()
                ll = ll.lstrip()
                if ll.startswith('#'):
                    continue
                if ll.lower().startswith('harmonic'):
                    if len(ll.split()) < 2:
                        raise IndexError('harmonic line definition has incorrect format:' + ll)
                    harmonic = int(ll.split()[1])

                    all_fields = []
                    for line in fi:
                        ll = line.rstrip()
                        ll = ll.lstrip()
                        if ll.startswith('#'):
                            continue
                        if ll == '':
                            break

                        fields = list(map(float, ll.split()))
                        if not (len(fields) == 3 or len(fields) == 2):
                            raise ValueError('Error in format of line:' + line)

                        # If table is without flux that's no problem, just add a zero
                        if len(fields) == 2:
                            fields.append(0)
                        all_fields.append(fields)

                    if len(all_fields) < 2:
                        continue

                    all_fields.sort(key=lambda x: x[0])

                    energy = [x[0] for x in all_fields]
                    gap = [x[1] for x in all_fields]
                    flux = [x[2] for x in all_fields]

                    #cs_gap_vs_energy = CubicSpline(energy, gap)
                    #cs_flux_vs_energy = CubicSpline(energy, flux)

                    # Forward and reverse splines
                    self.splines_gap_vs_energy[harmonic] = CubicSpline(energy, gap)
                    self.splines_flux_vs_energy[harmonic] = CubicSpline(energy, flux)
                    self.energy_range[harmonic] = [energy[0], energy[-1]]

                    self.splines_energy_vs_gap[harmonic] = CubicSpline(gap, energy)
                    self.splines_flux_vs_gap[harmonic] = CubicSpline(gap, flux)
                    self.gap_range[harmonic] = [gap[0], gap[-1]]
        return



    def summary (self, **kwargs):
        """Show the summary plots for this lut1d.  See get_gaps() for arguments"""

        kwargs['show'] = True
        self.get_gaps(**kwargs)

        return



    def get_gaps (self,
                  energy_eV=0,
                  gap=None,
                  show=False,
                  ofile=None,
                  name='',
                  harmonic_range=None,
                  odd=True,
                  even=False,
                  plots=['gap', 'flux'],
                  xscale='linear',
                  yscale='linear',
                  xlim=None,
                  gap_ylim=None,
                  flux_ylim=None,
                  grid=False
                 ):
        """Get the gaps at which to find the given energy

        Parameters
        ----------
        facility : str
            Name of the facility.  See below for currently supported.

        energy_eV : float
            Energy at which to find gaps and flux

        gap : float
            If you want to draw a line on the plot for a specific gap (does not change calculation)

        show : bool
            Set to true if you want to produce a plot

        ofile : str
            Name of the output file for the plot

        name : str
            Name for adding to plot

        harmonic_range : list[int, int]
            List containing the min and max harmonic to consider

        odd : bool
            Set to true if you want to see odd harmonics (default)

        even : bool
            Set to true if you want to see even harmonics (default is False)

        plots : list of str
            Can contain 'gap' and/or 'flux' depending on which (or both) plots you want to make.  both set as default.

        xscale : str
            Set the xscale to 'log' or 'linear'

        yscale : str
            Set the yscale to 'log' or 'linear'

        xlim : list[float, float]
            Range of the x-axis on the plot

        gap_ylim : list[float, float]
            Range of y-axis on the gap plot

        flux_ylim : list[float, float]
            Range of y-axis on the flux plot

        grid : bool
            Show gridlines or not (default not)

        """

        # Set name for plot
        if name == '' and self.name is not None:
            name = self.name
        
        # Grab the harmonics which cover the range
        harmonic_list = []
        for harmonic in self.energy_range:
            # Check for odd even
            if harmonic % 2 and not odd:
                continue
            if harmonic % 2 == 0 and not even:
                continue

            # Check input harmonic range and make sure this is within the range
            if energy_eV == 0 or (harmonic_range is not None and (harmonic < harmonic_range[0] or harmonic > harmonic_range[1])):
                continue
                
            if energy_eV >= self.energy_range[harmonic][0] and energy_eV <= self.energy_range[harmonic][1]:
                harmonic_list.append(harmonic)

        # Get results fro spline and add to list
        gap_list = []
        for harmonic in harmonic_list:
            gap_list.append([harmonic,
                             float(self.splines_gap_vs_energy[harmonic](energy_eV)),
                             float(self.splines_flux_vs_energy[harmonic](energy_eV))])
            
        # Sort list in order of decreasing flux (highest flux first)
        gap_list.sort(key=lambda x: -x[2])

        # If show, create the two plots
        if len(plots) != 0 and (show == True or ofile is not None):
            if len(plots) == 2:
                plt.figure(figsize=[12, 8])
            else:
                plt.figure()

           # Check for subplotting
            if len(plots) > 1 and 'gap' in plots:
                plt.subplot(221)
                
            # Do we plot the gap vs energy?
            if 'gap' in plots:
                plt.title(name + ' Gap vs. Photon Energy')
                plt.xlabel('Photon Energy [eV]')
                plt.ylabel('Gap [mm]')

                # Grab current axis
                ax = plt.gca()

                # Scale (for linear, log, etc)
                plt.xscale(xscale)
                plt.yscale(yscale)
 
                # Loop over harmonics
                for harmonic in self.energy_range:
                    # Check for odd even
                    if harmonic % 2 and not odd:
                        continue
                    if harmonic % 2 == 0 and not even:
                        continue

                    # Get a new color for this marmonic
                    color = next(ax._get_lines.prop_cycler)['color']

                    # Check input harmonic range and make sure this is within the range
                    if harmonic_range is not None and (harmonic < harmonic_range[0] or harmonic > harmonic_range[1]):
                        continue

                    # Plot the original points (as remembered by the spline)
                    xorig = self.splines_gap_vs_energy[harmonic].x
                    plt.plot(xorig, self.splines_gap_vs_energy[harmonic](xorig), '.', color=color)

                    # Plot the spline version of this harmonic
                    xnew = np.linspace(self.energy_range[harmonic][0],
                                       self.energy_range[harmonic][-1],
                                       1000,
                                       endpoint=True)
                    plt.plot(xnew, self.splines_gap_vs_energy[harmonic](xnew), color=color, label=str(harmonic), linestyle=['-', '--'][(harmonic+1)%2])

                # Legend for harmonics
                legend1 = plt.legend(title='Harmonic', ncol=2, prop={'size': 6}, loc='upper right')
                plt.gca().add_artist(legend1)

                # Plot a line repesenting the energy, and gap found
                if energy_eV != 0:
                    le = plt.axvline(x=energy_eV, linestyle='--', linewidth=1, label=str(energy_eV)+' [eV]')
                    lg = plt.axhline(y=gap_list[0][1], linestyle='-.', linewidth=1, label=str(gap_list[0][1]) + ' [mm]')
                    plt.legend([le, lg], [str(energy_eV)+' [eV]', str(round(gap_list[0][1], 3)) + ' [mm]'], loc='upper left')

                # If plot limits defined use it
                if xlim is not None:
                    plt.xlim(xlim)
                if gap_ylim is not None:
                    plt.ylim(gap_ylim)
                
                # Show grid?
                if grid:
                    plt.grid()

                # Draw gap-only line
                if gap is not None:
                    plt.axhline(y=gap, linestyle='--', linewidth=1, color='tab:green')

            
            # Second plot flux vs energy
            if len(plots) > 1 and 'flux' in plots:
                plt.subplot(222)
                
            # Plot flux vs energy if requested
            if 'flux' in plots:
                plt.title(name + ' Flux vs. Photon Energy')
                plt.xlabel('Photon Energy [eV]')
                plt.ylabel('Flux [a.u.]')

                # Grab current axis
                ax = plt.gca()

                # Scale (for linear, log, etc)
                plt.xscale(xscale)
                plt.yscale(yscale)

                # Loop over harmonics
                for harmonic in self.energy_range:
                    # Check for odd even
                    if harmonic % 2 and not odd:
                        continue
                    if harmonic % 2 == 0 and not even:
                        continue
                    # Get a new color for this marmonic
                    color = next(ax._get_lines.prop_cycler)['color']

                    # Check input harmonic range and make sure this is within the range
                    if harmonic_range is not None and (harmonic < harmonic_range[0] or harmonic > harmonic_range[1]):
                        continue

                    # Plot the original points (as remembered by the spline)
                    xorig = self.splines_flux_vs_energy[harmonic].x
                    plt.plot(xorig, self.splines_flux_vs_energy[harmonic](xorig), '.', color=color)

                    # Plot the spline version of this harmonic
                    xnew = np.linspace(self.energy_range[harmonic][0],
                                       self.energy_range[harmonic][-1],
                                       1000,
                                       endpoint=True)
                    plt.plot(xnew, self.splines_flux_vs_energy[harmonic](xnew), color=color, label=str(harmonic), linestyle=['-', '--'][(harmonic+1)%2])

                # Legend for harmonics
                legend1 = plt.legend(title='Harmonic', ncol=2, prop={'size': 6}, loc='upper right')
                plt.gca().add_artist(legend1)

                # Plot a line repesenting the energy, and flux found
                if energy_eV != 0:
                    le = plt.axvline(x=energy_eV, linestyle='--', linewidth=1)
                    lf = plt.axhline(y=gap_list[0][2], linestyle='-.', linewidth=1)
                    plt.legend([le, lf], [str(energy_eV)+' [eV]', "{:1.1e}".format(gap_list[0][2]) + ' [a.u.]'], loc='upper left')

                # If plot limits defined use it
                if xlim is not None:
                    plt.xlim(xlim)
                if flux_ylim is not None:
                    plt.ylim(flux_ylim)

                # Show grid?
                if grid:
                    plt.grid()

            
            # Save file if ofile name is given
            if ofile is not None:
                plt.savefig(ofile, bbox_inches='tight')

            # Show the plot
            if show == True:
                plt.show()
            plt.close()
            
        return gap_list
