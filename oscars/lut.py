import numpy as np
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt

class lut1d:
    """Class for 1D lookup tables for insertion devices as a function of gap"""

    def __init__(self, ifile=None, name=None):
        print('init called')
        self.name = name
        self.splines_gap_vs_energy = dict()
        self.splines_flux_vs_energy = dict()
        self.energy_range = dict()
        
        if ifile is not None:
            self.read_file(ifile)
            
        return

        
    def read_file(self, ifile):
        """read a file and setup data accordingly"""
        print('read called')
        
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
                        if len(fields) != 3:
                            raise ValueError('Error in format of line:' + line)
                        all_fields.append(fields)

                    if len(all_fields) < 2:
                        continue

                    all_fields.sort(key=lambda x: x[0])

                    energy = [x[0] for x in all_fields]
                    gap = [x[1] for x in all_fields]
                    flux = [x[2] for x in all_fields]

                    #cs_gap_vs_energy = CubicSpline(energy, gap)
                    #cs_flux_vs_energy = CubicSpline(energy, flux)

                    self.splines_gap_vs_energy[harmonic] = CubicSpline(energy, gap)
                    self.splines_flux_vs_energy[harmonic] = CubicSpline(energy, flux)
                    self.energy_range[harmonic] = [energy[0], energy[-1]]

        return

    def get_gaps (self,
                  energy=0,
                  show=False,
                  ofile=None,
                  name='',
                  harmonic_range=None,
                  plots=['gap', 'flux'],
                  xscale='linear',
                  yscale='linear',
                  xlim=None,
                  gap_ylim=None,
                  flux_ylim=None,
                  grid=False
                 ):
        """Get the gaps at which to find the given energy"""
        
        # Grab the harmonics which cover the range
        harmonic_list = []
        for harmonic in self.energy_range:
            # Check input harmonic range and make sure this is within the range
            if energy == 0 or (harmonic_range is not None and (harmonic < harmonic_range[0] or harmonic > harmonic_range[1])):
                continue
                
            if energy >= self.energy_range[harmonic][0] and energy <= self.energy_range[harmonic][1]:
                harmonic_list.append(harmonic)

        # Get results fro spline and add to list
        gap_list = []
        for harmonic in harmonic_list:
            gap_list.append([harmonic,
                             float(self.splines_gap_vs_energy[harmonic](energy)),
                             float(self.splines_flux_vs_energy[harmonic](energy))])
            
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
                    plt.plot(xnew, self.splines_gap_vs_energy[harmonic](xnew), color=color, label=str(harmonic))

                # Legend for harmonics
                legend1 = plt.legend(title='Harmonic', ncol=2, prop={'size': 6}, loc='upper right')
                plt.gca().add_artist(legend1)

                # Plot a line repesenting the energy, and gap found
                if energy != 0:
                    le = plt.axvline(x=energy, linestyle='--', linewidth=1, label=str(energy)+' [eV]')
                    lg = plt.axhline(y=gap_list[0][1], linestyle='-.', linewidth=1, label=str(gap_list[0][1]) + ' [mm]')
                    plt.legend([le, lg], [str(energy)+' [eV]', str(round(gap_list[0][1], 3)) + ' [mm]'], loc='upper left')

                # If plot limits defined use it
                if xlim is not None:
                    plt.xlim(xlim)
                if gap_ylim is not None:
                    plt.ylim(gap_ylim)
                
                # Show grid?
                if grid:
                    plt.grid()


            
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
                    plt.plot(xnew, self.splines_flux_vs_energy[harmonic](xnew), color=color, label=str(harmonic))

                # Legend for harmonics
                legend1 = plt.legend(title='Harmonic', ncol=2, prop={'size': 6}, loc='upper right')
                plt.gca().add_artist(legend1)

                # Plot a line repesenting the energy, and flux found
                if energy != 0:
                    le = plt.axvline(x=energy, linestyle='--', linewidth=1)
                    lf = plt.axhline(y=gap_list[0][2], linestyle='-.', linewidth=1)
                    plt.legend([le, lf], [str(energy)+' [eV]', "{:1.1e}".format(gap_list[0][2]) + ' [a.u.]'], loc='upper left')

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
