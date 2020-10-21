from __future__ import print_function

import oscars.sr
import oscars.th
import oscars.lut
import oscars.util
import oscars.fit
import oscars.plots_mpl

import warnings
import configparser
import os
import glob
import uuid


class bl(oscars.lut.lut1d):
    """Class for configuring beamlines with oscars

    Setup a known beamline.  Load magnetic field data tables if known, setup
    beam parameters, load summary tables.
    
    Parameters
    ----------
    facility : str
        Name of the facility.  See below for currently supported.
    
    beamline : str
        Name of beamline.  See below for currently supported.
    
    device : str
        Name of device.  See below for currently supported.

    current : float
        The current in [mA].  Not needed if specified in config.ini files

    base_path : str - optional
        Path to the oscars beamline data directory
    
    
    Returns
    -------
    None
    """

    def __init__ (self, facility=None, beamline=None, device=None, phase_mode=None, current=None, base_path=None, gpu=1, nthreads=8):
        oscars.lut.lut1d.__init__(self)

        # Set base path if defined, otherwise default
        if base_path is not None:
            self.base_path = base_path
        else:
            self.base_path = os.path.join(os.sep, 'GPFS', 'APC', 'dhidas', 'OSCARSDATA')

        if facility is None or beamline is None or device is None:
            print('You can select from the available facility, beamline, and device list:')
            self.list()
            return

        self.return_all = True
        self.show_all = True

        self.facility = facility
        self.beamline = beamline
        self.device = device

        self.name = None

        self.gap = None
        self.phase = 0

        self.phase_mode = 'planar'


        self.osr = oscars.sr.sr(gpu=gpu, nthreads=nthreads)
        self.oth = oscars.th.th(gpu=gpu, nthreads=nthreads)

        self.bfield_mapping_1d_filename = None
        self.bfield_mapping_2d_filename = None
        self.lut2d_filename = None

        self.has_lut1d = False
        self.has_lut2d = False

        # Read configuration file in order of precidence
        self.config = configparser.ConfigParser()
        self.config.read(os.path.join(self.base_path, 'Facilities', 'config.ini'))
        self.config.read(os.path.join(self.base_path, 'Facilities', facility, 'config.ini'))
        self.config.read(os.path.join(self.base_path, 'Facilities', facility, beamline, 'config.ini'))
        self.config.read(os.path.join(self.base_path, 'Facilities', facility, beamline, device, 'config.ini'))

        # Beamline and Device path
        self.beamline_path = os.path.join(self.base_path, 'Facilities', facility, beamline)
        self.device_path = os.path.join(self.beamline_path, device)
        self.bfield_path = os.path.join(self.device_path, 'bfield', self.phase_mode)

        # Setup beam
        if 'beam' in self.config:
            c = self.config['beam'] # Beam config
            a = dict()              # kwargs

            if 'type' in c: a['type'] = c['type']
            if 'name' in c: a['name'] = c['name']
            if 'energy_GeV' in c: a['energy_GeV'] = float(c['energy_GeV'])
            if 'd0' in c: a['d0'] = list(map(float, c['d0'].split()))
            if 'x0' in c: a['x0'] = list(map(float, c['x0'].split()))
            if 'beam' in c: a['beam'] = c['beam']
            if 'sigma_energy_GeV' in c: a['sigma_energy_GeV'] = float(c['sigma_energy_GeV'])
            if 't0' in c: a['t0'] = float(c['t0'])
            if current is None:
                if 'current' in c: a['current'] = float(c['current'])
            else:
                a['current'] = current
            if 'weight' in c: a['weight'] = float(c['weight'])
            if 'rotations' in c: a['rotations'] = list(map(float, c['rotations'].split()))
            if 'translation' in c: a['translation'] = list(map(float, c['translation'].split()))
            if 'horizontal_direction' in c: a['horizontal_direction'] = list(map(float, c['horizontal_direction'].split()))
            if 'beta' in c: a['beta'] = list(map(float, c['beta'].split()))
            if 'alpha' in c: a['alpha'] = list(map(float, c['alpha'].split()))
            if 'gamma' in c: a['gamma'] = list(map(float, c['gamma'].split()))
            if 'emittance' in c: a['emittance'] = list(map(float, c['emittance'].split()))
            if 'eta' in c: a['eta'] = list(map(float, c['eta'].split()))
            if 'lattice_reference' in c: a['lattice_reference'] = list(map(float, c['lattice_reference'].split()))
            if 'mass' in c: a['mass'] = float(c['mass'])
            if 'charge' in c: a['charge'] = float(c['charge'])

            self.osr.set_particle_beam(**a)

            if 'ctstartstop' in c:
                ctstartstop = list(map(float, c['ctstartstop'].split()))
                self.osr.set_ctstartstop(ctstartstop[0], ctstartstop[1])


        # Config bfield translation, etc
        self.bfield_kwargs = None
        if 'bfield' in self.config:
            c = self.config['bfield'] # bfield config
            a = dict()                # kwargs

            if 'translation' in c: a['translation'] = list(map(float, c['translation'].split()))
            if 'rotations' in c: raise ValueError('rotations keyword not allowed in config file, please remove')

            self.bfield_kwargs = a

        if 'general' in self.config:
            g = self.config['general']

            # Try to grab name
            if 'name' in g:
                self.name = g['name']
            else:
                self.name = self.device

            if phase_mode is not None:
                self.phase_mode = phase_mode
            elif 'phase_mode' in g:
                self.phase_mode = g['phase_mode']
        self.bfield_path = os.path.join(self.device_path, 'bfield', self.phase_mode)


        # Find all modes from directory
        self.modes = []
        for fname in os.listdir(os.path.join(self.device_path, 'bfield')):
            path = os.path.join(self.device_path, 'bfield', fname)
            if os.path.isdir(path):
                self.modes.append(fname)

        # Check default mode exists
        if self.phase_mode not in self.modes:
            warnings.warn('default phase mode is not in list of modes.  check default or directory structure: ' + self.phase_mode)
            self.phase_mode = None

        # Set bfield map and lookup table paths
        if self.phase_mode is not None:
            if os.path.isfile(os.path.join(self.bfield_path, 'file_list_1d.txt')):
                self.bfield_mapping_1d_filename = os.path.join(self.bfield_path, 'file_list_1d.txt')
            if os.path.isfile(os.path.join(self.bfield_path, 'file_list_2d.txt')):
                self.bfield_mapping_2d_filename = os.path.join(self.bfield_path, 'file_list_2d.txt')
            if os.path.isfile(os.path.join(self.bfield_path, 'lut1d.txt')):
                self.lut1d_filename = os.path.join(self.bfield_path, 'lut1d.txt')
                self.read_file_lut1d(self.lut1d_filename)
                self.has_lut1d = True
            if os.path.isfile(os.path.join(self.bfield_path, 'lut2d.txt')):
                self.lut2d_filename = os.path.join(self.bfield_path, 'lut2d.txt')
                #self.read_file_lut2d( self.lut2d_filename)
                self.has_lut2d = True
        else:
            warnings.warn('phase_mode is None.  no bfield or LUT setup')

        self.spectrum_path = os.path.join(self.device_path, 'spectrum', self.phase_mode)
        self.harmonics_path = os.path.join(self.device_path, 'harmonics', self.phase_mode)

        # Setup spectrum arguments
        self.spectrum_kwargs= None
        if 'spectrum' in self.config:
            c = self.config['spectrum'] # Spectrum config
            a = dict()                  # kwargs

            if 'obs' in c: a['obs'] = list(map(float, c['obs'].split()))
            if 'npoints' in c: a['npoints'] = int(c['npoints'])
            if 'energy_range_eV' in c: a['energy_range_eV'] = list(map(float, c['energy_range_eV'].split()))
            if 'points_eV' in c: a['points_eV'] = list(map(float, c['points_eV'].split()))
            if 'polarization' in c: a['polarization'] = c['polarization']
            if 'horizontal_direction' in c: a['horizontal_direction'] = list(map(float, c['horizontal_direction'].split()))
            if 'vertical_direction' in c: a['vertical_direction'] = list(map(float, c['vertical_direction'].split()))
            if 'precision' in c: a['precision'] = float(c['precision'])
            if 'max_level' in c: a['max_level'] = int(c['max_level'])
            if 'max_level_extended' in c: a['max_level_extended'] = int(c['max_level_extended'])
            if 'angle' in c: a['angle'] = float(c['angle'])
            if 'nparticles' in c: a['nparticles'] = int(c['nparticles'])
            if 'nthreads' in c: a['nthreads'] = int(c['nthreads'])
            if 'gpu' in c: a['gpu'] = int(c['gpu'])
            if 'ngpu' in c:
                a['ngpu'] = list(map(int, c['ngpu'].split()))
                if len(a[ngpu]) == 1: a['ngpu'] = a['ngpu'][0]
            if 'quantity' in c: a['quantity'] = c['quantity']
            if 'ofile' in c: a['ofile'] = c['ofile']
            if 'bofile' in c: a['bofile'] = c['bofile']

            self.spectrum_kwargs = a


        # Setup flux arguments
        self.flux_kwargs= None
        if 'flux' in self.config:
            c = self.config['flux'] # Flux config
            a = dict()              # kwargs

            if 'energy_eV' in c: a['energy_eV'] = float(c['energy_eV'])
            if 'npoints' in c: a['npoints'] = list(map(int, c['npoints'].split()))
            if 'plane' in c: a['plane'] = c['plane']
            if 'normal' in c: a['normal'] = int(c['normal'])
            if 'dim' in c: a['dim'] = int(c['dim'])
            if 'width' in c: a['width'] = list(map(float, c['width'].split()))
            if 'rotations' in c: a['rotations'] = list(map(float, c['rotations'].split()))
            if 'translation' in c: a['translation'] = list(map(float, c['translation'].split()))
            if 'x0x1x2' in c:
                v = c['x0x1x2'].split(',')
                if len(v) != 3:
                    warnings.warn('x0x1x2 does not contain 3 vectors.  Check input format')
                else:
                    x0x1x2 = [ list(map(float, v[0])), list(map(float, v[1])), list(map(float, v[2])) ]
                    a['x0x1x2'] = x0x1x2
            if 'polarization' in c: a['polarization'] = c['polarization']
            if 'angle' in c: a['angle'] = float(c['angle'])
            if 'horizontal_direction' in c: a['horizontal_direction'] = list(map(float, c['horizontal_direction'].split()))
            if 'vertical_direction' in c: a['vertical_direction'] = list(map(float, c['vertical_direction'].split()))
            if 'nparticles' in c: a['nparticles'] = int(c['nparticles'])
            if 'nthreads' in c: a['nthreads'] = int(c['nthreads'])
            if 'gpu' in c: a['gpu'] = int(c['gpu'])
            if 'ngpu' in c:
                a['ngpu'] = list(map(int, c['ngpu'].split()))
                if len(a[ngpu]) == 1: a['ngpu'] = a['ngpu'][0]
            if 'precision' in c: a['precision'] = float(c['precision'])
            if 'max_level' in c: a['max_level'] = int(c['max_level'])
            if 'max_level_extended' in c: a['max_level_extended'] = int(c['max_level_extended'])
            if 'quantity' in c: a['quantity'] = c['quantity']
            if 'ofile' in c: a['ofile'] = c['ofile']
            if 'bofile' in c: a['bofile'] = c['bofile']

            self.flux_kwargs = a

        # Setup flux arguments
        self.power_density_kwargs= None
        if 'power_density' in self.config:
            c = self.config['power_density']    # Flux config
            a = dict()                          # kwargs

            if 'npoints' in c: a['npoints'] = list(map(int, c['npoints'].split()))
            if 'plane' in c: a['plane'] = c['plane']
            if 'width' in c: a['width'] = list(map(float, c['width'].split()))
            if 'x0x1x2' in c:
                v = c['x0x1x2'].split(',')
                if len(v) != 3:
                    warnings.warn('x0x1x2 does not contain 3 vectors.  Check input format')
                else:
                    x0x1x2 = [ list(map(float, v[0])), list(map(float, v[1])), list(map(float, v[2])) ]
                    a['x0x1x2'] = x0x1x2
            if 'rotations' in c: a['rotations'] = list(map(float, c['rotations'].split()))
            if 'translation' in c: a['translation'] = list(map(float, c['translation'].split()))
            if 'ofile' in c: a['ofile'] = c['ofile']
            if 'bofile' in c: a['bofile'] = c['bofile']
            if 'normal' in c: a['normal'] = int(c['normal'])
            if 'nparticles' in c: a['nparticles'] = int(c['nparticles'])
            if 'gpu' in c: a['gpu'] = int(c['gpu'])
            if 'ngpu' in c:
                a['ngpu'] = list(map(int, c['ngpu'].split()))
                if len(a[ngpu]) == 1: a['ngpu'] = a['ngpu'][0]
            if 'nthreads' in c: a['nthreads'] = int(c['nthreads'])
            if 'precision' in c: a['precision'] = float(c['precision'])
            if 'max_level' in c: a['max_level'] = int(c['max_level'])
            if 'max_level_extended' in c: a['max_level_extended'] = int(c['max_level_extended'])
            if 'dim' in c: a['dim'] = int(c['dim'])
            if 'quantity' in c: a['quantity'] = c['quantity']

            self.power_density_kwargs = a

        return


    def load_lut1d (self, ifile):
        """Load a specific lut1d file"""


        if os.path.isfile(ifile):
            try:
                super(bl, self).clear()
                self.read_file_lut1d(ifile)
                self.lut1d_filename = ifile
                self.has_lut1d = True
            except:
                raise IOError('Unable to load file.  Probably incorrect format or permissions: ' + ifile)
        else:
            raise IOError('file does not exist: ' + ifile)

        return




    def list (self):
        """List facilities and beamlines"""

        fstring = '    {:12} {:12}     {}'

        print(fstring.format('Beamline', 'Device', 'Modes'))
        print(fstring.format('--------', '------', '-----'))
        print('')
        for facility_name in os.listdir(os.path.join(self.base_path, 'Facilities')):
            facility_path = os.path.join(self.base_path, 'Facilities', facility_name)
            if os.path.isdir(facility_path):
                print(facility_name)

                for beamline_name in os.listdir(facility_path):
                    beamline_path = os.path.join(facility_path, beamline_name)
                    if os.path.isdir(beamline_path):
                        #print('    Beamline:', beamline_name)

                        for device_name in os.listdir(beamline_path):
                            device_path = os.path.join(beamline_path, device_name)
                            if os.path.isdir(device_path):
                                #print('        Device:', device_name)

                                modes = []
                                try:
                                    for mode_name in os.listdir(os.path.join(device_path, 'bfield')):
                                        mode_path = os.path.join(device_path, 'bfield', mode_name)
                                        if os.path.isdir(mode_path):
                                            modes.append(mode_name)
                                except:
                                    warnings.warn('No bfields found for device: ' + device_name)
                                print(fstring.format(beamline_name, device_name, ' '.join(modes)))





        return






    def info (self):
        """Print info about this beamline setup"""

        print('***BEGIN CONFIG***')
        for section in self.config.sections():
            print(section)
            for key in self.config[section]: print(' ', key, self.config[section][key])
        print('***END CONFIG***')
        print('')
        print('***BEGIN OSCARS CONFIG***')
        self.osr.print_all()
        print('***END OSCARS CONFIG***')






    def return_return (self, ret):
        """
        Do I return the return value.  Based on self.return_all and ret input.  If nothing is specified
        the default is to return true.

        Parameters
        ----------
        ret : bool
            Input return

        Returns
        -------
        True of False depending on self.return_all and ret, defaulting to True
        """

        if ret is not None:
            return ret
        if self.return_all is not None:
            return self.return_all
        return True






    def show_plot (self, show):
        """
        Do I show the plot in a function call.  Based on self.show_all and show input.  If nothing is specified
        the default is to return true.

        Parameters
        ----------
        show : bool
            Input return

        Returns
        -------
        True of False depending on self.return_all and ret, defaulting to True
        """

        if show is not None:
            return show
        if self.show_all is not None:
            return self.show_all
        return True




    def plot_spectra (self, spectra, **kwargs):
        """
        Plot multiple spectra.  Forward call to oscars.plots_mpl.plot_spectra().  See there for documentation.
        """

        if 'title' not in kwargs:
            kwargs['title'] = self.facility + ' ' + self.beamline + ' ' + self.device + ' Spectra'
        if 'show' not in kwargs:
            kwargs['show'] = self.show_plot(show=None)

        return oscars.plots_mpl.plot_spectra(spectra, **kwargs)
        





    def bfield (self, gap=None, **kwargs):
        """Plot the magnetic field.

        Parameters
        ----------
        gap : str
            Gap of device - optional

        Returns
        -------
        None
        """

        if gap is not None:
            self.set_gap(gap)

        if self.gap is None:
            raise RuntimeError('gap is not set.  try set_gap() first')

        oscars.plots_mpl.plot_bfield(self.osr)

        return


    def trajectory (self, gap=None, **kwargs):
        """Plot the trajectory

        Parameters
        ----------
        gap : str
            Gap of device - optional

        Returns
        -------
        None
        """

        if gap is not None:
            self.set_gap(gap)

        if self.gap is None:
            raise RuntimeError('gap is not set.  try set_gap() first')

        self.osr.set_new_particle(particle='ideal')

        trajectory = self.osr.calculate_trajectory()
        oscars.plots_mpl.plot_trajectory_position(trajectory)

        return

    def spectrum (self, gap=None, phase=None, ret=None, show=None, **kwargs):
        """Calculate and plot spectrum
        
        Parameters
        ----------
        gap : float
            Gap of device - optional

        phase : float
            If used, default None

        All other parameters can be found in the function:
            oscars.sr.calculate_spectrum()

        Returns
        -------
        None
        """

        if gap is not None or phase is not None:
            self.set_gap(gap, phase)

        if self.gap is None:
            raise RuntimeError('gap is not set.  try set_gap() first')

        newargs = self.spectrum_kwargs.copy()

        for key in kwargs:
            newargs[key] = kwargs[key]

        if 'nparticles' in newargs:
            if newargs['nparticles'] <= 1:
                self.osr.set_new_particle(particle='ideal')
            else:
                self.osr.set_new_particle()
        else:
            self.osr.set_new_particle(particle='ideal')

        fofile = None
        if 'fofile' in newargs:
            fofile = newargs['fofile']
            del newargs['fofile']

        spectrum = self.osr.calculate_spectrum(**newargs)

        oscars.plots_mpl.plot_spectrum(spectrum, figsize=[16, 4], show=self.show_plot(show), ofile=fofile)

        # Check return request
        if self.return_return(ret):
            return spectrum
        return





    def flux (self, gap=None, phase=None, ret=None, show=None, **kwargs):
        """Calculate the flux
        
        Parameters
        ----------
        gap : float
            Gap of device - optional

        phase : float
            If used, default None

        All other parameters can be found in the function:
            oscars.sr.calculate_flux_rectangle()

        Returns
        -------
        None
        """

        if gap is not None or phase is not None:
            self.set_gap(gap)

        if self.gap is None:
            raise RuntimeError('gap is not set.  try set_gap() first')

        newargs = self.flux_kwargs.copy()

        for key in kwargs:
            newargs[key] = kwargs[key]

        if 'nparticles' in newargs:
            if newargs['nparticles'] <= 1:
                self.osr.set_new_particle(particle='ideal')
        else:
            self.osr.set_new_particle(particle='ideal')

        fofile = None
        if 'fofile' in newargs:
            fofile = newargs['fofile']
            del newargs['fofile']

        flux = self.osr.calculate_flux_rectangle(**newargs)

        oscars.plots_mpl.plot_flux(flux, show=self.show_plot(show), ofile=fofile)

        # If no return value return
        if self.return_return(ret):
            return flux
        return


    def set_phase (self, gap=None, phase=None):
        return set_gap(gap=gap, phase=phase)

    def set_gap (self, gap=None, phase=None):
        """Set the gap for a device.  Sets basic magnetic field data"""

        if gap is None and phase is None:
            raise RuntimeError('gap or phase must be something')

        if gap is None:
            gap = self.gap
        if phase is None:
            phase = self.phase

        if phase is None:
            phase = self.phase

        self.osr.clear_bfields()

        if self.bfield_mapping_1d_filename is not None:
            mapping = oscars.util.read_file_list_with_header(self.bfield_mapping_1d_filename)
        elif self.bfield_mapping_2d_filename is not None:
            if phase is None:
                phase = 0
            mapping = oscars.util.read_file_list_2d_with_header(self.bfield_mapping_2d_filename)
        else:
            raise RuntimeError('no mapping found.  make sure phase_mode is set correctly')

        translation = mapping.translation
        if 'translation' in self.bfield_kwargs:
            t = self.bfield_kwargs['translation']
            if len(t) != 3:
                print('t:', t)
                raise IndexError('translation input of incorrect length')
            for i in range(3):
                translation[i] += t[i]

        if self.bfield_mapping_1d_filename is not None:
            self.osr.add_bfield_interpolated(
                mapping=mapping.mapping,
                iformat=mapping.format,
                rotations=mapping.rotations,
                translation=translation,
                scale=mapping.scale,
                parameter=gap
            )
        elif self.bfield_mapping_2d_filename is not None:
            newfname = '.OSCARS_tmp_' + str(uuid.uuid4())
            newfile = oscars.util.interpolate_file_list_2d (mapping=mapping, gap=gap, phase=phase, ofile=newfname)
            self.osr.add_bfield_file(
                ifile=newfname,
                iformat=mapping.format,
                rotations=mapping.rotations,
                translation=translation,
                scale=mapping.scale
            )

            if os.path.exists(newfname): os.remove(newfname)


        self.gap = gap
        self.phase = phase

        return



    def total_power (self, gap=None):
        """Calculate the total power
        
        Parameters
        ----------
        gap : str
            Gap of device - optional

        Returns
        -------
        Estimate of the total power in [W]
        """

        if gap is not None:
            self.set_gap(gap)

        if self.gap is None:
            raise RuntimeError('gap is not set.  try set_gap() first')

        return self.osr.calculate_total_power()



    def power_density (self, gap=None, phase=None, ret=None, show=None, **kwargs):
        """Calculate the power_density
        
        Parameters
        ----------
        gap : foat
            Gap of device - optional

        phase : float
            Desired phase if phase mode allows

        All other parameters can be found in the function:
            oscars.sr.calculate_power_density_rectangle()

        Returns
        -------
        None
        """

        if gap is not None or phase is not None:
            self.set_gap(gap, phase)

        if self.gap is None:
            raise RuntimeError('gap is not set.  try set_gap() first')

        newargs = self.power_density_kwargs.copy()

        for key in kwargs:
            newargs[key] = kwargs[key]

        if 'nparticles' in newargs:
            if newargs['nparticles'] <= 1:
                self.osr.set_new_particle(particle='ideal')
        else:
            self.osr.set_new_particle(particle='ideal')

        fofile = None
        if 'fofile' in newargs:
            fofile = newargs['fofile']
            del newargs['fofile']

        power_density = self.osr.calculate_power_density_rectangle(**newargs)

        oscars.plots_mpl.plot_power_density(power_density, show=self.show_plot(show), ofile=fofile)

        # If no return value return
        if self.return_return(ret):
            return power_density
        return



    def summary (self, gap=None):
        """
        Print some summary information and plots
        
        Parameters
        ----------
        gap : str
            Gap of device - optional

        Returns
        -------
        None
        """

        if gap is not None:
            try:
                self.set_gap(gap)
            except:
                warnings.warn('could not set gap to that value.  likely it is outside of the interpolation range')

        if self.has_lut1d:
            self.get_gaps(show=True, gap=self.gap, name=self.name)

        if self.gap is not None:

            self.osr.set_new_particle(particle='ideal')
            spectrum = self.osr.calculate_spectrum(**self.spectrum_kwargs)

            tmp_flux_kwargs = self.flux_kwargs.copy()
            if 'energy_eV' not in self.flux_kwargs:
                tmp_flux_kwargs['energy_eV'] = oscars.fit.find_first_harmonic(spectrum)[1]

            oscars.plots_mpl.plot_spectrum(spectrum, figsize=[16, 4], axvlines=[tmp_flux_kwargs['energy_eV']])

            oscars.plots_mpl.plot_flux(self.osr.calculate_flux_rectangle(**tmp_flux_kwargs))

            oscars.plots_mpl.plot_power_density(self.osr.calculate_power_density_rectangle(**self.power_density_kwargs))

            print('Estimated total power radiated:', round(self.osr.calculate_total_power(), 3), '[W]')

        return


    def get_spectrum_files (self, pattern='spec*', directory=None):
        """
        Get the list of spectrum files that exists for this phase_mode.  Used in energy/harmonics calculations (Not so useful for general users)
        
        Parameters
        ----------
        pattern : str
            files will be matched to this pattern using glob

        directory : str
            Alternative directiry to look at.  Typically should use the default.

        Returns
        -------
        Sorted list of files
        """
        if directory is None:
            directory = self.spectrum_path

        files = glob.glob(os.path.join(directory, pattern))
        files.sort()

        return files


    def get_spectrum_file_gap (self, filename):
        """
        Get spectrum gap from filename (Not so useful for general users)
        
        Parameters
        ----------
        filename : str
            Name of file to parse


        Returns
        -------
        float - gap
        """

        return float(os.path.splitext(os.path.basename(filename).split('_')[-1])[0])


    def get_spectrum_file_phase (self, filename):
        """
        Get spectrum phase from filename (Not so useful for general users)
        
        Parameters
        ----------
        filename : str
            Name of file to parse


        Returns
        -------
        float - gap
        """

        return float(os.path.basename(f).split('_')[-3])


    def get_filename (self, gap, phase=0.0):
        """
        Get common output name (Not so useful for general users)
        
        Parameters
        ----------
        gap : float
            gap

        phase : float
            phase - default is 0.0

        Returns
        -------
        str : filename formated output
        """

        return '{}_{}_{}_phase_{:+08.3f}_gap_{:07.3f}'.format(self.beamline, self.device, self.phase_mode, phase, gap)



    def parse_filename (self, filename):
        """
        parse common output name (Not so useful for general users)
        
        Parameters
        ----------
        filename : str
            Name of file to parse

        Returns
        -------
        obj : with gap, phase, phase_mode, device, beamline
        """

        class fileinfo:

            def __init__ (self, filename):
                fields = os.path.basename(filename).split('_')[-7:]
                self.gap = float(os.path.splitext(fields[-1])[0])
                self.phase = float(fields[-3])
                self.phase_mode = fields[-5]
                self.device = fields[-6]
                self.beamline = fields[-7]

                del fields

        a = fileinfo(filename)

        return a
