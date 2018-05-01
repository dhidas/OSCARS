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

    def __init__ (self, facility=None, beamline=None, device=None, current=None, base_path=None, gpu=1, nthreads=8):
        oscars.lut.lut1d.__init__(self)

        # Set base path if defined, otherwise default
        if base_path is not None:
            self.base_path = base_path
        else:
            self.base_path = os.path.join(os.sep, 'Users', 'dhidas', 'OSCARSDATA')

        if facility is None or beamline is None or device is None:
            print('You can select from the available facility, beamline, and device list:')
            self.list()
            return

        self.facility = facility
        self.beamline = beamline
        self.device = device

        self.name = None

        self.gap = None
        self.phase = None
        self.phase_mode = None

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
        self.bfield_path = os.path.join(self.device_path, 'bfield')

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

        self.phase_mode = None
        if 'general' in self.config:
            g = self.config['general']

            # Try to grab name
            if 'name' in g:
                self.name = g['name']
            else:
                self.name = self.device

            if 'default_mode' in g: self.phase_mode = g['default_mode']
        else:
            self.phase_mode = 'planar'


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
            if os.path.isfile(os.path.join(self.bfield_path, self.phase_mode, 'file_list_1d.txt')):
                self.bfield_mapping_1d_filename = os.path.join(self.bfield_path, self.phase_mode, 'file_list_1d.txt')
            if os.path.isfile(os.path.join(self.bfield_path, self.phase_mode, 'file_list_2d.txt')):
                self.bfield_mapping_2d_filename = os.path.join(self.bfield_path, self.phase_mode, 'file_list_2d.txt')
            if os.path.isfile(os.path.join(self.bfield_path, self.phase_mode, 'lut1d.txt')):
                self.lut1d_filename = os.path.join(self.bfield_path, self.phase_mode, 'lut1d.txt')
                self.read_file_lut1d(self.lut1d_filename)
                self.has_lut1d = True
            if os.path.isfile(os.path.join(self.bfield_path, self.phase_mode, 'lut2d.txt')):
                self.lut2d_filename = os.path.join(self.bfield_path, self.phase_mode, 'lut2d.txt')
                #self.read_file_lut2d( self.lut2d_filename)
                self.has_lut2d = True
        else:
            warnings.warn('phase_mode is None.  no bfield or LUT setup')


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
                self.clear_lut1d()
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


    def spectrum (self, gap=None, **kwargs):
        """Calculate and plot spectrum
        
        Parameters
        ----------
        gap : str
            Gap of device - optional

        All other parameters can be found in the function:
            oscars.sr.calculate_spectrum()

        Returns
        -------
        None
        """

        if gap is not None:
            self.set_gap(gap)

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

        spectrum = self.osr.calculate_spectrum(**newargs)
        oscars.plots_mpl.plot_spectrum(spectrum, figsize=[16, 4])

        return





    def flux (self, gap=None, **kwargs):
        """Calculate the flux
        
        Parameters
        ----------
        gap : str
            Gap of device - optional

        All other parameters can be found in the function:
            oscars.sr.calculate_flux_rectangle()

        Returns
        -------
        None
        """

        if gap is not None:
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

        flux = self.osr.calculate_flux_rectangle(**newargs)
        oscars.plots_mpl.plot_flux(flux)

        return


    def set_gap (self, gap):
        """Set the gap for a device.  Sets basic magnetic field data"""

        self.gap = None

        self.osr.clear_bfields()
        mapping = oscars.util.read_file_list_with_header(self.bfield_mapping_1d_filename)

        translation = mapping[3]
        if 'translation' in self.bfield_kwargs:
            t = self.bfield_kwargs['translation']
            if len(t) != 3:
                print('t:', t)
                raise IndexError('translation input of incorrect length')
            for i in range(3):
                translation[i] += t[i]

        self.osr.add_bfield_interpolated(mapping=mapping[0], iformat=mapping[1], rotations=mapping[2], translation=translation, scale=mapping[4], parameter=gap)

        self.gap = gap

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



    def power_density (self, gap=None, **kwargs):
        """Calculate the power_density
        
        Parameters
        ----------
        gap : str
            Gap of device - optional

        All other parameters can be found in the function:
            oscars.sr.calculate_power_density_rectangle()

        Returns
        -------
        None
        """

        if gap is not None:
            self.set_gap(gap)

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

        power_density = self.osr.calculate_power_density_rectangle(**newargs)
        oscars.plots_mpl.plot_power_density(power_density)

        return



    def summary (self, gap=None):
        """Print some summary information and plots
        
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
