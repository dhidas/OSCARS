import oscars.sr
import oscars.th
import oscars.lut
import oscars.util
import oscars.plots_mpl

import warnings
import configparser
import os


class bl(oscars.lut.lut1d):
    """Class for configuring beamlines with oscars"""

    def __init__ (self, facility, beamline, device, base_path=None, gpu=1, nthreads=8):
        oscars.lut.lut1d.__init__(self)

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
        self.lut1d_filename = None
        self.lut2d_filename = None

        self.has_lut1d = False
        self.has_lut2d = False

        if base_path is not None:
            self.base_path = base_path
        else:
            self.base_path = '/Users/dhidas/OSCARSDATA'

        # Read configuration file in order of precidence
        self.config = configparser.ConfigParser()
        self.config.read(os.path.join(self.base_path, 'Facilities', 'config.ini'))
        self.config.read(os.path.join(self.base_path, 'Facilities', facility, 'config.ini'))
        self.config.read(os.path.join(self.base_path, 'Facilities', facility, beamline, 'config.ini'))
        self.config.read(os.path.join(self.base_path, 'Facilities', facility, beamline, device, 'config.ini'))

        for section in self.config.sections():
            print(section)
            for key in self.config[section]: print(' ', key, self.config[section][key])

        # Beamline and Device path
        self.beamline_path = os.path.join(self.base_path, 'Facilities', facility, beamline)
        self.device_path = os.path.join(self.beamline_path, device)
        self.bfield_path = os.path.join(self.device_path, 'bfield')

        # Setup beam
        if 'beam' in self.config:
            b = self.config['beam'] # Beam config
            a = dict()              # kwargs

            if 'type' in b: a['type'] = b['type']
            if 'name' in b: a['name'] = b['name']
            if 'energy_GeV' in b: a['energy_GeV'] = float(b['energy_GeV'])
            if 'd0' in b: a['d0'] = list(map(float, b['d0'].split()))
            if 'x0' in b: a['x0'] = list(map(float, b['x0'].split()))
            if 'beam' in b: a['beam'] = b['beam']
            if 'sigma_energy_GeV' in b: a['sigma_energy_GeV'] = float(b['sigma_energy_GeV'])
            if 't0' in b: a['t0'] = float(b['t0'])
            if 'current' in b: a['current'] = float(b['current'])
            if 'weight' in b: a['weight'] = float(b['weight'])
            if 'rotations' in b: a['rotations'] = list(map(float, b['rotations'].split()))
            if 'translation' in b: a['translation'] = list(map(float, b['translation'].split()))
            if 'horizontal_direction' in b: a['horizontal_direction'] = list(map(float, b['horizontal_direction'].split()))
            if 'beta' in b: a['beta'] = list(map(float, b['beta'].split()))
            if 'alpha' in b: a['alpha'] = list(map(float, b['alpha'].split()))
            if 'gamma' in b: a['gamma'] = list(map(float, b['gamma'].split()))
            if 'emittance' in b: a['emittance'] = list(map(float, b['emittance'].split()))
            if 'eta' in b: a['eta'] = list(map(float, b['eta'].split()))
            if 'lattice_reference' in b: a['lattice_reference'] = list(map(float, b['lattice_reference'].split()))
            if 'mass' in b: a['mass'] = float(b['mass'])
            if 'charge' in b: a['charge'] = float(b['charge'])

            self.osr.set_particle_beam(**a)

            if 'ctstartstop' in b:
                print('ctstartstop:', b['ctstartstop'])
                ctstartstop = list(map(float, b['ctstartstop'].split()))
                self.osr.set_ctstartstop(ctstartstop[0], ctstartstop[1])

        self.phase_mode = None
        if 'general' in self.config:
            g = self.config['general']
            if 'name' in g: self.name = g['name']
            if 'default_mode' in g: self.phase_mode = g['default_mode']

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
            b = self.config['spectrum'] # Spectrum config
            a = dict()                  # kwargs

            if 'obs' in b: a['obs'] = list(map(float, b['obs'].split()))
            print('obsobsobs', a['obs'])
            if 'npoints' in b: a['npoints'] = int(b['npoints'])
            if 'energy_range_eV' in b: a['energy_range_eV'] = list(map(float, b['energy_range_eV'].split()))
            if 'points_eV' in b: a['points_eV'] = list(map(float, b['points_eV'].split()))
            if 'polarization' in b: a['polarization'] = b['polarization']
            if 'horizontal_direction' in b: a['horizontal_direction'] = list(map(float, b['horizontal_direction'].split()))
            if 'vertical_direction' in b: a['vertical_direction'] = list(map(float, b['vertical_direction'].split()))
            if 'precision' in b: a['precision'] = float(b['precision'])
            if 'max_level' in b: a['max_level'] = int(b['max_level'])
            if 'max_level_extended' in b: a['max_level_extended'] = int(b['max_level_extended'])
            if 'angle' in b: a['angle'] = float(b['angle'])
            if 'nparticles' in b: a['nparticles'] = int(b['nparticles'])
            if 'nthreads' in b: a['nthreads'] = int(b['nthreads'])
            if 'gpu' in b: a['gpu'] = int(b['gpu'])
            if 'ngpu' in b:
                a['ngpu'] = list(map(int, b['ngpu'].split()))
                if len(a[ngpu]) == 1: a['ngpu'] = a['ngpu'][0]
            if 'quantity' in b: a['quantity'] = b['quantity']
            if 'ofile' in b: a['ofile'] = b['ofile']
            if 'bofile' in b: a['bofile'] = b['bofile']

            self.spectrum_kwargs = a


        return




    def setup_beamline (self, facility, beamline, device):
        """
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
    
        osr : oscars.sr object
            oscars.sr object to load
    
        oth : oscars.th object
            oscars.th object to load
    
        Currently Supported - facility beamline
            NSLSII SST U42
    
        Returns
        -------
        None
        """

        return   




    def set_gap (self, gap):
        """Set the gap for a device.  Sets basic magnetic field data"""

        self.gap = None

        self.osr.clear_bfields()
        mapping = oscars.util.read_file_list_with_header(self.bfield_mapping_1d_filename)
        self.osr.add_bfield_interpolated(mapping=mapping[0], iformat=mapping[1], rotations=mapping[2], translation=mapping[3], scale=mapping[4], parameter=gap)
        self.osr.set_new_particle()

        self.gap = gap

        return



    def summary (self, gap=None):
        """Print some summary information and plots"""

        if gap is not None:
            try:
                self.set_gap(gap)
            except:
                warnings.warn('could not set gap to that value.  likely it is outside of the interpolation range')

        if self.has_lut1d is not None:
            self.get_gaps(show=True, gap=self.gap)

        if self.gap is not None:
            oscars.plots_mpl.plot_bfield(self.osr)
            oscars.plots_mpl.plot_trajectory_position(self.osr.calculate_trajectory())
            oscars.plots_mpl.plot_spectrum(self.osr.calculate_spectrum(**self.spectrum_kwargs), figsize=[16, 4])

        return
