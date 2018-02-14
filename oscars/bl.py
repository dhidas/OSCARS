import oscars.sr
import oscars.th
import oscars.lut
import oscars.util


class bl:
    """Class for configuring beamlines with oscars"""

    def __init__ (self, facility, beamline, device, base_path=None, gpu=1, nthreads=8):
        self.facility = facility
        self.beamline = beamline
        self.device = device

        self.gap = None
        self.phase = None
        self.phase_mode = None

        self.osr = oscars.sr.sr(gpu=gpu, nthreads=nthreads)
        self.oth = oscars.th.th(gpu=gpu, nthreads=nthreads)

        self.lut1d = None
        self.base_path = '/Users/dhidas/OSCARSDATA'
        self.device_path = None
        self.bfield_lut1d_file = None

        self.spectrum_range = None

        self.setup_beamline(facility, beamline, device)

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
    
        self.device_path = self.base_path + '/Facilities/' + '/'.join([facility, beamline, device])
        
        if facility == 'NSLSII':
            if beamline == 'SST' and device == 'U42':
                print('Setting up SST U42')

                # Setup beam
                self.osr.set_particle_beam(beam='NSLSII', x0=[0, 0, -1])
                self.osr.set_ctstartstop(0, 2)

                self.phase_mode = 'planar'
                self.lut1d = oscars.lut.lut1d(self.device_path + '/lut/' + self.phase_mode + '/lut1d.txt')
                self.bfield_lut1d_file = self.device_path + '/bfield/' + self.phase_mode + '/file_list.txt'

                self.spectrum_range = [10, 8000]

                return
            elif beamline == 'SST' and device == 'EPU60':
                print('Setting up SST EPU60')

                # Setup beam
                self.osr.set_particle_beam(beam='NSLSII', x0=[0, 0, -1])
                self.osr.set_ctstartstop(0, 2)

                self.phase_mode = 'planar'
                self.lut1d = oscars.lut.lut1d(self.device_path + '/lut/' + self.phase_mode + '/lut1d.txt')
                self.bfield_lut1d_file = self.device_path + '/bfield/' + self.phase_mode + '/file_list.txt'

                self.spectrum_range = [10, 3000]

                pass
            else:
              raise  ValueError('beamline and/or device not recognized: ' + str(beamline) + ' ' + str(device))
        else:
          raise ValueError('facility not recognized: ' + str(facility))
    
        return




    def set_gap (self, gap):
        """Set the gap for a device.  Sets basic magnetic field data"""

        self.gap = gap

        self.osr.clear_bfields()
        mapping = oscars.util.read_file_list_with_header(self.bfield_lut1d_file)
        self.osr.add_bfield_interpolated(mapping=mapping[0], iformat=mapping[1], rotations=mapping[2], translation=mapping[3], scale=mapping[4], parameter=gap)
        self.osr.set_new_particle()

        return



    def summary (self):
        """Print some summary information and plots"""

        if self.lut1d is not None:
            self.lut1d.get_gaps(show=True)

        oscars.plots_mpl.plot_bfield(self.osr)
        oscars.plots_mpl.plot_trajectory_position(self.osr.calculate_trajectory())
        oscars.plots_mpl.plot_spectrum(self.osr.calculate_spectrum(obs=[0, 0, 30], energy_range_eV=self.spectrum_range), figsize=[16, 4])

        return
