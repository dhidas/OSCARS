# To test the sr module spectra functions and mpl plots

# Import the OSCARS SR module
import oscars.sr

# Import basic plot utilities
from oscars.plots_mpl import *

# Create a new OSCARS object
osr = oscars.sr.sr()

# Set default number of threads
osr.set_nthreads_global(8)

# Clear any existing fields (just good habit in notebook style) and add an undulator field
osr.clear_bfields()
osr.add_bfield_undulator(bfield=[0, 1, 0], period=[0, 0, 0.049], nperiods=21)

# Setup beam similar to NSLSII
osr.clear_particle_beams()
osr.set_particle_beam(type='electron',
                      name='beam_0',
                      x0=[0, 0, -1],
                      d0=[0, 0, 1],
                      energy_GeV=3,
                      current=0.500
                     )

# Set the start and stop times for the calculation
osr.set_ctstartstop(0, 2)

# Print settings
osr.print_all()


# Calculate spectrum
spectrum = osr.calculate_spectrum(obs=[0, 0, 30],
                                  energy_range_eV=[100, 200],
                                  npoints=200,
                                  ofile='sr_spectrum.txt',
                                  bofile='sr_spectrum.dat')

# Plot spectrum
plot_spectrum(spectrum, ofile='sr_spectrum.png', title='Spectrum: Calculated')



# Get new sr object, plot txt output
osr = oscars.sr.sr()
spectrum = osr.average_spectra(ifiles=['sr_spectrum.txt'])
plot_spectrum(spectrum, title='Spectrum: From text input')



# Get new sr object, plot binary output
osr = oscars.sr.sr()
spectrum = osr.average_spectra(bifiles=['sr_spectrum.dat'])
plot_spectrum(spectrum, title='Spectrum: From binary input')
