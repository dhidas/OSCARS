# To test the sr module flux functions and mpl plots

# Import the OSCARS SR module
import oscars.sr

# Import basic plot utilities (matplotlib).  You don't need these to run OSCARS, but it's used here for basic plots
from oscars.plots_mpl import *

# Create a new OSCARS object
osr = oscars.sr.sr()

# Set default number of threads
osr.set_nthreads_global(8)

# Clear any existing fields (just good habit in notebook style) and add an undulator field
osr.clear_bfields()
osr.add_bfield_undulator(bfield=[0, 1, 0], period=[0, 0, 0.049], nperiods=21)

# Just to check the field that we added seems visually correct

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

# Print info
osr.print_all()


# Calculate power density
power_density = osr.calculate_power_density_rectangle(plane='XY',
                                                      width=[0.04, 0.04],
                                                      npoints=[21, 21],
                                                      translation=[0, 0, 30],
                                                      ofile='sr_powerdensity.txt',
                                                      bofile='sr_powerdensity.dat')

# Plot power density
plot_power_density(power_density, ofile='sr_powerdensity.png', title='Power Density: Calculated')



# Get new sr object
osr = oscars.sr.sr()
power_density = osr.average_power_density(ifiles=['sr_powerdensity.txt'])
plot_power_density(power_density, title='Power Density: From text input')

# Get new sr object
osr = oscars.sr.sr()
power_density = osr.average_power_density(bifiles=['sr_powerdensity.dat'])
plot_power_density(power_density, title='Power Density: From binary input')

