# Import the OSCARS SR module
import oscars.sr

# Import basic plot utilities (matplotlib).  You don't need these to run OSCARS, but it's used here for basic plots
from oscars.plots_mpl import *

# Create a new OSCARS object
osr = oscars.sr.sr()

# Clear any existing fields (just good habit in notebook style) and add an undulator field
osr.clear_bfields()
osr.add_bfield_undulator(bfield=[0, 1, 0], period=[0, 0, 0.049], nperiods=31)

# Just to check the field that we added seems visually correct
plot_bfield(osr, -1, 1)

# Setup beam similar to NSLSII
osr.set_particle_beam(type='electron',
                      name='beam_0',
                      energy_GeV=3,
                      x0=[0, 0, -1],
                      d0=[0, 0, 1],
                      current=0.500,
                      sigma_energy_GeV=0.001*3,
                      beta=[1.5, 0.8],
                      emittance=[0.9e-9, 0.008e-9],
                      horizontal_direction=[1, 0, 0],
                      lattice_reference=[0, 0, 0])

# Set the start and stop times for the calculation
osr.set_ctstartstop(0, 2)

# Run the particle trajectory calculation on the ideal case
osr.set_new_particle(particle='ideal')
trajectory = osr.calculate_trajectory()

# Plot the trajectory position and velocity
plot_trajectory_position(trajectory)
plot_trajectory_velocity(trajectory)

# Single particle spectrum for the ideal case
osr.set_new_particle(particle='ideal')
spectrum_se = osr.calculate_spectrum(obs=[0, 0, 30], energy_range_eV=[100, 500], npoints=500)

# Multi-particle spectrum
spectrum_me = osr.calculate_spectrum(obs=[0, 0, 30], energy_range_eV=[100, 500], npoints=500, nparticles=100, nthreads=8)

# Plot the two results together
plot_spectra([spectrum_se, spectrum_me], ['single-electron', 'multi-electron'], ofile='UndulatorSpectrum.png')
