# Import the OSCARS SR module
import oscars.sr

# Import basic plot utilities (matplotlib).  You don't need these to run OSCARS, but it's used here for basic plots
from oscars.plots_mpl import *

# Create a new OSCARS object
osr = oscars.sr.sr()

# Phase difference between fields in [rad]
phase = osr.pi()/2.

# Clear any existing fields (just good habit in notebook style) and add an undulator field
osr.clear_bfields()
osr.add_bfield_undulator(bfield=[0, 0.7, 0], period=[0, 0, 0.049], nperiods=21, phase=-phase/2.)
osr.add_bfield_undulator(bfield=[0.7, 0, 0], period=[0, 0, 0.049], nperiods=21, phase=+phase/2.)

# Just to check the field that we added seems visually correct
plot_bfield(osr)

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

# Run the particle trajectory calculation
trajectory = osr.calculate_trajectory()

# Plot the trajectory position and velocity
plot_trajectory_position(trajectory)
plot_trajectory_velocity(trajectory)

# Calculate Flux
power_densty = osr.calculate_power_density_rectangle(plane='XY',
                                                     width=[0.05, 0.05],
                                                     npoints=[51, 51],
                                                     translation=[0, 0, 30],
                                                     nthreads=8)
plot_power_density(power_densty)

