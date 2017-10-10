# Import the OSCARS SR module
import oscars.sr

# Import basic plot utilities (matplotlib).  You don't need these to run OSCARS, but it's used here for basic plots
from oscars.plots_mpl import *


# Create a new OSCARS SR object
osr = oscars.sr.sr()


# Clear any existing fields (just good habit in notebook style) and add an undulator field
osr.clear_bfields()
osr.add_bfield_undulator(bfield=[0, 1, 0], period=[0, 0, 0.049], nperiods=21)


# Just to check the field that we added seems visually correct
plot_bfield(osr)


# Setup beam similar to NSLSII
osr.clear_particle_beams()
osr.set_particle_beam(x0=[0, 0, -1], energy_GeV=3, current=0.500)

# Set the start and stop times for the calculation
osr.set_ctstartstop(0, 2)


# Run the particle trajectory calculation
trajectory = osr.calculate_trajectory()

# Plot the trajectory position and velocity
plot_trajectory_position(trajectory)
plot_trajectory_velocity(trajectory)




