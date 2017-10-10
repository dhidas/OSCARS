# Import the OSCARS SR module
import oscars.sr

# Import basic plot utilities.  You don't need these to run OSCARS, but it's used here for basic plots
from oscars.plots_mpl import *


# Create a new OSCARS object.  Default to 8 threads and always use the GPU if available
osr = oscars.sr.sr(nthreads=8, gpu=1)


# Clear any existing fields (just good habit in notebook style) and add an undulator field
osr.clear_bfields()
osr.add_bfield_uniform(bfield=[0, -0.4, 0], width=[0, 0, 1])


# Just to check the field that we added seems visually correct
plot_bfield(osr)


# Setup beam similar to NSLSII
osr.clear_particle_beams()
osr.set_particle_beam(x0=[0, 0, -1], energy_GeV=3, current=0.500)

# Set the start and stop times for the calculation
osr.set_ctstartstop(0, 2)


# Verify input information - print all to screen
osr.print_all()


# Run the particle trajectory calculation
trajectory = osr.calculate_trajectory()

# Plot the trajectory position and velocity
plot_trajectory_position(trajectory)
plot_trajectory_velocity(trajectory)


# Setup beam similar to NSLSII
osr.clear_particle_beams()
osr.set_particle_beam(energy_GeV=3, current=0.500)

# Set the start and stop times for the calculation
osr.set_ctstartstop(-1, 1)


# Run the particle trajectory calculation
trajectory = osr.calculate_trajectory()

# Plot the trajectory position and velocity
plot_trajectory_position(trajectory)
plot_trajectory_velocity(trajectory)




