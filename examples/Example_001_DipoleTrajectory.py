# Import the OSCARS SR module
import oscars.sr

# Import basic plot utilities.  You don't need these to run OSCARS, but it's used here for basic plots
from oscars.plots_mpl import *

# Create a new OSCARS object
osr = oscars.sr.sr()

# Clear any existing fields (just good habit) and add an undulator field
osr.clear_bfields()
osr.add_bfield_uniform(bfield=[0, -0.4, 0], width=[0, 0, 1])

# Just to check the field that we added seems visually correct
plot_bfield(osr, -1, 1)

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


## Alternate Starting Position (and backward propogation)

# Setup beam similar to NSLSII
osr.clear_particle_beams()
osr.set_particle_beam(type='electron',
                      name='beam_0',
                      x0=[0, 0, 0],
                      d0=[0, 0, 1],
                      energy_GeV=3,
                      current=0.500
                     )

# Set the start and stop times for the calculation
osr.set_ctstartstop(-1, 1)

# Run the particle trajectory calculation
trajectory = osr.calculate_trajectory()

# Plot the trajectory position and velocity
plot_trajectory_position(trajectory)
plot_trajectory_velocity(trajectory)

