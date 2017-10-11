# Import the OSCARS SR module
import oscars.sr

# Import basic plot utilities (matplotlib).  You don't need these to run OSCARS, but it's used here for basic plots
from oscars.plots_mpl import *


# Create a new OSCARS object.  Default to 8 threads and always use the GPU if available
osr = oscars.sr.sr(nthreads=8, gpu=1)


# Clear any existing fields (just good habit in notebook style) and add an undulator field
osr.clear_bfields()
dist_between_centers = 0.049*15

osr.add_bfield_undulator(
    name='Undulator_1',
    bfield=[0, 1, 0],
    period=[0, 0, 0.049],
    nperiods=21,
    translation=[0, 0, -dist_between_centers]
)

osr.add_bfield_undulator(
    name='Undulator_2',
    bfield=[0, 1, 0],
    period=[0, 0, 0.049],
    nperiods=21,
    translation=[0, 0, +dist_between_centers]
)

osr.add_bfield_uniform(name='BM_1', bfield=[0, -0.4, 0], width=[0, 0, 0.5], translation=[0, 0, -2])
osr.add_bfield_uniform(name='BM_2', bfield=[0, -0.4, 0], width=[0, 0, 0.5], translation=[0, 0, +2])


# Just to check the field that we added seems visually correct
plot_bfield(osr, -2, 2)


osr.print_bfields()


# Setup beam similar to NSLSII
osr.clear_particle_beams()
osr.set_particle_beam(energy_GeV=3, current=0.500)

# Set the start and stop times for the calculation
osr.set_ctstartstop(-1.8, 1.8)


# Run the particle trajectory calculation
trajectory = osr.calculate_trajectory()

# Plot the trajectory position and velocity
plot_trajectory_position(trajectory)
plot_trajectory_velocity(trajectory)




