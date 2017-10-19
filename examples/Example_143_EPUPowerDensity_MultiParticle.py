# Import the OSCARS SR module
import oscars.sr

# Import basic plot utilities (matplotlib).  You don't need these to run OSCARS, but it's used here for basic plots
from oscars.plots_mpl import *


# Create a new OSCARS object.  Default to 8 threads and always use the GPU if available
osr = oscars.sr.sr(nthreads=8, gpu=1)


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
osr.set_particle_beam(
    energy_GeV=3,
    x0=[0, 0, -1],
    current=0.500,
    sigma_energy_GeV=0.001*3,
    beta=[1.5, 0.8],
    emittance=[0.9e-9, 0.008e-9]
)

# Set the start and stop times for the calculation
osr.set_ctstartstop(0, 2)


# Run the particle trajectory calculation for the ideal particle
# First must load the ideal particle initial conditions
osr.set_new_particle(particle='ideal')
trajectory = osr.calculate_trajectory()

# Plot the trajectory position and velocity
plot_trajectory_position(trajectory)
plot_trajectory_velocity(trajectory)


# Calculate the multi-particle Flux.  We'll pick a point around 154 [eV]
# Although 5 particles are used here, you shuold use a lot more.
pd = osr.calculate_power_density_rectangle(
    plane='XY',
    width=[0.05, 0.05],
    npoints=[101, 101],
    translation=[0, 0, 30],
    nparticles=5
)

plot_power_density(pd)

