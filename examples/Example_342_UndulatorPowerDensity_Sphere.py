# Import the OSCARS SR module
import oscars.sr

# Import basic plot utilities (matplotlib).  You don't need these to run OSCARS, but it's used here for
# basic plots.  Additionally import oscars parametric surfaces
from oscars.plots_mpl import *
from oscars.plots3d_mpl import *
from oscars.parametric_surfaces import *


# Create a new OSCARS object.  Default to 8 threads and always use the GPU if available
osr = oscars.sr.sr(nthreads=8, gpu=1)


# Clear any existing fields (just good habit in notebook style) and add an undulator field
osr.clear_bfields()
osr.add_bfield_undulator(bfield=[0, 1, 0], period=[0, 0, 0.049], nperiods=21)

# Just to check the field that we added seems visually correct
plot_bfield(osr)


# Setup beam similar to NSLSII
osr.clear_particle_beams()
osr.set_particle_beam(
    x0=[0, 0, -1],
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


# First create a sphere
sphere = PSSphere(R=0.01, nu=51, nv=51)

# Next run the calclation and plot
power_density_3d(osr, sphere, translation=[0, 0, 30])




