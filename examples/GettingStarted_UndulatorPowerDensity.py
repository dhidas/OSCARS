
# coding: utf-8

# Plots inline for notebook
get_ipython().run_line_magic('matplotlib', 'inline')

# Import the OSCARS SR module
import oscars.sr

# Optionally import the plotting tools (matplotlib)
from oscars.plots_mpl import *


# Create a new OSCARS object.  Default to 8 threads and always use the GPU if available
osr = oscars.sr.sr(nthreads=8, gpu=1)


# Clear all existing fields and create an undulator field
osr.clear_bfields()
osr.add_bfield_undulator(bfield=[0, 1, 0], period=[0, 0, 0.050], nperiods=31)


# Define simple electron beam
osr.set_particle_beam(energy_GeV=3, x0=[0, 0, -1], current=0.5)

# Define the start and stop times for the calculation
osr.set_ctstartstop(0, 2)


# Calculate spectrum at 30 [m].  Note use of the nthreads argument.
power_density = osr.calculate_power_density_rectangle(
    plane='XY',
    width=[0.05, 0.05],
    npoints=[101, 101],
    translation=[0, 0, 30]
)


# Plot power density
plot_power_density(power_density)

