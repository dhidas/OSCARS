# Import the OSCARS SR module
import oscars.sr

# Optionally import the plotting tools (matplotlib)
from oscars.plots_mpl import *

# Create an OSCARS SR object
osr = oscars.sr.sr()

# Clear all existing fields and create an undulator field
osr.clear_bfields()
osr.add_bfield_undulator(bfield=[0, 1, 0], period=[0, 0, 0.050], nperiods=31)

# Define simple electron beam
osr.set_particle_beam(type='electron', name='beam_0', energy_GeV=3, x0=[0, 0, -1], d0=[0, 0, 1], current=0.5)

# Define the start and stop times for the calculation
osr.set_ctstartstop(0, 2)

# If the GPU is available, let's use it from now on
if osr.check_gpu() > 0:
    osr.set_gpu_global(1)

# Calculate spectrum at 30 [m].  Note use of the nthreads argument.
power_density = osr.calculate_power_density_rectangle(plane='XY',
                                                      width=[0.05, 0.05],
                                                      npoints=[101, 101],
                                                      translation=[0, 0, 30], nthreads=8)

# Plot power density
plot_power_density(power_density)
