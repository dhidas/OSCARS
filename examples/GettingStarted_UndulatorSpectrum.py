
# coding: utf-8

# Plots inline for notebook
get_ipython().run_line_magic('matplotlib', 'inline')

# Import the OSCARS SR module
import oscars.sr


# Create a new OSCARS object.  Default to 8 threads and always use the GPU if available
osr = oscars.sr.sr(nthreads=8, gpu=1)


# Clear all existing fields and create an undulator field
osr.clear_bfields()
osr.add_bfield_undulator(bfield=[0, 1, 0], period=[0, 0, 0.050], nperiods=31)


# Define simple electron beam
osr.set_particle_beam(energy_GeV=3, x0=[0, 0, -1], current=0.5)

# Define the start and stop times for the calculation
osr.set_ctstartstop(0, 2)


# Calculate spectrum at 30 [m]
spectrum = osr.calculate_spectrum(obs=[0, 0, 30], energy_range_eV=[100, 2000])


# Optionally import the plotting tools (matplotlib)
from oscars.plots_mpl import *

# Plot spectrum
plot_spectrum(spectrum)

