
# coding: utf-8

# Plots inline for notebook
get_ipython().run_line_magic('matplotlib', 'inline')

# Import the OSCARS SR module
import oscars.sr

# Import the 3D and parametric surfaces utilities
from oscars.plots3d_mpl import *


# Create a new OSCARS object.  Default to 8 threads and always use the GPU if available
osr = oscars.sr.sr(nthreads=8, gpu=1)


# Clear all existing fields and create an undulator field
osr.clear_bfields()
osr.add_bfield_undulator(bfield=[0, 1, 0], period=[0, 0, 0.050], nperiods=31)


# Define simple electron beam
osr.set_particle_beam(energy_GeV=3, x0=[0, 0, -1], current=0.5)

# Define the start and stop times for the calculation
osr.set_ctstartstop(0, 2)


# Clculate power density on CAD surface from STL file.
# Here 'scale' is used since the input file units are in mm
# while OSCARS expects m (there is no STL standard so you must
# keep track of this yourself)
pd = osr.calculate_power_density_stl(
    ifile='sphere.stl',
    scale=0.001,
    translation=[0, 0, 5]
)

plot_power_density_stl(pd)

