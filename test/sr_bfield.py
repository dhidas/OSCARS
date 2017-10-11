# To test the sr module spectra functions and mpl plots

# Import the OSCARS SR module and plotting
import oscars.sr
from oscars.plots_mpl import *

# Import basic plot utilities
from oscars.plots_mpl import *

# Create a new OSCARS object
osr = oscars.sr.sr()

# Set default number of threads
osr.set_nthreads_global(8)

# Clear any existing fields (just good habit in notebook style) and add an undulator field
osr.clear_bfields()
osr.add_bfield_undulator(bfield=[0, 1, 0], period=[0, 0, 0.049], nperiods=21)

# Write binary bfield
osr.write_bfield(bofile='sr_bfield.dat', oformat='OSCARS', zlim=[-1, 1], nz=5000)

# Clear fields and import the binary file
osr.clear_bfields()
osr.add_bfield_file(bifile='sr_bfield.dat', iformat='OSCARS')

plot_bfield(osr)
