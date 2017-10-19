
# coding: utf-8

# matplotlib plots inline
get_ipython().run_line_magic('matplotlib', 'inline')

# Import the OSCARS SR module
import oscars.sr

# Import OSCARS plots (matplotlib)
from oscars.plots_mpl import *


# Create a new OSCARS object.  Default to 8 threads and always use the GPU if available
osr = oscars.sr.sr(nthreads=8, gpu=1)


# To illustrate the basic 1D format let's create a data file, then import it.
# It will be plotted before and after the import

# Create an undulator field
osr.clear_bfields()
osr.add_bfield_undulator(bfield=[0, 1, 0], period=[0, 0, 0.05], nperiods=31)
plot_bfield(osr)

# Now write the field to a file
osr.write_bfield(
    ofile='GettingStarted_OSCARS1D.dat',
    oformat='OSCARS1D Z By Bx Bz',
    zlim=[-1, 1],
    nz=5000
)   

# Clear fields and import the field from the file created above
osr.clear_bfields()
osr.add_bfield_file(ifile='GettingStarted_OSCARS1D.dat', iformat='OSCARS1D Z By Bx Bz')
plot_bfield(osr)


# To illustrate the basic 3D format let's create a data file, then import it.
# It will be plotted before and after the import

# Create an undulator field
osr.clear_bfields()
osr.add_bfield_undulator(bfield=[1, 0, 0], period=[0, 0, 0.050], nperiods=11)
plot_bfield(osr)

# Now write the field to a file
osr.write_bfield(
    ofile='GettingStarted_OSCARS.dat',
    oformat='OSCARS',
    xlim=[-1, 1],
    nx=2,
    ylim=[-1, 1],
    ny=2,
    zlim=[-1, 1],
    nz=5000
)

# Clear fields and import the field from the file created above
# You can also scale position, in case your input is not in [m]
osr.clear_bfields()
osr.add_bfield_file(ifile='GettingStarted_OSCARS.dat', iformat='OSCARS')

plot_bfield(osr)

