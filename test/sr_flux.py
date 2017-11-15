# To test the sr module flux functions and mpl plots

# Import the OSCARS SR module
import oscars.sr

# Import basic plot utilities (matplotlib).  You don't need these to run OSCARS, but it's used here for basic plots
from oscars.plots_mpl import *

# Create a new OSCARS object
osr = oscars.sr.sr(gpu=1, nthreads=16)

# Clear any existing fields (just good habit in notebook style) and add an undulator field
osr.clear_bfields()
osr.add_bfield_undulator(bfield=[0, 1, 0], period=[0, 0, 0.049], nperiods=21)

# Just to check the field that we added seems visually correct

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

# Print info
osr.print_all()


# Calculate spectrum zoom
flux = osr.calculate_flux_rectangle(plane='XY',
                                    energy_eV=152,
                                    width=[0.01, 0.01],
                                    npoints=[401, 401],
                                    translation=[0, 0, 30],
                                    ofile='sr_flux.txt',
                                    bofile='sr_flux.dat',
                                    quantity='flux')
plot_flux(flux, ofile='sr_flux.png', title='Flux: Calculated', show=False)
exit(0)



# Create a new OSCARS object
osr = oscars.sr.sr()
flux = osr.average_flux(ifiles=['sr_flux.txt'])
plot_flux(flux, title='Flux: From text input')


# Create a new OSCARS object
osr = oscars.sr.sr()
flux = osr.average_flux(bifiles=['sr_flux.dat'])
plot_flux(flux, title='Flux: From binary input')
