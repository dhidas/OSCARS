# Import the OSCARS SR module
import oscars.sr

# Import basic plot utilities (matplotlib).  You don't need these to run OSCARS, but it's used here for basic plots
from oscars.plots_mpl import *

# Create a new OSCARS object
osr = oscars.sr.sr()

# If you want to make this calculation go faster uncomment the following line
# osr.set_nthreads_global(8)

# Clear any existing fields (just good habit in notebook style) and add an undulator field
osr.clear_bfields()
osr.add_bfield_undulator(bfield=[0, 1, 0], period=[0, 0, 0.049], nperiods=21)

# Just to check the field that we added seems visually correct
plot_bfield(osr, -1, 1)

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

# Run the particle trajectory calculation
trajectory = osr.calculate_trajectory()

# Plot the trajectory position and velocity
plot_trajectory_position(trajectory)
plot_trajectory_velocity(trajectory)

# Calculate spectrum zoom
spectrum = osr.calculate_spectrum(obs=[0, 0, 30], energy_range_eV=[145, 160], npoints=200)
plot_spectrum(spectrum)

flux = osr.calculate_flux_rectangle(plane='XY',
                                    energy_eV=153,
                                    width=[0.01, 0.01],
                                    npoints=[51, 51],
                                    translation=[0, 0, 30]
                                   )
plot_flux(flux)
