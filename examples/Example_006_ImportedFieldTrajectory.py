# Import the OSCARS SR module
import oscars.sr

# Import basic plot utilities (matplotlib).  You don't need these to run OSCARS, but it's used here for basic plots
from oscars.plots_mpl import *

# Create a new OSCARS SR object
osr = oscars.sr.sr()

# Clear any existing fields (just good habit in notebook style) and add an undulator field
osr.clear_bfields()
osr.add_bfield_undulator(bfield=[0, 1, 0], period=[0, 0, 0.049], nperiods=41)

# Add some random imperfections
for i in range(10):
    osr.add_bfield_gaussian(bfield=[0, 0.001 * (osr.rand() - 0.5), 0],
                            sigma=[0, 0, 0.1*(osr.rand())],
                            translation=[0, 0, 2 * (osr.rand() - 0.5)])


# Export the field
osr.write_bfield(ofile='EX006.dat', oformat='OSCARS', zlim=[-3, 3], nz=50000)


# Clear all magnetic fields!
osr.clear_bfields()

# Import the field from the file created above
osr.add_bfield_file(ifile='EX006.dat', iformat='OSCARS')

# Plot imported field
plot_bfield(osr, -2, 2)

# Setup beam similar to NSLSII with different starting position from above
# (this makes more sense for some scenarios)
osr.clear_particle_beams()
osr.set_particle_beam(type='electron',
                      name='beam_0',
                      x0=[0, 0, -3],
                      d0=[0, 0, 1],
                      energy_GeV=3,
                      current=0.500
                     )

# Set the start and stop times for the calculation
osr.set_ctstartstop(0, 6)

# Run the particle trajectory calculation
trajectory = osr.calculate_trajectory()

# Plot the trajectory position and velocity
plot_trajectory_position(trajectory)
plot_trajectory_velocity(trajectory)


