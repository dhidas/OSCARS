# This is meant to be run on OSG or similar condor-type systems


# In this example MPI is used to calculate the spectrum in a multi-particle
# simulation.  The rank 0 process calculates the ideal single-particle
# spectrum and then waits for the other processes to return their data.
# When a process returns the results are added to the others.  The results
# are then plotted together and saved as a png.  If you want the plot to
# show you can remove the show=False

# Command line arguments are given by condor at runtime
import sys
sys.path.append('.')

Process = sys.argv[1]
NAME = sys.argv[2]

out_file_name = NAME + '_' + Process + '.dat'

# Import the OSCARS SR and helper modules
import oscars.sr


# Get an OSCARS SR object
osr = oscars.sr.sr()

# Let's just make sure that each process only uses 1 threads since
# we assume mpi is handeling this
osr.set_nthreads_global(1)

# Set a particle beam with non-zero emittance
osr.set_particle_beam(type='electron',
                      name='beam_0',
                      energy_GeV=3,
                      x0=[0, 0, -1],
                      d0=[0, 0, 1],
                      current=0.500,
                      sigma_energy_GeV=0.001*3,
                      beta=[1.5, 0.8],
                      emittance=[0.9e-9, 0.008e-9],
                      horizontal_direction=[1, 0, 0],
                      lattice_reference=[0, 0, 0])

# Must set the start and stop time for calculations
osr.set_ctstartstop(0, 2)

# Add an undulator field
osr.add_bfield_undulator(bfield=[0, 1, 0], period=[0, 0, 0.049], nperiods=31)

# Number of particles per node of rank > 1
particles_per_node = 2

# Observation point for spectrum
observation_point = [0, 0, 30]

# Number of points in the spectrum
npoints = 500

# Energy range for spectrum
range_eV = [100, 500]

# If not rank 0, calculate the desired spectrum
data = osr.calculate_spectrum(obs=observation_point,
                              energy_range_eV=range_eV,
                              npoints=npoints,
                              nparticles=particles_per_node,
                              ofile=out_file_name)





