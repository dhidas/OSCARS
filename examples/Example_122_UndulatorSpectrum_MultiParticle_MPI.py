# How to run this code:
#   mpiexec -n 5 python Example_122_UndulatorSpectrum_MultiParticle_MPI.py


# In this example MPI is used to calculate the spectrum in a multi-particle
# simulation.  The rank 0 process calculates the ideal single-particle
# spectrum and then waits for the other processes to return their data.
# When a process returns the results are added to the others.  The results
# are then plotted together and saved as a png.  If you want the plot to
# show you can remove the show=False


# Import the OSCARS SR and helper modules
import oscars.sr
from oscars.plots_mpl import *

# MPI imports
from mpi4py import MPI
from mpi4py.MPI import ANY_SOURCE

# Common MPI communication, rank, size
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

# Get an OSCARS SR object with nthreads=1 (meaning 1 thread per mpi execution)
osr = oscars.sr.sr(nthreads=1)

# Set a particle beam with non-zero emittance
osr.set_particle_beam(
    energy_GeV=3,
    x0=[0, 0, -1],
    current=0.500,
    sigma_energy_GeV=0.001*3,
    beta=[1.5, 0.8],
    emittance=[0.9e-9, 0.008e-9]
)

# Must set the start and stop time for calculations
osr.set_ctstartstop(0, 2)

# Add an undulator field
osr.add_bfield_undulator(bfield=[0, 1, 0], period=[0, 0, 0.049], nperiods=31)

# Number of particles per node of rank > 1
particles_per_node = 20

# Observation point for spectrum
observation_point = [0, 0, 30]

# Number of points in the spectrum
npoints = 500

# Energy range for spectrum
range_eV = [100, 500]




# Now the MPI fun
if rank == 0:
    # For rank 0 we calculate the ideal single-particle spectrum
    osr.set_new_particle(particle='ideal')
    spectrum_se = osr.calculate_spectrum(obs=observation_point, energy_range_eV=range_eV, npoints=npoints)

    # Weight for each spectrum in summation
    weight = 1. / size

    # Get a new OSCARS SR object to collect the data in
    osr_sum = oscars.sr.sr()

    # Now wait and collect data from all other processes when it comes in
    for i in range(1, size):
        # Get incoming data
        data = comm.recv(source=ANY_SOURCE)

        # Sum the spectra (internally this is a compensated sum)
        osr_sum.add_to_spectrum(spectrum=data, weight=weight)

    # Plot the single-particle and multi-particle data and save to file
    plot_spectra([spectrum_se, osr_sum.get_spectrum()], ['single-electron', 'multi-electron'], ofile='UndulatorSpectrum_MPI.png', show=False)



else:
    # If not rank 0, calculate the desired spectrum
    data = osr.calculate_spectrum(obs=observation_point, energy_range_eV=range_eV, npoints=npoints, nparticles=particles_per_node)

    # Send results back to rank 0
    comm.send(data, dest=0)





