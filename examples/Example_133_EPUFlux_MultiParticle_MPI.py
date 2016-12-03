# How to run this code:
#   mpirun -n 5 python Example_133_EPUFlux_MultiParticle_MPI.py


# In this example MPI is used to calculate the flux in a multi-particle
# simulation.  The rank 0 process calculates the ideal single-particle
# flux and then waits for the other processes to return their data.
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

# Phase difference between fields in [rad]
phase = osr.pi()/2.

# Clear any existing fields (just good habit in notebook style) and add an undulator field
osr.clear_bfields()
osr.add_bfield_undulator(bfield=[0, 0.7, 0], period=[0, 0, 0.049], nperiods=21, phase=-phase/2.)
osr.add_bfield_undulator(bfield=[0.7, 0, 0], period=[0, 0, 0.049], nperiods=21, phase=+phase/2.)

# Number of particles per node of rank > 1
particles_per_node = 1

# Translation for center of rectangle
rectangle_center = [0, 0, 30]

# Width of rectangle
width = [0.01, 0.01]

# Number of points in flux
npoints = [51, 51]

# Energy we are interested in
energy_eV = 154



# Now the MPI fun
if rank == 0:
    # For rank 0 we calculate the ideal single-particle flux
    osr.set_new_particle(particle='ideal')
    flux = osr.calculate_flux_rectangle(plane='XY',
                                        energy_eV=energy_eV,
                                        width=width,
                                        npoints=npoints,
                                        translation=rectangle_center
                                       )

    plot_flux(flux, ofile='EPUFlux_Ideal.png', show=False)

    # Weight for each flux in summation
    weight = 1. / size

    # Get a new OSCARS SR object to collect the data in
    osr_sum = oscars.sr.sr()

    # Now wait and collect data from all other processes when it comes in
    for i in range(1, size):
        # Get incoming data
        data = comm.recv(source=ANY_SOURCE)

        # Sum the flux (internally this is a compensated sum)
        osr_sum.add_to_flux(flux=data, weight=weight)

    # Plot the single-particle and multi-particle data and save to file
    plot_flux(osr_sum.get_flux(), ofile='EPUFlux_MultiParticle.png', show=False)



else:
    # If not rank 0, calculate the desired flux
    data = osr.calculate_flux_rectangle(plane='XY',
                                        energy_eV=energy_eV,
                                        width=width,
                                        npoints=npoints,
                                        translation=rectangle_center,
                                        nparticles=particles_per_node,
                                       )

    # Send results back to rank 0
    comm.send(data, dest=0)




