#This is meant to be run on OSG or similar condor-type systems

# In this example OSG is used to calculate a multi-particle
# simulation.  The rank 0 process calculates the ideal single-particle
# data and then waits for the other processes to return their data.
# When a process returns the results are added to the others.  The results
# are then outputted as data files, which can be plotted using Jupyter
# Notebook or something similar.

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

# Set particle beam with non-zeros emittance
osr.set_particle_beam(type='electron',
                      name='beam_0',
                      energy_GeV=3,
                      x0=[0,0,-1],
                      d0=[0,0,1],
                      current=0.500,
                      sigma_energy_GeV=0.001*3,
                      beta=[1.5,0.8],
                      emittance=[0.9e-9,0.008e-9],
                      horizontal_direction=[1,0,0],
                      lattice_reference=[0,0,0])

# Must set the start and stop time for calculations
osr.set_ctstartstop(0,2)

# Phase difference between fields in [rad]
phase = osr.pi()/2

# Clear any exsisting fields (just good habit in notebook style) and add
# an undulator field
osr.clear_bfields()
osr.add_bfield_undulator(bfield=[0,0.7,0],
                         period=[0,0,0.049],
                         nperiods=21,
                         phase=-phase/2)
osr.add_bfield_undulator(bfield=[0.7,0,0],
                         period=[0,0,0.049],
                         nperiods=21,
                         phase=+phase/2)

# Number of particles per node of rank > 1
particles_per_node = 1000

# Observation point for spectrum 
observation_point = [0,0,30]

# Number of points in the spectrum
npoints = 1000

# Energy range for spectrum
range_eV = [130, 180]

# Ideal single-particle data
if int(Process) == 0:
  osr.set_new_particle(particle='ideal')
  data = osr.calculate_spectrum(obs=observation_point,
                                energy_range_eV=range_eV,
                                npoints=npoints,
                                ofile=out_file_name)

# Multi-particle simulation
else:
  data = osr.calculate_spectrum(obs=observation_point,
                                energy_range_eV=range_eV,
                                npoints=npoints,
                                nparticles=particles_per_node,
                                ofile=out_file_name)
