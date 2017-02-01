import sys
sys.path.append('.')

Process = sys.argv[1]
NAME = sys.argv[2]

out_file_name = NAME + '_' + Process + '.dat'

import oscars.sr

osr = oscars.sr.sr()

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

osr.set_ctstartstop(0,2)

osr.clear_bfields()
osr.add_bfield_undulator(bfield=[0,1,0],period=[0,0,0.049],nperiods=31)

particles_per_node = 1000

rcenter = [0, 0, 30]

width = [0.01,0.01]

npoints = [101, 101]

energy_eV = 152

if int(Process) == 0:
  osr.set_new_particle(particle='ideal')
  flux = osr.calculate_flux_rectangle(plane='XY',
                                      energy_eV=energy_eV,
                                      width=width,
                                      npoints=npoints,
                                      translation=rcenter,
                                      ofile=out_file_name)

else:
  data = osr.calculate_flux_rectangle(plane='XY',
                                       energy_eV=energy_eV,
                                       width=width,
                                       npoints=npoints,
                                       translation=rcenter,
                                       nparticles=particles_per_node,
                                       ofile=out_file_name)
  
