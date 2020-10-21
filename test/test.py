import oscars.sr


osr = oscars.sr.sr(nthreads=1, gpu=0)

osr.set_particle_beam(beam='NSLSII', ctstartstop=[-1, 1])

osr.set_trajectory_calculation('RKAS')
osr.clear_bfields()
osr.add_bfield_uniform(bfield=[0, 1, 0])

pd_short = osr.calculate_power_density_rectangle(
    plane='XY',
    width=[0.05, 0.02],
    npoints=[2, 2],
    translation=[0, 0, 30],
    max_level=4,
)

print(pd_short)
