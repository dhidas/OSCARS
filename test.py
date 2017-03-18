import oscars.th

oth = oscars.th.th()

print(oth.undulator_K(1, 0.055))



#result = oth.dipole_spectrum(bfield=0.4, beam_energy_GeV=3, angle=0.0001, energy_eV=2390)


result = oth.undulator_flux(bfield=0.8, period=0.022, nperiods=41, beam_energy_GeV=3, angle_vertical=-0.0001, angle_horizontal=0.0001, energy_eV=2390)
