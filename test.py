import oscars.th

oth = oscars.th.th()

print(oth.undulator_K(1, 0.055))



oth.dipole_spectrum(bfield=0.4, beam_energy_GeV=3, angle=0.0001, energy_range_eV=[10, 1000])
