import oscars.th

oth = oscars.th.th()

print(oth.undulator_K(1, 0.055))



#result = oth.dipole_spectrum(bfield=0.4, beam_energy_GeV=3, angle=0.0001, energy_eV=2390)


result = oth.undulator_flux_k(k=1, period=0.033, nperiods=41, beam_energy_GeV=3, harmonic=1)
print(result)
