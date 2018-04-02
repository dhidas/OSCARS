from oscars.fit import *
import sys

class SpecRange:
    range_eV = [0, 0]
    done = False
    flux = 0
    energy = 0
    fwhm = 0
    spectrum = []
    nparticles = 0

    def __init__ (self, range_eV):
        self.range_eV=range_eV

    def get_points_eV (self):
        pp = []
        if not self.done:
            for i in range(self.range_eV[0], self.range_eV[1]+1):
                pp.append(i)
        return pp

    def spec_test (self, h):
        # Compare spectrum to previous
        flux_test = abs((self.flux - h[0]) / self.flux) < 0.01
        energy_test = abs(self.energy - h[1]) / self.energy < 0.01
        fwhm_test = True #abs((self.fwhm - h[2]) / self.fwhm) < 0.01

        return flux_test and energy_test and fwhm_test
    
    def save_spectrum_range(self, sin):
        self.spectrum = []
        for s in sin:
            if s[0] >= self.range_eV[0] and s[0] <= self.range_eV[1]:
                self.spectrum.append(s)
                
        return


def get_points(sr):
    pp = []
    for s in sr:
        if not s.done:
            for p in s.get_points_eV():
                pp.append(p)
    return pp






def calculate_harmonics_me (osr, obs=[0, 0, 30], nparticles=1000, niterations=10):

    spectrum_se = osr.calculate_spectrum(obs=[0, 0, 30], energy_range_eV=[10, 30000])
    harmonics = find_odd_harmonics(spectrum_se, show=False)
    
    spectrum_summed = []
    specr = []
    for h in harmonics:
        specr.append(SpecRange([int(h[1] - 6*h[2]), int(h[1] + 3*h[2])]))

    total_particles = nparticles

    spectrum_summed = osr.calculate_spectrum(obs=obs,
                                             energy_points_eV=get_points(specr),
                                             nparticles=nparticles)

    for sr in specr:
        #print(sr.range_eV)
        new_peak = find_peaks_parabola(spectrum_summed, [sr.range_eV])
        if len(new_peak) != 0:
            sr.flux = new_peak[0][0]
            sr.energy = new_peak[0][1]
            sr.fwhm = new_peak[0][2]
        else:
            print('failing to see spectrum in init')

    for i in range(niterations):
        # Number of particles this round and total after this round is calculated
        nparticles = total_particles
        total_particles += nparticles

        print('iteration', i, 'nparticles', nparticles, total_particles)
        sys.stdout.flush()


        # Calculate spectrum and harmonics
        spectrum_i = osr.calculate_spectrum(obs=obs,
                                            energy_points_eV=get_points(specr),
                                            nparticles=nparticles)
        print('average lens', len(spectrum_i), len(spectrum_summed))
        sys.stdout.flush()

        spectrum_summed = osr.average_spectra(spectra=[spectrum_i, spectrum_summed])

        for sr in specr:
            if not sr.done:
                peak_i = find_peaks_parabola(spectrum_i, [sr.range_eV])
                if len(peak_i) == 0:
                    print('did not see spectrum in ith spectrum')


                if len(peak_i) != 0 and sr.spec_test(peak_i[0]):
                    sr.done = True
                    sr.save_spectrum_range(spectrum_summed)
                    sr.nparticles = total_particles
                peak_summed = find_peaks_parabola(spectrum_summed, [sr.range_eV])
                if len(peak_summed) == 0:
                    print('did not see spectrum in summed spectrum')

                if len(peak_summed) != 0:
                    sr.flux = peak_summed[0][0]
                    sr.energy = peak_summed[0][1]
                    sr.fwhm = peak_summed[0][2]

        spectrum_new = []
        for p in spectrum_summed:
            for sr in specr:
                if not sr.done and p[0] >= sr.range_eV[0] and p[0] <= sr.range_eV[1]:
                    #print('adding', p)
                    spectrum_new.append(p)
                    break
        print('sizes', len(spectrum_summed), len(spectrum_new))
        sys.stdout.flush()

        spectrum_summed = spectrum_new

        for sr in specr:
            print('is done', sr.range_eV, sr.done)
        if sum([not x.done for x in specr]) == 0:
            print('done')
            break

    return specr


