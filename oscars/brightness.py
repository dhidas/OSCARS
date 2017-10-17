class Undulator:
    name = ''
    period = 0
    length = 0
    k_range = [0, 0]
    beta = [0, 0]
    alpha = [0, 0]
    minimum = 0
    npoints = 1000
    harmonic_range = [1, 51]
    eta = [0, 0]
    
    def __init__ (self, name, period, length, k_range, beta, minimum, npoints = 1000, harmonic_range = [1, 51], alpha=[0, 0], eta = [0, 0]):
        self.name = name
        self.period = period
        self.length = length
        self.k_range = k_range
        self.beta = beta
        self.alpha = alpha
        self.minimum = minimum
        self.npoints = npoints
        self.harmonic_range = harmonic_range
        self.eta = eta
        
class BendingMagnet:
    name = ''
    bfield = 0
    beta = [0, 0]
    alpha = [0, 0]
    eta = [0, 0]
    npoints = 1000

    def __init__ (self, name, bfield, beta, eta, alpha=[0, 0], npoints = 1000):
        self.name = name
        self.bfield = bfield
        self.beta = beta
        self.alpha = alpha
        self.eta = eta
        self.npoints = npoints


class Wiggler:
    name = ''
    period = 0
    length = 0
    bfield = 0
    beta = [0, 0]
    alpha = [0, 0]
    eta = [0, 0]
    npoints = 1000

    def __init__ (self, name, period, length, bfield, beta, eta, alpha=[0, 0], npoints = 1000):
        self.name = name
        self.period = period
        self.length = length
        self.bfield = bfield
        self.beta = beta
        self.alpha = alpha
        self.eta = eta
        self.npoints = npoints

        
from scipy.interpolate import interp1d

class Synchrotron:
    name = ''
    beam_energy_GeV = 0
    sigma_energy_GeV = 0
    current = 0
    emittance = [0, 0]
    energy_range_eV = [0, 0]
    devices = []
    curves = []
    
    def __init__ (self, name, beam_energy_GeV, sigma_energy_GeV, current, emittance, energy_range_eV, devices):
        self.name = name
        self.beam_energy_GeV = beam_energy_GeV
        self.sigma_energy_GeV = sigma_energy_GeV
        self.current = current
        self.emittance = emittance
        self.energy_range_eV = energy_range_eV
        self.devices = devices
        self.curves = []

        return
            
    def add_device(self, device):
        self.devices.append(device)
        return
    
    def add_devices(self, devices):
        self.devices += devices
        
    def get_brightness_curve(self, energy_range_eV=None):
        self.curves = []
        for d in self.devices:
            oth.set_particle_beam(
                energy_GeV=self.beam_energy_GeV,
                sigma_energy_GeV=self.sigma_energy_GeV,
                current=self.current,
                emittance=self.emittance,
                beta=d.beta,
                alpha=d.alpha,
                eta=d.eta
            )

            if isinstance(d, BendingMagnet):
                bm = oth.dipole_brightness(
                    bfield=d.bfield,
                    energy_range_eV=self.energy_range_eV,
                    npoints=d.npoints
                )
                
                x = [b[0] for b in bm]
                y = [b[1] for b in bm]
                f = interp1d(x, y, kind='cubic')
                
                self.curves.append([self.energy_range_eV, f])


            elif isinstance(d, Wiggler):
                br = oth.dipole_brightness(
                    bfield=d.bfield,
                    energy_range_eV=self.energy_range_eV,
                    npoints=d.npoints
                )
                
                x = [b[0] for b in br]
                y = [b[1] * 2. * int(d.length / d.period) for b in br]
                f = interp1d(x, y, kind='cubic')
                
                self.curves.append([self.energy_range_eV, f])


            elif isinstance(d, Undulator):
                nperiods = int(d.length / d.period)
                for i in range(d.harmonic_range[0], d.harmonic_range[1]+1):
                    if i % 2 == 0:
                        continue

                    br = oth.undulator_brightness(
                        K_range=d.k_range,
                        period=d.period,
                        nperiods=nperiods,
                        harmonic=i,
                        npoints=d.npoints,
                        minimum=d.minimum
                    )

                    x = [b[0] for b in br]
                    y = [b[1] for b in br]
                    f = interp1d(x, y, kind='cubic')
                    self.curves.append([[x[0], x[-1]], f])

            else:
                raise ValueError(d.__class__.__name__ + ' is an invalid type.  please check')

        return
    
    
def plot_brightness (experiments, energy_range_eV, figsize=[16, 12], ylim=None):

    plt.figure(1, figsize=figsize)
    
    for exp in experiments:
        print(exp.name)
        exp.get_brightness_curve(energy_range_eV)
        X = []
        Y = []
        
        for x in np.linspace(energy_range_eV[0], energy_range_eV[1], 5000):
            y = 0

            for curve in exp.curves:
                if x >= curve[0][0] and x <= curve[0][1] and curve[1](x) > y:
                    y = curve[1](x)
                    
            if y is not 0:
                X.append(x)
                Y.append(y)
                
        plt.plot(X, Y, label=exp.name)

    yl = plt.ylim()
    if ylim is not None:
        plt.ylim(ylim)
    plt.grid()
    plt.loglog()
    plt.legend()
    plt.show()
    return plt


