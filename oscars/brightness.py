from __future__ import print_function

import oscars.th
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from itertools import cycle



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
    color = None
    linestyle = None
    
    def __init__ (self, name, period, length, k_range, beta, minimum, npoints = 1000, harmonic_range = [1, 51], alpha=[0, 0], eta = [0, 0], color=None, linestyle=None):
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
        self.color = color
        self.linestyle = linestyle
        
class BendingMagnet:
    name = ''
    bfield = 0
    beta = [0, 0]
    alpha = [0, 0]
    eta = [0, 0]
    npoints = 1000
    color = None
    linestyle = None

    def __init__ (self, name, bfield, beta, eta, alpha=[0, 0], npoints = 1000, color=None, linestyle=None):
        self.name = name
        self.bfield = bfield
        self.beta = beta
        self.alpha = alpha
        self.eta = eta
        self.npoints = npoints
        self.color = color
        self.linestyle = color


class Wiggler:
    name = ''
    period = 0
    length = 0
    bfield = 0
    beta = [0, 0]
    alpha = [0, 0]
    eta = [0, 0]
    npoints = 1000
    color = None
    linestyle = None

    def __init__ (self, name, period, length, bfield, beta, eta, alpha=[0, 0], npoints = 1000, color=None, linestyle=None):
        self.name = name
        self.period = period
        self.length = length
        self.bfield = bfield
        self.beta = beta
        self.alpha = alpha
        self.eta = eta
        self.npoints = npoints
        self.color = color
        self.linestyle = linestyle

        

class Synchrotron:
    name = ''
    beam_energy_GeV = 0
    sigma_energy_GeV = 0
    current = 0
    emittance = [0, 0]
    energy_range_eV = [0, 0]
    devices = []
    curves = []
    color = None
    linestyle = None
    marker = None
    oth = None
    
    def __init__ (self, name, beam_energy_GeV, sigma_energy_GeV, current, emittance, energy_range_eV, devices, color=None, linestyle=None, marker=None):
        self.name = name
        self.beam_energy_GeV = beam_energy_GeV
        self.sigma_energy_GeV = sigma_energy_GeV
        self.current = current
        self.emittance = emittance
        self.energy_range_eV = energy_range_eV
        self.devices = devices
        self.curves = []
        self.oth = oscars.th.th()
        self.color = color
        self.linestyle = linestyle
        self.marker = marker

        return
            
    def add_device(self, device):
        self.devices.append(device)
        return
    
    def add_devices(self, devices):
        self.devices += devices
        
    def get_brightness_curves (self, energy_range_eV=None):
        self.curves = []
        for d in self.devices:
            self.oth.set_particle_beam(
                energy_GeV=self.beam_energy_GeV,
                sigma_energy_GeV=self.sigma_energy_GeV,
                current=self.current,
                emittance=self.emittance,
                beta=d.beta,
                alpha=d.alpha,
                eta=d.eta
            )

            if isinstance(d, BendingMagnet):
                bm = self.oth.dipole_brightness(
                    bfield=d.bfield,
                    energy_range_eV=self.energy_range_eV,
                    npoints=d.npoints
                )
                
                x = [b[0] for b in bm]
                y = [b[1] for b in bm]
                f = interp1d(x, y, kind='cubic')
                
                self.curves.append([d, self.energy_range_eV, f])


            elif isinstance(d, Wiggler):
                br = self.oth.dipole_brightness(
                    bfield=d.bfield,
                    energy_range_eV=self.energy_range_eV,
                    npoints=d.npoints
                )
                
                x = [b[0] for b in br]
                y = [b[1] * 2. * int(d.length / d.period) for b in br]
                f = interp1d(x, y, kind='cubic')
                
                self.curves.append([d, self.energy_range_eV, f])


            elif isinstance(d, Undulator):
                nperiods = int(d.length / d.period)
                for i in range(d.harmonic_range[0], d.harmonic_range[1]+1):
                    if i % 2 == 0:
                        continue

                    br = self.oth.undulator_brightness(
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
                    self.curves.append([d, [x[0], x[-1]], f])

            else:
                raise ValueError(d.__class__.__name__ + ' is an invalid type.  please check')

        return



    def get_flux_curves (self, energy_range_eV=None):
        self.curves = []
        for d in self.devices:
            self.oth.set_particle_beam(
                energy_GeV=self.beam_energy_GeV,
                sigma_energy_GeV=self.sigma_energy_GeV,
                current=self.current,
                emittance=self.emittance,
                beta=d.beta,
                alpha=d.alpha,
                eta=d.eta
            )

            if isinstance(d, BendingMagnet):
                bm = self.oth.dipole_spectrum(
                    bfield=d.bfield,
                    energy_range_eV=self.energy_range_eV,
                    npoints=d.npoints
                )
                
                x = [b[0] for b in bm]
                y = [b[1] for b in bm]
                f = interp1d(x, y, kind='cubic')
                
                self.curves.append([d, self.energy_range_eV, f])


            elif isinstance(d, Wiggler):
                br = self.oth.dipole_spectrum(
                    bfield=d.bfield,
                    energy_range_eV=self.energy_range_eV,
                    npoints=d.npoints
                )
                
                x = [b[0] for b in br]
                y = [b[1] * 2. * int(d.length / d.period) for b in br]
                f = interp1d(x, y, kind='cubic')
                
                self.curves.append([d, self.energy_range_eV, f])


            elif isinstance(d, Undulator):
                nperiods = int(d.length / d.period)
                for i in range(d.harmonic_range[0], d.harmonic_range[1]+1):
                    if i % 2 == 0:
                        continue

                    br = self.oth.undulator_flux(
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
                    self.curves.append([d, [x[0], x[-1]], f])

            else:
                raise ValueError(d.__class__.__name__ + ' is an invalid type.  please check')

        return




def plot_brightness (experiments, energy_range_eV, figsize=None, ylim=None, ofile=''):

    plt.figure(1, figsize=figsize)
    
    cc = cycle(['C0','C1','C2','C3','C4','C5','C6','C7','C8','C9'])

    for exp in experiments:
        labeldone = False
        color = next(cc) if exp.color is None else exp.color
        print(exp.name)
        exp.get_brightness_curves(energy_range_eV)
        X = []
        Y = []
        
        for xx in np.logspace(np.log(energy_range_eV[0]), np.log(energy_range_eV[1]), 5000):
            yy = 0

            for curve in exp.curves:
                if xx >= curve[1][0] and xx <= curve[1][1] and curve[2](xx) > yy:
                    yy = curve[2](xx)
                    
                
            if yy is not 0:
                X.append(xx)
                Y.append(yy)
            elif yy is 0 and len(X) is not 0:
                if labeldone:
                    plt.plot(X, Y, color=color, linestyle=exp.linestyle, marker=exp.marker)
                else:
                    plt.plot(X, Y, color=color, linestyle=exp.linestyle, marker=exp.marker, label=exp.name)
                    labeldone = True
                X=[]
                Y=[]

        if len(X) > 0:
            if labeldone:
                plt.plot(X, Y, color=color, linestyle=exp.linestyle, marker=exp.marker)
            else:
                plt.plot(X, Y, color=color, linestyle=exp.linestyle, marker=exp.marker, label=exp.name)
                labeldone = True

    yl = plt.ylim()
    if ylim is not None:
        plt.ylim(ylim)
    plt.grid()
    plt.loglog()
    plt.legend()
    plt.xlabel('Photon Energy [eV]')
    plt.ylabel('Brightness [$photons/s/0.1\%bw/mm^2/mrad^2$]')

    if ofile != '':
        plt.savefig(ofile, bbox_inches='tight')

    plt.show()
    return plt


def plot_flux (experiments, energy_range_eV, figsize=None, ylim=None, ofile=''):

    plt.figure(1, figsize=figsize)
    
    cc = cycle(['C0','C1','C2','C3','C4','C5','C6','C7','C8','C9'])

    for exp in experiments:
        labeldone = False
        print(exp.name)
        color = next(cc) if exp.color is None else exp.color
        exp.get_flux_curves(energy_range_eV)
        X = []
        Y = []
        
        for xx in np.linspace(energy_range_eV[0], energy_range_eV[1], 5000):
            yy = 0

            for curve in exp.curves:
                if xx >= curve[1][0] and xx <= curve[1][1] and curve[2](xx) > yy:
                    yy = curve[2](xx)
                    
                
            if yy is not 0:
                X.append(xx)
                Y.append(yy)
            elif yy is 0 and len(X) is not 0:
                if labeldone:
                    plt.plot(X, Y, color=color, linestyle=exp.linestyle)
                else:
                    plt.plot(X, Y, color=color, linestyle=exp.linestyle, label=exp.name)
                    labeldone = True
                X=[]
                Y=[]

        if len(X) > 0:
            if labeldone:
                plt.plot(X, Y, color=color, linestyle=exp.linestyle)
            else:
                plt.plot(X, Y, color=color, linestyle=exp.linestyle, label=exp.name)
                labeldone = True

    yl = plt.ylim()
    if ylim is not None:
        plt.ylim(ylim)
    plt.grid()
    plt.loglog()
    plt.legend()
    plt.xlabel('Photon Energy [eV]')
    plt.ylabel('Flux [$photons/s/0.1\%bw/mrad^2$]')

    if ofile != '':
        plt.savefig(ofile, bbox_inches='tight')

    plt.show()
    return plt


def plot_flux_all (experiments, energy_range_eV, figsize=None, ylim=None, ofile=''):

    plt.figure(1, figsize=figsize)
    
    cc = cycle(['C0','C1','C2','C3','C4','C5','C6','C7','C8','C9'])

    for exp in experiments:
        print(exp.name)
        for d in exp.devices:
            color = next(cc) if d.color is None else d.color
            exp.oth.set_particle_beam(
                energy_GeV=exp.beam_energy_GeV,
                sigma_energy_GeV=exp.sigma_energy_GeV,
                current=exp.current,
                emittance=exp.emittance,
                beta=d.beta,
                alpha=d.alpha,
                eta=d.eta
            )

            fl = []
            if isinstance(d, BendingMagnet):
                fl = exp.oth.dipole_spectrum(
                    bfield=d.bfield,
                    energy_range_eV=exp.energy_range_eV,
                    npoints=d.npoints,
                    #angle_integrated=True
                )
                
                x = [b[0] for b in fl]
                y = [b[1] for b in fl]

                plt.plot(x, y, color=color, linestyle=d.linestyle, label=d.name)

            elif isinstance(d, Wiggler):
                fl = exp.oth.dipole_spectrum(
                    bfield=d.bfield,
                    energy_range_eV=exp.energy_range_eV,
                    npoints=d.npoints
                )
                
                x = [b[0] for b in fl]
                y = [b[1] * 2. * int(d.length / d.period) for b in fl]
                f = interp1d(x, y, kind='cubic')
                
                plt.plot(x, y, color=color, linestyle=d.linestyle, label=d.name)


            elif isinstance(d, Undulator):
                labeldone = False
                nperiods = int(d.length / d.period)
                curves = []
                for i in range(d.harmonic_range[0], d.harmonic_range[1]+1):
                    if i % 2 == 0:
                        continue

                    fl = exp.oth.undulator_flux(
                        K_range=d.k_range,
                        period=d.period,
                        nperiods=nperiods,
                        harmonic=i,
                        npoints=d.npoints,
                        minimum=d.minimum
                    )

                    x = [b[0] for b in fl]
                    y = [b[1] for b in fl]
                    f = interp1d(x, y, kind='cubic')
                    curves.append([d, [x[0], x[-1]], f])

                X = []
                Y = []
                for xx in np.linspace(energy_range_eV[0], energy_range_eV[1], 5000):
                    yy = 0
                
                    for curve in curves:
                        if xx >= curve[1][0] and xx <= curve[1][1] and curve[2](xx) > yy:
                            yy = curve[2](xx)
                            
                    if yy is not 0:
                        X.append(xx)
                        Y.append(yy)
                    elif yy is 0 and len(X) is not 0:
                        if labeldone:
                            plt.plot(X, Y, color=color, linestyle=d.linestyle)
                        else:
                            plt.plot(X, Y, color=color, linestyle=d.linestyle, label=d.name)
                            labeldone = True
                        X=[]
                        Y=[]

                if len(X) > 0:
                    if labeldone:
                        plt.plot(X, Y, color=color, linestyle=d.linestyle)
                    else:
                        plt.plot(X, Y, color=color, linestyle=d.linestyle, label=d.name)
                        labeldone = True

            else:
                raise ValueError(d.__class__.__name__ + ' is an invalid type.  please check')

    yl = plt.ylim()
    if ylim is not None:
        plt.ylim(ylim)
    plt.grid()
    plt.loglog()
    plt.legend()
    plt.xlabel('Photon Energy [eV]')
    plt.ylabel('Flux [$photons/s/0.1\%bw/mrad^2$]')

    if ofile != '':
        plt.savefig(ofile, bbox_inches='tight')

    plt.show()
    return plt


