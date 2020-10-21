from __future__ import print_function

import oscars.th
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from itertools import cycle
import pickle
import os

from oscars.srwl_uti_brightness import srw_und_brightness, srw_epu_brightness, srw_und_flux, srw_epu_flux, srw_und_flux_onaxis, srw_epu_flux_onaxis

class Undulator:
    name = ''
    period = 0
    length = 0
    k_range = [0, 0]
    beta = [0, 0]
    alpha = [0, 0]
    minimum = 0.1
    npoints = 1000
    harmonic_range = [1, 51]
    eta = [0, 0]
    color = None
    linestyle = None
    calculation = 'SRW'
    
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
        
class EPU:
    name = ''
    period = 0
    length = 0
    kx_range = [0, 0]
    ky_range = [0, 0]
    beta = [0, 0]
    alpha = [0, 0]
    minimum = 0.1
    npoints = 1000
    eta = [0, 0]
    color = None
    linestyle = None
    calculation = 'SRW'
    
    def __init__ (self, name, period, length, kx_range, ky_range, beta, minimum, npoints = 1000, alpha=[0, 0], eta = [0, 0], color=None, linestyle=None):
        self.name = name
        self.period = period
        self.length = length
        self.kx_range = kx_range
        self.ky_range = ky_range
        self.beta = beta
        self.alpha = alpha
        self.minimum = minimum
        self.npoints = npoints
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
    linewidth = 1
    marker = None
    
    def __init__ (self, name, beam_energy_GeV, sigma_energy_GeV, current, emittance, devices, energy_range_eV=[10, 200000], color=None, linestyle=None, linewidth=1, marker=None):
        self.name = name
        self.beam_energy_GeV = beam_energy_GeV
        self.sigma_energy_GeV = sigma_energy_GeV
        self.current = current
        self.emittance = emittance
        self.energy_range_eV = energy_range_eV
        self.devices = devices
        self.curves = []
        self.color = color
        self.linestyle = linestyle
        self.linewidth = linewidth
        self.marker = marker

        return

    def __str__ (self):
        p =  'Synchrotron ' + self.name + '\n'
        p += '    beam_energy:     ' + str(self.beam_energy_GeV) + ' [GeV]\n'
        p += '    sigma_energy:    ' + str(self.sigma_energy_GeV) + ' [GeV]\n'
        p += '    current:         ' + str(self.current) + ' [A]\n'
        p += '    emittance:       ' + str(self.emittance) + '\n'
        return p
            
    def add_device(self, device):
        self.devices.append(device)
        return
    
    def add_devices(self, devices):
        self.devices += devices
        
    def get_brightness_curves (self, energy_range_eV=None, odir='data'):
        self.curves = []

        try:
            self.curves = pickle.load(open( os.path.join(odir, self.name + '_brightness.dat'), 'rb'))
            print('read from file', self.name, os.path.join(odir, self.name + '_brightness.dat'))
            return
        except:
            pass

        oth = oscars.th.th()
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

            if energy_range_eV is None:
                energy_range_eV = self.energy_range_eV

            if isinstance(d, BendingMagnet):
                bm = oth.dipole_brightness(
                    bfield=d.bfield,
                    energy_range_eV=energy_range_eV,
                    npoints=d.npoints
                )
                
                x = [b[0] for b in bm]
                y = [b[1] for b in bm]
                f = interp1d(x, y, kind='cubic')
                
                self.curves.append([d, energy_range_eV, f])


            elif isinstance(d, Wiggler):
                br = oth.dipole_brightness(
                    bfield=d.bfield,
                    energy_range_eV=energy_range_eV,
                    npoints=d.npoints
                )
                
                x = [b[0] for b in br]
                y = [b[1] * 2. * int(d.length / d.period) for b in br]
                f = interp1d(x, y, kind='cubic')
                
                self.curves.append([d, energy_range_eV, f])


            elif isinstance(d, Undulator):
                nperiods = int(d.length / d.period)
                for i in range(d.harmonic_range[0], d.harmonic_range[1]+1):
                    if i % 2 == 0:
                        continue


                    if d.calculation is 'SRW':
                        br = srw_und_brightness(
                            oth=oth,
                            Kx_range=d.k_range,
                            period=d.period,
                            length=d.length,
                            harmonic=i,
                            npoints=d.npoints)
                        x = [b[0]*1000. for b in br]
                        y = [b[1] for b in br]
                        x.reverse()
                        y.reverse()
                    else:
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
                    self.curves.append([d, [x[0], x[-1]], f])
            elif isinstance(d, EPU):
                if d.calculation is 'SRW':
                    br = srw_epu_brightness(
                        oth=oth,
                        Kx_range=d.kx_range,
                        Ky_range=d.ky_range,
                        period=d.period,
                        length=d.length,
                        npoints=d.npoints)
                    x = [b[0]*1000. for b in br]
                    y = [b[1] for b in br]
                    x.reverse()
                    y.reverse()
                else:
                    raise ValueError('Not implemented')

                f = interp1d(x, y, kind='cubic')
                self.curves.append([d, [x[0], x[-1]], f])



            else:
                raise ValueError(d.__class__.__name__ + ' is an invalid type.  please check')

        if odir is not None:
            pickle.dump(self.curves, open( os.path.join(odir, self.name + '_brightness.dat'), 'wb'))

        return



    def get_flux_onaxis_curves (self, energy_range_eV=None, odir='data'):
        self.curves = []

        try:
            self.curves = pickle.load(open( os.path.join(odir, self.name + '_fluxonaxis.dat'), 'rb'))
            print('read from file', self.name, os.path.join(odir, self.name + '_fluxonaxis.dat'))
            return
        except:
            pass

        oth = oscars.th.th()
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

            if energy_range_eV is None:
                energy_range_eV = self.energy_range_eV

            if isinstance(d, BendingMagnet):
                bm = oth.dipole_spectrum(
                    bfield=d.bfield,
                    energy_range_eV=energy_range_eV,
                    npoints=d.npoints
                )
                x = [b[0] for b in bm]
                y = [b[1] for b in bm]
                f = interp1d(x, y, kind='cubic')
                self.curves.append([d, energy_range_eV, f])


            elif isinstance(d, Wiggler):
                br = oth.dipole_spectrum(
                    bfield=d.bfield,
                    energy_range_eV=energy_range_eV,
                    npoints=d.npoints
                )
                x = [b[0] for b in br]
                y = [b[1] * 2. * int(d.length / d.period) for b in br]
                f = interp1d(x, y, kind='cubic')
                self.curves.append([d, energy_range_eV, f])


            elif isinstance(d, Undulator):
                nperiods = int(d.length / d.period)
                for i in range(d.harmonic_range[0], d.harmonic_range[1]+1):
                    if i % 2 == 0:
                        continue

                    #fl = oth.undulator_flux_onaxis(
                    #    K_range=d.k_range,
                    #    period=d.period,
                    #    nperiods=nperiods,
                    #    harmonic=i,
                    #    npoints=d.npoints,
                    #    minimum=d.minimum
                    #)
                    #x = [b[0] for b in fl]
                    #y = [b[1] for b in fl]

                    fl = srw_und_flux_onaxis(
                            oth,
                            d.k_range,
                            d.period,
                            d.length,
                            i,
                            d.npoints
                            )

                    x = [b[0]*1000 for b in fl]
                    y = [b[1] for b in fl]
                    x.reverse()
                    y.reverse()

                    f = interp1d(x, y, kind='cubic')
                    self.curves.append([d, [x[0], x[-1]], f])

            elif isinstance(d, EPU):
                fl = srw_epu_flux_onaxis(
                        oth,
                        d.kx_range,
                        d.ky_range,
                        d.period,
                        d.length,
                        d.npoints
                        )

                x = [b[0]*1000 for b in fl]
                y = [b[1] for b in fl]
                x.reverse()
                y.reverse()

                f = interp1d(x, y, kind='cubic')
                self.curves.append([d, [x[0], x[-1]], f])
                
            else:
                raise ValueError(d.__class__.__name__ + ' is an invalid type.  please check')

        if odir is not None:
            pickle.dump(self.curves, open( os.path.join(odir, self.name + '_fluxonaxis.dat'), 'wb'))


        return




    def get_flux_curves (self, energy_range_eV=None, odir='data'):
        self.curves = []

        try:
            self.curves = pickle.load(open( os.path.join(odir, self.name + '_flux.dat'), 'rb'))
            print('read from file', self.name, os.path.join(odir, self.name + '_flux.dat'))
            return
        except:
            pass

        oth = oscars.th.th()
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

            if energy_range_eV is None:
                energy_range_eV = self.energy_range_eV

            if isinstance(d, BendingMagnet):
                bm = oth.dipole_spectrum(
                    bfield=d.bfield,
                    energy_range_eV=energy_range_eV,
                    npoints=d.npoints,
                    angle_integrated=True
                )
                x = [b[0] for b in bm]
                y = [b[1] * 0.001 for b in bm]
                f = interp1d(x, y, kind='cubic')
                self.curves.append([d, [x[0], x[-1]], f])


            elif isinstance(d, Wiggler):
                fl = oth.dipole_spectrum(
                    bfield=d.bfield,
                    energy_range_eV=energy_range_eV,
                    npoints=d.npoints,
                )
                x = [b[0] for b in fl]
                y = [b[1] * 0.001 * 2. * int(d.length / d.period) for b in fl]
                f = interp1d(x, y, kind='cubic')
                self.curves.append([d, [x[0], x[-1]], f])


            elif isinstance(d, Undulator):
                nperiods = int(d.length / d.period)
                for i in range(d.harmonic_range[0], d.harmonic_range[1]+1):
                    if i % 2 == 0:
                        continue

                    #fl = oth.undulator_flux(
                    #    K_range=d.k_range,
                    #    period=d.period,
                    #    nperiods=nperiods,
                    #    harmonic=i,
                    #    npoints=d.npoints,
                    #    minimum=d.minimum
                    #)
                    #x = [b[0] for b in fl]
                    #y = [b[1] for b in fl]

                    fl = srw_und_flux(
                            oth,
                            d.k_range,
                            d.period,
                            d.length,
                            i,
                            d.npoints
                            )

                    x = [b[0]*1000 for b in fl]
                    y = [b[1] for b in fl]
                    x.reverse()
                    y.reverse()

                    f = interp1d(x, y, kind='cubic')
                    self.curves.append([d, [x[0], x[-1]], f])

            elif isinstance(d, EPU):
                fl = srw_epu_flux(
                        oth,
                        d.kx_range,
                        d.ky_range,
                        d.period,
                        d.length,
                        d.npoints
                        )

                x = [b[0]*1000 for b in fl]
                y = [b[1] for b in fl]
                x.reverse()
                y.reverse()

                f = interp1d(x, y, kind='cubic')
                self.curves.append([d, [x[0], x[-1]], f])
                
            else:
                raise ValueError(d.__class__.__name__ + ' is an invalid type.  please check')

        if odir is not None:
            pickle.dump(self.curves, open( os.path.join(odir, self.name + '_flux.dat'), 'wb'))

        return




    def get_coherentflux_curves (self, energy_range_eV=None, odir='data'):
        self.curves = []

        try:
            self.curves = pickle.load(open( os.path.join(odir, self.name + '_coherentflux.dat'), 'rb'))
            print('read from file', self.name, os.path.join(odir, self.name + '_coherentflux.dat'))
            return
        except:
            pass

        oth = oscars.th.th()
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

            if energy_range_eV is None:
                energy_range_eV = self.energy_range_eV

            if isinstance(d, BendingMagnet):
                pass
                #bm = oth.dipole_spectrum(
                #    bfield=d.bfield,
                #    energy_range_eV=energy_range_eV,
                #    npoints=d.npoints
                #)
                #
                #x = [b[0] for b in bm]
                #y = [b[1] for b in bm]
                #f = interp1d(x, y, kind='cubic')
                #
                #self.curves.append([d, energy_range_eV, f])


            elif isinstance(d, Wiggler):
                pass
                #br = oth.dipole_spectrum(
                #    bfield=d.bfield,
                #    energy_range_eV=energy_range_eV,
                #    npoints=d.npoints
                #)
                #
                #x = [b[0] for b in br]
                #y = [b[1] * 2. * int(d.length / d.period) for b in br]
                #f = interp1d(x, y, kind='cubic')
                #
                #self.curves.append([d, energy_range_eV, f])


            elif isinstance(d, Undulator):
                nperiods = int(d.length / d.period)
                for i in range(d.harmonic_range[0], d.harmonic_range[1]+1):
                    if i % 2 == 0:
                        continue

                    fl = srw_und_flux(
                            oth,
                            d.k_range,
                            d.period,
                            d.length,
                            i,
                            d.npoints
                            )

                    x = [b[0]*1000 for b in fl]
                    y = [b[1] for b in fl]
                    z = [b[2] for b in fl]
                    x.reverse()
                    y.reverse()
                    z.reverse()
                    cff = oth.undulator_coherentflux_fraction(
                        K_points=z,
                        period=d.period,
                        nperiods=nperiods,
                        harmonic=i,
                    )

                    fr = [cf[1] for cf in cff]
                    yy = [y[ik] * fr[ik] for ik in range(len(z))]
                    f = interp1d(x, yy, kind='cubic')
                    self.curves.append([d, [x[0], x[-1]], f])

            elif isinstance(d, EPU):
                nperiods = int(d.length / d.period)
                fl = srw_epu_flux(
                        oth,
                        d.kx_range,
                        d.ky_range,
                        d.period,
                        d.length,
                        d.npoints
                        )
                x = [b[0]*1000 for b in fl]
                y = [b[1] for b in fl]
                z = [b[2] for b in fl]
                x.reverse()
                y.reverse()
                z.reverse()
                cff = oth.undulator_coherentflux_fraction(
                    K_points=z,
                    period=d.period,
                    nperiods=nperiods,
                    harmonic=i,
                )

                fr = [cf[1] for cf in cff]
                yy = [y[ik] * fr[ik] for ik in range(len(z))]
                f = interp1d(x, yy, kind='cubic')
                self.curves.append([d, [x[0], x[-1]], f])

            else:
                raise ValueError(d.__class__.__name__ + ' is an invalid type.  please check')

        if odir is not None:
            pickle.dump(self.curves, open( os.path.join(odir, self.name + '_coherentflux.dat'), 'wb'))


        return









    def get_coherentflux_fraction_curves (self, energy_range_eV=None, odir='data'):
        self.curves = []

        try:
            self.curves = pickle.load(open( os.path.join(odir, self.name + '_coherentfluxfraction.dat'), 'rb'))
            print('read from file', self.name, os.path.join(odir, self.name + '_coherentfluxfraction.dat'))
            return
        except:
            pass

        oth = oscars.th.th()
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

            if energy_range_eV is None:
                energy_range_eV = self.energy_range_eV

            if isinstance(d, BendingMagnet):
                pass
                #bm = oth.dipole_spectrum(
                #    bfield=d.bfield,
                #    energy_range_eV=energy_range_eV,
                #    npoints=d.npoints
                #)
                #
                #x = [b[0] for b in bm]
                #y = [b[1] for b in bm]
                #f = interp1d(x, y, kind='cubic')
                #
                #self.curves.append([d, energy_range_eV, f])


            elif isinstance(d, Wiggler):
                pass
                #br = oth.dipole_spectrum(
                #    bfield=d.bfield,
                #    energy_range_eV=energy_range_eV,
                #    npoints=d.npoints
                #)
                #
                #x = [b[0] for b in br]
                #y = [b[1] * 2. * int(d.length / d.period) for b in br]
                #f = interp1d(x, y, kind='cubic')
                #
                #self.curves.append([d, energy_range_eV, f])


            elif isinstance(d, Undulator):
                nperiods = int(d.length / d.period)
                for i in range(d.harmonic_range[0], d.harmonic_range[1]+1):
                    if i % 2 == 0:
                        continue

                    br = oth.undulator_coherentflux_fraction(
                        K_range=d.k_range,
                        period=d.period,
                        nperiods=nperiods,
                        harmonic=i,
                        npoints=d.npoints,
                        minimum=0
                    )

                    x = [b[0] for b in br]
                    y = [b[1] for b in br]
                    f = interp1d(x, y, kind='cubic')
                    self.curves.append([d, [x[0], x[-1]], f])
            elif isinstance(d, EPU):
                pass
            else:
                raise ValueError(d.__class__.__name__ + ' is an invalid type.  please check')


        if odir is not None:
            pickle.dump(self.curves, open( os.path.join(odir, self.name + '_coherentfluxfraction.dat'), 'wb'))

        return











def plot_brightness (experiments, title='Brightness', energy_range_eV=None, figsize=None, xlim=None, ylim=None, ofile='', show=True, energy_eV=None, legendloc=''):

    plt.figure(figsize=figsize)

    cc = cycle(['C0','C1','C2','C3','C4','C5','C6','C7','C8','C9'])

    if energy_range_eV is None:
        energy_range_eV = [1e99, 1]
        for exp in experiments:
            if exp.energy_range_eV[0] < energy_range_eV[0]:
                energy_range_eV[0] = exp.energy_range_eV[0]
            if exp.energy_range_eV[1] > energy_range_eV[1]:
                energy_range_eV[1] = exp.energy_range_eV[1]


    for exp in experiments:

        max_br_at_eV = 0
        labeldone = False
        color = next(cc) if exp.color is None else exp.color
        print(exp.name)
        exp.get_brightness_curves(energy_range_eV)
        X = []
        Y = []


        for xx in np.logspace(np.log10(energy_range_eV[0]), np.log10(energy_range_eV[1]), 5000):
            yy = 0

            for curve in exp.curves:
                if xx >= curve[1][0] and xx <= curve[1][1] and curve[2](xx) > yy:
                    yy = curve[2](xx)
                    
                
            if yy is not 0:
                X.append(xx)
                Y.append(yy)
            elif yy is 0 and len(X) is not 0:
                if labeldone:
                    plt.plot(X, Y, color=color, linestyle=exp.linestyle, linewidth=exp.linewidth, marker=exp.marker)
                else:
                    plt.plot(X, Y, color=color, linestyle=exp.linestyle, linewidth=exp.linewidth, marker=exp.marker, label=exp.name)
                    labeldone = True
                X=[]
                Y=[]

        if energy_eV is not None:
            for curve in exp.curves:
                if energy_eV >= curve[1][0] and energy_eV <= curve[1][1] and curve[2](energy_eV) > max_br_at_eV:
                    max_br_at_eV = curve[2](energy_eV)
                    max_br_at_eV_name = curve[0].name
            print(f'Brightness max at {energy_eV} [eV] is {max_br_at_eV:1.2e} for {max_br_at_eV_name}')


        if len(X) > 0:
            if labeldone:
                plt.plot(X, Y, color=color, linestyle=exp.linestyle, linewidth=exp.linewidth, marker=exp.marker)
            else:
                plt.plot(X, Y, color=color, linestyle=exp.linestyle, linewidth=exp.linewidth, marker=exp.marker, label=exp.name)
                labeldone = True

    yl = plt.ylim()
    if xlim is not None:
        plt.xlim(xlim)
    if ylim is not None:
        plt.ylim(ylim)
    plt.grid()
    plt.loglog()
    plt.legend(loc=legendloc)
    plt.title(title)
    plt.xlabel('Photon Energy [eV]')
    plt.ylabel('Brightness [$photons/s/0.1\%bw/mm^2/mrad^2$]')

    if ofile != '':
        plt.savefig(ofile, bbox_inches='tight')

    if show:
        plt.show()

    return plt


def get_brightness_points (experiments, title='Brightness', energy_range_eV=None, figsize=None, xlim=None, ylim=None, ofile='', show=True, energy_eV=None, legendloc='', energy_points_eV=None):



    if energy_range_eV is None:
        energy_range_eV = [1e99, 1]
        for exp in experiments:
            if exp.energy_range_eV[0] < energy_range_eV[0]:
                energy_range_eV[0] = exp.energy_range_eV[0]
            if exp.energy_range_eV[1] > energy_range_eV[1]:
                energy_range_eV[1] = exp.energy_range_eV[1]


    points = []
    for exp in experiments:

        max_br_at_eV = 0
        print(exp.name)
        exp.get_brightness_curves(energy_range_eV)
        X = []
        Y = []


        xxpoints = energy_points_eV if energy_points_eV is not None else np.logspace(np.log10(energy_range_eV[0]), np.log10(energy_range_eV[1]), 5000)
        for xx in xxpoints:
            yy = 0

            for curve in exp.curves:
                if xx >= curve[1][0] and xx <= curve[1][1] and curve[2](xx) > yy:
                    yy = curve[2](xx)
                    
                
            X.append(xx)
            Y.append(yy)
            #elif yy is 0 and len(X) is not 0:
            #    X=[]
            #    Y=[]

        if energy_eV is not None:
            for curve in exp.curves:
                if energy_eV >= curve[1][0] and energy_eV <= curve[1][1] and curve[2](energy_eV) > max_br_at_eV:
                    max_br_at_eV = curve[2](energy_eV)
                    max_br_at_eV_name = curve[0].name
            print(f'Brightness max at {energy_eV} [eV] is {max_br_at_eV:1.2e} for {max_br_at_eV_name}')

        points.append([exp.name, X, Y])

    return points






def plot_flux (experiments, energy_range_eV=None, figsize=None, xlim=None, ylim=None, title='', ofile='', energy_eV=None):

    plt.figure(1, figsize=figsize)
    
    cc = cycle(['C0','C1','C2','C3','C4','C5','C6','C7','C8','C9'])

    if energy_range_eV is None:
        energy_range_eV = [1e99, 1]
        for exp in experiments:
            if exp.energy_range_eV[0] < energy_range_eV[0]:
                energy_range_eV[0] = exp.energy_range_eV[0]
            if exp.energy_range_eV[1] > energy_range_eV[1]:
                energy_range_eV[1] = exp.energy_range_eV[1]

    for exp in experiments:
        max_br_at_eV = 0
        labeldone = False
        print(exp.name)
        color = next(cc) if exp.color is None else exp.color
        exp.get_flux_curves(energy_range_eV)
        X = []
        Y = []
        
        for xx in np.logspace(np.log10(energy_range_eV[0]), np.log10(energy_range_eV[1]), 5000):
            yy = 0

            for curve in exp.curves:
                try:
                    if xx >= curve[1][0] and xx <= curve[1][1] and curve[2](xx) > yy:
                        yy = curve[2](xx)
                except:
                    print('ERROR in', exp.name, curve)
                
            if yy is not 0:
                X.append(xx)
                Y.append(yy)
            elif yy is 0 and len(X) is not 0:
                if labeldone:
                    plt.plot(X, Y, color=color, linestyle=exp.linestyle, linewidth=exp.linewidth)
                else:
                    plt.plot(X, Y, color=color, linestyle=exp.linestyle, linewidth=exp.linewidth, label=exp.name)
                    labeldone = True
                X=[]
                Y=[]
        if energy_eV is not None:
            for curve in exp.curves:
                if energy_eV >= curve[1][0] and energy_eV <= curve[1][1] and curve[2](energy_eV) > max_br_at_eV:
                    max_br_at_eV = curve[2](energy_eV)
                    max_br_at_eV_name = curve[0].name
            print(f'Flux max at {energy_eV} [eV] is {max_br_at_eV:1.2e} for {max_br_at_eV_name}')

        if len(X) > 0:
            if labeldone:
                plt.plot(X, Y, color=color, linestyle=exp.linestyle, linewidth=exp.linewidth)
            else:
                plt.plot(X, Y, color=color, linestyle=exp.linestyle, linewidth=exp.linewidth, label=exp.name)
                labeldone = True

    yl = plt.ylim()
    if xlim is not None:
        plt.xlim(xlim)
    if ylim is not None:
        plt.ylim(ylim)
    plt.grid()
    plt.loglog()
    plt.legend()
    plt.title(title)
    plt.xlabel('Photon Energy [eV]')
    plt.ylabel('Flux [$photons/s/0.1\%bw$]')

    if ofile != '':
        plt.savefig(ofile, bbox_inches='tight')

    plt.show()
    return plt



def get_flux_points (experiments, energy_range_eV=None, figsize=None, xlim=None, ylim=None, title='', ofile='', energy_eV=None, energy_points_eV=None):


    if energy_range_eV is None:
        energy_range_eV = [1e99, 1]
        for exp in experiments:
            if exp.energy_range_eV[0] < energy_range_eV[0]:
                energy_range_eV[0] = exp.energy_range_eV[0]
            if exp.energy_range_eV[1] > energy_range_eV[1]:
                energy_range_eV[1] = exp.energy_range_eV[1]

    points = []
    for exp in experiments:
        max_br_at_eV = 0
        print(exp.name)
        exp.get_flux_curves(energy_range_eV)
        X = []
        Y = []
        
        xxpoints = energy_points_eV if energy_points_eV is not None else np.logspace(np.log10(energy_range_eV[0]), np.log10(energy_range_eV[1]), 5000)
        for xx in xxpoints:
            yy = 0

            for curve in exp.curves:
                try:
                    if xx >= curve[1][0] and xx <= curve[1][1] and curve[2](xx) > yy:
                        yy = curve[2](xx)
                except:
                    print('ERROR in', exp.name, curve)
                
            X.append(xx)
            Y.append(yy)

        if energy_eV is not None:
            for curve in exp.curves:
                if energy_eV >= curve[1][0] and energy_eV <= curve[1][1] and curve[2](energy_eV) > max_br_at_eV:
                    max_br_at_eV = curve[2](energy_eV)
                    max_br_at_eV_name = curve[0].name
            print(f'Flux max at {energy_eV} [eV] is {max_br_at_eV:1.2e} for {max_br_at_eV_name}')
        
        points.append([exp.name, X, Y])
    return points





def plot_diff (experiment_sets, kind='flux_onaxis', energy_range_eV=None, figsize=None, xlim=None, ylim=None, title='', ofile='', energy_eV=None, ylabel=None):

    plt.figure(1, figsize=figsize)
    
    cc = cycle(['C0','C1','C2','C3','C4','C5','C6','C7','C8','C9'])

    for experiments in experiment_sets:
        A = experiments[0]
        B = experiments[1]
        color = next(cc) if A.color is None else A.color
        if energy_range_eV is None:
            energy_range_eV = [1e99, 1]
            for exp in experiments:
                if exp.energy_range_eV[0] < energy_range_eV[0]:
                    energy_range_eV[0] = exp.energy_range_eV[0]
                if exp.energy_range_eV[1] > energy_range_eV[1]:
                    energy_range_eV[1] = exp.energy_range_eV[1]
    
    
        X = np.logspace(np.log10(energy_range_eV[0]), np.log10(energy_range_eV[1]), 2000)
    
        for exp in experiments:
            if kind is 'flux_onaxis':
                exp.get_flux_onaxis_curves(energy_range_eV)
                if ylabel is None:
                    ylabel = 'Fractional On-Axis Flux Difference\n(A-B)/A'
            elif kind is 'flux':
                exp.get_flux_curves(energy_range_eV)
                if ylabel is None:
                    ylabel = 'Fractional Flux Difference\n(A-B)/A'
            elif kind is 'brightness':
                exp.get_brightness_curves(energy_range_eV)
                if ylabel is None:
                    ylabel = 'Fractional Brightness Difference\n(A-B)/A'
            elif kind is 'coherentflux':
                exp.get_coherentflux_curves(energy_range_eV)
                if ylabel is None:
                    ylabel = 'Fractional Coherent Flux Difference\n(A-B)/A'
            elif kind is 'coherentflux_fraction':
                exp.get_coherentflux_fraction_curves(energy_range_eV)
                if ylabel is None:
                    ylabel = 'Fractional Coherent Flux Fraction Difference\n(A-B)/A'
            else:
                raise ValueError('kind is incorrect', kind)
        # get points for A
        Y0 = []
        for xx in X:
            yy = 0
            for curve in A.curves:
                if xx >= curve[1][0] and xx <= curve[1][1] and curve[2](xx) > yy:
                    yy = curve[2](xx)
    
            Y0.append(yy)
    
        # get points for B
        Y1 = []
        for xx in X:
            yy = 0
            for curve in B.curves:
                if xx >= curve[1][0] and xx <= curve[1][1] and curve[2](xx) > yy:
                    yy = curve[2](xx)
    
            Y1.append(yy)
    
        # I am lazy right now so pop off the end
        X = X[0:-1]
        Y0 = Y0[0:-1]
        Y1 = Y1[0:-1]
    
        # check for zeros
        if 0 in Y0 or 0 in Y1:
            for i in range(len(X)):
                print(X[i], Y0[i], Y1[i])
            raise ValueError('zero found in Y0 or Y1')
    
        Y = [(a-b)/b for a, b in zip(Y0, Y1)]
    
    
        plt.plot(X, Y, color=color, linestyle=A.linestyle, linewidth=A.linewidth, label=A.name)

    yl = plt.ylim()
    if xlim is not None:
        plt.xlim(xlim)
    if ylim is not None:
        plt.ylim(ylim)
    plt.grid()
    plt.semilogx()
    #plt.loglog()
    plt.legend()
    plt.title(title)
    plt.xlabel('Photon Energy [eV]')
    plt.ylabel(ylabel)

    if ofile != '':
        plt.savefig(ofile, bbox_inches='tight')

    plt.show()
    return plt





def plot_flux_onaxis (experiments, energy_range_eV=None, figsize=None, xlim=None, ylim=None, title='', ofile='', energy_eV=None):

    plt.figure(1, figsize=figsize)
    
    cc = cycle(['C0','C1','C2','C3','C4','C5','C6','C7','C8','C9'])

    if energy_range_eV is None:
        energy_range_eV = [1e99, 1]
        for exp in experiments:
            if exp.energy_range_eV[0] < energy_range_eV[0]:
                energy_range_eV[0] = exp.energy_range_eV[0]
            if exp.energy_range_eV[1] > energy_range_eV[1]:
                energy_range_eV[1] = exp.energy_range_eV[1]

    for exp in experiments:
        max_br_at_eV = 0
        labeldone = False
        print(exp.name)
        color = next(cc) if exp.color is None else exp.color
        exp.get_flux_onaxis_curves(energy_range_eV)
        X = []
        Y = []
        
        for xx in np.logspace(np.log10(energy_range_eV[0]), np.log10(energy_range_eV[1]), 5000):
            yy = 0

            for curve in exp.curves:
                if xx >= curve[1][0] and xx <= curve[1][1] and curve[2](xx) > yy:
                    yy = curve[2](xx)
                    
                
            if yy is not 0:
                X.append(xx)
                Y.append(yy)
            elif yy is 0 and len(X) is not 0:
                if labeldone:
                    plt.plot(X, Y, color=color, linestyle=exp.linestyle, linewidth=exp.linewidth)
                else:
                    plt.plot(X, Y, color=color, linestyle=exp.linestyle, linewidth=exp.linewidth, label=exp.name)
                    labeldone = True
                X=[]
                Y=[]
        if energy_eV is not None:
            for curve in exp.curves:
                if energy_eV >= curve[1][0] and energy_eV <= curve[1][1] and curve[2](energy_eV) > max_br_at_eV:
                    max_br_at_eV = curve[2](energy_eV)
                    max_br_at_eV_name = curve[0].name
            print(f'Flux on axis max at {energy_eV} [eV] is {max_br_at_eV:1.2e} for {max_br_at_eV_name} ')

        if len(X) > 0:
            if labeldone:
                plt.plot(X, Y, color=color, linestyle=exp.linestyle, linewidth=exp.linewidth)
            else:
                plt.plot(X, Y, color=color, linestyle=exp.linestyle, linewidth=exp.linewidth, label=exp.name)
                labeldone = True

    yl = plt.ylim()
    if xlim is not None:
        plt.xlim(xlim)
    if ylim is not None:
        plt.ylim(ylim)
    plt.grid()
    plt.loglog()
    plt.legend()
    plt.title(title)
    plt.xlabel('Photon Energy [eV]')
    plt.ylabel('Flux [$photons/s/0.1\%bw/mrad^2$]')

    if ofile != '':
        plt.savefig(ofile, bbox_inches='tight')

    plt.show()
    return plt



def plot_coherentflux (experiments, energy_range_eV=None, figsize=None, xlim=None, ylim=None, title='', ofile='', energy_eV=None):

    plt.figure(1, figsize=figsize)
    
    cc = cycle(['C0','C1','C2','C3','C4','C5','C6','C7','C8','C9'])
    
    if energy_range_eV is None:
        energy_range_eV = [1e99, 1]
        for exp in experiments:
            if exp.energy_range_eV[0] < energy_range_eV[0]:
                energy_range_eV[0] = exp.energy_range_eV[0]
            if exp.energy_range_eV[1] > energy_range_eV[1]:
                energy_range_eV[1] = exp.energy_range_eV[1]

    for exp in experiments:
        max_br_at_eV = 0
        labeldone = False
        print(exp.name)
        color = next(cc) if exp.color is None else exp.color
        exp.get_coherentflux_curves(energy_range_eV)
        X = []
        Y = []
        
        for xx in np.logspace(np.log10(energy_range_eV[0]), np.log10(energy_range_eV[1]), 5000):
            yy = 0

            for curve in exp.curves:
                if xx >= curve[1][0] and xx <= curve[1][1] and curve[2](xx) > yy:
                    yy = curve[2](xx)
                    
                
            if yy is not 0:
                X.append(xx)
                Y.append(yy)
            elif yy is 0 and len(X) is not 0:
                if labeldone:
                    plt.plot(X, Y, color=color, linestyle=exp.linestyle, linewidth=exp.linewidth)
                else:
                    plt.plot(X, Y, color=color, linestyle=exp.linestyle, linewidth=exp.linewidth, label=exp.name)
                    labeldone = True
                X=[]
                Y=[]
        if energy_eV is not None:
            for curve in exp.curves:
                if energy_eV >= curve[1][0] and energy_eV <= curve[1][1] and curve[2](energy_eV) > max_br_at_eV:
                    max_br_at_eV = curve[2](energy_eV)
                    max_br_at_eV_name = curve[0].name
            print(f'Coherent Flux on axis max at {energy_eV} [eV] is {max_br_at_eV:1.2e} for {max_br_at_eV_name}')

        if len(X) > 0:
            if labeldone:
                plt.plot(X, Y, color=color, linestyle=exp.linestyle, linewidth=exp.linewidth)
            else:
                plt.plot(X, Y, color=color, linestyle=exp.linestyle, linewidth=exp.linewidth, label=exp.name)
                labeldone = True

    yl = plt.ylim()
    if xlim is not None:
        plt.xlim(xlim)
    if ylim is not None:
        plt.ylim(ylim)
    plt.grid()
    plt.loglog()
    plt.legend()
    plt.title(title)
    plt.xlabel('Photon Energy [eV]')
    plt.ylabel('Flux [$photons/s/0.1\%bw$]')

    if ofile != '':
        plt.savefig(ofile, bbox_inches='tight')

    plt.show()
    return plt











def plot_coherentflux_fraction (experiments, energy_range_eV=None, figsize=None, xlim=None, ylim=None, title='', ofile='', energy_eV=None):

    plt.figure(1, figsize=figsize)
    
    cc = cycle(['C0','C1','C2','C3','C4','C5','C6','C7','C8','C9'])

    if energy_range_eV is None:
        energy_range_eV = [1e99, 1]
        for exp in experiments:
            if exp.energy_range_eV[0] < energy_range_eV[0]:
                energy_range_eV[0] = exp.energy_range_eV[0]
            if exp.energy_range_eV[1] > energy_range_eV[1]:
                energy_range_eV[1] = exp.energy_range_eV[1]

    for exp in experiments:
        max_br_at_eV = 0
        labeldone = False
        print(exp.name)
        color = next(cc) if exp.color is None else exp.color
        exp.get_coherentflux_fraction_curves(energy_range_eV)
        X = []
        Y = []
        
        for xx in np.logspace(np.log10(energy_range_eV[0]), np.log10(energy_range_eV[1]), 5000):
            yy = 0

            for curve in exp.curves:
                if xx >= curve[1][0] and xx <= curve[1][1] and curve[2](xx) > yy:
                    yy = curve[2](xx)
                    
                
            if yy is not 0:
                X.append(xx)
                Y.append(yy)
            elif yy is 0 and len(X) is not 0:
                if labeldone:
                    plt.plot(X, Y, color=color, linestyle=exp.linestyle, linewidth=exp.linewidth)
                else:
                    plt.plot(X, Y, color=color, linestyle=exp.linestyle, linewidth=exp.linewidth, label=exp.name)
                    labeldone = True
                X=[]
                Y=[]
        if energy_eV is not None:
            for curve in exp.curves:
                if energy_eV >= curve[1][0] and energy_eV <= curve[1][1] and curve[2](energy_eV) > max_br_at_eV:
                    max_br_at_eV = curve[2](energy_eV)
                    max_br_at_eV_name = curve[0].name
            print(f'Coherent Flux Fraction on axis max at {energy_eV} [eV] is {max_br_at_eV} for {max_br_at_eV_name}')

        if len(X) > 0:
            if labeldone:
                plt.plot(X, Y, color=color, linestyle=exp.linestyle, linewidth=exp.linewidth)
            else:
                plt.plot(X, Y, color=color, linestyle=exp.linestyle, linewidth=exp.linewidth, label=exp.name)
                labeldone = True

    yl = plt.ylim()
    if xlim is not None:
        plt.xlim(xlim)
    if ylim is not None:
        plt.ylim(ylim)
    plt.grid()
    plt.loglog()
    plt.legend()
    plt.title(title)
    plt.xlabel('Photon Energy [eV]')
    plt.ylabel('Coherent Fraction')

    if ofile != '':
        plt.savefig(ofile, bbox_inches='tight')

    plt.show()
    return plt
















def plot_brightness_all (experiments, energy_range_eV, title='Brightness', figsize=None, xlim=None, ylim=None, ofile=''):

    plt.figure(1, figsize=figsize)
    
    cc = cycle(['C0','C1','C2','C3','C4','C5','C6','C7','C8','C9'])
    ll = cycle(['-', '--', ':'])

    for exp in experiments:
        print(exp.name)
        oth = oscars.th.th()
        for d in exp.devices:
            color = next(cc) if d.color is None else d.color
            linestyle = next(ll) if d.linestyle is None else d.linestyle
            oth.set_particle_beam(
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
                pass
                #print('BM')
                #fl = oth.dipole_spectrum(
                #    bfield=d.bfield,
                #    energy_range_eV=exp.energy_range_eV,
                #    npoints=d.npoints,
                #    angle_integrated=True
                #)
                # 
                #x = [b[0] for b in fl]
                #y = [b[1]*0.001 for b in fl]
                #
                #plt.plot(x, y, color=color, linestyle=d.linestyle, label=d.name)

            elif isinstance(d, Wiggler):
                pass
                #print('WG')
                #fl = oth.dipole_spectrum(
                #    bfield=d.bfield,
                #    energy_range_eV=exp.energy_range_eV,
                #    npoints=d.npoints,
                #    angle_integrated=True
                #)
                #
                #x = [b[0] for b in fl]
                #y = [b[1] *0.001 * 2. * int(d.length / d.period) for b in fl]
                #f = interp1d(x, y, kind='cubic')
                #
                #plt.plot(x, y, color=color, linestyle=d.linestyle, label=d.name)


            elif isinstance(d, Undulator):
                labeldone = False
                nperiods = int(d.length / d.period)
                curves = []
                for i in range(d.harmonic_range[0], d.harmonic_range[1]+1):
                    if i % 2 == 0:
                        continue

                    if d.calculation is 'SRW':
                        br = srw_und_brightness(
                            oth=oth,
                            Kx_range=d.k_range,
                            period=d.period,
                            length=d.length,
                            harmonic=i,
                            npoints=d.npoints)
                        x = [b[0]*1000. for b in br]
                        y = [b[1] for b in br]
                        x.reverse()
                        y.reverse()
                    else:
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
                            plt.plot(X, Y, color=color, linestyle=linestyle)
                        else:
                            plt.plot(X, Y, color=color, linestyle=linestyle, label=d.name)
                            labeldone = True
                        X=[]
                        Y=[]

                if len(X) > 0:
                    if labeldone:
                        plt.plot(X, Y, color=color, linestyle=linestyle)
                    else:
                        plt.plot(X, Y, color=color, linestyle=linestyle, label=d.name)
                        labeldone = True
            elif isinstance(d, EPU):
                labeldone = False
                curves = []

                if d.calculation is 'SRW':
                    br = srw_epu_brightness(
                        oth=oth,
                        Kx_range=d.kx_range,
                        Ky_range=d.ky_range,
                        period=d.period,
                        length=d.length,
                        npoints=d.npoints)
                    x = [b[0]*1000. for b in br]
                    y = [b[1] for b in br]
                    x.reverse()
                    y.reverse()
                else:
                    raise ValueError('Not implemented')

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
                            plt.plot(X, Y, color=color, linestyle=linestyle)
                        else:
                            plt.plot(X, Y, color=color, linestyle=linestyle, label=d.name)
                            labeldone = True
                        X=[]
                        Y=[]

                if len(X) > 0:
                    if labeldone:
                        plt.plot(X, Y, color=color, linestyle=linestyle)
                    else:
                        plt.plot(X, Y, color=color, linestyle=linestyle, label=d.name)
                        labeldone = True
            else:
                raise ValueError(d.__class__.__name__ + ' is an invalid type.  please check')

    yl = plt.ylim()
    if xlim is not None:
        plt.xlim(xlim)
    if ylim is not None:
        plt.ylim(ylim)
    plt.grid()
    plt.loglog()
    plt.legend()
    plt.title(title)
    plt.xlabel('Photon Energy [eV]')
    plt.ylabel('Brightness [$photons/s/0.1\%bw/mm^2/mrad^2$]')

    if ofile != '':
        plt.savefig(ofile, bbox_inches='tight')

    plt.show()
    return plt

























def plot_flux_all (experiments, energy_range_eV, title='Spectral Flux', figsize=None, xlim=None, ylim=None, ofile=''):

    plt.figure(1, figsize=figsize)
    
    cc = cycle(['C0','C1','C2','C3','C4','C5','C6','C7','C8','C9'])
    ll = cycle(['-', '--', ':'])

    for exp in experiments:
        print(exp.name)
        oth = oscars.th.th()
        for d in exp.devices:
            color = next(cc) if d.color is None else d.color
            linestyle = next(ll) if d.linestyle is None else d.linestyle
            oth.set_particle_beam(
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
                fl = oth.dipole_spectrum(
                    bfield=d.bfield,
                    energy_range_eV=exp.energy_range_eV,
                    npoints=d.npoints,
                    angle_integrated=True
                )
                
                x = [b[0] for b in fl]
                y = [b[1]*0.001 for b in fl]

                plt.plot(x, y, color=color, linestyle=linestyle, label=d.name)

            elif isinstance(d, Wiggler):
                fl = oth.dipole_spectrum(
                    bfield=d.bfield,
                    energy_range_eV=exp.energy_range_eV,
                    npoints=d.npoints,
                    angle_integrated=True
                )
                
                x = [b[0] for b in fl]
                y = [b[1] *0.001 * 2. * int(d.length / d.period) for b in fl]
                f = interp1d(x, y, kind='cubic')
                
                plt.plot(x, y, color=color, linestyle=linestyle, label=d.name)


            elif isinstance(d, Undulator):
                labeldone = False
                nperiods = int(d.length / d.period)
                curves = []
                for i in range(d.harmonic_range[0], d.harmonic_range[1]+1):
                    if i % 2 == 0:
                        continue

                    #fl = oth.undulator_flux(
                    #    K_range=d.k_range,
                    #    period=d.period,
                    #    nperiods=nperiods,
                    #    harmonic=i,
                    #    npoints=d.npoints,
                    #    minimum=d.minimum
                    #)
                    fl = srw_und_flux(
                            oth,
                            d.k_range,
                            d.period,
                            d.length,
                            i,
                            d.npoints
                            )

                    x = [b[0]*1000 for b in fl]
                    y = [b[1] for b in fl]
                    x.reverse()
                    y.reverse()
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
                            plt.plot(X, Y, color=color, linestyle=linestyle)
                        else:
                            plt.plot(X, Y, color=color, linestyle=linestyle, label=d.name)
                            labeldone = True
                        X=[]
                        Y=[]

                if len(X) > 0:
                    if labeldone:
                        plt.plot(X, Y, color=color, linestyle=linestyle)
                    else:
                        plt.plot(X, Y, color=color, linestyle=linestyle, label=d.name)
                        labeldone = True
            elif isinstance(d, EPU):
                labeldone = False
                nperiods = int(d.length / d.period)
                curves = []
                fl = srw_epu_flux(
                        oth,
                        d.kx_range,
                        d.ky_range,
                        d.period,
                        d.length,
                        d.npoints
                        )

                x = [b[0]*1000 for b in fl]
                y = [b[1] for b in fl]
                x.reverse()
                y.reverse()
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
                            plt.plot(X, Y, color=color, linestyle=linestyle)
                        else:
                            plt.plot(X, Y, color=color, linestyle=linestyle, label=d.name)
                            labeldone = True
                        X=[]
                        Y=[]

                if len(X) > 0:
                    if labeldone:
                        plt.plot(X, Y, color=color, linestyle=linestyle)
                    else:
                        plt.plot(X, Y, color=color, linestyle=linestyle, label=d.name)
                        labeldone = True
            else:
                raise ValueError(d.__class__.__name__ + ' is an invalid type.  please check')

    yl = plt.ylim()
    if xlim is not None:
        plt.xlim(xlim)
    if ylim is not None:
        plt.ylim(ylim)
    plt.grid()
    plt.loglog()
    plt.legend()
    plt.title(title)
    plt.xlabel('Photon Energy [eV]')
    plt.ylabel('Flux [$photons/s/0.1\%bw$]')

    if ofile != '':
        plt.savefig(ofile, bbox_inches='tight')

    plt.show()
    return plt












































def plot_flux_onaxis_all (experiments, energy_range_eV, title='Peak Flux / mrad$^2$', figsize=None, xlim=None, ylim=None, ofile=''):

    plt.figure(1, figsize=figsize)
    
    cc = cycle(['C0','C1','C2','C3','C4','C5','C6','C7','C8','C9'])
    ll = cycle(['-', '--', ':'])

    for exp in experiments:
        print(exp.name)
        oth = oscars.th.th()
        for d in exp.devices:
            color = next(cc) if d.color is None else d.color
            linestyle = next(ll) if d.linestyle is None else d.linestyle
            oth.set_particle_beam(
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
                fl = oth.dipole_spectrum(
                    bfield=d.bfield,
                    energy_range_eV=exp.energy_range_eV,
                    npoints=d.npoints,
                )
                
                x = [b[0] for b in fl]
                y = [b[1] for b in fl]

                plt.plot(x, y, color=color, linestyle=linestyle, label=d.name)

            elif isinstance(d, Wiggler):
                fl = oth.dipole_spectrum(
                    bfield=d.bfield,
                    energy_range_eV=exp.energy_range_eV,
                    npoints=d.npoints,
                )
                
                x = [b[0] for b in fl]
                y = [b[1] * 2. * int(d.length / d.period) for b in fl]
                f = interp1d(x, y, kind='cubic')
                
                plt.plot(x, y, color=color, linestyle=linestyle, label=d.name)


            elif isinstance(d, Undulator):
                labeldone = False
                nperiods = int(d.length / d.period)
                curves = []
                for i in range(d.harmonic_range[0], d.harmonic_range[1]+1):
                    if i % 2 == 0:
                        continue

                    #fl = oth.undulator_flux_onaxis(
                    #    K_range=d.k_range,
                    #    period=d.period,
                    #    nperiods=nperiods,
                    #    harmonic=i,
                    #    npoints=d.npoints,
                    #    minimum=d.minimum
                    #)
                    fl = srw_und_flux_onaxis(
                            oth,
                            d.k_range,
                            d.period,
                            d.length,
                            i,
                            d.npoints
                            )

                    x = [b[0]*1000 for b in fl]
                    y = [b[1] for b in fl]
                    x.reverse()
                    y.reverse()
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
                            plt.plot(X, Y, color=color, linestyle=linestyle)
                        else:
                            plt.plot(X, Y, color=color, linestyle=linestyle, label=d.name)
                            labeldone = True
                        X=[]
                        Y=[]

                if len(X) > 0:
                    if labeldone:
                        plt.plot(X, Y, color=color, linestyle=linestyle)
                    else:
                        plt.plot(X, Y, color=color, linestyle=linestyle, label=d.name)
                        labeldone = True
            elif isinstance(d, EPU):
                labeldone = False
                nperiods = int(d.length / d.period)
                curves = []
                fl = srw_epu_flux_onaxis(
                        oth,
                        d.kx_range,
                        d.ky_range,
                        d.period,
                        d.length,
                        d.npoints
                        )

                x = [b[0]*1000 for b in fl]
                y = [b[1] for b in fl]
                x.reverse()
                y.reverse()
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
                            plt.plot(X, Y, color=color, linestyle=linestyle)
                        else:
                            plt.plot(X, Y, color=color, linestyle=linestyle, label=d.name)
                            labeldone = True
                        X=[]
                        Y=[]

                if len(X) > 0:
                    if labeldone:
                        plt.plot(X, Y, color=color, linestyle=linestyle)
                    else:
                        plt.plot(X, Y, color=color, linestyle=linestyle, label=d.name)
                        labeldone = True
            else:
                raise ValueError(d.__class__.__name__ + ' is an invalid type.  please check')

    yl = plt.ylim()
    if xlim is not None:
        plt.xlim(xlim)
    if ylim is not None:
        plt.ylim(ylim)
    plt.grid()
    plt.loglog()
    plt.legend()
    plt.title(title)
    plt.xlabel('Photon Energy [eV]')
    plt.ylabel('Flux [$photons/s/0.1\%bw/mrad^2$]')

    if ofile != '':
        plt.savefig(ofile, bbox_inches='tight')

    plt.show()
    return plt


