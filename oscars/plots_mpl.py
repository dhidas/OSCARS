from __future__ import print_function

from matplotlib.colors import LogNorm
import matplotlib.ticker
import numpy as np
import matplotlib.pyplot as plt
from math import sqrt


def plot_trajectory_position(trajectory, show=True, ofile='', axis='Z', figsize=[18, 4.5], ret=False):
    """Plot the trajectory position.  You can optionally change the axis, output to a file, show or not, and return or not
       
       :param trajectory: Particle trajectory
       :type  trajectory: list
       :param show: to show the plot or not
       :type  show: bool
       :param ofile: output file name
       :type  ofile: str
       :param axis: which axis to plot trajectory on
       :type  axis: str
       :param figsize: dimension of the figure
       :type  figsize: list
       :param ret: to return the plot or not
       :type  ret: bool
       """

    # Get coordinate lists
    X  = [item[1][0] for item in trajectory]
    Y  = [item[1][1] for item in trajectory]
    Z  = [item[1][2] for item in trajectory]

    if axis is 'X':
        X1Label = 'X [m]'
        X2Label = 'Y [m]'
        X3Label = 'Z [m]'
        X1 = X
        X2 = Y
        X3 = Z
    elif axis is 'Y':
        X1Label = 'Y [m]'
        X2Label = 'Z [m]'
        X3Label = 'X [m]'
        X1 = Y
        X2 = Z
        X3 = X
    elif axis is 'Z':
        X1Label = 'Z [m]'
        X2Label = 'X [m]'
        X3Label = 'Y [m]'
        X1 = Z
        X2 = X
        X3 = Y

    # Plot X and Y vs. Z
    plt.figure(1, figsize=figsize)
    plt.subplot(131)
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.plot(X1, X2)
    plt.xlabel(X1Label)
    plt.ylabel(X2Label)
    plt.title('Particle Trajectory')

    plt.subplot(132)
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.plot(X1, X3)
    plt.xlabel(X1Label)
    plt.ylabel(X3Label)
    plt.title('Particle Trajectory')

    plt.subplot(133)
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.plot(X2, X3)
    plt.xlabel(X2Label)
    plt.ylabel(X3Label)
    plt.title('Particle Trajectory')


    if ofile != '':
        plt.savefig(ofile, bbox_inches='tight')

    if show == True:
        plt.show()

    if ret is True:
        return plt
    return


def plot_trajectory_velocity(trajectory, show=True, ofile='', figsize=[18, 4.5], ret=False):
    """Plot the trajectory velocity.  You can optionally change the axis, output to a file, show or not, and return or not
       
       :param trajectory: Particle trajectory
       :type  trajectory: list
       :param show: to show the plot or not
       :type  show: bool
       :param ofile: output file name
       :type  ofile: str
       :param axis: which axis to plot trajectory on
       :type  axis: str
       :param figsize: dimension of the figure
       :type  figsize: list
       :param ret: to return the plot or not
       :type  ret: bool
       """


    # Get coordinate lists
    T  = [item[0]    for item in trajectory]
    VX = [item[2][0] for item in trajectory]
    VY = [item[2][1] for item in trajectory]
    VZ = [item[2][2] for item in trajectory]

    # Plot VX, VY, VZ vs. T
    plt.figure(1, figsize=figsize)
    plt.subplot(131)
    plt.plot(T, VX)
    plt.xlabel('T [s]')
    plt.ylabel('$\\beta_x$')
    plt.title('Particle $\\beta_x$')

    plt.subplot(132)
    plt.plot(T, VY)
    plt.xlabel('T [s]')
    plt.ylabel('$\\beta_y$')
    plt.title('Particle $\\beta_y$')

    plt.subplot(133)
    plt.plot(T, VZ)
    plt.xlabel('T [s]')
    plt.ylabel('$\\beta_z$')
    plt.title('Particle $\\beta_z$')

    if ofile != '':
        plt.savefig(ofile, bbox_inches='tight')

    if show == True:
        plt.show()

    if ret:
        return plt
    return
    
    
def plot_trajectory_betaprime(trajectory, show=True, ofile='', figsize=[18, 4.5], ret=False):
    """Plot the trajectory beta^prime.  You can optionally change the axis, output to a file, show or not, and return or not
       
       :param trajectory: Particle trajectory
       :type  trajectory: list
       :param show: to show the plot or not
       :type  show: bool
       :param ofile: output file name
       :type  ofile: str
       :param axis: which axis to plot trajectory on
       :type  axis: str
       :param figsize: dimension of the figure
       :type  figsize: list
       :param ret: to return the plot or not
       :type  ret: bool
       """


    # Get coordinate lists
    T  = [item[0]    for item in trajectory]
    VX = [item[3][0] for item in trajectory]
    VY = [item[3][1] for item in trajectory]
    VZ = [item[3][2] for item in trajectory]

    # Plot VX, VY, VZ vs. T
    plt.figure(1, figsize=figsize)
    plt.subplot(131)
    plt.plot(T, VX)
    plt.xlabel('T [s]')
    plt.ylabel('$\\beta_x$')
    plt.title('Particle $\\beta^\\prime_x$')

    plt.subplot(132)
    plt.plot(T, VY)
    plt.xlabel('T [s]')
    plt.ylabel('$\\beta_y$')
    plt.title('Particle $\\beta^\\prime_y$')

    plt.subplot(133)
    plt.plot(T, VZ)
    plt.xlabel('T [s]')
    plt.ylabel('$\\beta_z$')
    plt.title('Particle $\\beta^\\prime_z$')

    if ofile != '':
        plt.savefig(ofile, bbox_inches='tight')

    if show == True:
        plt.show()

    if ret:
        return plt
    return
    
    
def plot_power_density(V, title=None, xlabel='X1 Axis [$m$]', ylabel='X2 Axis [$m$]', show=True, ofile='', figsize=None, ret=False, x1=None, x2=None):
    """Plot a 2D histogram with equal spacing"""
     
    if x1 is not None or x2 is not None:
        return plot_power_density_2d1d(V, x1=x1, x2=x2, title=title, show=show, ofile=ofile, figsize=figsize, ret=ret)

    X = [item[0][0] for item in V]
    Y = [item[0][1] for item in V]
    P = [item[1]    for item in V]

    NX = len(np.unique(X))
    NY = len(np.unique(Y))

    plt.figure(1, figsize=figsize)
    plt.hist2d(X, Y, bins=[NX, NY], weights=P)
    plt.colorbar(format='%.0e')
    #cb.formatter.set_scientific(True)
    #cb.formatter.set_powerlimits((0, 0))
    #cb.ax.yaxis.set_offset_position('right')
    #cb.update_ticks()
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

    if title is None:
        title = 'Power Density [$W / mm^2$]'

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)

    if ofile != '':
        plt.savefig(ofile, bbox_inches='tight')

    if show == True:
        plt.show()

    if ret:
        return plt

    return


def write_power_density_csv(V, ofile=''):
    """write histogram with equal spacing"""
 

    if ofile is '':
        raise ValueError('Must specify an output file')

    X = [item[0][0] for item in V]
    Y = [item[0][1] for item in V]
    P = [item[1]    for item in V]

    NX = len(np.unique(X))
    NY = len(np.unique(Y))

    with open(ofile, 'w') as f:
        f.write('            , ')
        for ix in range(NX):
            f.write('{: 11.5e}'.format(X[ix * NY]))
            if ix < NX:
                f.write(', ')
        f.write('\n')

        for j in range(NY):
            f.write('{: 11.5e}'.format(Y[j]))
            f.write(', ')
            for i in range(NX):
                ip = i * NY + j
 
                f.write('{: 11.5e}'.format(P[ip]))
                if i < NX:
                    f.write(', ')
            f.write('\n')


    return


def plot_power_density_2d1d(V, x1=None, x2=None, title=None, xlabel='[$m$]', ylabel='[$W / mm^2$]', xlim=None, ylim=None, show=True, ofile='', figsize=None, ret=False):
    """Plot a 2D histogram with equal spacing"""

    if x1 is not None and x2 is not None:
        raise ValueError('x1 and x2 cannot both be defined')
    if x1 is None and x2 is None:
        raise ValueError('x1 or x2 must be defined')

    X = [item[0][0] for item in V]
    Y = [item[0][1] for item in V]
    P = [item[1]    for item in V]

    XValues = np.unique(X).tolist()
    YValues = np.unique(Y).tolist()

    NX = len(XValues)
    NY = len(YValues)

    x1_index0 = -1
    x1_index1 = -1
    x1_frac0 = 0.5
    x2_index0 = -1
    x2_index1 = -1
    x2_frac0 = 0.5

    XC = []
    YP = []

    title_position = '' 

    if x1 is not None:
        if xlabel is None:
            xlabel = 'x1 [m]'
        title_position = 'x1 = ' + str(round(x1, 4))
        if x1 in XValues:
            x1_index0 = XValues.index(x1)
            x1_index1 = x1_index0
        else:
            x1_index1 = next(x[0] for x in enumerate(XValues) if x[1] > x1)
            x1_index0 = x1_index1 - 1
            if x1_index0 < 0:
                raise ValueError('x1 is not in range')
        for iy in range(NY):
            v0 = P[NY*x1_index0 + iy] * x1_frac0
            v1 = P[NY*x1_index1 + iy] * (1 - x1_frac0)
            YP.append(Y[NY*x1_index0 + iy])
            XC.append(v0 + v1)

            
    elif x2 is not None:
        if xlabel is None:
            xlabel = 'x2 [m]'
        title_position = 'x2 = ' + str(round(x2, 4))
        if x2 in YValues:
            x2_index0 = YValues.index(x2)
            x2_index1 = x2_index0
        else:
            x2_index1 = next(x[0] for x in enumerate(YValues) if x[1] > x2)
            x2_index0 = x2_index1 - 1
            if x2_index0 < 0:
                raise ValueError('x2 is not in range')
        for ix in range(NX):
            v0 = P[NY*ix + x2_index0] * x2_frac0
            v1 = P[NY*ix + x2_index1] * (1 - x2_frac0)
            YP.append(X[NY*ix + x2_index0])
            XC.append(v0 + v1)





    plt.figure(1, figsize=figsize)
    if ylim is not None: plt.ylim(ylim[0], ylim[1])
    if xlim is not None: plt.xlim(xlim[0], xlim[1])
    plt.plot(YP, XC)

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    if title is None:
        title = 'Power Density [$W / mm^2$]'
        plt.title(title + ' at ' + title_position)
    else:
        plt.title(title)

    if ofile != '':
        plt.savefig(ofile, bbox_inches='tight')

    if show == True:
        plt.show()

    if ret:
        return plt

    return


def plot_power_density_1d(V, title='Power Density [$W / mm^2$]', xlabel='[$m$]', ylabel='[$W / mm^2$]', xlim=None, ylim=None, show=True, ofile='', figsize=None, ret=False):
    """Plot a 1D power density"""

    X = [item[0][0] for item in V]
    P = [item[1]    for item in V]

    plt.figure(1, figsize=figsize)
    plt.plot(X, P)
    if ylim is not None: plt.ylim(ylim[0], ylim[1])
    if xlim is not None: plt.xlim(xlim[0], xlim[1])

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)

    if ofile != '':
        plt.savefig(ofile, bbox_inches='tight')

    if show == True:
        plt.show()

    if ret:
        return plt

    return


def plot_flux(V, title='Flux [$\gamma / mm^2 / 0.1\%bw / s$]', xlabel='X1 Axis [$m$]', ylabel='X2 Axis [$m$]', clim=None, show=True, ofile='', figsize=None, ylim=None, xlim=None, colorbar=True, ret=False, nticks_cb=None):
    """Plot a 2D histogram with equal spacing"""
        
    X = [item[0][0] for item in V]
    Y = [item[0][1] for item in V]
    P = [item[1]    for item in V]

    NX = len(np.unique(X))
    NY = len(np.unique(Y))

    # Size and limits
    plt.figure(1, figsize=figsize)
    if ylim is not None: plt.ylim(ylim[0], ylim[1])
    if xlim is not None: plt.xlim(xlim[0], xlim[1])

    plt.hist2d(X, Y, bins=[NX, NY], weights=P)
    #cb.formatter.set_scientific(True)
    #cb.formatter.set_powerlimits((0, 0))
    #cb.ax.yaxis.set_offset_position('right')
    #cb.update_ticks()
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    if colorbar is True:
        cb = plt.colorbar(format='%.2e')
        if nticks_cb is not None:
            tick_locator = matplotlib.ticker.MaxNLocator(nbins=nticks_cb)
            cb.locator = tick_locator
            cb.update_ticks()

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)


    if clim is not None:
      plt.clim(clim)

    if ofile != '':
        plt.savefig(ofile, bbox_inches='tight')

    if show == True:
        plt.show()

    if ret:
        return plt

    plt.close()
    return


def plot_spectrum(S, log=False, show=True, ofile='', title='Spectrum', xlabel='Energy [eV]', ylabel='$\gamma / mm^2 / 0.1\%bw / s$', figsize=None, ylim=None, xlim=None, transparent=True, ret=False, xticks=None, axhlines=None, axvlines=None, **kwargs):
    """Plot the spectrum"""

    # Size and limits
    plt.figure(1, figsize=figsize)
    if ylim is not None: plt.ylim(ylim[0], ylim[1])
    if xlim is not None: plt.xlim(xlim[0], xlim[1])

    X = [item[0] for item in S]
    Y = [item[1] for item in S]
    plt.plot(X, Y, **kwargs)
    if log:
        plt.yscale('log')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)

    if xticks is not None:
        plt.xticks(xticks)

    if axvlines is not None:
        for line in axvlines:
            plt.axvline(x=line, color='red', linestyle='--')
    if axhlines is not None:
        for line in axhlines:
            plt.axhline(y=line, color='green', linestyle='--')

    if ofile != '':
        plt.savefig(ofile, bbox_inches='tight', transparent=transparent)

    if show == True:
        plt.show()

    if ret is True:
        return plt

    return



def plot_spectra(spectra, label=None, legend=None, colors=None, show=True, ofile='', title='', loc=None, log=False, loglog=False, xlabel='Energy [eV]', ylabel='[$\gamma / mm^2 / 0.1\%bw / s$]', figsize=None, ylim=None, xlim=None, ret=False, axis=None, transparent=True, xticks=None, xvlines=None, **kwargs):


    # Size and limits
    plt.figure(1, figsize=figsize)
    if loglog:
        plt.loglog()
    if ylim is not None: plt.ylim(ylim[0], ylim[1])
    if xlim is not None: plt.xlim(xlim[0], xlim[1])

    for i in range(len(spectra)):
        s = spectra[i]
        color = colors[i] if colors is not None else None
        X = [item[0] for item in s]
        Y = [item[1] for item in s]
        if label is not None:
            if label[i] is not None:
                plt.plot(X, Y, label=label[i], c=color, **kwargs)
            else:
                plt.plot(X, Y, c=color, **kwargs)
        else:
            plt.plot(X, Y, c=color, **kwargs)

    if xvlines is not None:
        for xvline in xvlines:
            plt.axvline(x=xvline, color='red', linestyle='dashed', linewidth=1)

    if log:
        plt.yscale('log')

    plt.legend(title=legend, loc=loc)
    if xlabel is not None:
        plt.xlabel(xlabel)
    if ylabel is not None:
        plt.ylabel(ylabel)
    if title == None:
        title='Spectrum'
    plt.title(title)
    
    if axis is not None:
        plt.axis(axis)

    if xticks is not None:
        plt.xticks(xticks)

    if ofile != '':
        plt.savefig(ofile, bbox_inches='tight', transparent=transparent)

    if show == True:
        plt.show()

    if ret:
        return plt
    return





def plot_bfield(osr, mymin=-1, mymax=1, t=0, name='', ylim=None, show=True, ofile='', axis='Z', npoints=20000, between_two_points=None, ret=False):
    """Plot the magnetic field as a function of Z"""


    P = []
    Bx = []
    By = []
    Bz = []


    if between_two_points is not None:
        p0 = between_two_points[0]
        p1 = between_two_points[1]
        step = [ (p1[0] - p0[0]) / float(npoints - 1), (p1[1] - p0[1]) / float(npoints - 1), (p1[2] - p0[2]) / float(npoints - 1)]
        distance = sqrt( pow(p1[0]-p0[0], 2) + pow(p1[1]-p0[1], 2) + pow(p1[2]-p0[2], 2) )
        P = np.linspace(0, distance, npoints)


        for i in range(npoints):
            x = p0[0] + step[0] * float(i)
            y = p0[1] + step[1] * float(i)
            z = p0[2] + step[2] * float(i)

            B = osr.get_bfield([x, y, z], t, name)

            Bx.append(B[0])
            By.append(B[1])
            Bz.append(B[2])
            axis = 'Position'
    else:
        P = np.linspace(mymin, mymax, npoints)
        B = []
        for p in P:
            if axis is 'X':
                B = osr.get_bfield([p, 0, 0], t, name)
            elif axis is 'Y':
                B = osr.get_bfield([0, p, 0], t, name)
            elif axis is 'Z':
                B = osr.get_bfield([0, 0, p], t, name)
            else:
                raise
            Bx.append(B[0])
            By.append(B[1])
            Bz.append(B[2])

    plt.figure(1, figsize=(18, 4.5))
    plt.subplot(131)
    plt.plot(P, Bx)
    plt.xlabel(axis + ' [m]')
    plt.ylabel('Bx [T]')
    plt.ylim(ylim)

    plt.subplot(132)
    plt.plot(P, By)
    plt.xlabel(axis + ' [m]')
    plt.ylabel('By [T]')
    plt.ylim(ylim)

    plt.subplot(133)
    plt.plot(P, Bz)
    plt.xlabel(axis + ' [m]')
    plt.ylabel('Bz [T]')
    plt.ylim(ylim)

    if ofile != '':
        plt.savefig(ofile, bbox_inches='tight')

    if show == True:
        plt.show()

    if ret:
        return plt
    return







def plot_efield(osr, mymin=-1, mymax=1, t=0, name='', ylim=None, show=True, ofile='', axis='Z', npoints=20000, between_two_points=None, ret=False):
    """Plot the electric field as a function of Z"""


    P = []
    Bx = []
    By = []
    Bz = []


    if between_two_points is not None:
        p0 = between_two_points[0]
        p1 = between_two_points[1]
        step = [ (p1[0] - p0[0]) / float(npoints - 1), (p1[1] - p0[1]) / float(npoints - 1), (p1[2] - p0[2]) / float(npoints - 1)]
        distance = sqrt( pow(p1[0]-p0[0], 2) + pow(p1[1]-p0[1], 2) + pow(p1[2]-p0[2], 2) )
        P = np.linspace(0, distance, npoints)


        for i in range(npoints):
            x = p0[0] + step[0] * float(i)
            y = p0[1] + step[1] * float(i)
            z = p0[2] + step[2] * float(i)

            B = osr.get_bfield([x, y, z], t, name)

            Bx.append(B[0])
            By.append(B[1])
            Bz.append(B[2])
            axis = 'Position'
    else:
        P = np.linspace(mymin, mymax, npoints)
        B = []
        for p in P:
            if axis is 'X':
                B = osr.get_bfield([p, 0, 0], t, name)
            elif axis is 'Y':
                B = osr.get_bfield([0, p, 0], t, name)
            elif axis is 'Z':
                B = osr.get_bfield([0, 0, p], t, name)
            else:
                raise
            Bx.append(B[0])
            By.append(B[1])
            Bz.append(B[2])

    plt.figure(1, figsize=(18, 4.5))
    plt.subplot(131)
    plt.plot(P, Bx)
    plt.xlabel(axis + ' [m]')
    plt.ylabel('Ex [V/m]')
    plt.ylim(ylim)

    plt.subplot(132)
    plt.plot(P, By)
    plt.xlabel(axis + ' [m]')
    plt.ylabel('Ey [V/m]')
    plt.ylim(ylim)

    plt.subplot(133)
    plt.plot(P, Bz)
    plt.xlabel(axis + ' [m]')
    plt.ylabel('Ez [V/m]')
    plt.ylim(ylim)

    if ofile != '':
        plt.savefig(ofile, bbox_inches='tight')

    if show == True:
        plt.show()

    if ret:
        return plt
    return





def total_power(pd):
    """Calculate the total power in a uniform grid.
    
    This will not work for a non-uniform grid.  Different NX and NY are ok."""
    
    X = [item[0][0] for item in pd]
    Y = [item[0][1] for item in pd]
    P = [item[1]    for item in pd]

    NX = len(np.unique(X))
    NY = len(np.unique(Y))
    
    dx = (max(X) - min(X)) / (NX - 1)
    dy = (max(Y) - min(Y)) / (NY - 1)
    
    return dx * dy * sum(P) * 1e6




def plot_electric_field_vs_time(efield, show=True, ofile='', ret=False):
    """Plot the electric field as a function of time"""

    T  = [item[0]    for item in efield]
    Ex = [item[1][0] for item in efield]
    Ey = [item[1][1] for item in efield]
    Ez = [item[1][2] for item in efield]

    plt.figure(1, figsize=(18, 9))
    plt.subplot(131)
    plt.plot(T, Ex)
    plt.xlabel('T [s]')
    plt.ylabel('Ex')
    plt.title('Electric Field (Ex)')

    plt.figure(1, figsize=(18, 9))
    plt.subplot(132)
    plt.plot(T, Ey)
    plt.xlabel('T [s]')
    plt.ylabel('Ey')
    plt.title('Electric Field (Ey)')

    plt.figure(1, figsize=(18, 9))
    plt.subplot(133)
    plt.plot(T, Ez)
    plt.xlabel('T [s]')
    plt.ylabel('Ez')
    plt.title('Electric Field (Ez)')

    if ret:
        return plt
    return
    
    
    



def plot_undulator_flux_onaxis(oth, period, nperiods, harmonics, minimum=0, bfield_range=None, k_range=None, show=True, ret=False, title='Flux On-Axis', npoints=500, figsize=None, ofile=None):
    '''Plot the on-axis flux of an undulator.  More docstring needed'''

    if bfield_range is None and k_range is None:
        raise ValueError('bfield_range or k_range must be defined')

    if bfield_range is not None and k_range is not None:
        raise ValueError('bfield_range and k_range cannot both be defined.  pick one or the other')

    if bfield_range is not None:
        if len(bfield_range) is not 2:
            raise ValueError('bfield_range must be: [min, max]')
        else:
            k_range = [oth.undulator_K(bfield_range[0], period), oth.undulator_K(bfield_range[1], period)]

    # Size and limits
    plt.figure(1, figsize=figsize)
    plt.loglog()


    # Loop over all harmonics
    for i in harmonics:
        if i % 2 != 1:
            continue

        fl = oth.undulator_flux_onaxis(
            K_range=k_range,
            period=period,
            nperiods=nperiods,
            harmonic=i,
            npoints=npoints
            )

        X = [f[0] for f in fl]
        Y = [f[1] for f in fl]

        plt.plot(X, Y, label=str(i))

    plt.legend(title='Harmonics')

    plt.xlabel('Energy [eV]')
    plt.ylabel('[$\gamma / mrad^2 / 0.1\%bw / s$]')
    plt.title(title)
    
    if ofile is not None:
        plt.savefig(ofile, bbox_inches='tight', transparent=transparent)

    if show == True:
        plt.show()

    if ret:
        return plt
    return



def plot_undulator_brightness(oth, period, nperiods, harmonics, minimum=0, bfield=None, K=None, show=True, ret=False, title='Brightness', figsize=None, ofile=None):
    '''Plot the brightness of an undulator.'''

    if bfield is None and K is None:
        raise ValueError('bfield or K must be defined')

    if bfield is not None and K is not None:
        raise ValueError('bfield and K cannot both be defined.  pick one or the other')

    if bfield is not None:
        if len(bfield) is not 2:
            raise ValueError('bfield must be: [min, max]')
        else:
            K = [oth.undulator_K(bfield[0], period), oth.undulator_K(bfield[1], period)]

    # Size and limits
    plt.figure(1, figsize=figsize)
    plt.loglog()


    # Loop over all harmonics
    for i in harmonics:
        
        R = []
        for k in np.linspace(K[1], K[0], 300):
            ev_brightness = oth.undulator_brightness(K=k,
                                                     period=period,
                                                     nperiods=nperiods,
                                                     harmonic=i
                                                    )
            if ev_brightness[1] >= minimum:
                R.append(ev_brightness)

        X = []
        Y = []
        for j in range(len(R)):
            X.append(R[j][0])
            Y.append(R[j][1])
        plt.plot(X, Y, label=str(i))

    plt.legend(title='Harmonics')

    plt.xlabel('Energy [eV]')
    plt.ylabel('[$\gamma / mm^2 / mrad^2 / 0.1\%bw / s$]')
    plt.title(title)
    
    if ofile is not None:
        plt.savefig(ofile, bbox_inches='tight', transparent=True)

    if show == True:
        plt.show()

    if ret:
        return plt
    return





def plot_flux_spectrum(F, S, energy=None, title='Flux [$\gamma / mm^2 / 0.1\%bw / s]$', xlabel='X1 Axis [$m$]', ylabel='X2 Axis [$m$]', show=True, ofile='', figsize=[10, 3], ylim=None, xlim=None, colorbar=True, ret=False):
    """Plot a 2D histogram with equal spacing"""
        
    X = [item[0][0] for item in F]
    Y = [item[0][1] for item in F]
    P = [item[1]    for item in F]

    NX = len(np.unique(X))
    NY = len(np.unique(Y))

    # Size and limits
    plt.figure(1, figsize=figsize)
    plt.subplot(121)
    if ylim is not None: plt.ylim(ylim[0], ylim[1])
    if xlim is not None: plt.xlim(xlim[0], xlim[1])

    plt.hist2d(X, Y, bins=[NX, NY], weights=P)
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    if colorbar is True:
        plt.colorbar(format='%.0e')

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)



    plt.subplot(122)

    X2 = [item[0] for item in S]
    Y2 = [item[1] for item in S]
    plt.plot(X2, Y2)
    plt.xlabel('Energy [eV]')
    plt.ylabel('$\gamma / mm^2 / 0.1\%bw / s$')
    plt.title('Spectrum')

    if energy is not None:
        plt.axvline(x=energy, color='red')

    plt.figure(1)
    if ofile != '':
        plt.savefig(ofile, bbox_inches='tight')

    if show == True:
        plt.show()

    if ret:
        return plt

    return


