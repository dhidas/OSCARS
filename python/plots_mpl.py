from matplotlib.colors import LogNorm
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
    X  = [item[0][0] for item in trajectory]
    Y  = [item[0][1] for item in trajectory]
    Z  = [item[0][2] for item in trajectory]

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
    VX = [item[1][0] for item in trajectory]
    VY = [item[1][1] for item in trajectory]
    VZ = [item[1][2] for item in trajectory]
    T = range(len(VX))

    # Plot VX, VY, VZ vs. T
    plt.figure(1, figsize=figsize)
    plt.subplot(131)
    plt.plot(T, VX)
    plt.xlabel('T [step]')
    plt.ylabel('BX []')
    plt.title('Particle Beta X')

    plt.subplot(132)
    plt.plot(T, VY)
    plt.xlabel('T [step]')
    plt.ylabel('BY []')
    plt.title('Particle Beta Y')

    plt.subplot(133)
    plt.plot(T, VZ)
    plt.xlabel('T [step]')
    plt.ylabel('BZ []')
    plt.title('Particle Beta Z')

    if ofile != '':
        plt.savefig(ofile, bbox_inches='tight')

    if show == True:
        plt.show()

    if ret:
        return plt
    return
    
    
def plot_power_density(V, title='Power Density [$W / mm^2$]', xlabel='X1 Axis [$m$]', ylabel='X2 Axis [$m$]', show=True, ofile='', figsize=None, ret=False):
    """Plot a 2D histogram with equal spacing"""
        
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


def plot_flux(V, title='Flux [$\gamma / mm^2 / 0.1\%bw / s]$', xlabel='X1 Axis [$m$]', ylabel='X2 Axis [$m$]', show=True, ofile='', figsize=None, ylim=None, xlim=None, colorbar=True, ret=False):
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

    plt.figure(1, figsize=figsize)
    plt.hist2d(X, Y, bins=[NX, NY], weights=P)
    #cb.formatter.set_scientific(True)
    #cb.formatter.set_powerlimits((0, 0))
    #cb.ax.yaxis.set_offset_position('right')
    #cb.update_ticks()
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    if colorbar is True:
        plt.colorbar(format='%.0e')

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


def plot_spectrum(S, log=False, show=True, ofile='', title='Spectrum', figsize=None, ylim=None, xlim=None, transparent=True, ret=False, **kwargs):
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
    plt.xlabel('Energy [eV]')
    plt.ylabel('$\gamma / mm^2 / 0.1\%bw / s$')
    plt.title(title)

    if ofile != '':
        plt.savefig(ofile, bbox_inches='tight', transparent=transparent)

    if show == True:
        plt.show()

    if ret is True:
        return plt
    return



def plot_spectra(spectra, label=None, show=True, ofile='', title='', loc='upper left', log=False, xlabel='Energy [eV]', ylabel='[$\gamma / mm^2 / 0.1\%bw / s$]', figsize=None, ylim=None, xlim=None, ret=False, axis=None, transparent=True, **kwargs):


    # Size and limits
    plt.figure(1, figsize=figsize)
    if ylim is not None: plt.ylim(ylim[0], ylim[1])
    if xlim is not None: plt.xlim(xlim[0], xlim[1])

    for i in range(len(spectra)):
        s = spectra[i]
        
        X = [item[0] for item in s]
        Y = [item[1] for item in s]
        if label is not None:
            if label[i] is not None:
                plt.plot(X, Y, label=label[i], **kwargs)
            else:
                plt.plot(X, Y, **kwargs)
        else:
            plt.plot(X, Y, **kwargs)


    if log:
        plt.yscale('log')
        
    plt.legend(loc=loc)
    if xlabel is not None:
        plt.xlabel(xlabel)
    if ylabel is not None:
        plt.ylabel(ylabel)
    if title == None:
        title='Spectrum'
    plt.title(title)
    
    if axis is not None:
        plt.axis(axis)

    if ofile != '':
        plt.savefig(ofile, bbox_inches='tight', transparent=transparent)

    if show == True:
        plt.show()

    if ret:
        return plt
    return





def plot_bfield(osc, mymin=-1, mymax=1, ylim=None, show=True, ofile='', axis='Z', npoints=20000, between_two_points=None, ret=False):
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

            Bx.append(osc.get_bfield([x, y, z])[0])
            By.append(osc.get_bfield([x, y, z])[1])
            Bz.append(osc.get_bfield([x, y, z])[2])
            axis = 'Position'
    else:
        P = np.linspace(mymin, mymax, npoints)
        if axis is 'X':
            Bx = [osc.get_bfield([p, 0, 0])[0] for p in P]
            By = [osc.get_bfield([p, 0, 0])[1] for p in P]
            Bz = [osc.get_bfield([p, 0, 0])[2] for p in P]
        elif axis is 'Y':
            Bx = [osc.get_bfield([0, p, 0])[0] for p in P]
            By = [osc.get_bfield([0, p, 0])[1] for p in P]
            Bz = [osc.get_bfield([0, p, 0])[2] for p in P]
        elif axis is 'Z':
            Bx = [osc.get_bfield([0, 0, p])[0] for p in P]
            By = [osc.get_bfield([0, 0, p])[1] for p in P]
            Bz = [osc.get_bfield([0, 0, p])[2] for p in P]
        else:
            raise

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







def plot_efield(osc, mymin=-1, mymax=1, ylim=None, show=True, ofile='', axis='Z', npoints=20000, between_two_points=None, ret=False):
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

            Bx.append(osc.get_efield([x, y, z])[0])
            By.append(osc.get_efield([x, y, z])[1])
            Bz.append(osc.get_efield([x, y, z])[2])
            axis = 'Position'
    else:
        P = np.linspace(mymin, mymax, npoints)
        if axis is 'X':
            Bx = [osc.get_efield([p, 0, 0])[0] for p in P]
            By = [osc.get_efield([p, 0, 0])[1] for p in P]
            Bz = [osc.get_efield([p, 0, 0])[2] for p in P]
        elif axis is 'Y':
            Bx = [osc.get_efield([0, p, 0])[0] for p in P]
            By = [osc.get_efield([0, p, 0])[1] for p in P]
            Bz = [osc.get_efield([0, p, 0])[2] for p in P]
        elif axis is 'Z':
            Bx = [osc.get_efield([0, 0, p])[0] for p in P]
            By = [osc.get_efield([0, 0, p])[1] for p in P]
            Bz = [osc.get_efield([0, 0, p])[2] for p in P]
        else:
            raise

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
    
    
    



