from __future__ import print_function

import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from math import sqrt
import copy

def power_density_3d(srs, surface,
                    normal=1, rotations=[0, 0, 0], translation=[0, 0, 0], nparticles=0, gpu=0, nthreads=0,
                    title='Power Density [$W / mm^2$]', xlim=None, ylim=None, zlim=None, colorbar=True, figsize=None,
                    alpha=0.4, ofile=None, show=True, view_init=None, axis=None, transparent=True,
                    xticks=None, yticks=None, zticks=None, bbox_inches='tight', max_level=24, quantity='power density'):
    """calculate power density for and plot a parametric surface in 3d"""

    points = []


    for u in np.linspace(surface.ustart, surface.ustop, surface.nu):
        for v in np.linspace(surface.vstart, surface.vstop, surface.nv):
            points.append([surface.position(u, v), surface.normal(u, v)])


    power_density = srs.calculate_power_density(points=points, normal=normal, rotations=rotations, translation=translation, nparticles=nparticles, gpu=gpu, nthreads=nthreads, max_level=max_level, quantity=quantity)
    P = [item[1] for item in power_density]

    X2 = []
    Y2 = []
    Z2 = []
    for i in range(surface.nu):
        tX = []
        tY = []
        tZ = []
        for j in range(surface.nv):
            tX.append(power_density[i * surface.nv + j][0][0])
            tY.append(power_density[i * surface.nv + j][0][1])
            tZ.append(power_density[i * surface.nv + j][0][2])
        X2.append(tX)
        Y2.append(tY)
        Z2.append(tZ)

    colors =[]
    MAXP = max(P)
    MINP = min(P)
    if MAXP == 0:
        MAXP=1
    PP = []
    for i in range(surface.nu):
        tmpP = []
        tmpPP = []
        for j in range(surface.nv):
            tmpP.append(P[i * surface.nv + j] / MAXP)
            tmpPP.append(P[i * surface.nv + j])

        colors.append(tmpP)
        PP.append(tmpPP)
  


    fig = plt.figure(1, figsize=figsize)
    ax = fig.gca(projection = '3d')
    if view_init is not None:
        ax.view_init(view_init[0], view_init[1])


    if xticks is not None:
        ax.set_xticks(xticks)
    if zticks is not None:
        ax.set_yticks(zticks)
    if yticks is not None:
        ax.set_zticks(yticks)

    #ax.view_init(30, -40)
    ax.set_xlabel('X [m]')
    ax.set_ylabel('Z [m]')
    ax.set_zlabel('Y [m]')
    #ax.set_ylim(translation_[2] - 0.1, translation_[2] + 0.1)
    #ax.set_axis_off()

    if xlim is not None:
        ax.set_xlim(xlim[0], xlim[1])
    if ylim is not None:
        ax.set_zlim(ylim[0], ylim[1])
    if zlim is not None:
        ax.set_ylim(zlim[0], zlim[1])


    ax.plot_surface(np.asarray(X2), np.asarray(Z2), np.asarray(Y2), facecolors=cm.viridis(colors), rstride=1, cstride=1, alpha=alpha, linewidth=0)
    ax.invert_xaxis()

    if axis is not None:
        plt.axis(axis)


    m = cm.ScalarMappable()
    m.set_array(PP)
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.ticklabel_format(style='sci', axis='z', scilimits=(0,0))
    if colorbar is True:
        sm = plt.cm.ScalarMappable(cmap=cm.viridis, norm=plt.Normalize(vmin=MINP, vmax=MAXP))
        sm._A = []
        plt.colorbar(sm, format='%.0e')


    plt.title(title)

    if ofile is not None:
        plt.savefig(ofile, bbox_inches=bbox_inches, transparent=transparent)

    if show is True:
        plt.show()

    return power_density













def power_density_3ds(srs, surface, ax,
                    normal=1, rotations=[0, 0, 0], translation=[0, 0, 0], nparticles=0, gpu=0, nthreads=0,
                    alpha=0.4, transparent=True, max_level=-2, ofile=''):
    """calculate power density for and plot a parametric surface in 3d"""

    points = []


    for u in np.linspace(surface.ustart, surface.ustop, surface.nu):
        for v in np.linspace(surface.vstart, surface.vstop, surface.nv):
            points.append([surface.position(u, v), surface.normal(u, v)])


    power_density = srs.calculate_power_density(points=points, normal=normal, rotations=rotations, translation=translation, nparticles=nparticles, gpu=gpu, nthreads=nthreads, max_level=max_level, ofile=ofile)
    P = [item[1] for item in power_density]

    X2 = []
    Y2 = []
    Z2 = []
    for i in range(surface.nu):
        tX = []
        tY = []
        tZ = []
        for j in range(surface.nv):
            tX.append(power_density[i * surface.nv + j][0][0])
            tY.append(power_density[i * surface.nv + j][0][1])
            tZ.append(power_density[i * surface.nv + j][0][2])
        X2.append(tX)
        Y2.append(tY)
        Z2.append(tZ)

    colors =[]
    MAXP = max(P)
    PP = []
    for i in range(surface.nu):
        tmpP = []
        tmpPP = []
        for j in range(surface.nv):
            tmpP.append(P[i * surface.nv + j] / MAXP)
            tmpPP.append(P[i * surface.nv + j])

        colors.append(tmpP)
        PP.append(tmpPP)
  



    #ax.invert_xaxis()
    return ax.plot_surface(np.asarray(X2), np.asarray(Z2), np.asarray(Y2), facecolors=cm.viridis(colors), rstride=1, cstride=1, alpha=alpha)













def plot_bfield2D (srs, xlim=[-0.001, 0.001], zlim=[-1, 1], nx=20, nz=50):
    """plot the bfield in 3D vector form"""

    xx = []
    zz = []
    uu = []
    ww = []
    cc = []
    BMax = 1.



#    for i in np.linspace(xlim[0], xlim[1], nx):
#        for k in np.linspace(zlim[0], zlim[1], nz):
#            B = srs.get_bfield([i, 0, k])
#            BMag = sqrt(B[0]*B[0] + B[1]*B[1] + B[2]*B[2])
#            if BMag > BMax:
#                BMax = BMag
    
    for i in np.linspace(xlim[0], xlim[1], nx):
        xt = []
        zt = []
        ut = []
        wt = []
        ct = []
        for k in np.linspace(zlim[0], zlim[1], nz):
            xt.append(i)
            zt.append(k)

            B = srs.get_bfield([i, 0, k])
            BMag = sqrt(B[0]*B[0] + B[1]*B[1] + B[2]*B[2])

            ut.append(B[0])
            wt.append(B[2])
            ct.append(BMag / BMax)

        xx.append(xt)
        zz.append(zt)
        uu.append(ut)
        ww.append(wt)
        cc.append(ct)




    plt.figure()

    #plt.quiver(xx, zz, uu, ww, cc, cmap=cm.jet)
    plt.quiver(xx, zz, uu, ww)


    plt.show()

    return


def plot_bfield3D (srs, xlim=[-0.02, 0.02], ylim=[-0.02, 0.02], zlim=[-0.2, 0.02], nx=10, ny=10, nz=10):
    """plot the bfield in 3D vector form"""

    BMax = 0.

    fig = plt.figure()
    ax = fig.gca(projection='3d')


    dx = (xlim[1] - xlim[0]) / (nx - 1)
    dy = (ylim[1] - ylim[0]) / (ny - 1)
    dz = (zlim[1] - zlim[0]) / (nz - 1)

    dl = min([dx, dy, dz])

    for i in np.linspace(xlim[0], xlim[1], nx):
        for j in np.linspace(ylim[0], ylim[1], ny):
            for k in np.linspace(zlim[0], zlim[1], nz):
                B = srs.get_bfield([i, j, k])
                BMag = sqrt(B[0]*B[0] + B[1]*B[1] + B[2]*B[2])
                if BMag > BMax:
                    BMax = BMag


    length_scale = dl / BMax

    for i in np.linspace(xlim[0], xlim[1], nx):
        for j in np.linspace(ylim[0], ylim[1], ny):
            for k in np.linspace(zlim[0], zlim[1], nz):
                B = srs.get_bfield([i, j, k])
                BMag = sqrt(B[0]*B[0] + B[1]*B[1] + B[2]*B[2])
                ax.quiver([[[i]]], [[[j]]], [[[k]]], [[[B[0]]]], [[[B[1]]]], [[[B[2]]]], length=BMag*length_scale, cmap='Reds')




    ax.set_title('3D Vector Field')             # title
    ax.view_init(elev=18, azim=30)              # camera elevation and angle
    ax.dist=8                                   # camera distance

    plt.show()

    return






def plot_surface(surface, xlim=None, ylim=None, zlim=None, alpha=0.5, **kwargs):
    """plot a parametric surface in 3d"""


    X2 = []
    Y2 = []
    Z2 = []
    for u in np.linspace(surface.ustart, surface.ustop, surface.nu):
        tX = []
        tY = []
        tZ = []
        for v in np.linspace(surface.vstart, surface.vstop, surface.nv):
            point = surface.position(u, v)
            tX.append(point[0])
            tY.append(point[1])
            tZ.append(point[2])
        X2.append(tX)
        Y2.append(tY)
        Z2.append(tZ)


    fig = plt.figure(1)
    ax = fig.gca(projection = '3d')
    ax.set_xlabel('X [m]')
    ax.set_ylabel('Z [m]')
    ax.set_zlabel('Y [m]')

    if xlim is not None:
        ax.set_xlim(xlim[0], xlim[1])
    if ylim is not None:
        ax.set_zlim(ylim[0], ylim[1])
    if zlim is not None:
        ax.set_ylim(zlim[0], zlim[1])


    ax.plot_surface(np.asarray(X2), np.asarray(Z2), np.asarray(Y2), rstride=1, cstride=1, alpha=0.5, **kwargs)
    ax.invert_xaxis()

    plt.show()

    return




def plot_trajectory3d(trajectory, figsize=None):
    """Plot the trajectory in 3D"""

    #mpl.rcParams['legend.fontsize'] = 10

    fig = plt.figure(1, figsize=figsize)
    ax = fig.gca(projection='3d')

    # Get coordinate lists
    X  = [item[0][0] for item in trajectory]
    Y  = [item[0][1] for item in trajectory]
    Z  = [item[0][2] for item in trajectory]

    ax.set_xlabel('X [m]')
    ax.set_ylabel('Z [m]')
    ax.set_zlabel('Y [m]')

    plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.ticklabel_format(style='sci', axis='z', scilimits=(0,0))

    ax.plot(X, Z, Y, label='Trajectory')
    ax.invert_xaxis()

    ax.view_init(elev=30, azim=5)              # camera elevation and angle

    plt.show()






def get_surface_points (surface=None, surfaces=None, translation=None):

    points = []
    if surface is not None and surfaces is not None:
        raise ValueError('pick one: surface or surfaces')

    if surface is not None:
        surfaces = [surface]

    for surface in surfaces:
        for u in np.linspace(surface.ustart, surface.ustop, surface.nu):
            for v in np.linspace(surface.vstart, surface.vstop, surface.nv):
                p = surface.position(u, v)
                if translation is not None:
                    p[0] += translation[0]
                    p[1] += translation[1]
                    p[2] += translation[2]
                points.append([p, surface.normal(u, v)])

    return points






def plot_power_density_scatter (V, s=10, title='Power Density [$W/mm^2$]', figsize=None, bbox_inches='tight', alpha=1, transparent=True, ofile=None, colorbar=True):

    if len(V) == 0:
        return

    X = [item[0][0] for item in V]
    Y = [item[0][1] for item in V]
    Z = [item[0][2] for item in V]
    P = [item[1]    for item in V]

    pmax = max(P)
    pmin = min(P)

    C = []
    if pmax == pmin:
        C=[0] * len(P)
    else:
        C = [p / pmax for p in P]

    Cen3D = plt.figure(figsize=figsize)
    ax = Cen3D.add_subplot(111, projection='3d')
    plt.title(title)

    ax.scatter(X, Z, Y, c=C, s=s, alpha=alpha)
    ax.invert_xaxis()

    ax.set_xlabel('X [m]')
    ax.set_ylabel('Z [m]')
    ax.set_zlabel('Y [m]')

    m = cm.ScalarMappable()
    m.set_array(P)
    if colorbar is True:
        sm = plt.cm.ScalarMappable(cmap=cm.viridis, norm=plt.Normalize(vmin=pmin, vmax=pmax))
        sm._A = []
        plt.colorbar(sm, format='%.0e')

    if ofile is not None:
        plt.savefig(ofile, bbox_inches=bbox_inches, transparent=transparent)

    plt.show()

    return









def plot_power_density_stl (P, title='Power Density [$W/mm^2$]', elev=30, azim=30, alpha=0.8, bbox_inches='tight', transparent=True, ofile=None, colorbar=True):

    pmax = max([p[1] for p in P])
    pmin = min([p[1] for p in P])

    m = cm.ScalarMappable(cm.viridis)
#m.set_array([p[1]/pmax for p in P])

    xmax = ymax = zmax = -99999
    xmin = ymin = zmin = +99999
    for p in P:
        xmax = max([xmax, p[0][0][0], p[0][1][0], p[0][2][0]])
        xmin = min([xmin, p[0][0][0], p[0][1][0], p[0][2][0]])
        zmax = max([zmax, p[0][0][1], p[0][1][1], p[0][2][1]])
        zmin = min([zmin, p[0][0][1], p[0][1][1], p[0][2][1]])
        ymax = max([ymax, p[0][0][2], p[0][1][2], p[0][2][2]])
        ymin = min([ymin, p[0][0][2], p[0][1][2], p[0][2][2]])

    fig = plt.figure()
    ax = Axes3D(fig)
    for p in P:
        pp = copy.deepcopy(p)
        pp[0][0][1] = p[0][0][2]
        pp[0][0][2] = p[0][0][1]
        pp[0][1][1] = p[0][1][2]
        pp[0][1][2] = p[0][1][1]
        pp[0][2][1] = p[0][2][2]
        pp[0][2][2] = p[0][2][1]
        triangle = Poly3DCollection([pp[0]], alpha=alpha, edgecolors='b', linewidths=0.05)
        triangle.set_facecolor(m.to_rgba(pp[1]/pmax))
        #triangle.set_facecolor([1, 1-p[1]/pmax, 1-p[1]/pmax])
        ax.add_collection3d(triangle)
 
    #ax.view_init(elev=elev, azim=azim)
    plt.title(title)
    ax.invert_xaxis()
    ax.set_xlim([xmin, xmax])
    ax.set_ylim([ymin, ymax])
    ax.set_zlim([zmin, zmax])
    ax.set_xlabel('X [m]')
    ax.set_ylabel('Z [m]')
    ax.set_zlabel('Y [m]')

    if colorbar is True:
        sm = plt.cm.ScalarMappable(cmap=cm.viridis, norm=plt.Normalize(vmin=pmin, vmax=pmax))
        sm._A = []
        plt.colorbar(sm, format='%.0e')

    if ofile is not None:
        plt.savefig(ofile, bbox_inches=bbox_inches, transparent=transparent)


    plt.show()

    return







