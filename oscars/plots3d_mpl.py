import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from math import sqrt

def power_density_3d(srs, surface,
                    normal=1, rotations=[0, 0, 0], translation=[0, 0, 0], nparticles=0, gpu=0, nthreads=0, ret=False,
                    title='Power Density [$W / mm^2$]', xlim=None, ylim=None, zlim=None, colorbar=True, figsize=None, alpha=0.4, ofile=None, show=True, view_init=None, axis=None, transparent=True, xticks=None, yticks=None, zticks=None, bbox_inches='tight'):
    """calculate power density for and plot a parametric surface in 3d"""

    points = []


    for u in np.linspace(surface.ustart, surface.ustop, surface.nu):
        for v in np.linspace(surface.vstart, surface.vstop, surface.nv):
            points.append([surface.position(u, v), surface.normal(u, v)])


    power_density = srs.calculate_power_density(points=points, normal=normal, rotations=rotations, translation=translation, nparticles=nparticles, gpu=gpu, nthreads=nthreads)
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


    ax.plot_surface(X2, Z2, Y2, facecolors=cm.viridis(colors), rstride=1, cstride=1, alpha=alpha)
    ax.invert_xaxis()

    if axis is not None:
        plt.axis(axis)


    m = cm.ScalarMappable()
    m.set_array(PP)
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.ticklabel_format(style='sci', axis='z', scilimits=(0,0))
    if colorbar is True:
        plt.colorbar(m, format='%.0e')
        #plt.colorbar(m, format=r'%3.1f')


    plt.title(title)

    if ofile is not None:
        plt.savefig(ofile, bbox_inches=bbox_inches, transparent=transparent)

    if show is True:
        plt.show()

    if ret:
        return plt
    return













def power_density_3ds(srs, surface, ax,
                    normal=1, rotations=[0, 0, 0], translation=[0, 0, 0], nparticles=0, gpu=0, nthreads=0,
                    alpha=0.4, transparent=True):
    """calculate power density for and plot a parametric surface in 3d"""

    points = []


    for u in np.linspace(surface.ustart, surface.ustop, surface.nu):
        for v in np.linspace(surface.vstart, surface.vstop, surface.nv):
            points.append([surface.position(u, v), surface.normal(u, v)])


    power_density = srs.calculate_power_density(points=points, normal=normal, rotations=rotations, translation=translation, nparticles=nparticles, gpu=gpu, nthreads=nthreads)
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
    return ax.plot_surface(X2, Z2, Y2, facecolors=cm.viridis(colors), rstride=1, cstride=1, alpha=alpha)













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






def plot_surface(surface, xlim=None, ylim=None, zlim=None, **kwargs):
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


    ax.plot_surface(X2, Z2, Y2, rstride=1, cstride=1, alpha=0.5, **kwargs)
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






























