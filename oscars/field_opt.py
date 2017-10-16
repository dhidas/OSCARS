import math
import copy
import numpy as np
import scipy.optimize as op
import scipy.integrate as integrate
import scipy.spatial.distance as sd
import oscars.sr
from oscars.plots_mpl import *
from oscars.fit import *
from oscars.util import *

def ideal_undulator(osr, field, length, period, pieces, height, gap, 
                    packing, t_distance, t_width, t_field):
    """Assign the field of an ideal undulator to an sr object.

    Keyword arguments:
    osr -- oscars sr object to which the undulator field will be added.
    field -- absolute value of maximum field strength of each magnet
    length -- length of undulator centered at 0 without term. magnets
    period -- length of one period in undulator's  B-field
    pieces -- number of discrete magnets forming oscillating field
    packing -- factor describing closeness of magnets
    t_distance -- distance from origin to each terminating magnet
    t_width -- width of each terminating magnet
    t_field -- |maximum field strength| of each terminating magnet
    """
    br, l, p, m, h, g, e, td, tw, tbr = (field, length, period, pieces,
                                         height, gap, packing, t_distance,
                                         t_width, t_field)
    w = e * p / m
    d = l / 2.0
    constants = (-2 * br * math.sin(e*math.pi/m) * math.e**(-math.pi*g/p)
                 *(1 - math.e**(-2 * math.pi * h / p)) * m / math.pi)
    
    def f_und(z):
        if -l/2 <= z <= l/2:
            return [0, constants * math.cos(2 * math.pi * z / p), 0]
        else:
            return [0, 0, 0]
    
    osr.add_bfield_function(f_und) 
    return osr


def b_y_pre_op_ints(osr=None):
    """Integrate B field of unoptimized undulator.

    Keyword arguments:
    pmu -- object of class Undulator describing physical parameters
    osr -- oscars sr object with a b field added to it
    """
    lower_bound = -1.33
    upper_bound = 1.34
    a = integrate.quad(lambda z: osr.get_bfield([0, 0, z])[1],
                       lower_bound, upper_bound)
    b = integrate.dblquad(lambda y, x: osr.get_bfield([0, 0, z][1]),
                         lower_bound, upper_bound,
                         lambda x: lower_bound, lambda x: x)
    return [a[0], b[0]]
    

def b_y_pre_op_traj(osr=None):
    """Return unoptimized B field function for trajectory.
    
    Keyword arguments:
    pmu -- object of class Undulator describing physical parameters
    osr -- oscars sr object with a b field added to it
    """
    lowest_bound = -1.35
    lower_bound = -1.33
    upper_bound = 1.34
    uppest_bound = 1.36   
    
    # Return unoptimized B field function
    def field(x, y, z, t):
        if z >= lower_bound and z <= upper_bound:
            return [0, osr.get_bfield([0, 0, z])[1], 0]
        elif z >= lowest_bound and z < lower_bound:
            return [0, 0, 0]
        elif z > upper_bound and z <= uppest_bound:
            return [0, 0, 0]
        else:
            return [0, 0, 0]
    
    return field
    

def b_y_pre_op_plot(z_list, osr=None):
    """Return unoptimized B field function for plotting.
    
    Keyword arguments:
    z_list -- sequence of all z values of interest
    pmu -- object of class Undulator describing physical parameters
    osr -- oscars sr object with a b field added to it
    """
    lower_bound = -1.33
    upper_bound = 1.34
    
    return [osr.get_bfield([0, 0, ]z)[1]  if 
            (z >= lower_bound and z <= upper_bound) else 0) for z in z_list]
    

def b_y(osr=None, sr_info=None):
    """Return optimized B field function.
    
    Keyword arguments:
    osr -- oscars sr object with a b field added to it
    sr_info -- if using sr,  array containing the following info:
        [length, t_distance, t_width, t_field]
        length -- length of undulator centered at 0 without t mags
        t_distance -- distance from origin to each terminating magnet
        t_width -- width of each terminating magnet
        t_field -- max field strength of each terminating magnet
    """
    # Perform single integral on B_y
    def b_y_integral_1(beta, gamma):
        return (integrate.quad(lambda z: beta[1], 
                               lowest_bound, lower_bound)[0]
                + integrate.quad(lambda z: osr.get_bfield([0, 0, z])[1],  
                                 lower_bound, upper_bound)[0]
                + integrate.quad(lambda z: gamma[1], 
                                 upper_bound, uppest_bound)[0])
    
    # Perform double integral on B_y
    def b_y_integral_2(beta, gamma):
        return  (integrate.dblquad(lambda y, x: beta[1], 
                                   lowest_bound, lower_bound, 
                                   lambda x: lowest_bound, lambda x: x)[0]
                 + integrate.dblquad(lambda y, x: osr.get_bfield([0, 0, y])[1], 
                                     lower_bound, upper_bound, 
                                     lambda x: lower_bound, lambda x: x)[0]
                 + integrate.dblquad(lambda y, x: gamma[1], 
                                     upper_bound, uppest_bound,
                                     lambda x: upper_bound, lambda x: x)[0])
    
    # Perform single integral on B_x
    def b_x_integral_1(beta, gamma):
        return (integrate.quad(lambda z: beta[0],
                               lowest_bound, lower_bound)[0]
                + integrate.quad(lambda z: osr.get_bfield([0, 0, z])[0], 
                                 lower_bound, upper_bound)[0]
                + integrate.quad(lambda z: gamma[0],
                                 upper_bound, uppest_bound)[0])
    
    # Perform double integral on B_x
    def b_x_integral_2(beta, gamma):
        return (integrate.dblquad(lambda y, x: beta[0],
                                  lowest_bound, lower_bound,
                                  lambda x: lowest_bound, lambda x: x)[0]
                + integrate.dblquad(lambda y, x: osr.get_bfield([0, 0, y])[1],
                                    lower_bound, upper_bound,
                                    lambda x: lower_bound, lambda x: x)[0]
                + integrate.dblquad(lambda y, x: gamma[0],
                                    upper_bound, uppest_bound,
                                    lambda x: upper_bound, lambda x: x)[0])
    
    # Return weighted sum of the first and second integral of B_y
    def field_ints_norm(beta, gamma):
        w1, w2, w3, w4  = (0.25, 0.25, 0.25, 0.25) 
        a = b_y_integral_1(beta, gamma)
        b = b_y_integral_2(beta, gamma)
        c = b_x_integral_1(beta, gamma)
        d = b_x_integral_2(beta, gamma)
        return (w1*math.fabs(a) + w2*math.fabs(b) + w3*math.fabs(c) 
                + w4*math.fabs(d))
    
    # Return the minimum distance between a point and the trajectory
    def min_dist(beta, gamma, osr=None, point=[0, 0, 0]):
        osr.remove_bfield('term')
        osr.add_bfield_function(b_y_term(beta, gamma), name='term')
        global trajectory 
        trajectory = osr.calculate_trajectory()
        trajec = []
        for i in range(len(trajectory)):
            trajec.append(trajectory[i][1])
        point = [point,]
        distances = sd.cdist(point, trajec) 
        return (np.argmin(distances), np.amin(distances))
    
    # Return the weighted sum of the position and direction differences
    def traj_norm(beta, gamma, osr, point):
        (min_index, min_distance) = min_dist(beta, gamma, osr, point)
        vel = trajectory[min_index][2]
        vel_ideal = np.array([0, 0, 1])
        vel = np.array(vel)
        cross_product = np.cross(vel, vel_ideal)
        cross_product_norm = np.linalg.norm(cross_product)
        return math.fabs(min_distance) + math.fabs(cross_product_norm)
    
    # Return function for only the terminating fields.
    def b_y_term(beta, gamma):
        def field(x, y, z, t):
            if z >= lowest_bound and z < lower_bound:
                return [beta[0], beta[1], 0]
            elif z > upper_bound and z <= uppest_bound:
                return [gamma[0], gamma[1], 0]
            else:
                return [0, 0, 0]
        return field
    
    # Return callable for weighted sum function depending on pmu or osr
    def min_fun(u):
        return traj_norm(point=[0,0,5], osr=osr, beta=[u[0], u[1]], 
                             gamma=[u[2], u[3]])
    
    # Return field function in compatible format for oscars
    def field(x, y, z, t):    
        if z >= lower_bound and z <= upper_bound:
            return [osr.get_bfield([0, 0, z])[0], osr.get_bfield([0, 0, z])[1], 0]
        elif z >= lowest_bound and z < lower_bound:
            return [beta[0], beta[1], 0]
        elif z > upper_bound and z <= uppest_bound:
            return [gamma[0], gamma[1], 0]
        else:
            return [0, 0, 0]

    lowest_bound = -sr_info[1] - sr_info[2] # -1.35
    lower_bound = -sr_info[1] # -1.33
    upper_bound = sr_info[1] # 1.34
    uppest_bound = sr_info[1] + sr_info[2] # 1.36
    field_bound = sr_info[3]
    
    cons = ({'type' : 'ineq', 'fun' : lambda u: field_bound - math.fabs(u[0])},
            {'type' : 'ineq', 'fun' : lambda u: field_bound - math.fabs(u[1])},
            {'type' : 'ineq', 'fun' : lambda u: field_bound - math.fabs(u[2])},
            {'type' : 'ineq', 'fun' : lambda u: field_bound - math.fabs(u[3])})
    op_f = op.minimize(fun=lambda u: min_fun(u), 
                       x0=(0, 0, 0, 0),
                       constraints=cons,
                       tol=1e-10,
                       options={'disp' : True})
    beta  = [op_f.x[0], op_f.x[1]]
    gamma = [op_f.x[2], op_f.x[3]]
    print('\nSolution array: ' + str(op_f.x))
    plot_trajectory_position(trajectory)
    plot_trajectory_velocity(trajectory)
    plot_bfield(osr, -sr_info[1]-2*sr_info[2], sr_info[1]+2*sr_info[2])
    
    return field
