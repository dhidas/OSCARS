from __future__ import print_function

import oscars.plots_mpl

from scipy.optimize import curve_fit, minimize
import numpy as np
import warnings

try:
    import matplotlib.pyplot as plt
except ImportError:
    warnings.warn('matplotlib cannot be imported')




def find_harmonics (spectrum, first=None, xwidth=50, parity='all', figsize=None, quiet=True, ofile=None, show=True):

    if first is None:
        first = find_first_harmonic(spectrum, quiet=True)

    xvalues = []
    if parity == 'all':
        xvalues = np.arange(first[1], spectrum[-1][0], first[1], dtype=float)
    elif parity == 'odd':
        xvalues = np.arange(first[1], spectrum[-1][0], 2*first[1], dtype=float)
    elif parity == 'even':
        xvalues = np.arange(2*first[1], spectrum[-1][0], 2*first[1], dtype=float)


    xranges=[]
    for x in xvalues:
        start = x - xwidth
        stop = x + xwidth
        if start < spectrum[0][0]:
            start = spectrum[0][0]
        if stop > spectrum[-1][0]:
            stop = spectrum[-1][0]
        xranges.append([start, stop])

    #fits = fit_spectrum_gaussian(spectrum, xranges=xranges, figsize=figsize, quiet=quiet, show=show, ofile=ofile)
    fits = find_peaks_parabola(spectrum, xranges)

    make_plot = show or ofile

    if make_plot:
        plt = oscars.plots_mpl.plot_spectrum(spectrum, figsize=figsize, show=False, ret=True)
        for fit in fits:
            plt.axvline(x=fit[1], linestyle='dashed')
        plt.show()

        


    return fits





def find_odd_harmonics (spectrum, first=None, xwidth=50, figsize=None, quiet=True, ofile=None, show=True):

    return find_harmonics(spectrum=spectrum, first=first, xwidth=xwidth, parity='odd', figsize=figsize, quiet=quiet, ofile=ofile, show=show)


def find_even_harmonics (spectrum, first=None, xwidth=50, figsize=None, quiet=True, ofile=None, show=True):

    return find_harmonics(spectrum=spectrum, first=first, xwidth=xwidth, parity='even', figsize=figsize, quiet=quiet, ofile=ofile, show=show)

def find_all_harmonics (spectrum, first=None, xwidth=50, figsize=None, quiet=True, ofile=None, show=True):

    return find_harmonics(spectrum=spectrum, first=first, xwidth=xwidth, parity='all', figsize=figsize, quiet=quiet, ofile=ofile, show=show)










def find_first_harmonic (spectrum, last_fit=None, quiet=True):
    """Find the first harmonic of a spectrum"""
    
    if len(spectrum) == 0:
        return last_fit

    # Just try finding one!
    try:
        #this_fit = fit_spectrum_gaussian(spectrum, n=1, quiet=quiet, show=False)[0]
        returned_fits = find_peaks_parabola(spectrum)
        if len(returned_fits) == 1:
            this_fit = returned_fits[0]
        else:
            raise ValueError('did not get a fit returned')
    except ValueError:
        Y = [s[1] for s in spectrum]
        X = [s[0] for s in spectrum]
        ymax = max(Y)
        ind = Y.index(ymax)
        if last_fit is not None:
            return last_fit
        return [ymax, X[ind], 0]

    if last_fit is not None:
        if this_fit[0] < last_fit[0] * 0.10:
            return last_fit

    cutoff = this_fit[1] - 5*this_fit[2]
    new_spectrum = [s for s in spectrum if s[0] < cutoff]
    return find_first_harmonic(new_spectrum, this_fit)





def find_first_harmonic_old (spectrum, quiet=True):
    """Find the first harmonic of a spectrum"""
    
    # Just try finding one!
    try:
        fit_1 = fit_spectrum_gaussian(spectrum, n=1, quiet=quiet, show=False)[0]
    except ValueError:
        Y = [s[1] for s in spectrum]
        X = [s[0] for s in spectrum]
        ymax = max(Y)
        ind = Y.index(ymax)
        return [ymax, X[ind], 0]
 
    # I guess we found one good peak, so let's try dividing the energy and finding lower ones
    # until we can't reasonably do it anymore
    last_amplitude = fit_1[0]
    last_energy = fit_1[1]
    last_sigma = fit_1[2]

    is_fit_bad = False
    last_fit_bad = False
    myfit = []

    for ifit in range(1, 51):
        try:
            myfit = fit_spectrum_gaussian(spectrum, xranges=[[fit_1[1] * 0.9 / ifit, fit_1[1] * 1.1 / ifit]], quiet=quiet, show=False)
        except:
            myfit = []
            if_fit_bad = True

        if not quiet:
            print('myfit for ifit=', ifit, myfit)

        if is_fit_bad is False and len(myfit) != 0 and myfit[0][0] > last_amplitude * 0.01 and myfit[0][2] < last_sigma * 1.5:
            last_amplitude = myfit[0][0]
            last_energy = myfit[0][1]
            last_sigma = myfit[0][2]
            last_fit_bad = False
        else:
            if last_fit_bad:
                return [last_amplitude, last_energy, last_sigma]

            last_fit_bad = True

    raise RuntimeError('Could not find first harmonic')

    return





def find_peaks_parabola (spectrum, xranges=None):
    """
    Fit the highest point and it's left and right neighbors with a parabola
    to estimate the actual peak.
    
    Parameters
    ----------
    spectrum : list of lists
        oscars.sr spectrum-like object: [[en0, flux0], [en1, flux1], ...]

    xranges : list of lists
        ranges to search in, eg: [[90, 110], [290, 310], ...].  Can be None.

    Returns
    -------
    peaks - a list of lists, each list as such: [flux, energy, 0] (where the
    last is typically sigma in gaussian fits)
    """
    harmonics = []
    if xranges is None:
        xranges = [[spectrum[0][0], spectrum[-1][0]]]
    for i in range(len(xranges)):
        xr = xranges[i]
        XP = [s[0] for s in spectrum if s[0] >= xr[0] and s[0] <= xr[1]]
        YP = [s[1] for s in spectrum if s[0] >= xr[0] and s[0] <= xr[1]]
    
        ymax = max(YP)
        ind = YP.index(ymax)
        if ind == 0 or ind == len(YP)-1:
            #raise IndexError('Maximum Y found at edge of range:', xranges[i])
            continue
        x1 = XP[ind-1]
        x2 = XP[ind]
        x3 = XP[ind+1]
        y1 = YP[ind-1]
        y2 = YP[ind]
        y3 = YP[ind+1]
        
        denom = (x1 - x2)*(x1 - x3)*(x2 - x3)
        A = (x3 * (y2 - y1) + x2 * (y1 - y3) + x1 * (y3 - y2)) / denom
        B = (x3*x3 * (y1 - y2) + x2*x2 * (y3 - y1) + x1*x1 * (y2 - y3)) / denom
        C = (x2 * x3 * (x2 - x3) * y1 + x3 * x1 * (x3 - x1) * y2 + x1 * x2 * (x1 - x2) * y3) / denom

        xx = -B / (2 * A)
        yy = A * xx * xx + B * xx + C
        

        x_up = -1
        x_down = -1
        for ip in range(ind, len(XP)):
            if YP[ip] < ymax / 2 and x_up < 0:
                x_up = (ymax / 2 - YP[ip-1]) * (XP[ip] - XP[ip-1]) / (YP[ip] - YP[ip-1]) + XP[ip-1]

        for ip in range(ind, 0, -1):
            if YP[ip] < ymax / 2 and x_down < 0:
                x_down = (ymax / 2 - YP[ip]) * (XP[ip+1] - XP[ip]) / (YP[ip+1] - YP[ip]) + XP[ip]

        fwhm = 0 if x_up < 0 or x_down < 0 else x_up - x_down

        harmonics.append([yy, xx, fwhm])
        
    return harmonics




def fit_spectrum_parabolic_gaussian (spectrum, xranges=[], n=None, figsize=None, quiet=True, show=True, ofile=None):

    # If no number is specified and xranges is empty use the whole range
    if n is None and len(xranges) == 0:
        xranges=[]
        xranges.append([spectrum[0][0], spectrum[-1][0]])

    # Sigma to remove from spectrum when peak is found, initial width guess
    nsigma_rm = 10.
    sigma_guess = 10

    # Not quiet
    nq = not quiet


    # X and Y data from spectrum
    X = [s[0] for s in spectrum]
    Y = [s[1] for s in spectrum]

    # Gaussian function
    def func(x, a, x0, sigma):
        return a*np.exp(-(x-x0)**2/(2*sigma**2))

    # Make Plot or not
    mp = show or (ofile is not None)

    # If make plot configure it
    if mp:
        plt.figure(1, figsize=figsize)
        plt.plot(X, Y, marker='.')
        plt.ylim(0, plt.ylim()[1])
        plt.xlim(X[0], X[-1])

    # Fit results we will return
    fit_results = []

    # If n is not specified use the xranges ranges ranges ranges (just kidding about the last two)
    if n is None:
        # Loop over the ranges in xranges
        for xr in xranges:

            # X and Y points in this range
            XP = [s[0] for s in spectrum if s[0] >= xr[0] and s[0] <= xr[1]]
            YP = [s[1] for s in spectrum if s[0] >= xr[0] and s[0] <= xr[1]]

            # Guess that the peak is at the maximum point
            amplitude_guess = max(YP)
            x_guess = XP[YP.index(amplitude_guess)]

            # Find peak using parabola
            #fit = find_peaks_parabola(spectrum, 
            try:
                popt, pcov = curve_fit(func, XP, YP, p0=[amplitude_guess, x_guess, sigma_guess])
                if popt[1] >= xr[0] and popt[1] <= xr[1]:
                    fit_results.append(list(popt))
            except RuntimeError:
                pass

    else:
        pass







    return






def fit_spectrum_gaussian (spectrum, xranges=[], n=None, figsize=None, quiet=True, show=True, ofile=None):
    """Fit multiple gaussians in ranges given to spectrum

    """

    if n is None and len(xranges) == 0:
        xranges=[]
        xranges.append([spectrum[0][0], spectrum[-1][0]])

    nsigma_rm = 10.
    sigma_guess = 10

    nq = not quiet

    # X and Y data from spectrum
    X = [s[0] for s in spectrum]
    Y = [s[1] for s in spectrum]

    # Gaussian function
    def func(x, a, x0, sigma):
        return a*np.exp(-(x-x0)**2/(2*sigma**2))

    mp = show or (ofile is not None)

    if mp: plt.figure(1, figsize=figsize)
    if mp: plt.plot(X, Y, marker='.')
    if mp: plt.ylim(0, plt.ylim()[1])
    if mp: plt.xlim(X[0], X[-1])


    fit_results = []

    if n is None:
        for i in range(len(xranges)):
            xr = xranges[i]
            XP = [s[0] for s in spectrum if s[0] >= xr[0] and s[0] <= xr[1]]
            YP = [s[1] for s in spectrum if s[0] >= xr[0] and s[0] <= xr[1]]


            amplitude_guess = max(YP)
            x_guess = XP[YP.index(amplitude_guess)]

            try:
                popt, pcov = curve_fit(func, XP, YP, p0=[amplitude_guess, x_guess, sigma_guess])
                if popt[1] >= xr[0] and popt[1] <= xr[1]:
                    fit_results.append(list(popt))
            except RuntimeError:
                pass

    else:
        XP = [s[0] for s in spectrum]
        YP = [s[1] for s in spectrum]

        for i in range(n):
            amplitude_guess = max(YP)
            x_guess = XP[YP.index(amplitude_guess)]

            try:
                popt, pcov = curve_fit(func, XP, YP, p0=[amplitude_guess, x_guess, sigma_guess])
                fit_results.append(list(popt))
            except RuntimeError:
                popt = [amplitude_guess, x_guess, 10]
                fit_results.append(list(popt))

            # Remove these points from X and Y
            XNEW = []
            YNEW = []
            for j in range(len(XP)):
                if XP[j] <= popt[1] - nsigma_rm * popt[2] or XP[j] >= popt[1] + nsigma_rm * popt[2]:
                    XNEW.append(XP[j])
                    YNEW.append(YP[j])
            XP = XNEW
            YP = YNEW

    if n is not None:
        fit_results.sort(key=lambda x: x[1])

    for i in range(len(fit_results)):
        popt = fit_results[i]
        width = (xranges[i][1] - xranges[i][0])/2. if n is None else nsigma_rm * popt[2]
        nppx = 20. / popt[2]
        npx = int(nppx * 2. * width)
        xp = np.linspace(popt[1] - width, popt[1] + width, npx)
        ym = func(xp, popt[0], popt[1], popt[2])

        label = str(round(popt[1], 1)) + ' eV, $\\sigma = $' + str(round(popt[2], 2))
        if mp:
            p = plt.plot(xp, ym, linestyle='dashed', label=label)
            y0 = abs(plt.ylim()[0]) / (plt.ylim()[1] - plt.ylim()[0])
            plt.axvline(x=popt[1], ymin=y0, ymax=y0 + popt[0] / (plt.ylim()[1] - plt.ylim()[0]), linestyle='dashed', color=p[0].get_color())

    if mp: plt.legend()
    if ofile is not None:
        plt.savefig(ofile, bbox_inches='tight')
    if show: plt.show()
    
    return fit_results







def correct_trajectory(osr, position=[0, 0, 1], beta=[0, 0, 1], bfields=[], tol=1e-18):
    """
    Correct the trajectory using bfield kicks
    Parameters
    ----------
    osr : oscars.sr object
        oscars.sr object where you wish the trajectory to be corrected
        
    position : [float, float, float]
        Position to minimize distance of trajectory to
        
    beta : [float, float, float]
        The desired beta value (velocity / speed of light) at 'position'
        
    bfields : list
        List of bfield peak values, widths, positions, and names
        [[[BxMax, ByMax, BzMax], [SigmaX, SigmaY, SigmaZ], [x, y, z], 'name'], [...]]
        
    Examples
    --------
    correct_trajectory(osr, position=[0, 0, 0.8], beta=[0, 0, 1],
                   bfields=[[[0, 0.1, 0], [0, 0, 0.05], [0, 0, -0.7], 'kick_entry_x'],
                            [[0.1, 0, 0], [0, 0, 0.05], [0, 0, -0.7], 'kick_entry_y'],
                            [[0, 0.1, 0], [0, 0, 0.05], [0, 0, +0.7], 'kick_exit_x'],
                            [[0.1, 0, 0], [0, 0, 0.05], [0, 0, +0.7], 'kick_exit_y']
                           ])
    """
    
    # Number of available fields
    nfields = len(bfields)
    if nfields < 1:
        raise ValueError('number of bfields must be > 0')

    # Ranges
    bounds = [[-1, 1]] * nfields
    
    # Dimension and initiol guess
    x0 = np.array([0] * nfields)
    
    # Remove all fields with same name
    for bf in bfields:
        osr.remove_bfield(bf[3])

    # Function to add and remove kick for testing
    def bfield_kicks(x):
        
        for ib in range(len(bfields)):
            osr.add_bfield_gaussian(bfield=[bf * x[ib] for bf in bfields[ib][0]],
                                    sigma=bfields[ib][1],
                                    translation=bfields[ib][2],
                                    name=bfields[ib][3])


        # Set a new particle and calculate trajectory.  After remove kick
        osr.set_new_particle(particle='ideal')
        trj = osr.calculate_trajectory()

        for ib in range(len(bfields)):
            osr.remove_bfield(bfields[ib][3])
            
        ms = 999999
        mb = [0, 0, 0]
        for tpt in trj:
            tp = tpt[1]
            tb = tpt[2]
            
            sq = (tp[0] - position[0])**2 + (tp[1] - position[1])**2 + (tp[2] - position[2])**2
            if sq < ms:
                ms = sq
                mp = tp
                mb = tb

        # Weighted sum of position and beta offset (Beta_z is close enough to 1)
        return (  ms + (mb[0] - beta[0])**2 + (mb[1] - beta[1])**2 + (mb[2] - beta[2])**2 )

    #res = minimize(bfield_kicks, x0, method='nelder-mead', options={'xtol': 1e-8}, bounds=[[1, -1], [1, -1]])
    res = minimize(bfield_kicks, x0, tol=tol, bounds=bounds)

    print('Minimization bfield factors:', res.x)

    for ib in range(len(bfields)):
        osr.add_bfield_gaussian(bfield=[bf * res.x[ib] for bf in bfields[ib][0]],
                                sigma=bfields[ib][1],
                                translation=bfields[ib][2],
                                name=bfields[ib][3])

    return
