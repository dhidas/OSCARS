from scipy.optimize import curve_fit, minimize
import numpy as np
import warnings

try:
    import matplotlib.pyplot as plt
except ImportError:
    warnings.warn('matplotlib cannot be imported')



def fit_spectrum_gaussian (spectrum, xranges=[], n=None, figsize=None, quiet=True):
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


    if nq: plt.figure(1, figsize=figsize)
    if nq: plt.plot(X, Y, marker='.')
    if nq: plt.ylim(0, plt.ylim()[1])
    if nq: plt.xlim(X[0], X[-1])


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
                #popt = [amplitude_guess, x_guess, 10]
                #max_s = 0
                #max_x = 0
                #for ix in range(len(XP)):
                #    if XP[i] >= xr[0] and XP[i] <= xr[1] and YP[i] > max_s:
                #        max_s = YP[i]
                #        max_x = XP[i]
                #if max_s != 0:
                #    fit_results.append([max_s, max_x, 10])

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
        if nq:
            p = plt.plot(xp, ym, linestyle='dashed', label=label)
            y0 = abs(plt.ylim()[0]) / (plt.ylim()[1] - plt.ylim()[0])
            plt.axvline(x=popt[1], ymin=y0, ymax=y0 + popt[0] / (plt.ylim()[1] - plt.ylim()[0]), linestyle='dashed', color=p[0].get_color())

    if nq: plt.legend()
    if nq: plt.show()
    
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
