from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import numpy as np
def fit_spectrum_gaussian (spectrum, xranges=[], n=None, figsize=None, quiet=False):
    """Fit multiple gaussians in ranges given to spectrum

    asdd
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
            XP = [s[0] for s in spectrum if s[0] >= xranges[i][0] and s[0] <= xranges[i][1]]
            YP = [s[1] for s in spectrum if s[0] >= xranges[i][0] and s[0] <= xranges[i][1]]


            amplitude_guess = max(YP)
            x_guess = XP[YP.index(amplitude_guess)]

            try:
                popt, pcov = curve_fit(func, XP, YP, p0=[amplitude_guess, x_guess, sigma_guess])
                fit_results.append(list(popt))
            except RuntimeError:
                popt = [amplitude_guess, x_guess, 10]
                fit_results.append(list(popt))

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
