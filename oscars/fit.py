from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import numpy as np
import math
def fit_spectrum_quadratic (spectrum, n, erange=None):
    """Quadratic fit for spectral peaks.

    Keyword Arguments:
    spectrum - energy spectrum data as list of lists
    n - number of peaks to search for
    erange - range of energy values to search for peaks within
    """
    
    erange = ([spectrum[0][0], spectrum[-1][0]] if erange==None 
                    else erange)
    
    def parabola(x1, y1, x2, y2, x3, y3):
        denom = (x1 - x2) * (x1 - x3) * (x2 - x3)
        a = (x3 * (y2 - y1) + x2 * (y1 - y3) + x1 * (y3 - y2)) / denom
        b = (x3 * x3 * (y1 - y2) 
             + x2 * x2 * (y3 - y1) 
             + x1 * x1 * (y2 - y3)
            ) / denom
        c = (x2 * x3 * (x2 - x3) * y1
             + x3 * x1 * (x3 - x1) * y2 
             + x1 * x2 * (x1 - x2) * y3
            ) / denom
        inter = [(-b + math.sqrt(b**2 - 4*a*c))/(2*a), 
                 (-b - math.sqrt(b**2 - 4*a*c))/(2*a)]
        if inter[0] > inter[1]:
            intercepts = [inter[1], inter[0]]
        else:
            intercepts = inter
        vertex = [-b/(2*a), c-b**2/(4*a)]
        return [intercepts, vertex]

    # X and Y data from spectrum
    X = [s[0] for s in spectrum]
    Y = [s[1] for s in spectrum]

    fit_results = []
  
    # X and Y data within specified range of X (energy) values
    XP = [s[0] for s in spectrum if (s[0] > erange[0] and s[0] < erange[1])]
    YP = [s[1] for s in spectrum if (s[0] > erange[0] and s[0] < erange[1])]

    for i in range(n):
        amplitude_guess = max(YP)
        index_guess = YP.index(amplitude_guess)
        x_guess = XP[index_guess]
        amplitude_before = YP[index_guess - 1]
        x_before = XP[index_guess - 1]
        amplitude_after = YP[index_guess + 1]
        x_after = XP[index_guess + 1]

        [intercepts, vertex] = parabola(x_before, amplitude_before, x_guess, 
                                        amplitude_guess, x_after, amplitude_after)
        fit_results.append([intercepts, vertex])

        XNEW = []
        YNEW = []
        for j in range(len(XP)):
            if XP[j] <= intercepts[0] or XP[j] >= intercepts[1]:
                XNEW.append(XP[j])
                YNEW.append(YP[j])
        XP = XNEW
        YP = YNEW

    if n is not None:
        fit_results.sort(key=lambda x: x[1][0])
       
    for i in range(len(fit_results)):
        intercepts = fit_results[i][1]
        vertex = fit_results[i][0] 
        #print('#' + str(i) + ': intercepts: ' + str(intercepts) + ', vertex: ' + str(vertex))
       
        #popt = fit_results[i]
        #width = (xranges[i][1] - xranges[i][0])/2. if n is None else nsigma_rm * popt[2]
        #nppx = 20. / popt[2]
        #npx = int(nppx * 2. * width)
        #xp = np.linspace(popt[1] - width, popt[1] + width, npx)
        #ym = gaussian(xp, popt[0], popt[1], popt[2])
        #
        #label = str(round(popt[1], 1)) + ' eV, $\\sigma = $' + str(round(popt[2], 2))
        #if nq:
        #    p = plt.plot(xp, ym, linestyle='dashed', label=label)
        #    y0 = abs(plt.ylim()[0]) / (plt.ylim()[1] - plt.ylim()[0])
        #    plt.axvline(x=popt[1], ymin=y0, ymax=y0 + popt[0] / (plt.ylim()[1] - plt.ylim()[0]), 
        #                linestyle='dashed', color=p[0].get_color())

    #if nq: plt.legend()
    #if nq: plt.show()
    
    return fit_results

def show_fit(spectrum, fit_results):
    X = [s[0] for s in spectrum]
    Y = [s[1] for s in spectrum]
    
    plt.plot(X, Y)
    for peak in fit_results:
        plt.axvline(x=peak[1][0], ymin=0, ymax=1, linestyle='dashed')
    plt.show()

