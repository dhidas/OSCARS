from __future__ import print_function

import numpy as np
import matplotlib.pyplot as plt
import random

from math import pi, sin, cos, atan, sqrt, log

class twiss:
    """Class for twiss parameters and transformations with some basic monte carlo and plotting"""

    def __init__ (self, beta=None, alpha=None, gamma=None, emittance=1e-9):

        self.beta = beta
        self.alpha = alpha
        self.gamma = gamma
        self.emittance = emittance

        self.set(beta, alpha, gamma)
            
        return

    def set (self, beta=None, alpha=None, gamma=None, emittance=None):
        """Set the twiss parameters given 2 of the three parameters.  If one parameter given, alpha is assumed to be zero"""

        if emittance is not None:
            self.emittance = emittance

        if beta is None and gamma is None:
            raise ValueError('Must provide some arguments')

        if beta is not None and alpha is not None and gamma is not None:
            gamma_check = self.calculate_gamma(beta, alpha)
            if gamma_check != gamma:
                raise ValueError('given gamma does not match calculated gamma from beta and alpha.  These are not proper twiss parameters')

        if alpha is None and gamma is None:
            alpha = 0
            gamma = self.calculate_gamma(beta, alpha)
        elif alpha is None and beta is None:
            alpha = 0
            beta = self.calculate_beta(alpha, gamma)
        elif alpha is None:
            alpha = self.calculate_alpha(beta, gamma)
        elif beta is None:
            beta = self.calculate_beta(alpha, gamma)
        elif gamma is None:
            gamma = self.calculate_gamma(beta, alpha)

        self.beta = beta
        self.alpha = alpha
        self.gamma = gamma

        return

    def get (self):
        """Get beta alpha gamma"""

        return [self.beta, self.alpha, self.gamma]

    def calculate_beta (self, alpha, gamma):
        """Calculate beta given beta and alpha"""

        return (1 + alpha*alpha) / gamma

    def calculate_alpha (self, beta, gamma):
        """Calculate alpha given beta and alpha"""

        return sqrt(beta*gamma - 1)

    def calculate_gamma (self, beta, alpha):
        """Calculate gamma given beta and alpha"""

        return (1 + alpha*alpha) / beta


    def self_drift (self, L):
        """Drift for a length of L"""

        self.beta = self.beta - 2*L*self.alpha + L*L*self.gamma
        self.alpha = self.alpha - L*self.gamma

        return

    def get_drift (self, L):
        """Drift for a length of L"""

        beta = self.beta - 2*L*self.alpha + L*L*self.gamma
        alpha = self.alpha - L*self.gamma

        return [beta, alpha, self.gamma]


    def plot_drift (self, start=-1, stop=1):
        """Plot twiss parameters in a drift between start and stop"""

        XB=[]
        XA=[]
        XG=[]
        Z=[]
        for z in np.linspace(start, stop, 201, endpoint=True):
            bag = self.get_drift(z)
            XB.append(bag[0])
            XA.append(bag[1])
            XG.append(bag[2])
            Z.append(z)

        plt.figure()
        plt.plot(Z, XB, label='$\\beta$')
        plt.plot(Z, XA, label='$\\alpha$')
        plt.plot(Z, XG, label='$\\gamma$')
        plt.title('Twiss parameters in drift')
        plt.xlabel('s [m]')
        plt.legend()
        plt.show()
        plt.close()

    def get_ellipse_points (self, n=101):
        """Get n points that outline the ellipse for the internal twiss parameters in format [[x0, x1, ..], [y0, y1, ..]]"""

        p = [[], []]
        for theta in np.linspace(0, 2*pi, n):
            x = sqrt(self.emittance * self.beta) * cos(theta)
            y = -sqrt(self.emittance / self.beta) * (self.alpha * cos(theta) + sin(theta))
            p[0].append(x)
            p[1].append(y)

        return p



    def random (self, distribution='kv'):
        """Get a random point for x and x' from the distribution specified"""

        if distribution == 'gaussian':
            sigx = sqrt(self.emittance * self.beta)
            sigxp = sqrt(self.emittance * self.gamma)
            if self.alpha == 0:
                rho = 0
            rho = sqrt(1 - pow(self.emittance / (sigx * sigxp), 2))

            u = random.random()
            v = random.random()

            conv = 1 if self.alpha <= 0 else -1

            myx = sigx * sqrt(-2 * log(u)) * (sqrt(1 - rho*rho) * cos(2 * pi * v) + rho * sin(2 * pi * v))
            myy = conv * sigxp * sqrt(-2 * log(u)) * sin(2 * pi * v)

            return [myx, myy]


        elif distribution == 'kv':
            theta = random.random() * 2 * pi

            x = sqrt(self.emittance * self.beta) * cos(theta)
            y = -sqrt(self.emittance / self.beta) * (self.alpha * cos(theta) + sin(theta))
 
            r = sqrt(x*x + y*y)
 
            real_theta = atan(y/x)
            if x < 0 and y > 0:
                real_theta += pi
            elif x < 0 and y < 0:
                real_theta -= pi
            my_r = sqrt(abs(random.random())) * r
 
            myx = my_r * cos(real_theta)
            myy = my_r * sin(real_theta)
 
            return [myx, myy]

        else:
            raise ValueError('incorrect distribution type specified')

        return





