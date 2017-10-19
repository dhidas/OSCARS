#!/usr/bin/env python

import sys

import oscars.sr
from oscars.plots_mpl import plot_flux

osr = oscars.sr.sr()

plot_flux( osr.average_flux(ifiles=[sys.argv[1]]))

