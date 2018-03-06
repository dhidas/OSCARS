#!/usr/bin/env python

import sys

import oscars.sr
from oscars.plots_mpl import plot_spectra

osr = oscars.sr.sr()


spectra = []
labels  = []

for i in range(1, len(sys.argv)):
    spectra.append(osr.average_spectra(ifiles=[sys.argv[i]]))
    labels.append(i)


plot_spectra(spectra, labels, show=True, ofile='Spectrum.pdf', figsize=[12, 3])

