Introduction
============

OSCARS (Open Source Code for Advanced Radiation Simulation) is a general purpose code to compute radiation of charged particles moving in magnetic and electric fields.  Other then numerical discritisation, no approximations are made and the calculations herein are valid in both the non-relativistic and utra-relativistic regimes.  Radiation spectra and flux densities are calculated in the so-called near field.

OSCARS is written in modern c++ for speed.  The main interface to OSCARS is a simple but powerful python interface (technically written as a c extension).  OSCARS works with both python 2 and 3.


Advantages of OSCARS
--------------------

* 100% Open source
* Valid in the *near-field* regime
* Valid for ultra-relativistic and non-relativistic particles/beams
* Import arbitrary field data in several formats
* Use built-in fields
* Easily specify any field as a python function
* Power density on arbitrary shaped surfaces in 3D
* Easily mix different particle beam types


Practical Applications
----------------------

* Calculate radiation from moving charged particles (electrons, protons, pions, muons, etc)
* Calculate power density on sensative equiptment for synchrotrons, linear accelerators, and particle and photon beamlines
