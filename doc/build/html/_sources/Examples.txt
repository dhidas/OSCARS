Examples
========

Only the most basic examples are given here.  An extensive set of examples can be found on the OSCARS website at http://oscars.bnl.gov/examples.php

Please note that there is no preferred direction in OSCARS.  Any object (points, surfaces, beams, etc) can easily be rotated and translated in space.

Common Keyword Arguments
^^^^^^^^^^^^^^^^^^^^^^^^

An attempt has been made to use common keyword arguments wherever possible to simplify the interaction with OSCARS.  A list of common keyword arguments to functions is given below.

.. csv-table:: Common Keywords
   :header: "Keyword", "Format", "Description"
   :widths: 20, 50, 100

   "ifile", "string", "Full or relative path to the input file.  Text format."
   "ofile", "string", "Full or relative path to the output file.  Text format.  If outputting to a file the function will not return data (will return an empty list or None).  This is to allow you access to very large scale processing processing where memory may otherwise be an issue."
   "iformat", "string", "Format of the input text file.  e.g. 'iformat=[X Y Z Bx]'"
   "oformat", "string", "Format of the output text file."
   "rotations", "List :math:`[\theta_X, \theta_Y, \theta_Z]`", "A list of angles in [rad] used to perform a rotation on an object or field.  These rotations are performed in the order rotate about X, Y, and finally the Z axis."
   "translation", "List [X, Y, Z]", "An object is translated in space by the given amount in X, Y, and Z."
   "bfield", "List :math:`[B_X, B_Y, B_Z]`", "The direction and magnitude of a magnetic field given in [T].  This is the *peak* field if used in a distribution."
   "efield", "List :math:`[E_X, E_Y, E_Z]`", "The direction and magnitude of a electric field given in [V/m].  This is the *peak* field if used in a distribution."
   "name", "string", "The internal name which you give to an element.  Do not give the same name twice unless the previous object with this name has already been removed."
   "current", "float", "Beam current given in [A]."
   "mass", "float", "Mass given in [kg]"
   "charge", "float", "Charge given in [C]"
   "normal", "int", "Allowed values are -1, 0, 1.  -1 is for reversing the direction of the normal vertor of a surfafe internally to a calculation.  0 is for disregarding the sign of the scalar product of the surface with regard to the direction of radiation.  The default is always 0.  Setting this to 1 will allow you to see the direction with respect to the surface."
   "dim", "int", "Number of dimensions you want in the return data (either 2, or 3).  The default is 2 which gives you coordinates in the local surface defined by that surface.  dim=3 will return the 3D coordinates in the lab frame."


Example Calculations
^^^^^^^^^^^^^^^^^^^^

The example calculations include non-zero emittance beams, multi-particle simulation, multi-threaded operation, and use of the GPU.  They are in order of increasing complexity, although any of the calculations can be run in this way.

Spectrum Calculation
--------------------

Assume we have an NSLS2-like beam (3 GeV electron beam) and an undulator with a peak field of 1 [T], period of 49 [mm], with 31 periods.  We wish to calculate the on-axis spectrum 30 [m] downstream from the center of the undulator.

.. code-block:: py

   # Import the OSARS SR module
   import oscars.sr

   # Create an OSCARS SR object
   osr = oscars.sr.sr()

   # Add an undulator with peak field of 1 [T], period of 49 [mm], nperiods of 31
   osr.add_bfield_undulator(bfield=[0, 1, 0], period=[0, 0, 0.049], nperiods=31)

   # Add an electron beam with initial position X=Y=Z=0, and initially in the Z direction.
   osr.set_particle_beam(type='electron', name='beam_0', x0=[0, 0, -1], d0=[0, 0, 1], energy_GeV=3, current=0.500)

   # Set the initial and final times for trajectory calculation
   osr.set_ctstartstop(0, 2)

   # Calculate the spectrum
   spectrum = osr.calculate_spectrum(obs=[0, 0, 30], energy_range_eV=[100, 2000], npoints=500)






Power Density Calculation
-------------------------

* Multi-threaded

Assume we have an NSLS2-like beam (3 GeV electron beam) and an undulator with a peak field of 1 [T], period of 49 [mm], with 31 periods.  We wish to calculate the power density on two different surfaces, one flat surface downstream from the device which is almost perpendicular to the photon beam, and another which models the top inner surface of the beampipe passing through the undulator.  Here nthreads=8 is added to speed up the calculation.

.. code-block:: py

   # Import the OSARS SR module
   import oscars.sr

   # Create an OSCARS SR object
   osr = oscars.sr.sr()

   # Add an undulator with peak field of 1 [T], period of 49 [mm], nperiods of 31
   osr.add_bfield_undulator(bfield=[0, 1, 0], period=[0, 0, 0.049], nperiods=31)

   # Add an electron beam with initial position X=Y=Z=0, and initially in the Z direction.
   osr.set_particle_beam(type='electron', name='beam_0', x0=[0, 0, -1], d0=[0, 0, 1], energy_GeV=3, current=0.500)

   # Set the initial and final times for trajectory calculation
   osr.set_ctstartstop(0, 2)

   # Calculate the power density in a simple XY plane rotated 15 degrees at 30 [m] downstream.
   power_density_0 = osr.calculate_power_density_rectangle(plane='XY',
                                                           width=[0.06, 0.06],
                                                           npoints=[101, 101],
                                                           rotations=[0, 15. * osr.pi() / 180., 0],
                                                           translation=[0, 0, 30],
                                                           nthreads=8)

   # Calculate the power density on the top surface of beampipe inside of undulator
   power_density_1 = osr.calculate_power_density_rectangle(plane='XZ',
                                                           width=[0.02, 2.00],
                                                           npoints=[31, 201],
                                                           translation=[0, 0.005, 0],
                                                           nthreads=8)




Flux Density Calculation
------------------------

* Multi-particle
* Non-zero emittance
* Use of GPU

Assume we have an NSLS2-like beam (3 GeV electron beam) and an undulator with a peak field of 1 [T], period of 49 [mm], with 31 periods.  We wish to calculate the flux density on a plane 30 [m] downstream from the undulator.  Here we check for a GPU and if you have one, we use it.  The GPU takes priority, but if it doesn't exist the global nthreads will be used.

.. code-block:: py

   # Import the OSARS SR module
   import oscars.sr

   # Create an OSCARS SR object
   osr = oscars.sr.sr()

   # Set the global number of threads to use
   osr.set_nthreads_global(8)

   # If we have a GPU, use it!
   if osr.check_gpu() >= 1:
       osr.set_gpu_global(1)

   # Add an undulator with peak field of 1 [T], period of 49 [mm], nperiods of 31
   osr.add_bfield_undulator(bfield=[0, 1, 0], period=[0, 0, 0.049], nperiods=31)

   # Add an electron beam with non-zero emittance
   osr.add_particle_beam(type='electron',
                         name='beam_0',
                         energy_GeV=3,
                         x0=[0, 0, -1],
                         d0=[0, 0, 1],
                         current=0.500,
                         sigma_energy_GeV=0.001*3,
                         beta=[1.5, 0.8],
                         emittance=[0.9e-9, 0.008e-9],
                         horizontal_direction=[1, 0, 0],
                         lattice_reference=[0, 0, 0])

   # Set the initial and final times for trajectory calculation
   osr.set_ctstartstop(0, 2)

   # Calculate the flux in a simple XY plane for a multi-particle non-zero emittance beam
   flux = osr.calculate_flux_rectangle(plane='XY',
                                       energy_eV=455,
                                       width=[0.06, 0.06],
                                       npoints=[101, 101],
                                       translation=[0, 0, 30],
                                       nparticles=1000)

On Plotting Results
^^^^^^^^^^^^^^^^^^^

OSCARS comes with additional plotting tools which use matplotlib.  We are also considering implementing the ROOT framework for plotting.  Currently one can plot results from the above examples simply as follows.

.. code-block:: py

   # Import the OSCARS MLP tools
   from oscars.plots_mpl import *

   plot_spectrum(spectrum)
   plot_power_density(power_density_0)
   plot_power_density(power_density_1)
   plot_flux(flux)



