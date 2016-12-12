Member Functions
================

Primer for the Examples
-----------------------

All examples here assume that you have the OSCARS SR module in your path and have an OSCARS SR object called osr, perhaps as in the following example

.. code-block:: py

   # Import the OSCARS SR module
   import oscars.sr

   # Creat an OSCARS SR object
   osr = oscars.sr.sr()



oscars.sr
---------


.. py:method:: oscars.sr.pi()

   Get the value of pi: 3.14159265358979323846

   :returns: float




.. py:method:: oscars.sr.qe()

   Get the value of elementary charge: +1.602176462e-19 [C]

   :returns: float




.. py:method:: oscars.sr.rand()

   Random float between 0 and 1

   :returns: float




.. py:method:: oscars.sr.norm()

   Random float from the normal distribution

   :returns: float




.. py:method:: oscars.sr.set_seed(seed)

   Set the internal random seed

   :param seed: Random seed
   :type  seed: int



.. py:method:: oscars.sr.set_gpu_global(gpu)

   If set to 1, OSCARS will try to use the GPU for all calculations

   :param gpu: 1 or 0
   :type  gpu: int



.. py:method:: oscars.sr.check_gpu()

   Will return the number of GPUs available, or a negative number for any other error, for instance if your distribution was not compiled with GPU support.

   :return: integer


.. py:method:: oscars.sr.set_nthreads_global(nthreads)

   Set the number of threads you wish to use for all calculations.  If the GPU is requested it will take precedence over multi-threading.







.. py:method:: oscars.sr.get_ctstart()

   Gets the start time for particle propogation.  The *time* is in units of [m].

   :returns: float





.. py:method:: oscars.sr.get_ctstop()

   Gets the stop time for particle propogation.  The *time* is in units of [m].

   :returns: float




.. py:method:: oscars.sr.set_ctstartstop(time_start, time_stop)

   Set the start and stop times for the trajectory calculation.  Start time must be less than or equal to the T0 defined in the initial conditions for the particle beam.  Similarly, stop time must be greater than that T0.  This also sets a default value for the number of points used in the trajectory calculation.

   :param time_start: Start time
   :type time_start: float
   :param time_stop: Stop time
   :type time_stop: float
   :returns: None





.. py:method:: set_npoints_trajectory(npoints)

   Sets the number of points to be used in the trajectory calculation.  Only call this function **after** set_ctstartstop().

   :param npoints: Number of points
   :type npoints: int
   :returns: None
   :example:

   .. code-block:: py

      # First set the start and stop times
      osr.set_ctstartstop(0, 2)

      # Use 12345 points in the trajectory calculations
      osr.set_npoints_trajectory(12345)





.. py:method:: oscars.sr.get_npoints_trajectory

   Gets the number of points to be used in the trajectory calculation

   :returns: int - Number of points





.. py:method:: oscars.sr.add_bfield_file(ifile, iformat, [rotation, translation, scaling])

   Add a magnetic field from a text file *ifile* according to the format *iformat*.
   
   Currently OSCARS accepts the following file formats for iformat
      * 'OSCARS'
      * 'OSCARS1D [plus format string]'
      * 'SPECTRA'
      * 'SRW'

   For the OSCARS1D format you must also include the order of the data columns, which would typically look something like: 'OSCARS1D Z Bx By Bz'.  You may use X or Y instead of Z and the order and number of B[xyz] does not matter.  This mode also accepts non-uniformly distributed data.

   Optionally you can rotate and translate the field in space.  You can use an input *scaling* list to scale the input (which must be in SI units of [m] for distances/positions and [T] for magnetic field values.

   The rotation is performed first in the order: :math:`\theta_x, \theta_y, \theta_z`

   *scaling* is a list of less than or equal length to the number of elements in *iformat* when OSCARS1D is selected.  This will scale the input values of the i-th column before any rotation or translation.  This is useful if your data file is not in [T] and [m]

   :param ifile: Name of input file
   :type  ifile: str
   :param iformat: Format of the input file
   :type  iformat: str
   :param rotation: [:math:`\theta_x, \theta_y, \theta_z`]
   :type  rotation: list[3]
   :param translation: Translation in space [x, y, z]
   :type  filename: list[3]
   :param scaling: Scaling of input parameters
   :type  filename: list
   :returns: None

   :example:

   .. code-block:: py

      # Add a magnetic field from a file where the columns are in the order Z, Bx, By, Bz where Z is in [m] and Bx, By, Bz are in [T].
      osr.add_bfield_file(ifile='file.txt', iformat='OSCARS1D Z Bx By Bz')







.. py:method:: oscars.sr.add_bfield_function(function)

   Adds a magnetic field in the form of a user defined python function.  The input for this function must be (x, y, z, t).

   :param function: Python function
   :type  function: func
   :returns: None

   :example:

   .. code-block:: py

      # Create a function in python and use it as a magnetic field in OSCARS

      def myfunc(x, y, z, t):
          "Do not forget to write a docstring"
          if (z > 0):
            return 1
          return 0

      osr.add_bfield_function(myfunc)






.. py:method:: oscars.sr.add_bfield_gaussian(bfield, sigma, [rotations, translation])

   Add a gaussian magnetic field in 3D with the peak field magnitude and direction given by *bfield*, centered about a point with a given sigma in each coordinate.  If any component of *sigma* is less than or equal to zero it is ignored (ie. spatial extent is infinite).


   The functional form for this is the following:
   :math:`\\exp(-(x - x_0)^2 / \\sigma_x^2)`

   :param bfield: A list representing the peak field [Bx, By, Bz]
   :type  bfield: list[3]
   :param center: A list representing the coordinates of the center point [X0, Y0, Z0]
   :type  center: list[3]
   :param sigma: A list representing the sigma in each coordinate [S0, S0, S0]
   :type  sigma: list[3]
   :returns: None

   :Example: Add a magnetic field of 0.5 [T] in the X-direction centered at X=0, Y=0, Z=0 [in meters], with a sigma of 0.1 [m] in the Z-direction

   .. code-block:: py

      # Will add a magnetic field of 1 Tesla in the Y-direction centered
      # at X=0, Y=0, Z=0 [in meters], with a sigma of 0.1 in the Z-direction
      osr.add_bfield_gaussian(bfield=[1, 0, 0], sigma=[0, 0, 0.10])

   :Example: Add a magnetic field of 1 Tesla in the Y-direction centered at X=1, Y=1, Z=1 [in meters], with a sigma of 0.05 in the Z-direction

   .. code-block:: py

      # Will add a magnetic field of 1 Tesla in the Y-direction centered
      # at X=1, Y=1, Z=1 [in meters], with a sigma of 0.05 in the Z-direction
      osr.add_bfield_gaussian(bfield=[0, 1, 0], sigma=[0, 0, 0.05], translation=[1, 1, 1])






.. py:method:: oscars.sr.add_bfield_uniform(bfield, [width, rotations, translation])

   Add a uniform magnetic field in a given range or for all space.  The *bfield* is given as a 3D vector representing the field magnitude and direction.  *width* is an optional parameters, if not present the field permeates all space.  If a component of the 3D list *width* is less than or equal to zero, that coordinate will be ignored when calculating the field.


   :param bfield: A list representing the magnetic field [Bx, By, Bz]
   :type  bfield: list[3]
   :param width: A list representing the spetial extent of the field [Wx, Wy, Wz]
   :type  width: list[3]
   :param rotations: A list representing the rotations of this object about the X, Y amd Z axis (in that order)
   :type  rotations: list[3]
   :param translation: A list representing the translation of the center of this object
   :type  translation: list[3]
   :returns: None

   :Example: Add a magnetic field of 0.0001 [T] in the X-direction for all space

   .. code-block:: py

      # Will add a magnetic field of 0.0001 [T] in the X-direction over all space
      osr.add_bfield_uniform(bfield=[0.0001, 0, 0])

   :Example: Add a magnetic field of 0.0005 [T] in the Y-direction with a width in the Z-direction of 1.5 [m] (the other directions are ignored) centered at X=0, Y=0, Z=0.75 [m].

   .. code-block:: py

      # Will add a magnetic field of 1 Tesla in the Y-direction centered
      # at X=0, Y=0, Z=0.75 [in meters], extending in Z a width of 1.5 [m]
      osr.add_bfield_uniform(bfield=[0, 1, 0], width=[0, 0, 1.5], translation=[0, 0, 0.75])






.. py:method:: oscars.sr.add_bfield_undulator(bfield, period, nperiods, [phase, rotations, translation])

   Adds an ideal sinusoidal undulator field with a given maximum b-field amplitude, period, and number of periods.  Optionally one can specify the phase offset (in [rad]), rotations and translation.

   :param bfield: A list representing the peak field [Bx, By, Bz] in [T]
   :type  bfield: list[3]
   :param period: Length of one period
   :type  period: float
   :param nperiods: Number of periods
   :type  nperiods: int
   :param phase: Phase offset in [rad].
   :type  phase: float
   :param rotations: A list representing the rotations of this object about the X, Y amd Z axis (in that order)
   :type  rotations: list[3]
   :param translation: A list representing the translation of the center of this object
   :type  translation: list[3]
   :returns: None

   :Example: Add the magnetic field for an undulator of length 3 [m] with a period of 0.050 [m] (60 periods) having a peak magnetic field in the Y direction of 1.1 [T]

   .. code-block:: py

      # Will add magnetic field 3 [m] undulaotor having 1.1 [T] peak field
      # in the Y-direction with a period of 0.050 [m]
      osr.add_bfield_undulator(bfield=[0, 1.1, 0], period=0.050, nperiods=60)




.. py:method:: oscars.sr.get_bfield(X)

   Get the 3D field at any point in space.  This is the sum of all fields added to this OSCARS object.

   :param X: A 3D list representing a point in space [x, y, z]
   :type  X: list[3]
   :returns: [float, float, float] - list representing the B-field [Bx, By, Bz]






.. py:method:: oscars.sr.clear_bfields()

   Remove all of the existing magnetic fields.

   :returns: None













.. py:method:: oscars.sr.add_efield_file(ifile, iformat, [rotation, translation, scaling])

   Add an electric field from a text file *ifile* according to the format *iformat*.
   
   Currently OSCARS accepts the following file formats for iformat
      * 'OSCARS'
      * 'OSCARS1D [plus format string]'
      * 'SPECTRA'
      * 'SRW'

   For the OSCARS1D format you must also include the order of the data columns, which would typically look something like: 'OSCARS1D Z Bx By Bz'.  You may use X or Y instead of Z and the order and number of B[xyz] does not matter.  This mode also accepts non-uniformly distributed data.

   Optionally you can rotate and translate the field in space.  You can use an input *scaling* list to scale the input (which must be in SI units of [m] for distances/positions and [T] for electric field values.

   The rotation is performed first in the order: :math:`\theta_x, \theta_y, \theta_z`

   *scaling* is a list of less than or equal length to the number of elements in *iformat* when OSCARS1D is selected.  This will scale the input values of the i-th column before any rotation or translation.  This is useful if your data file is not in [T] and [m]

   :param ifile: Name of input file
   :type  ifile: str
   :param iformat: Format of the input file
   :type  iformat: str
   :param rotation: [:math:`\theta_x, \theta_y, \theta_z`]
   :type  rotation: list[3]
   :param translation: Translation in space [x, y, z]
   :type  filename: list[3]
   :param scaling: Scaling of input parameters
   :type  filename: list
   :returns: None

   :example:

   .. code-block:: py

      # Add electric field from a file where the columns are in the order Z, Bx, By, Bz where Z is in [m] and Bx, By, Bz are in [T].
      osr.add_efield_file(ifile='file.txt', iformat='OSCARS1D Z Bx By Bz')







.. py:method:: oscars.sr.add_efield_function(function)

   Adds electric field in the form of a user defined python function.  The input for this function must be (x, y, z, t).

   :param function: Python function
   :type  function: func
   :returns: None

   :example:

   .. code-block:: py

      # Create a function in python and use it as a electric field in OSCARS

      def myfunc(x, y, z, t):
          "Do not forget to write a docstring"
          if (z > 0):
            return 1
          return 0

      osr.add_efield_function(myfunc)






.. py:method:: oscars.sr.add_efield_gaussian(efield, sigma, [rotations, translation])

   Add a gaussian electric field in 3D with the peak field magnitude and direction given by *efield*, centered about a point with a given sigma in each coordinate.  If any component of *sigma* is less than or equal to zero it is ignored (ie. spatial extent is infinite).


   The functional form for this is the following:
   :math:`\\exp(-(x - x_0)^2 / \\sigma_x^2)`

   :param efield: A list representing the peak field [Bx, By, Bz]
   :type  efield: list[3]
   :param center: A list representing the coordinates of the center point [X0, Y0, Z0]
   :type  center: list[3]
   :param sigma: A list representing the sigma in each coordinate [S0, S0, S0]
   :type  sigma: list[3]
   :returns: None

   :Example: Add a electric field of 0.5 [T] in the X-direction centered at X=0, Y=0, Z=0 [in meters], with a sigma of 0.1 [m] in the Z-direction

   .. code-block:: py

      # Will add a electric field of 1 Tesla in the Y-direction centered
      # at X=0, Y=0, Z=0 [in meters], with a sigma of 0.1 in the Z-direction
      osr.add_efield_gaussian(efield=[1, 0, 0], sigma=[0, 0, 0.10])

   :Example: Add a electric field of 1 Tesla in the Y-direction centered at X=1, Y=1, Z=1 [in meters], with a sigma of 0.05 in the Z-direction

   .. code-block:: py

      # Will add a electric field of 1 Tesla in the Y-direction centered
      # at X=1, Y=1, Z=1 [in meters], with a sigma of 0.05 in the Z-direction
      osr.add_efield_gaussian(efield=[0, 1, 0], sigma=[0, 0, 0.05], translation=[1, 1, 1])






.. py:method:: oscars.sr.add_efield_uniform(efield, [width, rotations, translation])

   Add a uniform electric field in a given range or for all space.  The *efield* is given as a 3D vector representing the field magnitude and direction.  *width* is an optional parameters, if not present the field permeates all space.  If a component of the 3D list *width* is less than or equal to zero, that coordinate will be ignored when calculating the field.


   :param efield: A list representing the electric field [Bx, By, Bz]
   :type  efield: list[3]
   :param width: A list representing the spetial extent of the field [Wx, Wy, Wz]
   :type  width: list[3]
   :param rotations: A list representing the rotations of this object about the X, Y amd Z axis (in that order)
   :type  rotations: list[3]
   :param translation: A list representing the translation of the center of this object
   :type  translation: list[3]
   :returns: None

   :Example: Add a electric field of 0.0001 [T] in the X-direction for all space

   .. code-block:: py

      # Will add a electric field of 0.0001 [T] in the X-direction over all space
      osr.add_efield_uniform(efield=[0.0001, 0, 0])

   :Example: Add a electric field of 0.0005 [T] in the Y-direction with a width in the Z-direction of 1.5 [m] (the other directions are ignored) centered at X=0, Y=0, Z=0.75 [m].

   .. code-block:: py

      # Will add a electric field of 1 Tesla in the Y-direction centered
      # at X=0, Y=0, Z=0.75 [in meters], extending in Z a width of 1.5 [m]
      osr.add_efield_uniform(efield=[0, 1, 0], width=[0, 0, 1.5], translation=[0, 0, 0.75])






.. py:method:: oscars.sr.add_efield_undulator(efield, period, nperiods, [phase, rotations, translation])

   Adds an ideal sinusoidal undulator field with a given maximum b-field amplitude, period, and number of periods.  Optionally one can specify the phase offset (in [rad]), rotations and translation.

   :param efield: A list representing the peak field [Bx, By, Bz] in [T]
   :type  efield: list[3]
   :param period: Length of one period
   :type  period: float
   :param nperiods: Number of periods
   :type  nperiods: int
   :param phase: Phase offset in [rad].
   :type  phase: float
   :param rotations: A list representing the rotations of this object about the X, Y amd Z axis (in that order)
   :type  rotations: list[3]
   :param translation: A list representing the translation of the center of this object
   :type  translation: list[3]
   :returns: None

   :Example: Add the electric field for an undulator of length 3 [m] with a period of 0.050 [m] (60 periods) having a peak electric field in the Y direction of 1.1 [T]

   .. code-block:: py

      # Will add electric field 3 [m] undulaotor having 1.1 [T] peak field
      # in the Y-direction with a period of 0.050 [m]
      osr.add_efield_undulator(efield=[0, 1.1, 0], period=0.050, nperiods=60)




.. py:method:: oscars.sr.get_efield(X)

   Get the 3D field at any point in space.  This is the sum of all fields added to this OSCARS object.

   :param X: A 3D list representing a point in space [x, y, z]
   :type  X: list[3]
   :returns: [float, float, float] - list representing the B-field [Bx, By, Bz]






.. py:method:: oscars.sr.clear_efields()

   Remove all of the existing electric fields.

   :returns: None








.. py:method:: oscars.sr.write_bfield(ofile, oformat [, xlim, nx, ylim, ny, zlim, nz, comment])

   Write the magnetic field to ofile in the format described in oformat.  You must specify at least one limits and one number of points (e.g. zlim and nz).  All output is in SI units.

   The formats (*oformat*) available are:
   * OSCARS
   * OSCARS1D
   * SPECTRA
   * SRW

   For OSCARS1D you also must specify what you want in the output.  One spatial dimension must be specified along with at least one field dimension (in any order you like).  For examples theses are all valid: 'OSCARS1D Z Bx By Bz', 'OSCARS1D By Bx Z Bz'.

   :param ofile: Name of output file
   :type  ofile: str
   :param oformat: format of output file
   :type  oformat: str
   :param xlim: min and max for x dimension
   :type  xlim: list[min, max]
   :param ylim: min and max for y dimension
   :type  ylim: list[min, max]
   :param zlim: min and max for z dimension
   :type  zlim: list[min, max]
   :param comment: comment string to be added to file header.  LF and CR are removed.
   :type  comment: str
   :returns: None








.. py:method:: oscars.sr.write_efield(ofile, oformat [, xlim, nx, ylim, ny, zlim, nz, comment])

   Same as write_bfield, but writes the electric field.























.. py:method:: oscars.sr.set_particle_beam(type, name, energy_GeV, d0, x0, [sigma_energy_GeV, t0, current, weight, rotations, translation, horizontal_direction, beta, emittance, lattice_center, mass, charge])

   This function is the same as oscars.sr.add_particle_beam(), but it clears all particle beams before the 'add'.




.. py:method:: oscars.sr.add_particle_beam(type, name, energy_GeV, x0, d0, [sigma_energy_GeV, t0, current, weight, rotations, translation, horizontal_direction, beta, emittance, lattice_reference, mass, charge])

   Add a particle beam to the OSCARS object with a name given by *name*.  There is no limit to the number of different particle beams one can add.  They are added with a *weight* which is by default 1.  The weight is used in random sampling when asking for a new particle, for example in oscars.sr.set_new_particle().

   Supported particle types for *type* are:
      * electron
      * positron
      * muon
      * anti-muon
      * proton
      * anti-proton
      * pi+
      * pi-

   :param type: one of the built-in types of particle beams, or 'custom'.  If you use custom you must also specify *mass* and *charge*.
   :type  type: str
   :param name: User identified of this beam
   :type  name: str
   :param energy_GeV: Beam energy in [GeV]
   :type  energy_GeV: float
   :param d0: Vector representing the default direction of the beam at the initial position.  The normalization of this vector does not matter.
   :type  d0: [float, float, float]
   :param x0: Coordinates of the initial position in [m]
   :type  x0: [float, float, float]
   :param sigma_energy_GeV: Beam energy in [GeV]
   :type  sigma_energy_GeV: float
   :param t0: Initial time in [m] at the initial_position.  Time here is in terms of ct.
   :type  t0: float
   :param current: Beam current in [A].  If this parameter is 0, the current defaults to the single particle charge
   :type  current: float
   :param weight: Weight to give this beam for random sampling when there are multiple beams
   :type  weight: float
   :param rotations: A list representing the rotations of this beam about the X, Y amd Z axis (in that order)
   :type  rotations: list[3]
   :param translation: A list representing the translation of the x0 of this object
   :type  translation: list[3]
   :param horizontal_direction: A list representing the *horizontal* beam direction.  This must be perpendicular to the beam *direction*.  The vertical direction is defined internally as :math:`[direction] \times [horizontal\_direction]`.
   :type  horizontal_direction: list[3]
   :param beta: values of the horizontal and vertical beta funtion at the point *lattice_center* [beta_x, beta_y]
   :type  beat: [float, float]
   :param emittance: values of the horizontal and vertical emittance [emittance_x, emittance_y]
   :type  emittance: [float, float]
   :param lattice_reference: Coordinates of the lattice center (must be on-axis with respect to the beam)
   :type  lattice_reference: list[3]
   :param mass: mass of a *custom* particle.  This is only used if *type* = 'custom'.  Must be non-zero.
   :type  mass: float
   :param charge: Charge of a *custom* particle.  This is only used if *type* = 'custom'.  Must be non-zero.
   :type  charge: float
   :returns: None

   :Example: Add an electron beam with 0.500 [A] current at an initial position of [0, 0, 0] in the Y direction with energy of 3 [GeV]

   .. code-block:: py

      # Will add an electron beam called beam_0 starting at [0, 0, 0] headed in the
      # direction of +Y [0, 1, 0], with an energy of 3 [GeV] at time T0 = 0 [m] having a
      # current of 0.500 [A] and a weight 1
      osr.add_particle_beam(type='electron', name='beam_0', x0=[0, 0, 0], d0=[0, 1, 0], energy_GeV=3, current=0.500)

   :Example: Add a positron beam with 0.500 [A] current at an initial position of [-2, 0, 0] in the direction given by theta in the X-Y plane with energy of 3 [GeV]

   .. code-block:: py

      # Create a positron beam in the X-Y plane at an angle theta
      from math import sin, cos
      theta = 0.25 * osr.pi()
      osr.add_particle_beam(type='positron', name='beam_0', x0=[-2, 0, 0], d0=[sin(theta), cos(theta), 0], energy_GeV=3, current=0.500)






.. py:method:: oscars.sr.clear_particle_beams()

   Remove all of the existing particle beams

   :returns: None






.. py:method:: oscars.sr.set_new_particle(,[beam, particle])

   If no arguments are given sets the current internal particle to a random new particle.  The randomization is based on the weights given for each beam.  This also sets the initial conditions for the particle used in trajectory calculations based on the beam parameters within the randomly sepected beam.  You can specify which beam you want a random particle from using the *beam* parameter.  The *particle* parameter can be 'ideal' if you want the ideal initial conditions for a particle without randomization.

   :param beam: The name of the beam from which to get a particle from
   :type  beam: str
   :param particle: 'ideal' or 'random', for random, you may omit this.
   :type  particle: str
   :returns: None




.. py:method:: oscars.sr.get_particle_x0()

   Get the initial position for the current particle in [m]

   :returns: [float, float, float]






.. py:method:: oscars.sr.get_particle_beta0()

   Get the initial :math:`\vec \beta` for the current particle

   :returns: [float, float, float]






.. py:method:: oscars.sr.get_particle_e0()

   Get the initial energy for the current particle in [GeV]

   :returns: float






.. py:method:: oscars.sr.calculate_trajectory()

   Calculates the trajectory for the current internal particle.  This calculates the trajectory in 3D from the time set by oscars.sr.set_ctstart() to the time set by oscars.sr.set_ctstop() beginning at the *t0* given by the particle beam from which this particle comes from.  It first does a forward propogation to the stop time, then a backward propogation to the start time.

   It is not necessary to call this method before other calculations such as spectrum, power density, or flux calculation methods.

   If you have a current particle loaded using *SetNewParticle* this method will calculate the trajectory for that particle.  If no particle is defined one will be randomly selected based on the beam weights and beam parameters.

   :returns: A list of points of the form [[[x, y, z], [Beta_x, Beta_y, Beta_z]], ...]




.. py:method:: oscars.sr.get_trajectory()

   Get the current trajectory.  If the trajectory has not been calculated this will return an empty list.  The format of the returned list consists of a list of lists giving you the position and beta (v/c) of the particle at each position.  For a trajectory returnned is of the form: [[[x, y, z], [Beta_x, Beta_y, Beta_z]], ...]

   :returns: A list of points of the form [[[x, y, z], [Beta_x, Beta_y, Beta_z]], ...]





.. py:method:: oscars.sr.calculate_spectrum(obs, [[npoints, energy_range_eV], points_eV], [nparticles, ofile, nthreads, gpu])

   Calculate the spectrum given a point in space, the range in energy, and the number of points.  The calculation uses the current particle and its initial conditions.  If the trajectory has not been calculated it is calculated first.  The units of this calculation are [:math:`photons / mm^2 / 0.1% bw / s`]

   Previously to calling this function you must define a particle beam and define start and stop times at very minimum.

   You **must** provide either (*npoints* and *energy_range_eV*) or *points_eV*.

   :param obs: Point where you wish to calculate the spectrum
   :type  obs: [float, float, float]
   :param npoints: Number of points to calculate in the given energy range
   :type  npoints: int
   :param energy_range_eV: energy range in eV as a list of length 2
   :type  energy_range_eV: [float, float]
   :param points_eV: A list of points to calculate the flux at
   :type  points_eV: [float, ...]
   :param nparticles: The number of particles you wish to run for a multi-particle simulation
   :type  nparticles: int
   :param ofile: File name you would like to output this data to
   :type  ofile: str
   :param nthreads: Number of threads to use
   :type  nthreads: int
   :param gpu: Use the gpu or not
   :type  gpu: int
   :returns: A list of 2D lists, each of which is a pair representing the energy [eV] and flux [:math:`photons / mm^2 / 0.1% bw / s`] at that energy.  eg [[energy_0, flux_0], [energy_1, flux_1], ...]

   :Example: Calculate the spectrum at a point 30 [m] downstream (assuming beam direction is in the Z direction) in the energy range 100 to 1000 [eV] with 900 points.

   .. code-block:: py

      # 0.4 [T] 2 [m] long dipole centered at Z=0 [m]
      osr.add_bfield_uniform(bfield=[0, 0.4, 0], width=[0, 0, 2])

      # NSLS2 electron beam
      osr.add_particle_beam(type='electron', name='beam_0', x0=[0, 0, 0], direction=[0, 0, 1], t0=0, energy_GeV=3, current=0.500)

      # Set start and stop time for calculation
      osr.set_ctstartstop(-0.5, 0.5)


      # Set number of points for trajectory calculation (hight for dipole)
      osr.set_npoints_trajectory(20000)

      # Calculate spectrum zt Z=30 [m] in the energy range of 100 to
      # 1000 [eV] at 900 equally spaced points
      osr.calculate_spectrum(obs=[0, 0, 30], energy_range_eV=[100, 1000], npoints=900)









.. py:method:: oscars.sr.calculate_total_power()

   Calculate the total radiated power based on the current particle and that particle's beam current.

   See the :doc:`MathematicalNotes` section for the expression used in this calculation.

   :returns: float - Total power in [W]









.. py:method:: oscars.sr.calculate_power_density_rectangle(npoints, [[plane, width], x0x1x2], [rotations, translation, ofile, normal, dim, nthreads, gpu])

   Calculate the power density in a rectangle either defined by three points, or by defining the plane the rectangle is in and the width, and then rotating and translating it to where it needs be.  The simplest is outlined in the first example below.  By default (dim=2) this returns a list whose position coordinates are in the local coordinate space x1 and x2 (*ie* they do not include the rotations and translation).  if dim=3 the coordinates in the return list are in absolute 3D space.

   See the :doc:`MathematicalNotes` section for the expression used in this calculation.

   You **must** specify either both (*plane* and *width*) or *x0x1x2*

   :param npoints: number of in each dimension for surface
   :type  npoints: [int, int]
   :param plane: The plane to start in (XY, XZ, YZ, YX, ZX, ZY).  The normal to the surface is defined using the right handed cross product (ie the last three have opposite normal vectors from the first three)
   :type  plane: str
   :param width: Width if rectangle in X1 and X2
   :type  width: [float, float]
   :param x0x1x2: List of three points [[x0, y0, z0], [x1, y1, z1], [x2, y2, z2]]
   :type  x0x1x2: [[float, float, float], [float, float, float], [float, float, float]]
   :param rotations: A list representing the rotations of this beam about the X, Y amd Z axis (in that order)
   :type  rotations: list[3]
   :param translation: A list representing the translation of the x0 of this object
   :type  translation: list[3]
   :param ofile: Output file name
   :type  ofile: str
   :param normal: -1 if you wish to reverse the normal vector, 0 if you wish to ignore the +/- direction in computations, 1 if you with to use the direction of the normal vector as given. 
   :type  normal: int
   :param dim: Defaults to 2 where output is in the local plane coordinates X1 and X2.  If you want the return to be given in 3D set dim=3 which will return with X, Y, and Z in absolute coordinates.
   :type  dim: int
   :param nthreads: Number of threads to use
   :type  nthreads: int
   :param gpu: Use the gpu or not
   :type  gpu: int
   :returns: A list, each element of which is a pair representing the position (2D relative or 3D absolute) and power density [:math:`W / mm^2`] at that position.  eg [[x1_0, x2_0], pd_0, [x1_1, x2_1], pd_1],  ...]

   :Example: Calculate the power density within a simple rectangle 1 [cm] x 1 [cm], 30 [m] downstream

   .. code-block:: py

      # Calculate power density in rectangle in the XY plane 30 [m] downstream
      # from a photon beam in the +Z direction
      power_density = osr.calculate_power_density(plane='XY', width=[0.01, 0.01], npoints=[51, 51], translation=[0, 0, 30])


   :Example: Calculate the power density within a simple rectangle 1 [cm] x 1 [cm], 30 [m] downstream defined using the three-point method

   .. code-block:: py

      # Calculate power density in rectangle in the XY plane 30 [m] downstream
      # from a photon beam in the +Z direction
      power_density = osr.calculate_power_density(x0x1x2=[[-0.005, -0.005, 30], [+0.005, -0.005, 30], [-0.005, +0.005, 30]], npoints=[51, 51])


   :Example: Calculate the power density on a flat surface close to and parallel to the beam direction

   .. code-block:: py

      # Calculate power density in rectangle on the upper flat portion of a beampipe in the middle
      # of an undulator.  Here the 'normal' is reversed to get the correct sign
      osr.add_undulator(bfield=[0, 1, 0], period=0.050, nperiods=40)
      power_density = osr.calculate_power_density(plane='XZ', width=[0.008, 2], npoints=[51, 101], normal=-1)





.. py:method:: oscars.sr.calculate_power_density(points, [normal, rotations, translation, ofile, nthreads, gpu])

   Calculate the power density for each point in the list *points*.

   See the :doc:`MathematicalNotes` section for the expression used in this calculation.

   :param points: A list of points, each point containing a position in 3D (as a list) and a normal vector at that position (also as a 3D list).
   :type  points: list[ list[list[float, float, float], list[float, float, float]], ...]
   :param normal: -1 if you wish to reverse the normal vector, 0 if you wish to ignore the +/- direction in computations, 1 if you with to use the direction of the normal vector as given. 
   :type  normal: int
   :param rotations: A list representing the rotations of this beam about the X, Y amd Z axis (in that order)
   :type  rotations: list[3]
   :param translation: A list representing the translation of the x0 of this object
   :type  translation: list[3]
   :param ofile: Output file name
   :type  ofile: str
   :param nthreads: Number of threads to use
   :type  nthreads: int
   :param gpu: Use the gpu or not
   :type  gpu: int
   :returns: A list, each element of which is a pair representing the position (2D relative or 3D absolute) and power density [:math:`W / mm^2`] at that position.  eg [[x1_0, x2_0], pd_0, [x1_1, x2_1], pd_1],  ...]




.. py:method:: oscars.sr.calculate_flux_rectangle(energy_eV, npoints, [[plane, width], x0x1x2], [rotations, translation, ofile, normal, dim, nthreads, gpu])

   Calculate the flux density in a rectangle either defined by three points, or by defining the plane the rectangle is in and the width, and then rotating and translating it to where it needs be.  The simplest is outlined in the first example below.  By default (dim=2) this returns a list whose position coordinates are in the local coordinate space x1 and x2 (*ie* they do not include the rotations and translation).  if dim=3 the coordinates in the return list are in absolute 3D space.

   You **must** specify either both (*plane* and *width*) or *x0x1x2*

   See the :doc:`MathematicalNotes` section for the expression used in this calculation.


   :param energy_eV: The photon energy you are interested in
   :type  energy_eV: float
   :param npoints: number of in each dimension for surface
   :type  npoints: [int, int]
   :param plane: The plane to start in (XY, XZ, YZ, YX, ZX, ZY).  The normal to the surface is defined using the right handed cross product (ie the last three have opposite normal vectors from the first three)
   :type  plane: str
   :param width: Width if rectangle in X1 and X2
   :type  width: [float, float]
   :param x0x1x2: List of three points [[x0, y0, z0], [x1, y1, z1], [x2, y2, z2]]
   :type  x0x1x2: [[float, float, float], [float, float, float], [float, float, float]]
   :param rotations: A list representing the rotations of this beam about the X, Y amd Z axis (in that order)
   :type  rotations: list[3]
   :param translation: A list representing the translation of the x0 of this object
   :type  translation: list[3]
   :param ofile: Output file name
   :type  ofile: str
   :param normal: -1 if you wish to reverse the normal vector, 0 if you wish to ignore the +/- direction in computations, 1 if you with to use the direction of the normal vector as given. 
   :type  normal: int
   :param dim: Defaults to 2 where output is in the local plane coordinates X1 and X2.  If you want the return to be given in 3D set dim=3 which will return with X, Y, and Z in absolute coordinates.
   :type  dim: int
   :param nthreads: Number of threads to use
   :type  nthreads: int
   :param gpu: Use the gpu or not
   :type  gpu: int
   :returns: A list, each element of which is a pair representing the position (2D relative or 3D absolute) and power density [:math:`W / mm^2`] at that position.  eg [[x1_0, x2_0], pd_0, [x1_1, x2_1], pd_1],  ...]






.. py:method:: oscars.sr.average_spectra([, ifiles, bifiles, ofile, bofile])

   Average spectra from different files and output to specified file.  The input files must have the same format.  Binary operation will be supported in a future release.

   You **must** specify an input and an output

   :param ifiles: The input file names
   :type  ifiles: list
   :param ofile: The output file name
   :type  ofile: str
   :param bifiles: The binary input file names
   :type  bifiles: list
   :param bofile: The binary output file name
   :type  bofile: str
   :returns: None


.. py:method:: oscars.sr.average_flux([, ifiles, bifiles, ofile, bofile, dim])

   Average flux from different files and output to specified file.  The input files must have the same format.  Binary operation will be supported in a future release.

   You **must** specify an input and an output

   :param ifiles: The input file names
   :type  ifiles: list
   :param ofile: The output file name
   :type  ofile: str
   :param bifiles: The binary input file names
   :type  bifiles: list
   :param bofile: The binary output file name
   :type  bofile: str
   :param dim: in 2 or 3 dimensions (default is 2)
   :type  dim: integer
   :returns: None




.. py:method:: oscars.sr.average_power_density([, ifiles, bifiles, ofile, bofile, dim])

   Average power densities from different files and output to specified file.  The input files must have the same format.  Binary operation will be supported in a future release.

   You **must** specify an input and an output

   :param ifiles: The input file names
   :type  ifiles: list
   :param ofile: The output file name
   :type  ofile: str
   :param bifiles: The binary input file names
   :type  bifiles: list
   :param bofile: The binary output file name
   :type  bofile: str
   :param dim: in 2 or 3 dimensions (default is 2)
   :type  dim: integer
   :returns: None




.. py:method:: oscars.sr.add_to_spectrum(spectrum [, weight])

   Add a spectrum to the current spectrum with a given weight.

   :param spectrum: a list of pairs of numbers (spectrum format)
   :type  spectrum: list
   :param weight: Weight for *this* spectrum
   :type  weight: float
   :returns: None


.. py:method:: oscars.sr.get_spectrum()

   Get the current spectrum stored in memory

   :returns: list





.. py:method:: oscars.sr.add_to_flux(flux [, weight])

   Add a flux to the current flux with a given weight.

   :param spectrum: a list of [[x, y, z], flux] pairs (flux format)
   :type  spectrum: list
   :param weight: Weight for *this* flux
   :type  weight: float
   :returns: None


.. py:method:: oscars.sr.get_flux()

   Get the current flux map stored in memory

   :returns: list of [[x, y, z], flux]







.. py:method:: oscars.sr.add_to_power_density(power_density [, weight])

   Add a power density to the current power density with a given weight.

   :param spectrum: a list of [[x, y, z], power_density] pairs (power density format)
   :type  spectrum: list
   :param weight: Weight for *this* power density
   :type  weight: float
   :returns: None


.. py:method:: oscars.sr.get_power_density()

   Get the current power_density map stored in memory

   :returns: list of [[x, y, z], power_density]




























.. py:method:: oscars.sr.calculate_efield_vs_time(obs, [ofile])

   Calculate the electric field in the time domain for a single particle

   See the :doc:`MathematicalNotes` section for the expression used in this calculation.

   :param obs: Point where you wish to calculate the electric field
   :type  obs: [float, float, float]
   :param ofile: Output file name
   :type  ofile: str
   :returns: A list, each element of which has a time (in [s]) and a 3-dimensional list representing the x, y, and z componemts of the electric field at that time: [[t, [Ex, Ey, Ez]], ...]









oscars.plots_mpl
----------------

.. automodule:: oscars.plots_mpl
   :members:







oscars.plots3d_mpl
------------------

.. automodule:: oscars.plots3d_mpl
   :members:







oscars.parametric_surfaces
--------------------------

.. automodule:: oscars.parametric_surfaces
   :members:

