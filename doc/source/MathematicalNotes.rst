Mathematical Notes
==================


Coordinate System
-----------------

A *right-handed* cartesian coordinate system is used throughout.  To be clear:

.. math:: \hat z = \hat x \times \hat y



Undulator Parameters
--------------------

The *K* parameter is defined as follows

.. math::
   K = \frac{e B_{max} \lambda_0}{2 \pi m c^2}

or in practical units

.. math::
    K = 0.0943 B_{max} [T] \lambda_0 [mm]

where :math:`B_{max}` is the maximum magnetic field and :math:`\lambda_0` is the undulator period.


Electric Field in Frequency Domain
----------------------------------

Electric field calculation are done including the *near-field* terms.  For calculations of the flux and spectra, the electric field is calculated from the Liénard–Wiechert potentials in the frequency domain.  The derivation of this can be found in many texts ([Hoffman]_).  Some approximations remove the :math:`\frac{1}{r}` term from the integral, treating the source as a fixed distance from the observer.  No such approximation is made in this code and the calculation is valid for relativistic and non-relavistic, near and far particles.

.. math::
   \vec E(\vec r, \omega) = \frac{ie\omega}{4 \pi c \epsilon_0} \int_{-\infty}^{+\infty} \frac{\vec \beta - \hat n (1 + \frac{i c}{ \omega D})}{D} e^{i \omega ( \tau + \frac{D}{c} )} d\tau


Electric Field in Time Domain
-----------------------------

The electric field in the time domain is calculated including the *near-field* terms.  The observer time for these calculations is calculated as the retarded time from the particle position.  Currently the differential equation solver is a Runge-Kutta 4th order method with a constant time-step in the particle frame.  This means the that the time increments in the observer frame are not necessairly uniform.  It is tempting to employ a discrete fourier transform on this distribution, but if you do so please understand that the data is not uniformly distributed in the observer frame.


.. math::
   \vec E(x, y, z, t) = {\Huge\{} \frac{q}{4 \pi \epsilon_0} (\frac{1}{1 - \vec n \cdot \vec \beta})^3 {\Large [} \frac{1 - |\vec \beta|^2(\vec n - \vec \beta)}{ | \vec r |^2} + \frac{1}{ | \vec r |} \frac{1}{c} ( \vec n \times(( \vec n - \vec \beta ) \times \frac{\vec a}{c}) ) {\Large ]} {\Huge\}}_{ret}

where *ret* denotes where the retarded time :math:`\tau` is used.  The time *t* is calculated relative to :math:`\tau = 0` at the time *ctstart* as

.. math::
   t = \tau + \frac{| \vec r(\tau) |}{c}




Power Density
-------------

The power density calculated is that of the radiated power and is given by the formula

.. math::
   P(\vec x) = \frac{q I}{16 \pi^2 \epsilon_0 c} \int_{-\infty}^{\infty} \frac{\vec n \times ((\vec n - \vec \beta) \times \frac{\vec a}{c})}{(1 - \vec \beta \cdot \vec n)^5} \frac{1}{|\vec r|^2} (\hat n \cdot \hat S) \; \textrm{d}t

On a surface the power is given as :math:`[W/{mm}^{2}]`.  The value is reported is the scalar product of the power density at a given point in space with the unit normal to the surface


.. [Hoffman] Synchrotron Radiation, 1998


Trajectory Calculation
----------------------

At the moment the trajectory is a fourth order Runge-Kutta method.  If you are interested in very long distances, alternatives may exist, and are possible to implement.

.. math::
   \frac{\textrm{d} \vec p}{\textrm{d}\tau} = q (\vec E + \vec v \times \vec B)

Energy loss due to emission is not considered.  In most applications at accelerators this is rightfully neglected.  However, if this is of interest for you it can be implemented.

