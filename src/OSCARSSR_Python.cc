////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Fri Jun 17 10:39:12 EDT 2016
//
// This is the python-C extension which allows access to the c++
// class OSCARSSR.
//
////////////////////////////////////////////////////////////////////

// Include Python.h first!
#include <Python.h>

#include "OSCARSSR_Python.h"

#include "OSCARSPY.h"
#include "OSCARSSR.h"
#include "Version.h"

#include "TSurfacePoints_Rectangle.h"
#include "TSurfacePoints_3D.h"
#include "T3DScalarContainer.h"
#include "TFieldPythonFunction.h"
#include "TField3D_Gaussian.h"
#include "TField3D_UniformBox.h"
#include "TField3D_IdealUndulator.h"
#include "TField3D_Halbach.h"
#include "TField3D_Quadrupole.h"
#include "TDriftBox.h"
#include "TRandomA.h"
#include "TTriangle3DContainer.h"

#include <iostream>
#include <vector>
#include <algorithm>
#include <cstring>
#include <stdexcept>
#include <sstream>


// External global random generator
extern TRandomA* gRandomA;





static void OSCARSSR_dealloc(OSCARSSRObject* self)
{
  // Python needs to know how to deallocate things in the struct

  delete self->obj;
  //self->ob_type->tp_free((PyObject*) self);
  Py_TYPE(self)->tp_free((PyObject*) self);
}




static PyObject* OSCARSSR_new (PyTypeObject* type, PyObject* args, PyObject* keywds)
{
  // Python needs to know how to create things in this struct

  // Grab the values
  int NThreads = 0;
  int GPU = 0;

  // Input variables and parsing
  static const char *kwlist[] = {"nthreads",
                                 "gpu",
                                 NULL};

  if (!PyArg_ParseTupleAndKeywords(args, keywds, "|ii",
                                   const_cast<char **>(kwlist),
                                   &NThreads,
                                   &GPU)) {
    PyErr_SetString(PyExc_ValueError, "allowed inputs are currentl: 'nthreads', 'gpu'");
    return NULL;
  }

  // Allocate objects
  OSCARSSRObject* self = (OSCARSSRObject*) type->tp_alloc(type, 0);
  if (self != NULL) {

    // Create the new object for self
    self->obj = new OSCARSSR();
  }


  if (NThreads > 0) {
    self->obj->SetNThreadsGlobal(NThreads);
  }
  if (GPU != 0 && GPU != 1) {
    PyErr_SetString(PyExc_ValueError, "global gpu settign must be 0 or 1");
    return NULL;
  }

  // If it was not successful print an error message
  if (!self->obj->SetUseGPUGlobal(GPU)) {
    OSCARSPY::PyPrint_stderr("GPU is not available: Setting gpu global setting to 0.\n");
  }


  // Return myself
  return (PyObject*) self;
}





//static int sr_init(OSCARSSRObject* self, PyObject* args, PyObject* kwds)
//{
//  return 0;
//}














const char* DOC_OSCARSSR_Version = R"docstring(
version()

Version ID

Returns
-------
version : str
)docstring";
static PyObject* OSCARSSR_Version (OSCARSSRObject* self, PyObject* arg)
{
  return Py_BuildValue("s", OSCARSPY::GetVersionString().c_str());
}






const char* DOC_OSCARSSR_Pi = R"docstring(
pi()

The value of pi

Returns
-------
pi : float
)docstring";
static PyObject* OSCARSSR_Pi (OSCARSSRObject* self, PyObject* arg)
{
  // Return the internal OSCARSSR number constant pi
  return Py_BuildValue("d", TOSCARSSR::Pi());
}






const char* DOC_OSCARSSR_Qe = R"docstring(
qe()

Get the value of elementary charge: +1.602176462e-19 [C]

Returns
-------
charge of electron : float
)docstring";
static PyObject* OSCARSSR_Qe (OSCARSSRObject* self, PyObject* arg)
{
  // Return the internal OSCARSSR number for elementary charge
  return Py_BuildValue("d", TOSCARSSR::Qe());
}




const char* DOC_OSCARSSR_Me = R"docstring(
me()

Get the value of the mass of the electron 9.10938356e-31 [kg]

Returns
-------
electron mass : float
)docstring";
static PyObject* OSCARSSR_Me (OSCARSSRObject* self, PyObject* arg)
{
  // Return the internal OSCARSSR number for mass of the electron [kg]
  return Py_BuildValue("d", TOSCARSSR::Me());
}




const char* DOC_OSCARSSR_Random = R"docstring(
rand()

Uniformally distributed random float in the range [0, 1)

Returns
-------
float
)docstring";
static PyObject* OSCARSSR_Random (OSCARSSRObject* self, PyObject* arg)
{
  return Py_BuildValue("d", gRandomA->Uniform());
}




const char* DOC_OSCARSSR_RandomNormal = R"docstring(
norm()

Random float from the normal distribution centered at 0 with a sigma of 1

Returns
-------
float
)docstring";
static PyObject* OSCARSSR_RandomNormal (OSCARSSRObject* self, PyObject* arg)
{
  return Py_BuildValue("d", gRandomA->Normal());
}




const char* DOC_OSCARSSR_SetSeed = R"docstring(
set_seed(n)

Set the internal random seed

Parameters
----------
n : float
    Seed number you wish to use

Returns
-------
None
)docstring";
static PyObject* OSCARSSR_SetSeed (OSCARSSRObject* self, PyObject* arg)
{
  // Grab the value from input
  double Seed = PyFloat_AsDouble(arg);

  gRandomA->SetSeed(Seed);

  // Must return python object None in a special way
  Py_INCREF(Py_None);
  return Py_None;
}




const char* DOC_OSCARSSR_SetGPUGlobal = R"docstring(
set_gpu_global(gpu)

If set to 1, OSCARS will try to use the GPU for all calculations

Parameters
----------
gpu : int
    This will set the gpu active

Returns
-------
None
)docstring";
static PyObject* OSCARSSR_SetGPUGlobal (OSCARSSRObject* self, PyObject* arg)
{
  // Grab the value from input
  int const GPU = (int) PyLong_AsLong(arg);

  if (GPU != 0 && GPU != 1) {
    PyErr_SetString(PyExc_ValueError, "global gpu settign must be 0 or 1");
    return NULL;
  }

  // If it was not successful print an error message
  if (!self->obj->SetUseGPUGlobal(GPU)) {
    OSCARSPY::PyPrint_stderr("GPU is not available: Setting gpu global setting to 0.\n");
  }

  // Must return python object None in a special way
  Py_INCREF(Py_None);
  return Py_None;
}




const char* DOC_OSCARSSR_CheckGPU = R"docstring(
check_gpu()

Will return the number of GPUs available, or a negative number for any other error, for instance if your distribution was not compiled with GPU support.

Returns
-------
ngpu : int
)docstring";
static PyObject* OSCARSSR_CheckGPU (OSCARSSRObject* self, PyObject* arg)
{

  int const NGPUStatus = self->obj->CheckGPU();

  if (NGPUStatus == -1) {
    OSCARSPY::PyPrint_stderr("It appears this binary version of OSCARSSR was not compiled with GPU capability enabled.\n");
  }

  return PyLong_FromLong((long) NGPUStatus);
}





const char* DOC_OSCARSSR_SetNThreadsGlobal = R"docstring(
set_nthreads_global(nthreads)

Set the number of threads you wish to use for all calculations.  If the GPU is requested it will take precedence over this multi-threading.

Parameters
----------
nthreads : int
    Default number of threads to use for calculations

Returns
-------
None
)docstring";
static PyObject* OSCARSSR_SetNThreadsGlobal (OSCARSSRObject* self, PyObject* arg)
{
  // Grab the value from input
  int const NThreads = (int) PyLong_AsLong(arg);

  self->obj->SetNThreadsGlobal(NThreads);

  // Must return python object None in a special way
  Py_INCREF(Py_None);
  return Py_None;
}













const char* DOC_OSCARSSR_GetCTStart = R"docstring(
get_ctstart()

Gets the start time for particle propogation.  The *time* is in units of [m].

Returns
-------
time, float
)docstring";
static PyObject* OSCARSSR_GetCTStart (OSCARSSRObject* self)
{
  // Get the start time in [m] for calculation
  return Py_BuildValue("d", self->obj->GetCTStart());
}



const char* DOC_OSCARSSR_GetCTStop = R"docstring(
get_ctstop()

Gets the stop time for particle propogation.  The *time* is in units of [m].

Returns
-------
time : float
)docstring";
static PyObject* OSCARSSR_GetCTStop (OSCARSSRObject* self)
{
  // Get the CTStop variable from OSCARSSR
  return Py_BuildValue("d", self->obj->GetCTStop());
}



const char* DOC_OSCARSSR_SetCTStartStop = R"docstring(
set_ctstartstop(start, stop)

Set the start and stop times for the trajectory calculation.  Start time must be less than or equal to the T0 defined in the initial conditions for the particle beam.  Similarly, stop time must be greater than that T0.  This also sets a default value for the number of points used in the trajectory calculation.

Parameters
----------
start : float
    Start time

stop : float
    Stop time

Returns
-------
None
)docstring";
static PyObject* OSCARSSR_SetCTStartStop (OSCARSSRObject* self, PyObject* args)
{
  // Set the start and stop times for OSCARSSR in [m]

  // Grab the values
  double Start, Stop;
  if (!PyArg_ParseTuple(args, "dd", &Start, &Stop)) {
    PyErr_SetString(PyExc_ValueError, "Incorrect format");
    return NULL;
  }

  // Set the object variable
  self->obj->SetCTStartStop(Start, Stop);

  // Must return python object None in a special way
  Py_INCREF(Py_None);
  return Py_None;
}





const char* DOC_OSCARSSR_GetNPointsTrajectory = R"docstring(
get_npoints_trajectory()

Get the number of points to be used for the trajectory calculation

Returns
-------
npoints : int
)docstring";
static PyObject* OSCARSSR_GetNPointsTrajectory (OSCARSSRObject* self)
{
  // Get the numper of points for trajectory calculaton
  return PyLong_FromSize_t(self->obj->GetNPointsTrajectory());
}



const char* DOC_OSCARSSR_SetNPointsTrajectory = R"docstring(
set_npoints_trajectory(npoints)

Sets the number of points to be used in the trajectory calculation.  Only call this function **after** set_ctstartstop().

Parameters
----------
npoints : int
    Number of points

Returns
-------
None

example:

   .. code-block:: py

      # First set the start and stop times
      osr.set_ctstartstop(0, 2)

      # Use 12345 points in the trajectory calculations
      osr.set_npoints_trajectory(12345)
)docstring";
static PyObject* OSCARSSR_SetNPointsTrajectory (OSCARSSRObject* self, PyObject* arg)
{
  // Set the number of points for trajectory calculation

  // Grab the value from input
  size_t N = PyLong_AsSsize_t(arg);

  // Set the object variable
  self->obj->SetNPointsTrajectory(N);

  // Must return python object None in a special way
  Py_INCREF(Py_None);
  return Py_None;
}







const char* DOC_OSCARSSR_SetNPointsPerMeterTrajectory = R"docstring(
set_npoints_per_meter_trajectory(npoints)

Sets the number of points per meter to be used in the trajectory calculation.

Parameters
----------
npoints : int
    Number of points per neter

Returns
-------
None
)docstring";
static PyObject* OSCARSSR_SetNPointsPerMeterTrajectory (OSCARSSRObject* self, PyObject* arg)
{
  // Set the number of points for trajectory calculation

  // Grab the value from input
  size_t N = PyLong_AsSsize_t(arg);

  // Set the object variable
  self->obj->SetNPointsPerMeterTrajectory(N);

  // Must return python object None in a special way
  Py_INCREF(Py_None);
  return Py_None;
}







const char* DOC_OSCARSSR_AddMagneticField = R"docstring(
add_bfield_file([, ifile, bifile, iformat, rotations, translation, scale, name])

Add a magnetic field from a text file *ifile* according to the format *iformat*.

Currently OSCARS accepts the following file formats for iformat
    * 'OSCARS'
    * 'OSCARS1D [plus format string]'
    * 'SPECTRA'
    * 'SRW'

For the OSCARS1D format you must also include the order of the data columns, which would typically look something like: 'OSCARS1D Z Bx By Bz'.  You may use X or Y instead of Z and the order and number of B[xyz] does not matter.  This mode also accepts non-uniformly distributed data.

Optionally you can rotate and translate the field in space.  You can use an input *scale* list to scale the input (which must be in SI units of [m] for distances/positions and [T] for magnetic field values.

The rotation is performed first in the order: :math:`\theta_x, \theta_y, \theta_z`

*scale* is a list of less than or equal length to the number of elements in *iformat* when OSCARS1D is selected.  This will scale the input values of the i-th column before any rotation or translation.  This is useful if your data file is not in [T] and [m]

Parameters
----------
ifile : str
    Name of input file

bifile : str
    Name of binary input file

iformat : str
    Format of the input file.  Must be specified with ifile, but not with bifile

rotations : list, optional
    3-element list representing rotations around x, y, and z axes: [:math:`\theta_x, \theta_y, \theta_z`]

translation : list, optional
    3-element list representing a translation in space [x, y, z]

scale : list, optional
    List of scale factors to be used for multiplying the inputs in the order of iformat (equal in length or less than the number of input parameters)

name : str
    Name of this magnetic field

Returns
-------
None

Examples
--------
Add a magnetic field from a data file where the data is in columns in the order Z Bx By Bz

    >>> osr.add_bfield_file(ifile='file.txt', iformat='OSCARS1D Z Bx By Bz')
)docstring";
static PyObject* OSCARSSR_AddMagneticField (OSCARSSRObject* self, PyObject* args, PyObject* keywds)
{
  // Add a magnetic field from a file.
  // UPDATE: add binary file reading


  // Grab the values
  char const* FileNameText     = "";
  char const* FileNameBinary   = "";
  char const* FileFormat       = "";
  PyObject*   List_Rotations   = PyList_New(0);
  PyObject*   List_Translation = PyList_New(0);
  PyObject*   List_Scaling     = PyList_New(0);
  char const* Name             = "";

  TVector3D Rotations(0, 0, 0);
  TVector3D Translation(0, 0, 0);
  std::vector<double> Scaling;


  // Input variables and parsing
  static const char *kwlist[] = {"ifile",
                                 "bifile",
                                 "iformat",
                                 "rotations",
                                 "translation",
                                 "scale",
                                 "name",
                                 NULL};

  if (!PyArg_ParseTupleAndKeywords(args, keywds, "|sssOOOs",
                                   const_cast<char **>(kwlist),
                                   &FileNameText,
                                   &FileNameBinary,
                                   &FileFormat,
                                   &List_Rotations,
                                   &List_Translation,
                                   &List_Scaling,
                                   &Name)) {
    return NULL;
  }

  // Check that ifile and bifile both not input
  if (std::strlen(FileNameText) != 0 && std::strlen(FileNameBinary) != 0) {
    PyErr_SetString(PyExc_ValueError, "cannot specify both 'ifile' and 'bifile'");
    return NULL;
  }

  // Check that filename and format exist
  if (std::strlen(FileNameText) != 0 && std::strlen(FileFormat) == 0) {
    PyErr_SetString(PyExc_ValueError, "'iformat' is blank");
    return NULL;
  }

  // Check for Rotations in the input
  if (PyList_Size(List_Rotations) != 0) {
    try {
      Rotations = OSCARSPY::ListAsTVector3D(List_Rotations);
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'rotations'");
      return NULL;
    }
  }

  // Check for Translation in the input
  if (PyList_Size(List_Translation) != 0) {
    try {
      Translation = OSCARSPY::ListAsTVector3D(List_Translation);
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'translation'");
      return NULL;
    }
  }

  // Get any scaling factors
  // UPDATE: Check against fileformat number of strings
  for (int i = 0; i < PyList_Size(List_Scaling); ++i) {
    Scaling.push_back(PyFloat_AsDouble(PyList_GetItem(List_Scaling, i)));
  }

  // Name check
  if (std::string(Name).size() > 0 && Name[0] == '_') {
    PyErr_SetString(PyExc_ValueError, "'name' cannot begin with '_'.  This is reserved for internal use.  Please pick a different name");
    return NULL;
  }

  // Add the magnetic field to the OSCARSSR object
  try {
    if (std::strlen(FileNameBinary) != 0) {
      self->obj->AddMagneticField(FileNameBinary, "BINARY", Rotations, Translation, Scaling, Name);
    } else {
      self->obj->AddMagneticField(FileNameText, FileFormat, Rotations, Translation, Scaling, Name);
    }
  } catch (...) {
    PyErr_SetString(PyExc_ValueError, "Could not import magnetic field.  Check 'ifile' and 'iformat' are correct");
    return NULL;
  }

  // Must return python object None in a special way
  Py_INCREF(Py_None);
  return Py_None;
}







const char* DOC_OSCARSSR_AddMagneticFieldInterpolated = R"docstring(
add_bfield_interpolated(mapping, iformat, parameter [, rotations, translation, scale, name])

Add field given a paramater and a mapping between known parameters and field data files where the field is interpolated from the known data points.

Parameters
----------
mapping : list [[float, str], [float, str], ...]
    List of parameters and associated filenames [[p1, file1], [p2, file2], ...]

iformat : str
    Which input format to use (see add_bfield_file() for formats)

parameter : float
    Value of parameter you are interested in

rotations : list, optional
    3-element list representing rotations around x, y, and z axes: [:math:`\theta_x, \theta_y, \theta_z`]

translation : list, optional
    3-element list representing a translation in space [x, y, z]

scale : list, optional
    List of scale factors to be used for multiplying the inputs in the order of iformat (equal in length or less than the number of input parameters)

name : str
    Name of this magnetic field

Returns
-------
None

Examples
--------
Undulator gap interpolation between known magnetic measurement points

    >>> file_list = [
    ...     [10.9, 'file_10.9.dat'],
    ...     [11.0, 'file_11.0.dat'],
    ...     [13.2, 'file_13.2.dat'],
    ...     [16.9, 'file_16.9.dat']
    ... ]
    >>> osr.add_bfield_interpolated(mapping=file_list, iformat='OSCARS1D Z Bx By Bz', parameter=12.123)
)docstring";
static PyObject* OSCARSSR_AddMagneticFieldInterpolated (OSCARSSRObject* self, PyObject* args, PyObject* keywds)
{
  // Add a magnetic field from a file.
  // UPDATE: add binary file reading

  // Grab the values
  PyObject*   List_Mapping     = PyList_New(0);
  char const* FileFormat       = "";
  double      Parameter        = 0;
  PyObject*   List_Rotations   = PyList_New(0);
  PyObject*   List_Translation = PyList_New(0);
  PyObject*   List_Scaling     = PyList_New(0);
  char const* Name             = "";

  TVector3D Rotations(0, 0, 0);
  TVector3D Translation(0, 0, 0);
  std::vector<double> Scaling;

  // Mapping that is passed in for field interpolation
  std::vector<std::pair<double, std::string> > Mapping;

  // Input variables and parsing
  static const char *kwlist[] = {"mapping",
                                 "iformat",
                                 "parameter",
                                 "rotations",
                                 "translation",
                                 "scale",
                                 "name",
                                 NULL};

  if (!PyArg_ParseTupleAndKeywords(args, keywds, "Osd|OOOs",
                                   const_cast<char **>(kwlist),
                                   &List_Mapping,
                                   &FileFormat,
                                   &Parameter,
                                   &List_Rotations,
                                   &List_Translation,
                                   &List_Scaling,
                                   &Name)) {
    return NULL;
  }

  // Grab mapping from input
  if (PyList_Size(List_Mapping) != 0) {
    double ParameterValue;
    for (int i = 0; i != PyList_Size(List_Mapping); ++i) {
      PyObject* ThisPair = PyList_GetItem(List_Mapping, i);
      if (PyList_Size(ThisPair) != 2) {
        PyErr_SetString(PyExc_ValueError, "Incorrect format in 'mapping'");
        return NULL;
      }

      ParameterValue = PyFloat_AsDouble(PyList_GetItem(ThisPair, 0));
      std::string const FileName = OSCARSPY::GetAsString(PyList_GetItem(ThisPair, 1));

      Mapping.push_back(std::make_pair(ParameterValue, FileName));
    }
  }


  // Check that filename and format exist
  if (std::strlen(FileFormat) == 0) {
    PyErr_SetString(PyExc_ValueError, "'iformat' is blank");
    return NULL;
  }

  // Parameter value is what it is...

  // Check for Rotations in the input
  if (PyList_Size(List_Rotations) != 0) {
    try {
      Rotations = OSCARSPY::ListAsTVector3D(List_Rotations);
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'rotations'");
      return NULL;
    }
  }


  // Check for Translation in the input
  if (PyList_Size(List_Translation) != 0) {
    try {
      Translation = OSCARSPY::ListAsTVector3D(List_Translation);
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'translation'");
      return NULL;
    }
  }

  // Get any scaling factors
  // UPDATE: Check against fileformat number of strings
  for (int i = 0; i < PyList_Size(List_Scaling); ++i) {
    Scaling.push_back(PyFloat_AsDouble(PyList_GetItem(List_Scaling, i)));
  }

  // Name check
  if (std::string(Name).size() > 0 && Name[0] == '_') {
    PyErr_SetString(PyExc_ValueError, "'name' cannot begin with '_'.  This is reserved for internal use.  Please pick a different name");
    return NULL;
  }

  // Add the magnetic field to the OSCARSSR object
  try {
    self->obj->AddMagneticFieldInterpolated(Mapping, FileFormat, Parameter, Rotations, Translation, Scaling, Name);
  } catch (...) {
    PyErr_SetString(PyExc_ValueError, "Could not import magnetic field.  Check filenames and 'iformat' are correct");
    return NULL;
  }

  // Must return python object None in a special way
  Py_INCREF(Py_None);
  return Py_None;
}







const char* DOC_OSCARSSR_AddMagneticFieldFunction = R"docstring(
add_bfield_function(func [, name])

Adds field in the form of a user defined python function.  The input for this function must be (x, y, z, t) and return [Fx, Fy, Fz]

Parameters
----------

function : func
    A python function with input [x, y, z, t] and return of [Fx, Fy, Fz]

name : str
    Name of this field

Returns
-------
None

Examples
--------
Create a function in python and use it as a field in OSCARS

    >>> def myfunc(x, y, z, t):
    ...     "Do not forget to write a docstring"
    ...     if z > 0:
    ...         return 1
    ...     return 0
    >>> osr.add_bfield_function(myfunc)
)docstring";
static PyObject* OSCARSSR_AddMagneticFieldFunction (OSCARSSRObject* self, PyObject* args, PyObject* keywds)
{
  // Add a python function as a magnetic field object

  // Variables for function and name
  PyObject* Function;
  char const* Name = "";

  // Input variables and parsing
  static const char *kwlist[] = {"function",
                                 "name",
                                 NULL};

  if (!PyArg_ParseTupleAndKeywords(args, keywds, "O|s",
                                   const_cast<char **>(kwlist),
                                   &Function,
                                   &Name)) {
    return NULL;
  }

  // Increment ref to function for python
  Py_INCREF(Function);

  // Add the function as a field to the OSCARSSR object
  try {
    self->obj->AddMagneticField( (TField*) new TFieldPythonFunction(Function, Name));
  } catch (std::invalid_argument e) {
    PyErr_SetString(PyExc_ValueError, e.what());
    return NULL;
  }

  // Decrement reference to function for python
  Py_DECREF(Function);

  // Must return python object None in a special way
  Py_INCREF(Py_None);
  return Py_None;
}







const char* DOC_OSCARSSR_AddMagneticFieldGaussian = R"docstring(
add_bfield_gaussian(bfield, sigma [, rotations, translation, name])

Add a gaussian field in 3D with the peak field magnitude and direction given by *bfield*, centered about a point with a given sigma in each coordinate.  If any component of *sigma* is less than or equal to zero it is ignored (ie. spatial extent is infinite).

The functional form for this is :math:`\exp(-(x - x_0)^2 / \sigma_x^2)`


Parameters
----------
bfield : list
    A list representing the peak field [Fx, Fy, Fz]

sigma : list
    A list representing the sigma in each coordinate :math:`[\sigma_x, \sigma_y, \sigma_z]`

rotations : list, optional
    3-element list representing rotations around x, y, and z axes: [:math:`\theta_x, \theta_y, \theta_z`]

translation : list, optional
    3-element list representing a translation in space [x, y, z]

name : str
    Name of this field

Returns
-------
None

Examples
--------
Add a field of 0.5 [T] in the X-direction centered at X=0, Y=0, Z=0 [m], with a sigma of 0.1 [m] in the Z-direction

    >>> osr.add_bfield_gaussian(bfield=[1, 0, 0], sigma=[0, 0, 0.10])

Add a field of 1 [T] in the Y-direction centered at X=1, Y=1, Z=1 [m], with a sigma of 0.05 in the Z-direction

    >>> osr.add_bfield_gaussian(bfield=[0, 1, 0], sigma=[0, 0, 0.05], translation=[1, 1, 1])
)docstring";
static PyObject* OSCARSSR_AddMagneticFieldGaussian (OSCARSSRObject* self, PyObject* args, PyObject* keywds)
{
  // Add a magnetic field that is a gaussian

  // Lists and variables
  PyObject* List_BField       = PyList_New(0);
  PyObject* List_Translation  = PyList_New(0);
  PyObject* List_Rotations    = PyList_New(0);
  PyObject* List_Sigma        = PyList_New(0);
  const char* Name            = "";

  TVector3D BField(0, 0, 0);
  TVector3D Sigma(0, 0, 0);
  TVector3D Rotations(0, 0, 0);
  TVector3D Translation(0, 0, 0);


  // Input variables and parsing
  static const char *kwlist[] = {"bfield",
                                 "sigma",
                                 "rotations",
                                 "translation",
                                 "name",
                                 NULL};

  if (!PyArg_ParseTupleAndKeywords(args, keywds, "OO|OOs",
                                   const_cast<char **>(kwlist),
                                   &List_BField,
                                   &List_Sigma,
                                   &List_Rotations,
                                   &List_Translation,
                                   &Name)) {
    return NULL;
  }



  // Check BField
  try {
    BField = OSCARSPY::ListAsTVector3D(List_BField);
  } catch (std::length_error e) {
    PyErr_SetString(PyExc_ValueError, "Incorrect format in 'bfield'");
    return NULL;
  }

  // Check Width
  try {
    Sigma = OSCARSPY::ListAsTVector3D(List_Sigma);
  } catch (std::length_error e) {
    PyErr_SetString(PyExc_ValueError, "Incorrect format in 'sigma'");
    return NULL;
  }

  // Check for Rotations in the input
  if (PyList_Size(List_Rotations) != 0) {
    try {
      Rotations = OSCARSPY::ListAsTVector3D(List_Rotations);
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'rotations'");
      return NULL;
    }
  }

  // Check for Translation in the input
  if (PyList_Size(List_Translation) != 0) {
    try {
      Translation = OSCARSPY::ListAsTVector3D(List_Translation);
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'translation'");
      return NULL;
    }
  }

  // Name check
  if (std::string(Name).size() > 0 && Name[0] == '_') {
    PyErr_SetString(PyExc_ValueError, "'name' cannot begin with '_'.  This is reserved for internal use.  Please pick a different name");
    return NULL;
  }

  // Add field
  self->obj->AddMagneticField( (TField*) new TField3D_Gaussian(BField, Translation, Sigma, Rotations, Name));

  // Must return python object None in a special way
  Py_INCREF(Py_None);
  return Py_None;
}







const char* DOC_OSCARSSR_AddMagneticFieldUniform = R"docstring(
add_bfield_uniform(bfield [, width, rotations, translation, name])

Add a uniform field in a given range or for all space.  The *bfield* is given as a 3D vector representing the field magnitude and direction.  *width* is an optional parameters, if not present the field permeates all space.  If a component of the 3D list *width* is less than or equal to zero, that coordinate will be ignored when calculating the field.

Parameters
----------
bfield : list
    A list representing the field [Fx, Fy, Fz]

width : list
    A list representing the spetial extent of the field [Wx, Wy, Wz]

rotations : list, optional
    3-element list representing rotations around x, y, and z axes: [:math:`\theta_x, \theta_y, \theta_z`]

translation : list, optional
    3-element list representing a translation in space [x, y, z]

name : str
    Name of this field

Returns
-------
None

Examples
--------
Add a field of 0.0001 [T] in the X-direction for all space

    >>> osr.add_bfield_uniform(bfield=[0.0001, 0, 0])

Add a field of 0.0005 [T] in the Y-direction with a width in the Z-direction of 1.5 [m] (the other directions are ignored) centered at X=0, Y=0, Z=0.75 [m].

    >>> osr.add_bfield_uniform(bfield=[0, 1, 0], width=[0, 0, 1.5], translation=[0, 0, 0.75])

Add a field of 1 [T] in the X-direction in a volume of dimensions 1m x 1m x 1m.

    >>> osr.add_bfield_uniform(bfield=[1, 0, 0], width=[1, 1, 1])
)docstring";
static PyObject* OSCARSSR_AddMagneticFieldUniform (OSCARSSRObject* self, PyObject* args, PyObject* keywds)
{
  // Add a uniform field with a given width in a given direction, or for all space

  // Lists and vectors
  PyObject*   List_Field       = PyList_New(0);
  PyObject*   List_Translation = PyList_New(0);
  PyObject*   List_Rotations   = PyList_New(0);
  PyObject*   List_Width       = PyList_New(0);
  char const* Name             = "";

  TVector3D Field(0, 0, 0);
  TVector3D Width (0, 0, 0);
  TVector3D Rotations(0, 0, 0);
  TVector3D Translation(0, 0, 0);


  // Input variables and parsing
  static const char *kwlist[] = {"bfield",
                                 "width",
                                 "rotations",
                                 "translation",
                                 "name",
                                 NULL};

  if (!PyArg_ParseTupleAndKeywords(args, keywds, "O|OOOs",
                                   const_cast<char **>(kwlist),
                                   &List_Field,
                                   &List_Width,
                                   &List_Rotations,
                                   &List_Translation,
                                   &Name
                                   )) {
    return NULL;
  }


  // Check Field
  try {
    Field = OSCARSPY::ListAsTVector3D(List_Field);
  } catch (std::length_error e) {
    PyErr_SetString(PyExc_ValueError, "Incorrect format in 'bfield'");
    return NULL;
  }

  // Check Width
  if (PyList_Size(List_Width) != 0) {
    try {
      Width = OSCARSPY::ListAsTVector3D(List_Width);
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'width'");
      return NULL;
    }
  }

  // Check for Rotations in the input
  if (PyList_Size(List_Rotations) != 0) {
    try {
      Rotations = OSCARSPY::ListAsTVector3D(List_Rotations);
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'rotations'");
      return NULL;
    }
  }

  // Check for Translation in the input
  if (PyList_Size(List_Translation) != 0) {
    try {
      Translation = OSCARSPY::ListAsTVector3D(List_Translation);
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'translation'");
      return NULL;
    }
  }

  // Name check
  if (std::string(Name).size() > 0 && Name[0] == '_') {
    PyErr_SetString(PyExc_ValueError, "'name' cannot begin with '_'.  This is reserved for internal use.  Please pick a different name");
    return NULL;
  }

  // Add the field
  self->obj->AddMagneticField((TField*) new TField3D_UniformBox(Field, Width, Translation, Rotations, Name));

  // Must return python object None in a special way
  Py_INCREF(Py_None);
  return Py_None;
}




const char* DOC_OSCARSSR_AddMagneticFieldIdealUndulator = R"docstring(
add_bfield_undulator(bfield, period, nperiods [, phase, rotations, translation, taper, name])

Adds an ideal sinusoidal undulator field with a given maximum bfield amplitude, period, and number of periods.  Optionally one can specify the phase offset (in [rad]), rotations and translation.  The number of periods given is the full number of fields not counting the terminating fields.

Parameters
----------
bfield : list
    A list representing the peak field [Bx, By, Bz] in [T]

period : list
    Length of one period

nperiods : int
    Number of periods

phase : float
    Phase offset in [rad]

rotations : list, optional
    3-element list representing rotations around x, y, and z axes: [:math:`\theta_x, \theta_y, \theta_z`]

translation : list, optional
    3-element list representing a translation in space [x, y, z]

taper : float
    The fractional change of the bfield per meter along the magnetic axis

name : str
    Name of this field

Returns
-------
None

Examples
--------
Add an idealized undulator with 41 periods having a period of 0.050 [m] with a maximum field of 1 [T] in the y-direction where the magnetic axis is along the z-axis

    >>> osr.add_bfield_undulator(bfield=[0, 1, 0], period=[0, 0, 0.050], nperiods=41)
)docstring";
static PyObject* OSCARSSR_AddMagneticFieldIdealUndulator (OSCARSSRObject* self, PyObject* args, PyObject* keywds)
{
  // Add a magnetic field for undulator

  // Lists and variables
  PyObject*   List_Field       = PyList_New(0);
  PyObject*   List_Period      = PyList_New(0);
  PyObject*   List_Rotations   = PyList_New(0);
  PyObject*   List_Translation = PyList_New(0);
  int         NPeriods         = 0;
  double      Phase            = 0;
  double      Taper            = 0;
  char const* Name             = "";

  TVector3D Field(0, 0, 0);
  TVector3D Period(0, 0, 0);
  TVector3D Rotations(0, 0, 0);
  TVector3D Translation(0, 0, 0);

  // Input variables and parsing
  static const char *kwlist[] = {"bfield",
                                 "period",
                                 "nperiods",
                                 "phase",
                                 "rotations",
                                 "translation",
                                 "taper",
                                 "name",
                                 NULL};

  if (!PyArg_ParseTupleAndKeywords(args, keywds, "OOi|dOOds",
                                   const_cast<char **>(kwlist),
                                   &List_Field,
                                   &List_Period,
                                   &NPeriods,
                                   &Phase,
                                   &List_Rotations,
                                   &List_Translation,
                                   &Taper,
                                   &Name
                                   )) {
    return NULL;
  }


  // Check Field
  try {
    Field = OSCARSPY::ListAsTVector3D(List_Field);
  } catch (std::length_error e) {
    PyErr_SetString(PyExc_ValueError, "Incorrect format in 'bfield'");
    return NULL;
  }

  // Check Period
  try {
    Period = OSCARSPY::ListAsTVector3D(List_Period);
  } catch (std::length_error e) {
    PyErr_SetString(PyExc_ValueError, "Incorrect format in 'period'");
    return NULL;
  }

  // Check for Rotations in the input
  if (PyList_Size(List_Rotations) != 0) {
    try {
      Rotations = OSCARSPY::ListAsTVector3D(List_Rotations);
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'rotations'");
      return NULL;
    }
  }

  // Check for Translation in the input
  if (PyList_Size(List_Translation) != 0) {
    try {
      Translation = OSCARSPY::ListAsTVector3D(List_Translation);
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'translation'");
      return NULL;
    }
  }

  // Name check
  if (std::string(Name).size() > 0 && Name[0] == '_') {
    PyErr_SetString(PyExc_ValueError, "'name' cannot begin with '_'.  This is reserved for internal use.  Please pick a different name");
    return NULL;
  }


  // Rotate field and sigma
  // UPDATE: check this
  Field.RotateSelfXYZ(Rotations);
  Period.RotateSelfXYZ(Rotations);


  // Add field
  self->obj->AddMagneticField( (TField*) new TField3D_IdealUndulator(Field, Period, NPeriods, Translation, Phase, Taper, Name));

  // Must return python object None in a special way
  Py_INCREF(Py_None);
  return Py_None;
}







const char* DOC_OSCARSSR_AddMagneticFieldHalbach = R"docstring(
add_bfield_halbach()


Parameters
----------

Returns
-------
None

Examples
--------
Add an idealized undulator with 41 periods having a period of 0.050 [m] with a maximum field of 1 [T] in the y-direction where the magnetic axis is along the z-axis

    >>> osr.add_
)docstring";
static PyObject* OSCARSSR_AddMagneticFieldHalbach (OSCARSSRObject* self, PyObject* args, PyObject* keywds)
{
  // Add a magnetic field for undulator

  // Lists and variables
  double      BField           = 0;
  double      Period           = 0;
  double      MagnetWidth      = 0;
  int         NPerPeriod       = 4;
  double      Gap              = 0;
  double      MagnetHeight     = 0;
  int         NPeriods         = 0;
  char const* Name             = "";

  TVector3D Rotations(0, 0, 0);
  TVector3D Translation(0, 0, 0);

  // Input variables and parsing
  static const char *kwlist[] = {"bfield",
                                 "period",
                                 "nperiods",
                                 "gap",
                                 "magnet_height",
                                 "magnet_width",
                                 "n_per_period",
                                 "name",
                                 NULL};

  if (!PyArg_ParseTupleAndKeywords(args, keywds, "ddidd|dis",
                                   const_cast<char **>(kwlist),
                                   &BField,
                                   &Period,
                                   &NPeriods,
                                   &Gap,
                                   &MagnetHeight,
                                   &MagnetWidth,
                                   &NPerPeriod,
                                   &Name
                                   )) {
    return NULL;
  }


  // Check period
  if (Period <= 0) {
    PyErr_SetString(PyExc_ValueError, "'period' must be > 0");
    return NULL;
  }

  // Check nperiods
  if (NPeriods < 1) {
    PyErr_SetString(PyExc_ValueError, "'nperiods' must be > 0");
    return NULL;
  }

  // Check gap
  if (Gap <= 0) {
    PyErr_SetString(PyExc_ValueError, "'gap' must be > 0");
    return NULL;
  }

  // Check Magnet Height
  if (MagnetHeight <= 0) {
    PyErr_SetString(PyExc_ValueError, "'magnet_height' must be > 0");
    return NULL;
  }

  // Check NPerPeriod (sometimes called M)
  if (NPerPeriod < 2) {
    PyErr_SetString(PyExc_ValueError, "'n_per_period' must be >= 2");
    return NULL;
  }

  // Check magnet width
  if (MagnetWidth == 0) {
    MagnetWidth = Period / (double) NPerPeriod;
  } else if (MagnetWidth < 0) {
    PyErr_SetString(PyExc_ValueError, "'magnet_width' must be > 0");
    return NULL;
  } else if (MagnetWidth > Period / (double) NPerPeriod) {
    PyErr_SetString(PyExc_ValueError, "'magnet_width' * 'n_per_period' must be < 'period'");
    return NULL;
  }

  // Rotate field and sigma
  // UPDATE: check this
  //Field.RotateSelfXYZ(Rotations);
  //Period.RotateSelfXYZ(Rotations);


  // Add field
  self->obj->AddMagneticField( (TField*) new TField3D_Halbach(BField, Period, NPeriods, Gap, MagnetHeight, MagnetWidth, NPerPeriod, Name));

  // Must return python object None in a special way
  Py_INCREF(Py_None);
  return Py_None;
}







const char* DOC_OSCARSSR_AddMagneticFieldQuadrupole = R"docstring(
add_bfield_quadrupole(K, width [, rotations, translation, name])

Adds a quadrupole field in a given volume according to:

.. math::

    Bx = K * y

    By = K * x

    Bz = 0

For other directions or rotated skew quad fields one can use the *rotations* parameter

Parameters
----------
K : float
    Quadrupole strength given in [T/m]

width : float
    length in z direction (quadrupole axis direction)

rotations : list, optional
    3-element list representing rotations around x, y, and z axes: [:math:`\theta_x, \theta_y, \theta_z`]

translation : list, optional
    3-element list representing a translation in space [x, y, z]

name : str
    Name of this field

Returns
-------
None
)docstring";
static PyObject* OSCARSSR_AddMagneticFieldQuadrupole (OSCARSSRObject* self, PyObject* args, PyObject* keywds)
{
  // Add a magnetic field for undulator

  // UPDATE: Add axis [x, y, z] to be like others
  // Lists and variables
  double K = 0;
  double Width = 0;
  PyObject*   List_Rotations   = PyList_New(0);
  PyObject*   List_Translation = PyList_New(0);
  char const* Name             = "";

  TVector3D Rotations(0, 0, 0);
  TVector3D Translation(0, 0, 0);

  // Input variables and parsing
  static const char *kwlist[] = {"K",
                                 "width",
                                 "rotations",
                                 "translation",
                                 "name",
                                 NULL};

  if (!PyArg_ParseTupleAndKeywords(args, keywds, "dd|OOs",
                                   const_cast<char **>(kwlist),
                                   &K,
                                   &Width,
                                   &List_Rotations,
                                   &List_Translation,
                                   &Name
                                   )) {
    return NULL;
  }


  // Check for Rotations in the input
  if (PyList_Size(List_Rotations) != 0) {
    try {
      Rotations = OSCARSPY::ListAsTVector3D(List_Rotations);
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'rotations'");
      return NULL;
    }
  }

  // Check for Translation in the input
  if (PyList_Size(List_Translation) != 0) {
    try {
      Translation = OSCARSPY::ListAsTVector3D(List_Translation);
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'translation'");
      return NULL;
    }
  }

  // Name check
  if (std::string(Name).size() > 0 && Name[0] == '_') {
    PyErr_SetString(PyExc_ValueError, "'name' cannot begin with '_'.  This is reserved for internal use.  Please pick a different name");
    return NULL;
  }

  // Add field
  self->obj->AddMagneticField( (TField*) new TField3D_Quadrupole(K, Width, Rotations, Translation, Name));

  // Must return python object None in a special way
  Py_INCREF(Py_None);
  return Py_None;
}



const char* DOC_OSCARSSR_RemoveMagneticField = R"docstring(
remove_bfield(name)

Remove all fields with the given name

Parameters
----------
name : str
    Name of field to remove

Returns
-------
None
)docstring";
static PyObject* OSCARSSR_RemoveMagneticField (OSCARSSRObject* self, PyObject* args, PyObject* keywds)
{
  // Remove magnetic field

  char const* Name = "";

  // Input variables and parsing
  static const char *kwlist[] = {"name",
                                 NULL};

  if (!PyArg_ParseTupleAndKeywords(args, keywds, "s",
                                   const_cast<char **>(kwlist),
                                   &Name
                                   )) {
    return NULL;
  }

  // Remove fields with name name
  self->obj->RemoveMagneticField(Name);

  // Must return python object None in a special way
  Py_INCREF(Py_None);
  return Py_None;
}




const char* DOC_OSCARSSR_GetBField = R"docstring(
get_bfield(x [, t, name])

Get the 3D field at any point in space.  This is the sum of all fields added to this OSCARS object.

Parameters
----------
x : list
    A 3D list representing a point in space [x, y, z]

t : float
    time at which you wish to know the field.  optional, default is 0

name : str
    Name of the field you want.  All fields with this name will be summed.  Optional.

Returns
-------
bfield : [float, float, float]
    list representing the field [Fx, Fy, Fz]
)docstring";
static PyObject* OSCARSSR_GetBField (OSCARSSRObject* self, PyObject* args, PyObject* keywds)
{
  // Get the magnetic field at a point as a 3D list [Bx, By, Bz]

  // Python list object
  PyObject* List = 0x0;//PyList_New(0);
  double T = 0;
  char const* Name = "";

  static const char *kwlist[] = {"x",
                                 "t",
                                 "name",
                                 NULL};

  if (!PyArg_ParseTupleAndKeywords(args, keywds, "O|ds",
                                   const_cast<char **>(kwlist),
                                   &List,
                                   &T,
                                   &Name)) {
    return NULL;
  }

  // Grab the values
  TVector3D X;
  try {
    X = OSCARSPY::ListAsTVector3D(List);
  } catch (std::length_error e) {
    PyErr_SetString(PyExc_ValueError, "Incorrect format in input");
    return NULL;
  }

  // Set the object variable
  TVector3D const B = self->obj->GetB(X, T, Name);

  // Create a python list
  PyObject *PList = OSCARSPY::TVector3DAsList(B);

  // Return the python list
  return PList;
}




const char* DOC_OSCARSSR_ClearMagneticFields = R"docstring(
clear_bfields()

Remove all of the existing magnetic fields

Parameters
----------
None

Returns
-------
None
)docstring";
static PyObject* OSCARSSR_ClearMagneticFields (OSCARSSRObject* self)
{
  // Clear all magnetic fields in the OSCARSSR object
  self->obj->ClearMagneticFields();

  // Must return python object None in a special way
  Py_INCREF(Py_None);
  return Py_None;
}




const char* DOC_OSCARSSR_PrintMagneticFields = R"docstring(
print_bfields()

Print information about all magnetic fields to standard out

Parameters
----------
None

Returns
-------
None
)docstring";
static PyObject* OSCARSSR_PrintMagneticFields (OSCARSSRObject* self)
{
  // Print all magnetic stored in OSCARSSR

  // Out string stream for printing beam information
  std::ostringstream ostream;
  ostream << "*Magnetic Fields*\n";
  ostream << self->obj->GetBFieldContainer() << std::endl;

  OSCARSPY::PyPrint_stdout(ostream.str());

  // Must return python object None in a special way
  Py_INCREF(Py_None);
  return Py_None;
}













const char* DOC_OSCARSSR_AddElectricField = R"docstring(
add_efield_file(ifile, iformat [, rotations, translation, scale, name])

Add a electric field from a text file *ifile* according to the format *iformat*.

Currently OSCARS accepts the following file formats for iformat
    * 'OSCARS'
    * 'OSCARS1D [plus format string]'
    * 'SPECTRA'
    * 'SRW'

For the OSCARS1D format you must also include the order of the data columns, which would typically look something like: 'OSCARS1D Z Ex Ey Ez'.  You may use X or Y instead of Z and the order and number of E[xyz] does not matter.  This mode also accepts non-uniformly distributed data.

Optionally you can rotate and translate the field in space.  You can use an input *scale* list to scale the input (which must be in SI units of [m] for distances/positions and [V/m] for electric field values.

The rotation is performed first in the order: :math:`\theta_x, \theta_y, \theta_z`

*scale* is a list of less than or equal length to the number of elements in *iformat* when OSCARS1D is selected.  This will scale the input values of the i-th column before any rotation or translation.  This is useful if your data file is not in [V/m] and [m]

Parameters
----------
ifile : str
    Name of input file

iformat : str
    Format of the input file

rotations : list, optional
    3-element list representing rotations around x, y, and z axes: [:math:`\theta_x, \theta_y, \theta_z`]

translation : list, optional
    3-element list representing a translation in space [x, y, z]

scale : list, optional
    List of scale factors to be used for multiplying the inputs in the order of iformat (equal in length or less than the number of input parameters)

name : str
    Name of this field

Returns
-------
None

Examples
--------
Add a field from a data file where the data is in columns in the order Z Ex Ey Ez

    >>> osr.add_bfield_file(ifile='file.txt', iformat='OSCARS1D Z Bx By Bz')
)docstring";
static PyObject* OSCARSSR_AddElectricField (OSCARSSRObject* self, PyObject* args, PyObject* keywds)
{
  // Add a field from a file.
  // UPDATE: add binary file reading


  // Grab the values
  char const* FileName = "";
  char const* FileFormat = "";
  PyObject*   List_Rotations   = PyList_New(0);
  PyObject*   List_Translation = PyList_New(0);
  PyObject*   List_Scaling     = PyList_New(0);
  char const* Name             = "";

  TVector3D Rotations(0, 0, 0);
  TVector3D Translation(0, 0, 0);
  std::vector<double> Scaling;


  // Input variables and parsing
  static const char *kwlist[] = {"ifile",
                                 "iformat",
                                 "rotations",
                                 "translation",
                                 "scale",
                                 "name",
                                 NULL};

  if (!PyArg_ParseTupleAndKeywords(args, keywds, "ss|OOOs",
                                   const_cast<char **>(kwlist),
                                   &FileName,
                                   &FileFormat,
                                   &List_Rotations,
                                   &List_Translation,
                                   &List_Scaling,
                                   &Name)) {
    return NULL;
  }

  // Check that filename and format exist
  if (std::strlen(FileName) == 0 || std::strlen(FileFormat) == 0) {
    PyErr_SetString(PyExc_ValueError, "'ifile' or 'iformat' is blank");
    return NULL;
  }

  // Check for Rotations in the input
  if (PyList_Size(List_Rotations) != 0) {
    try {
      Rotations = OSCARSPY::ListAsTVector3D(List_Rotations);
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'rotations'");
      return NULL;
    }
  }


  // Check for Translation in the input
  if (PyList_Size(List_Translation) != 0) {
    try {
      Translation = OSCARSPY::ListAsTVector3D(List_Translation);
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'translation'");
      return NULL;
    }
  }

  // Get any scaling factors
  // UPDATE: Check against fileformat number of strings
  for (int i = 0; i < PyList_Size(List_Scaling); ++i) {
    Scaling.push_back(PyFloat_AsDouble(PyList_GetItem(List_Scaling, i)));
  }

  // Name check
  if (std::string(Name).size() > 0 && Name[0] == '_') {
    PyErr_SetString(PyExc_ValueError, "'name' cannot begin with '_'.  This is reserved for internal use.  Please pick a different name");
    return NULL;
  }

  // Add the magnetic field to the OSCARSSR object
  try {
    self->obj->AddElectricField(FileName, FileFormat, Rotations, Translation, Scaling, Name);
  } catch (...) {
    PyErr_SetString(PyExc_ValueError, "Could not import electric field.  Check 'ifile' and 'iformat' are correct");
    return NULL;
  }

  // Must return python object None in a special way
  Py_INCREF(Py_None);
  return Py_None;
}







const char* DOC_OSCARSSR_AddElectricFieldInterpolated = R"docstring(
add_efield_interpolated(mapping, iformat, parameter [, rotations, translation, scale, name])

Add field given a paramater and a mapping between known parameters and field data files where the field is interpolated from the known data points.

Parameters
----------
mapping : list [[float, str], [float, str], ...]
    List of parameters and associated filenames [[p1, file1], [p2, file2], ...]

iformat : str
    Which input format to use (see add_efield_file() for formats)

parameter : float
    Value of parameter you are interested in

rotations : list, optional
    3-element list representing rotations around x, y, and z axes: [:math:`\theta_x, \theta_y, \theta_z`]

translation : list, optional
    3-element list representing a translation in space [x, y, z]

scale : list, optional
    List of scale factors to be used for multiplying the inputs in the order of iformat (equal in length or less than the number of input parameters)

name : str
    Name of this field

Returns
-------
None

Examples
--------
See add_bfield_interpolated() for examples
)docstring";
static PyObject* OSCARSSR_AddElectricFieldInterpolated (OSCARSSRObject* self, PyObject* args, PyObject* keywds)
{
  // Add a magnetic field from a file.
  // UPDATE: add binary file reading

  // Grab the values
  PyObject*   List_Mapping     = PyList_New(0);
  char const* FileFormat       = "";
  double      Parameter        = 0;
  PyObject*   List_Rotations   = PyList_New(0);
  PyObject*   List_Translation = PyList_New(0);
  PyObject*   List_Scaling     = PyList_New(0);
  char const* Name             = "";

  TVector3D Rotations(0, 0, 0);
  TVector3D Translation(0, 0, 0);
  std::vector<double> Scaling;

  // Mapping that is passed in for field interpolation
  std::vector<std::pair<double, std::string> > Mapping;

  // Input variables and parsing
  static const char *kwlist[] = {"mapping",
                                 "iformat",
                                 "parameter",
                                 "rotations",
                                 "translation",
                                 "scale",
                                 "name",
                                 NULL};

  if (!PyArg_ParseTupleAndKeywords(args, keywds, "Osd|OOOs",
                                   const_cast<char **>(kwlist),
                                   &List_Mapping,
                                   &FileFormat,
                                   &Parameter,
                                   &List_Rotations,
                                   &List_Translation,
                                   &List_Scaling,
                                   &Name)) {
    return NULL;
  }

  // Grab mapping from input
  if (PyList_Size(List_Mapping) != 0) {
    double ParameterValue;
    for (int i = 0; i != PyList_Size(List_Mapping); ++i) {
      PyObject* ThisPair = PyList_GetItem(List_Mapping, i);
      if (PyList_Size(ThisPair) != 2) {
        PyErr_SetString(PyExc_ValueError, "Incorrect format in 'mapping'");
        return NULL;
      }

      ParameterValue = PyFloat_AsDouble(PyList_GetItem(ThisPair, 0));
      std::string const FileName = OSCARSPY::GetAsString(PyList_GetItem(ThisPair, 1));

      Mapping.push_back(std::make_pair(ParameterValue, FileName));
    }
  }


  // Check that filename and format exist
  if (std::strlen(FileFormat) == 0) {
    PyErr_SetString(PyExc_ValueError, "'iformat' is blank");
    return NULL;
  }

  // Parameter value is what it is...

  // Check for Rotations in the input
  if (PyList_Size(List_Rotations) != 0) {
    try {
      Rotations = OSCARSPY::ListAsTVector3D(List_Rotations);
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'rotations'");
      return NULL;
    }
  }


  // Check for Translation in the input
  if (PyList_Size(List_Translation) != 0) {
    try {
      Translation = OSCARSPY::ListAsTVector3D(List_Translation);
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'translation'");
      return NULL;
    }
  }

  // Get any scaling factors
  // UPDATE: Check against fileformat number of strings
  for (int i = 0; i < PyList_Size(List_Scaling); ++i) {
    Scaling.push_back(PyFloat_AsDouble(PyList_GetItem(List_Scaling, i)));
  }

  // Name check
  if (std::string(Name).size() > 0 && Name[0] == '_') {
    PyErr_SetString(PyExc_ValueError, "'name' cannot begin with '_'.  This is reserved for internal use.  Please pick a different name");
    return NULL;
  }

  // Add the magnetic field to the OSCARSSR object
  try {
    self->obj->AddElectricFieldInterpolated(Mapping, FileFormat, Parameter, Rotations, Translation, Scaling, Name);
  } catch (...) {
    PyErr_SetString(PyExc_ValueError, "Could not import magnetic field.  Check filenames and 'iformat' are correct");
    return NULL;
  }

  // Must return python object None in a special way
  Py_INCREF(Py_None);
  return Py_None;
}







const char* DOC_OSCARSSR_AddElectricFieldFunction = R"docstring(
add_efield_function(func [, name])

Adds field in the form of a user defined python function.  The input for this function must be (x, y, z, t) and return [Fx, Fy, Fz]

Parameters
----------

function : func
    A python function with input [x, y, z, t] and return of [Fx, Fy, Fz]

name : str
    Name of this field

Returns
-------
None

Examples
--------
Create a function in python and use it as a field in OSCARS

    >>> def myfunc(x, y, z, t):
    ...     "Do not forget to write a docstring"
    ...     if z > 0:
    ...         return 1
    ...     return 0
    >>> osr.add_efield_function(myfunc)
)docstring";
static PyObject* OSCARSSR_AddElectricFieldFunction (OSCARSSRObject* self, PyObject* args, PyObject* keywds)
{
  // Add a python function as an electric field object

  // Variables for function and name
  PyObject* Function;
  char const* Name = "";

  // Input variables and parsing
  static const char *kwlist[] = {"function",
                                 "name",
                                 NULL};

  if (!PyArg_ParseTupleAndKeywords(args, keywds, "O|s",
                                   const_cast<char **>(kwlist),
                                   &Function,
                                   &Name)) {
    return NULL;
  }


  // Increment ref to function for python
  Py_INCREF(Function);

  // Add the function as a field to the OSCARSSR object
  try {
    self->obj->AddElectricField( (TField*) new TFieldPythonFunction(Function, Name));
  } catch (std::invalid_argument e) {
    PyErr_SetString(PyExc_ValueError, e.what());
    return NULL;
  }

  // Decrement reference to function for python
  Py_DECREF(Function);

  // Must return python object None in a special way
  Py_INCREF(Py_None);
  return Py_None;
}







const char* DOC_OSCARSSR_AddElectricFieldGaussian = R"docstring(
add_efield_gaussian(efield, sigma [, rotations, translation])

Add a gaussian field in 3D with the peak field magnitude and direction given by *efield*, centered about a point with a given sigma in each coordinate.  If any component of *sigma* is less than or equal to zero it is ignored (ie. spatial extent is infinite).

The functional form for this is :math:`\exp(-(x - x_0)^2 / \sigma_x^2)`

Parameters
----------
efield : list
    A list representing the peak field [Fx, Fy, Fz]

sigma : list
    A list representing the sigma in each coordinate :math:`[\sigma_x, \sigma_y, \sigma_z]`

rotations : list, optional
    3-element list representing rotations around x, y, and z axes: [:math:`\theta_x, \theta_y, \theta_z`]

translation : list, optional
    3-element list representing a translation in space [x, y, z]

Returns
-------
None

Examples
--------
Add a field of 0.5 [V/m] in the X-direction centered at X=0, Y=0, Z=0 [m], with a sigma of 0.1 [m] in the Z-direction

    >>> osr.add_efield_gaussian(efield=[1, 0, 0], sigma=[0, 0, 0.10])

Add a field of 1 [V/m] in the Y-direction centered at X=1, Y=1, Z=1 [m], with a sigma of 0.05 in the Z-direction

    >>> osr.add_efield_gaussian(efield=[0, 1, 0], sigma=[0, 0, 0.05], translation=[1, 1, 1])
)docstring";
static PyObject* OSCARSSR_AddElectricFieldGaussian (OSCARSSRObject* self, PyObject* args, PyObject* keywds)
{
  // Add an electric field that is a gaussian

  // Lists and variables
  PyObject* List_Field       = PyList_New(0);
  PyObject* List_Translation  = PyList_New(0);
  PyObject* List_Rotations    = PyList_New(0);
  PyObject* List_Sigma        = PyList_New(0);
  const char* Name            = "";

  TVector3D Field(0, 0, 0);
  TVector3D Sigma(0, 0, 0);
  TVector3D Rotations(0, 0, 0);
  TVector3D Translation(0, 0, 0);


  // Input variables and parsing
  static const char *kwlist[] = {"efield",
                                 "sigma",
                                 "rotations",
                                 "translation",
                                 "name",
                                 NULL};

  if (!PyArg_ParseTupleAndKeywords(args, keywds, "OO|OOs",
                                   const_cast<char **>(kwlist),
                                   &List_Field,
                                   &List_Sigma,
                                   &List_Rotations,
                                   &List_Translation,
                                   &Name)) {
    return NULL;
  }



  // Check Field
  try {
    Field = OSCARSPY::ListAsTVector3D(List_Field);
  } catch (std::length_error e) {
    PyErr_SetString(PyExc_ValueError, "Incorrect format in 'efield'");
    return NULL;
  }

  // Check Width
  try {
    Sigma = OSCARSPY::ListAsTVector3D(List_Sigma);
  } catch (std::length_error e) {
    PyErr_SetString(PyExc_ValueError, "Incorrect format in 'sigma'");
    return NULL;
  }

  // Check for Rotations in the input
  if (PyList_Size(List_Rotations) != 0) {
    try {
      Rotations = OSCARSPY::ListAsTVector3D(List_Rotations);
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'rotations'");
      return NULL;
    }
  }

  // Check for Translation in the input
  if (PyList_Size(List_Translation) != 0) {
    try {
      Translation = OSCARSPY::ListAsTVector3D(List_Translation);
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'translation'");
      return NULL;
    }
  }

  // Name check
  if (std::string(Name).size() > 0 && Name[0] == '_') {
    PyErr_SetString(PyExc_ValueError, "'name' cannot begin with '_'.  This is reserved for internal use.  Please pick a different name");
    return NULL;
  }

  // Add field
  self->obj->AddElectricField( (TField*) new TField3D_Gaussian(Field, Translation, Sigma, Rotations, Name));

  // Must return python object None in a special way
  Py_INCREF(Py_None);
  return Py_None;
}







const char* DOC_OSCARSSR_AddElectricFieldUniform = R"docstring(
add_efield_uniform(efield [, width, rotations, translation, name])

Add a uniform field in a given range or for all space.  The *efield* is given as a 3D vector representing the field magnitude and direction.  *width* is an optional parameters, if not present the field permeates all space.  If a component of the 3D list *width* is less than or equal to zero, that coordinate will be ignored when calculating the field.

Parameters
----------
efield : list
    A list representing the field [Fx, Fy, Fz]

width : list
    A list representing the spetial extent of the field [Wx, Wy, Wz]

rotations : list, optional
    3-element list representing rotations around x, y, and z axes: [:math:`\theta_x, \theta_y, \theta_z`]

translation : list, optional
    3-element list representing a translation in space [x, y, z]

name : str
    Name of this field

Returns
-------
None

Examples
--------
Add a field of 0.0001 [V/m] in the X-direction for all space

    >>> osr.add_efield_uniform(efield=[0.0001, 0, 0])

Add a field of 0.0005 [V/m] in the Y-direction with a width in the Z-direction of 1.5 [m] (the other directions are ignored) centered at X=0, Y=0, Z=0.75 [m].

    >>> osr.add_efield_uniform(efield=[0, 1, 0], width=[0, 0, 1.5], translation=[0, 0, 0.75])

Add a field of 1 [V/m] in the X-direction in a volume of dimensions 1m x 1m x 1m.

    >>> osr.add_efield_uniform(efield=[1, 0, 0], width=[1, 1, 1])
)docstring";
static PyObject* OSCARSSR_AddElectricFieldUniform (OSCARSSRObject* self, PyObject* args, PyObject* keywds)
{
  // Add a uniform field with a given width in a given direction, or for all space

  // Lists and vectors
  PyObject*   List_Field       = PyList_New(0);
  PyObject*   List_Translation = PyList_New(0);
  PyObject*   List_Rotations   = PyList_New(0);
  PyObject*   List_Width       = PyList_New(0);
  char const* Name             = "";

  TVector3D Field(0, 0, 0);
  TVector3D Width (0, 0, 0);
  TVector3D Rotations(0, 0, 0);
  TVector3D Translation(0, 0, 0);


  // Input variables and parsing
  static const char *kwlist[] = {"efield",
                                 "width",
                                 "rotations",
                                 "translation",
                                 "name",
                                 NULL};

  if (!PyArg_ParseTupleAndKeywords(args, keywds, "O|OOOs",
                                   const_cast<char **>(kwlist),
                                   &List_Field,
                                   &List_Width,
                                   &List_Rotations,
                                   &List_Translation,
                                   &Name
                                   )) {
    return NULL;
  }


  // Check Field
  try {
    Field = OSCARSPY::ListAsTVector3D(List_Field);
  } catch (std::length_error e) {
    PyErr_SetString(PyExc_ValueError, "Incorrect format in 'efield'");
    return NULL;
  }

  // Check Width
  if (PyList_Size(List_Width) != 0) {
    try {
      Width = OSCARSPY::ListAsTVector3D(List_Width);
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'width'");
      return NULL;
    }
  }

  // Check for Rotations in the input
  if (PyList_Size(List_Rotations) != 0) {
    try {
      Rotations = OSCARSPY::ListAsTVector3D(List_Rotations);
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'rotations'");
      return NULL;
    }
  }

  // Check for Translation in the input
  if (PyList_Size(List_Translation) != 0) {
    try {
      Translation = OSCARSPY::ListAsTVector3D(List_Translation);
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'translation'");
      return NULL;
    }
  }

  // Name check
  if (std::string(Name).size() > 0 && Name[0] == '_') {
    PyErr_SetString(PyExc_ValueError, "'name' cannot begin with '_'.  This is reserved for internal use.  Please pick a different name");
    return NULL;
  }

  // Add the field
  self->obj->AddElectricField((TField*) new TField3D_UniformBox(Field, Width, Translation, Rotations, Name));

  // Must return python object None in a special way
  Py_INCREF(Py_None);
  return Py_None;
}










const char* DOC_OSCARSSR_AddElectricFieldIdealUndulator = R"docstring(
add_efield_undulator(efield, period, nperiods [, phase, rotations, translation, taper])

Adds an ideal sinusoidal undulator field with a given maximum efield amplitude, period, and number of periods.  Optionally one can specify the phase offset (in [rad]), rotations and translation.  The number of periods given is the full number of fields not counting the terminating fields.

Parameters
----------
efield : list
    A list representing the peak field [Fx, Fy, Fz] in [V/m]

period : list
    Length of one period

nperiods : int
    Number of periods

phase : float
    Phase offset in [rad]

rotations : list, optional
    3-element list representing rotations around x, y, and z axes: [:math:`\theta_x, \theta_y, \theta_z`]

translation : list, optional
    3-element list representing a translation in space [x, y, z]

taper : float
    The fractional change of the efield per meter along the magnetic axis

name : str
    Name of this field

Returns
-------
None

Examples
--------
Add an idealized undulator with 41 periods having a period of 0.050 [m] with a maximum field of 1 [T] in the y-direction where the magnetic axis is along the z-axis

    >>> osr.add_efield_undulator(efield=[0, 1, 0], period=[0, 0, 0.050], nperiods=41)
)docstring";
static PyObject* OSCARSSR_AddElectricFieldIdealUndulator (OSCARSSRObject* self, PyObject* args, PyObject* keywds)
{
  // Add an electric field undulator to OSCARSSR

  // Lists and variables
  PyObject*   List_Field       = PyList_New(0);
  PyObject*   List_Period      = PyList_New(0);
  PyObject*   List_Rotations   = PyList_New(0);
  PyObject*   List_Translation = PyList_New(0);
  int         NPeriods         = 0;
  double      Phase            = 0;
  double      Taper            = 0;
  char const* Name             = "";

  TVector3D Field(0, 0, 0);
  TVector3D Period(0, 0, 0);
  TVector3D Rotations(0, 0, 0);
  TVector3D Translation(0, 0, 0);

  // Input variables and parsing
  static const char *kwlist[] = {"efield",
                                 "period",
                                 "nperiods",
                                 "phase",
                                 "rotations",
                                 "translation",
                                 "taper",
                                 "name",
                                 NULL};

  if (!PyArg_ParseTupleAndKeywords(args, keywds, "OOi|dOOds",
                                   const_cast<char **>(kwlist),
                                   &List_Field,
                                   &List_Period,
                                   &NPeriods,
                                   &Phase,
                                   &List_Rotations,
                                   &List_Translation,
                                   &Taper,
                                   &Name
                                   )) {
    return NULL;
  }


  // Check Field
  try {
    Field = OSCARSPY::ListAsTVector3D(List_Field);
  } catch (std::length_error e) {
    PyErr_SetString(PyExc_ValueError, "Incorrect format in 'bfield'");
    return NULL;
  }

  // Check Period
  try {
    Period = OSCARSPY::ListAsTVector3D(List_Period);
  } catch (std::length_error e) {
    PyErr_SetString(PyExc_ValueError, "Incorrect format in 'period'");
    return NULL;
  }

  // Check for Rotations in the input
  if (PyList_Size(List_Rotations) != 0) {
    try {
      Rotations = OSCARSPY::ListAsTVector3D(List_Rotations);
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'rotations'");
      return NULL;
    }
  }

  // Check for Translation in the input
  if (PyList_Size(List_Translation) != 0) {
    try {
      Translation = OSCARSPY::ListAsTVector3D(List_Translation);
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'translation'");
      return NULL;
    }
  }

  // Name check
  if (std::string(Name).size() > 0 && Name[0] == '_') {
    PyErr_SetString(PyExc_ValueError, "'name' cannot begin with '_'.  This is reserved for internal use.  Please pick a different name");
    return NULL;
  }


  // Rotate field and sigma
  // UPDATE: check this
  Field.RotateSelfXYZ(Rotations);
  Period.RotateSelfXYZ(Rotations);


  // Add field
  self->obj->AddElectricField( (TField*) new TField3D_IdealUndulator(Field, Period, NPeriods, Translation, Phase, Taper, Name));

  // Must return python object None in a special way
  Py_INCREF(Py_None);
  return Py_None;
}







const char* DOC_OSCARSSR_RemoveElectricField = R"docstring(
remove_efield(name)

Remove all fields with the given name

Parameters
----------
name : str
    Name of field to remove

Returns
-------
None
)docstring";
static PyObject* OSCARSSR_RemoveElectricField (OSCARSSRObject* self, PyObject* args, PyObject* keywds)
{
  // Remove magnetic field

  char const* Name             = "";

  // Input variables and parsing
  static const char *kwlist[] = {"name",
                                 NULL};

  if (!PyArg_ParseTupleAndKeywords(args, keywds, "s",
                                   const_cast<char **>(kwlist),
                                   &Name
                                   )) {
    return NULL;
  }

  // Remove fields with name name
  self->obj->RemoveElectricField(Name);

  // Must return python object None in a special way
  Py_INCREF(Py_None);
  return Py_None;
}




const char* DOC_OSCARSSR_GetEField = R"docstring(
get_efield(x [, t, name])

Get the 3D field at any point in space.  This is the sum of all fields added to this OSCARS object.

Parameters
----------
x : list
    A 3D list representing a point in space [x, y, z]

t : float
    time at which you wish to know the field.  optional, default is 0

name : str
    Name of the field you want.  All fields with this name will be summed.  Optional.


Returns
-------
bfield : [float, float, float]
    list representing the field [Fx, Fy, Fz]
)docstring";
static PyObject* OSCARSSR_GetEField (OSCARSSRObject* self, PyObject* args, PyObject* keywds)
{
  // Get the magnetic field at a point as a 3D list [Ex, Ey, Ez]

  // Python list object
  PyObject* List = 0x0;
  double T = 0;
  char const* Name = "";

  static const char *kwlist[] = {"x",
                                 "t",
                                 "name",
                                 NULL};

  if (!PyArg_ParseTupleAndKeywords(args, keywds, "O|ds",
                                   const_cast<char **>(kwlist),
                                   &List,
                                   &T,
                                   &Name)) {
    return NULL;
  }

  // Grab the values
  TVector3D X;
  try {
    X = OSCARSPY::ListAsTVector3D(List);
  } catch (std::length_error e) {
    PyErr_SetString(PyExc_ValueError, "Incorrect format in input");
    return NULL;
  }

  // Set the object variable
  TVector3D const F = self->obj->GetE(X, T, Name);

  // Create a python list
  PyObject *PList = OSCARSPY::TVector3DAsList(F);

  // Return the python list
  return PList;
}




const char* DOC_OSCARSSR_ClearElectricFields = R"docstring(
clear_efields()

Remove all of the existing electric fields.

Parameters
----------
None

Returns
-------
None
)docstring";
static PyObject* OSCARSSR_ClearElectricFields (OSCARSSRObject* self)
{
  // Clear all magnetic fields in the OSCARSSR object
  self->obj->ClearElectricFields();

  // Must return python object None in a special way
  Py_INCREF(Py_None);
  return Py_None;
}




const char* DOC_OSCARSSR_PrintElectricFields = R"docstring(
print_efields()

Print information about all electric fields to the standard out

Parameters
----------
None

Returns
-------
None
)docstring";
static PyObject* OSCARSSR_PrintElectricFields (OSCARSSRObject* self)
{
  // Print all magnetic stored in OSCARSSR

  // Out string stream for printing beam information
  std::ostringstream ostream;
  ostream << "*Electric Fields*\n";
  ostream << self->obj->GetEFieldContainer() << std::endl;

  OSCARSPY::PyPrint_stdout(ostream.str());

  // Must return python object None in a special way
  Py_INCREF(Py_None);
  return Py_None;
}
















const char* DOC_OSCARSSR_WriteMagneticField = R"docstring(
write_bfield(oformat [, ofile, bofile, xlim, nx, ylim, ny, zlim, nz, comment, version])

Write the field to ofile in the format described in oformat.  You must specify at least one limits and one number of points (e.g. zlim and nz).  All output is in SI units.

The formats (*oformat*) available are:
* OSCARS
* OSCARS1D
* SPECTRA
* SRW

For OSCARS1D you also must specify what you want in the output.  One spatial dimension must be specified along with at least one field dimension (in any order you like).  For examples theses are all valid: 'OSCARS1D Z Fx Fy Fz', 'OSCARS1D Fy Fx Z Fz'.

Parameters
----------

oformat : str
    Format of output file

ofile : str
    Name of output file

bofile : str
    Name of binary output file

xlim : list
    [min, max] for x dimension

ylim : list
    [min, max] for y dimension

zlim : list
    [min, max] for z dimension

comment : str
    Comment string to be added to file header.  LF and CR are removed.

version : int
    Which version of output format (you should use the default unless you have good reason)

Returns
-------
None
)docstring";
static PyObject* OSCARSSR_WriteMagneticField (OSCARSSRObject* self, PyObject* args, PyObject* keywds)
{
  // Add a magnetic field that is a gaussian

  const char* OutFileName       = "";
  const char* OutFileNameBinary = "";
  const char* OutFormat         = "";
  const char* Comment           = "";

  // Lists and variables
  PyObject* List_XLim = PyList_New(0);
  PyObject* List_YLim = PyList_New(0);
  PyObject* List_ZLim = PyList_New(0);

  int NX = 0;
  int NY = 0;
  int NZ = 0;

  TVector2D XLim;
  TVector2D YLim;
  TVector2D ZLim;

  int Version = 0;


  // Input variables and parsing
  static const char *kwlist[] = {"oformat",
                                 "ofile",
                                 "bofile",
                                 "xlim",
                                 "nx",
                                 "ylim",
                                 "ny",
                                 "zlim",
                                 "nz",
                                 "comment",
                                 "version",
                                 NULL};

  if (!PyArg_ParseTupleAndKeywords(args, keywds, "s|ssOiOiOisi",
                                   const_cast<char **>(kwlist),
                                   &OutFormat,
                                   &OutFileName,
                                   &OutFileNameBinary,
                                   &List_XLim,
                                   &NX,
                                   &List_YLim,
                                   &NY,
                                   &List_ZLim,
                                   &NZ,
                                   &Comment,
                                   &Version)) {
    return NULL;
  }

  // Initial values for limits are all 0
  XLim.SetXY(0, 0);
  YLim.SetXY(0, 0);
  ZLim.SetXY(0, 0);


  // Check that filename and format exist
  if (std::strlen(OutFormat) == 0) {
    PyErr_SetString(PyExc_ValueError, "'oformat' is blank");
    return NULL;
  }

  // Check for XLim in the input
  if (PyList_Size(List_XLim) != 0) {
    try {
      XLim = OSCARSPY::ListAsTVector2D(List_XLim);
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'xlim'");
      return NULL;
    }
  }

  // Check for YLim in the input
  if (PyList_Size(List_YLim) != 0) {
    try {
      YLim = OSCARSPY::ListAsTVector2D(List_YLim);
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'ylim'");
      return NULL;
    }
  }

  // Check for ZLim in the input
  if (PyList_Size(List_ZLim) != 0) {
    try {
      ZLim = OSCARSPY::ListAsTVector2D(List_ZLim);
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'zlim'");
      return NULL;
    }
  }

  if (std::strlen(OutFileName) > 0) {
    try {
      self->obj->WriteField("B", OutFileName, OutFormat, XLim, NX, YLim, NY, ZLim, NZ, Comment);
    } catch (...) {
      PyErr_SetString(PyExc_ValueError, "could not write output file");
      return NULL;
    }
  }

  if (std::strlen(OutFileNameBinary) > 0) {
    try {
      self->obj->WriteFieldBinary("B", OutFileNameBinary, OutFormat, XLim, NX, YLim, NY, ZLim, NZ, Comment, Version);
    } catch (...) {
      PyErr_SetString(PyExc_ValueError, "could not write output file");
      return NULL;
    }
  }

  // Must return python object None in a special way
  Py_INCREF(Py_None);
  return Py_None;
}




const char* DOC_OSCARSSR_WriteElectricField = R"docstring(
write_efield(oformat [, ofile, bofile, xlim, nx, ylim, ny, zlim, nz, comment, version])

Write the field to ofile in the format described in oformat.  You must specify at least one limits and one number of points (e.g. zlim and nz).  All output is in SI units.

The formats (*oformat*) available are:
* OSCARS
* OSCARS1D
* SPECTRA
* SRW

For OSCARS1D you also must specify what you want in the output.  One spatial dimension must be specified along with at least one field dimension (in any order you like).  For examples theses are all valid: 'OSCARS1D Z Fx Fy Fz', 'OSCARS1D Fy Fx Z Fz'.

Parameters
----------

oformat : str
    Format of output file

ofile : str
    Name of output file

bofile : str
    Name of binary output file

xlim : list
    [min, max] for x dimension

ylim : list
    [min, max] for y dimension

zlim : list
    [min, max] for z dimension

comment : str
    Comment string to be added to file header.  LF and CR are removed.

version : int
    Which version of output format (you should use the default unless you have good reason)

Returns
-------
None
)docstring";
static PyObject* OSCARSSR_WriteElectricField (OSCARSSRObject* self, PyObject* args, PyObject* keywds)
{
  // Add a magnetic field that is a gaussian

  const char* OutFileName       = "";
  const char* OutFileNameBinary = "";
  const char* OutFormat         = "";
  const char* Comment           = "";

  // Lists and variables
  PyObject* List_XLim = PyList_New(0);
  PyObject* List_YLim = PyList_New(0);
  PyObject* List_ZLim = PyList_New(0);

  int NX = 0;
  int NY = 0;
  int NZ = 0;

  TVector2D XLim;
  TVector2D YLim;
  TVector2D ZLim;

  int Version = 0;


  // Input variables and parsing
  static const char *kwlist[] = {"oformat",
                                 "ofile",
                                 "bofile",
                                 "xlim",
                                 "nx",
                                 "ylim",
                                 "ny",
                                 "zlim",
                                 "nz",
                                 "comment",
                                 "version",
                                 NULL};

  if (!PyArg_ParseTupleAndKeywords(args, keywds, "s|ssOiOiOisi",
                                   const_cast<char **>(kwlist),
                                   &OutFormat,
                                   &OutFileName,
                                   &OutFileNameBinary,
                                   &List_XLim,
                                   &NX,
                                   &List_YLim,
                                   &NY,
                                   &List_ZLim,
                                   &NZ,
                                   &Comment,
                                   &Version)) {
    return NULL;
  }

  // Initial values for limits are all 0
  XLim.SetXY(0, 0);
  YLim.SetXY(0, 0);
  ZLim.SetXY(0, 0);


  // Check that filename and format exist
  if (std::strlen(OutFormat) == 0) {
    PyErr_SetString(PyExc_ValueError, "'oformat' is blank");
    return NULL;
  }

  // Check for XLim in the input
  if (PyList_Size(List_XLim) != 0) {
    try {
      XLim = OSCARSPY::ListAsTVector2D(List_XLim);
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'xlim'");
      return NULL;
    }
  }

  // Check for YLim in the input
  if (PyList_Size(List_YLim) != 0) {
    try {
      YLim = OSCARSPY::ListAsTVector2D(List_YLim);
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'ylim'");
      return NULL;
    }
  }

  // Check for ZLim in the input
  if (PyList_Size(List_ZLim) != 0) {
    try {
      ZLim = OSCARSPY::ListAsTVector2D(List_ZLim);
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'zlim'");
      return NULL;
    }
  }

  if (std::strlen(OutFileName) > 0) {
    try {
      self->obj->WriteField("E", OutFileName, OutFormat, XLim, NX, YLim, NY, ZLim, NZ, Comment);
    } catch (...) {
      PyErr_SetString(PyExc_ValueError, "could not write output file");
      return NULL;
    }
  }

  if (std::strlen(OutFileNameBinary) > 0) {
    try {
      self->obj->WriteFieldBinary("E", OutFileNameBinary, OutFormat, XLim, NX, YLim, NY, ZLim, NZ, Comment, Version);
    } catch (...) {
      PyErr_SetString(PyExc_ValueError, "could not write output file");
      return NULL;
    }
  }

  // Must return python object None in a special way
  Py_INCREF(Py_None);
  return Py_None;
}
























































const char* DOC_OSCARSSR_SetParticleBeam = R"docstring(
set_particle_beam([, type, name, energy_GeV, d0, x0, beam, sigma_energy_GeV, t0, current, weight, rotations, translation, horizontal_direction, beta, alpha, gamma, emittance, eta, lattice_reference, mass, charge])

This function is the same as add_particle_beam(), but it clears all particle beams before the 'add'.
)docstring";
static PyObject* OSCARSSR_SetParticleBeam (OSCARSSRObject* self, PyObject* args, PyObject* keywds)
{
  // Clear all particle beams, add this beam, and set a new particle

  self->obj->ClearParticleBeams();

  PyObject* ret = OSCARSSR_AddParticleBeam(self, args, keywds);
  if (ret == NULL) {
    return ret;
  }

  self->obj->SetNewParticle("", "ideal");

  return ret;
}




const char* DOC_OSCARSSR_AddParticleBeam = R"docstring(
add_particle_beam([, type, name, energy_GeV, d0, x0, beam, sigma_energy_GeV, t0, current, weight, rotations, translation, horizontal_direction, beta, alpha, gamma, emittance, eta, lattice_reference, mass, charge])

Add a particle beam to the OSCARS object with a name given by *name*.  There is no limit to the number of different particle beams one can add.  They are added with a *weight* which is by default 1.  The weight is used in random sampling when asking for a new particle, for example in oscars.sr.set_new_particle().  If the *beam* parameter is given you only need to specify *name* and *x0*.

Supported particle types for *type* are:
    * electron
    * positron
    * muon
    * anti-muon
    * proton
    * anti-proton
    * pi+
    * pi-

Parameters
----------

type : str
    One of the built-in types of particle beams, or 'custom'.  If you use custom you must also specify *mass* and *charge*.

name : str
    User identified of this beam

energy_GeV : float
    Beam energy in [GeV]

d0 : list
    Vector [float, float, float] representing the default direction of the beam at the initial position.  The normalization of this vector does not matter.

x0 : list
    Coordinates of the initial position in [m] [x, y, z]

beam : str
    Name of predefined beam

sigma_energy_GeV : float
    Beam energy spread in [GeV]

t0 : float
    Initial time in [m] at the initial_position.  Time here is in terms of ct.

current : float
    Beam current in [A].  If this parameter is 0, the current defaults to the single particle charge

weight : float
    Weight to give this beam for random sampling when there are multiple beams

rotations : list, optional
    3-element list representing rotations around x, y, and z axes: [:math:`\theta_x, \theta_y, \theta_z`]

translation : list, optional
    3-element list representing a translation in space [x, y, z]

horizontal_direction : list
    A list representing the *horizontal* beam direction.  This must be perpendicular to the beam *direction*.  The vertical direction is defined internally as :math:`[direction] \times [horizontal\_direction]`.

beta : list
    Values of the horizontal and vertical beta funtion at the point *lattice_center* [beta_x, beta_y]

alpha : list
    Values of the horizontal and vertical alpha funtion at the point *lattice_center* [alpha_x, alpha_y]

gamma : list
    Values of the horizontal and vertical gamma funtion at the point *lattice_center* [gamma_x, gamma_y]

emittance : list
    values of the horizontal and vertical emittance [emittance_x, emittance_y]

eta : list
    Values of the horizontal and vertical dispersion at the point *lattice_center* [eta_x, eta_y].  Currently eta_y is ignored and dispersion is only used in oscars.th calculations

lattice_reference : list
    Coordinates of the lattice center [x, y, z] (must be on-axis with respect to the beam)

mass : float
    mass of a *custom* particle.  This is only used if *type* = 'custom'.  Must be non-zero.

charge : float
    Charge of a *custom* particle.  This is only used if *type* = 'custom'.  Must be non-zero.

Returns
-------
None

Currently the predefined beam types are:
* NSLSII
* NSLSII-ShortStraight
* NSLSII-LongStraight

Examples
--------
Add an electron beam with 0.500 [A] current at an initial position of [0, 0, 0] in the Y direction with energy of 3 [GeV]

    >>> osr.add_particle_beam(type='electron', name='beam_0', x0=[0, 0, 0], d0=[0, 1, 0], energy_GeV=3, current=0.500)

Add a positron beam with 0.500 [A] current at an initial position of [-2, 0, 0] in the direction given by theta in the X-Y plane with energy of 3 [GeV]

    >>> from math import sin, cos
    >>> theta = 0.25 * osr.pi()
    >>> osr.add_particle_beam(type='positron', name='beam_0', x0=[-2, 0, 0], d0=[sin(theta), cos(theta), 0], energy_GeV=3, current=0.500)

Add a predefined beam for the NSLSII short straight section

    >>> osr.add_particle_beam(beam='NSLSII-ShortStraight', name='beam_0', x0=[-2, 0, 0])
)docstring";
static PyObject* OSCARSSR_AddParticleBeam (OSCARSSRObject* self, PyObject* args, PyObject* keywds)
{
  // Add a particle beam to the experiment

  // Lists and variables some with initial values
  char const* Type                       = "electron";
  char const* Name                       = "";
  double      Energy_GeV                 = -1;
  double      Sigma_Energy_GeV           = 0;
  double      T0                         = 0;
  double      Current                    = 0;
  double      Weight                     = 1;
  double      Mass                       = 0;
  double      Charge                     = 0;
  char const* Beam                       = "";
  PyObject*   List_Position              = 0x0;
  PyObject*   List_Direction             = 0x0;
  PyObject*   List_Rotations             = 0x0;
  PyObject*   List_Translation           = 0x0;
  PyObject*   List_Horizontal_Direction  = 0x0;
  PyObject*   List_Beta                  = 0x0;
  PyObject*   List_Alpha                 = 0x0;
  PyObject*   List_Gamma                 = 0x0;
  PyObject*   List_Emittance             = 0x0;
  PyObject*   List_Eta                   = 0x0;
  PyObject*   List_Lattice_Reference     = 0x0;

  TVector3D Position(0, 0, 0);
  TVector3D Direction(0, 0, 1);
  TVector3D Rotations(0, 0, 0);
  TVector3D Translation(0, 0, 0);
  TVector3D Horizontal_Direction;
  TVector2D Beta(0, 0);
  TVector2D Alpha(0, 0);
  TVector2D Gamma(0, 0);
  TVector2D Emittance(0, 0);
  TVector3D Lattice_Reference(0, 0, 0);
  TVector2D Eta(0, 0);


  // Input variables and parsing
  static const char *kwlist[] = {"type",
                                 "name",
                                 "energy_GeV",
                                 "d0",
                                 "x0",
                                 "beam",
                                 "sigma_energy_GeV",
                                 "t0",
                                 "current",
                                 "weight",
                                 "rotations",
                                 "translation",
                                 "horizontal_direction",
                                 "beta",
                                 "alpha",
                                 "gamma",
                                 "emittance",
                                 "eta",
                                 "lattice_reference",
                                 "mass",
                                 "charge",
                                 NULL};

  if (!PyArg_ParseTupleAndKeywords(args, keywds, "|ssdOOsddddOOOOOOOOOdd",
                                   const_cast<char **>(kwlist),
                                   &Type,
                                   &Name,
                                   &Energy_GeV,
                                   &List_Direction,
                                   &List_Position,
                                   &Beam,
                                   &Sigma_Energy_GeV,
                                   &T0,
                                   &Current,
                                   &Weight,
                                   &List_Rotations,
                                   &List_Translation,
                                   &List_Horizontal_Direction,
                                   &List_Beta,
                                   &List_Alpha,
                                   &List_Gamma,
                                   &List_Emittance,
                                   &List_Eta,
                                   &List_Lattice_Reference,
                                   &Mass,
                                   &Charge)) {
    return NULL;
  }


  // Are you asking for one of the predefined beams?
  bool const HasPredefinedBeam = std::strlen(Beam) != 0 ? true : false;

  TParticleBeam* ThisBeam = 0x0;

  // Check if beam is defined (for predefined beams)
  if (HasPredefinedBeam) {
    try {
      ThisBeam = &(self->obj->AddParticleBeam(Beam, Name, Weight));
    } catch (...) {
      PyErr_SetString(PyExc_ValueError, "Error in predefined beam name / definition");
      return NULL;
    }

    if (Energy_GeV >= 0) {
      ThisBeam->SetE0(Energy_GeV);
    }

  }


  // default is 0
  if (Energy_GeV == -1) {
    Energy_GeV = 0;
  }

  // Initial position
  if (List_Position != 0x0) {
    try {
      Position = OSCARSPY::ListAsTVector3D(List_Position);
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'x0'");
      return NULL;
    }

  }

  // Energy check
  if (Energy_GeV < 0) {
    PyErr_SetString(PyExc_ValueError, "We do not currently support negative energy beams.  Please change 'energy_GeV' to >= 0");
    return NULL;
  }

  // Initial direction
  if (List_Direction != 0x0) {
    try {
      Direction = OSCARSPY::ListAsTVector3D(List_Direction);
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'd0'");
      return NULL;
    }
  }

  // Check for Rotations in the input
  if (List_Rotations != 0x0) {
    try {
      Rotations = OSCARSPY::ListAsTVector3D(List_Rotations);
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'rotations'");
      return NULL;
    }
  }


  // Check for Translation in the input
  if (List_Translation != 0x0) {
    try {
      Translation = OSCARSPY::ListAsTVector3D(List_Translation);
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'translation'");
      return NULL;
    }
  }


  // Check for Horizontal_Direction in the input
  if (List_Horizontal_Direction != 0x0) {
    try {
      Horizontal_Direction = OSCARSPY::ListAsTVector3D(List_Horizontal_Direction);
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'horizontal_direction'");
      return NULL;
    }
  } else {
    Horizontal_Direction = -Direction.Orthogonal();
  }
  Horizontal_Direction = Horizontal_Direction.UnitVector();


  // Rotate beam parameters
  Position.RotateSelfXYZ(Rotations);
  Direction.RotateSelfXYZ(Rotations);
  Position += Translation;

  // Add the particle beam
  if (std::strlen(Beam) == 0) {
    try {
      if (std::string(Type) == "custom") {
        if (Mass == 0 || Charge == 0) {
          PyErr_SetString(PyExc_ValueError, "'mass' or 'charge' is zero");
          return NULL;
        }
        // UPDATE: for custom beams
        ThisBeam = &(self->obj->AddParticleBeam(Type, Name, Position, Direction, Energy_GeV, T0, Current, Weight, Charge, Mass));
      } else {
        ThisBeam = &(self->obj->AddParticleBeam(Type, Name, Position, Direction, Energy_GeV, T0, Current, Weight));
      }
    } catch (std::invalid_argument e) {
      PyErr_SetString(PyExc_ValueError, "invalid argument in adding particle beam.  possibly 'name' already exists");
      return NULL;
    }
  }

  // Change predefined beam accordingly
  if (HasPredefinedBeam) {
    if (List_Direction != 0x0) {
      ThisBeam->SetU0(Direction);
    }
    ThisBeam->SetX0(Position);
  }

  // Set horizontal direction
  ThisBeam->SetHorizontalDirection(Horizontal_Direction);

  // UPDATE
  // Check for Emittance in the input
  if (List_Emittance != 0x0) {
    try {
      Emittance = OSCARSPY::ListAsTVector2D(List_Emittance);
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'emittance'");
      return NULL;
    }

    // Set emittance and beam distribution
    ThisBeam->SetEmittance(Emittance);
    ThisBeam->SetBeamDistribution(TParticleBeam::kBeamDistribution_Gaussian);
  } else {
    // Beam distribution to filament
    if (std::strlen(Beam) == 0) {
      ThisBeam->SetBeamDistribution(TParticleBeam::kBeamDistribution_Filament);
    }
  }

  // Check for no beam
  if (ThisBeam == 0x0) {
    // UPDATE: add specific throw
    std::cerr << "ERROR: No beam at checkpoint" << std::endl;
    throw;
  }

  if (Sigma_Energy_GeV == 0) {
    // Do nothing.  zero energy diff is alright
  } else if (Sigma_Energy_GeV < 0) {
    PyErr_SetString(PyExc_ValueError, "'sigma_energy_GeV' cannot be less than zero");
    return NULL;
  } else {
    ThisBeam->SetSigmaEnergyGeV(Sigma_Energy_GeV);
  }

  // Check for beta, alpha, gammain the input
  int HasBAG = 0x0;
  if (List_Beta != 0x0) {
    try {
      Beta = OSCARSPY::ListAsTVector2D(List_Beta);
      HasBAG |= 0x4;
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'beta'");
      return NULL;
    }
  }
  if (List_Alpha != 0x0) {
    try {
      Alpha = OSCARSPY::ListAsTVector2D(List_Alpha);
      HasBAG |= 0x2;
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'alpha'");
      return NULL;
    }
  }
  if (List_Gamma != 0x0) {
    try {
      Gamma = OSCARSPY::ListAsTVector2D(List_Gamma);
      HasBAG |= 0x1;
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'gamma'");
      return NULL;
    }
  }

  // Check for Lattice reference in the input
  bool HasReferencePoint = false;
  if (List_Lattice_Reference != 0x0) {
    try {
      Lattice_Reference = OSCARSPY::ListAsTVector3D(List_Lattice_Reference);
      HasReferencePoint = true;
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'lattice_reference'");
      return NULL;
    }
  }
  ThisBeam->SetTwissLatticeReference(Lattice_Reference);

  // Set correct twiss parameters
  switch (HasBAG) {
    case 0x7:
      ThisBeam->SetTwissParameters(Beta, Alpha, Gamma, Lattice_Reference, HasReferencePoint);
      break;
    case 0x6:
      ThisBeam->SetTwissBetaAlpha(Beta, Alpha, Lattice_Reference, HasReferencePoint);
      break;
    case 0x5:
      ThisBeam->SetTwissBetaGamma(Beta, Gamma, Lattice_Reference, HasReferencePoint);
      break;
    case 0x3:
      ThisBeam->SetTwissAlphaGamma(Alpha, Gamma, Lattice_Reference, HasReferencePoint);
      break;
    case 0x4:
      ThisBeam->SetTwissBetaAlpha(Beta, TVector2D(0, 0), Lattice_Reference, HasReferencePoint);
      break;
    default:
      break;
  }

  if (List_Eta != 0x0) {
    try {
      Eta = OSCARSPY::ListAsTVector2D(List_Eta);
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'gamma'");
      return NULL;
    }
  }
  ThisBeam->SetEta(Eta);


  if (T0 != 0) {
    ThisBeam->SetT0(T0);
  }

  // Should set weight on input
  //self->obj->GetParticleBeam(Name).SetWeight(Weight);

  // Must return python object None in a special way
  Py_INCREF(Py_None);
  return Py_None;
}





const char* DOC_OSCARSSR_SetParticleBeamSize = R"docstring(
add_particle_beam_size(sigma, sigmap [, beam, lattice_reference, eta=0])

Set the beamsize of a particular beam (if the name is given) or all beams (if no name is given) at a reference point called lattice_reference.  If lattice_reference is not given the reference point is assumed to be 'x0' from the beam in question.  This function translates the sizes into twiss parameters.  If eta (dispersion) is not given it is assumed to be zero.

Parameters
----------

sigma : list[2]
    One sigma for the spatial distribution in the horizontal and vertical directions [sigma_h, sigma_v].

sigmap : list[2]
    One sigma for the angular distribution in the horizontal and vertical directions [sigmap_h, sigmap_v].

beam : str
    Name of the beam for which these paramters apply.

lattice_reference : list[3]
    3D coordinate of there reference point for the known beam dimensions.

eta : float
    The dispersion term.  Assumed to be zero if not given.

Returns
-------
None

)docstring";
static PyObject* OSCARSSR_SetParticleBeamSize (OSCARSSRObject* self, PyObject* args, PyObject* keywds)
{
  // Add a particle beam to the experiment

  // Lists and variables some with initial values
  PyObject*   List_Sigma                 = PyList_New(0);
  PyObject*   List_SigmaP                = PyList_New(0);
  char const* Beam                       = "";
  PyObject*   List_Lattice_Reference     = PyList_New(0);
  double      Eta                        = 0;

  TVector2D Sigma(0, 0);
  TVector2D SigmaP(0, 0);
  TVector3D Lattice_Reference(0, 0, 0);


  // Input variables and parsing
  static const char *kwlist[] = {"sigma",
                                 "sigmap",
                                 "beam",
                                 "lattice_reference",
                                 "eta",
                                 NULL};

  if (!PyArg_ParseTupleAndKeywords(args, keywds, "OO|sOd",
                                   const_cast<char **>(kwlist),
                                   &List_Sigma,
                                   &List_SigmaP,
                                   &Beam,
                                   &List_Lattice_Reference,
                                   &Eta)) {
    return NULL;
  }


  // Check for Lattice reference in the input
  if (PyList_Size(List_Sigma) != 0) {
    try {
      Sigma = OSCARSPY::ListAsTVector2D(List_Sigma);
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'sigma'");
      return NULL;
    }
  }

  // Check for Lattice reference in the input
  if (PyList_Size(List_SigmaP) != 0) {
    try {
      SigmaP = OSCARSPY::ListAsTVector2D(List_SigmaP);
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'sigmap'");
      return NULL;
    }
  }

  // Check for Lattice reference in the input
  if (PyList_Size(List_Lattice_Reference) != 0) {
    try {
      Lattice_Reference = OSCARSPY::ListAsTVector3D(List_Lattice_Reference);
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'lattice_reference'");
      return NULL;
    }
  }

  //self->obj->SetBeamSize(Beam, Sigma, SigmaP, Lattice_Reference, Eta);

  // Must return python object None in a special way
  Py_INCREF(Py_None);
  return Py_None;
}






const char* DOC_OSCARSSR_ClearParticleBeams = R"docstring(
clear_particle_beams()

Remove all of the existing particle beams

Parameters
----------
None

Returns
-------
None
)docstring";
static PyObject* OSCARSSR_ClearParticleBeams (OSCARSSRObject* self)
{
  // Clear the contents of the particle beam container in OSCARSSR

  self->obj->ClearParticleBeams();

  // Must return python object None in a special way
  Py_INCREF(Py_None);
  return Py_None;
}




const char* DOC_OSCARSSR_PrintParticleBeams = R"docstring(
print_particle_beams()

Print information for all internal particle beams

Parameters
----------
None

Returns
-------
None
)docstring";
static PyObject* OSCARSSR_PrintParticleBeams (OSCARSSRObject* self)
{
  // Print all particle beams stored in OSCARSSR

  // Out string stream for printing beam information
  std::ostringstream ostream;
  ostream << self->obj->GetParticleBeamContainer() << std::endl;

  OSCARSPY::PyPrint_stdout(ostream.str());

  // Must return python object None in a special way
  Py_INCREF(Py_None);
  return Py_None;
}







const char* DOC_OSCARSSR_SetTwissParameters = R"docstring(
set_twiss_parameters(, [beam, beta, alpha, gamma, lattice_reference])

Set the twiss parameters for a given beam (if a beam name is given) or for all beams (if no name is given).  You should only set 2 of the three [beta, alpha, gamma] as they are not independent.

Parameters
----------
beam : str
    Name of the beam for which to set the twiss parameters

beta : list[2]
    [horizontal, vertical]

alpha : list[2]
    [horizontal, vertical]

gamma : list[2]
    [horizontal, vertical]

lattice_reference : list[3]
    Location of the values of the twiss parameters.  These must be places on the beam axis (though it is not checked).

Returns
-------
None
)docstring";
static PyObject* OSCARSSR_SetTwissParameters (OSCARSSRObject* self, PyObject* args, PyObject* keywds)
{
  // Set the twiss parameters for a given beam or for all beams

  // Lists and variables some with initial values
  char const* Beam                       = "";
  PyObject*   List_Beta                  = PyList_New(0);
  PyObject*   List_Alpha                 = PyList_New(0);
  PyObject*   List_Gamma                 = PyList_New(0);
  PyObject*   List_Lattice_Reference     = PyList_New(0);

  TVector2D Beta(0, 0);
  TVector2D Alpha(0, 0);
  TVector2D Gamma(0, 0);
  TVector3D Lattice_Reference(0, 0, 0);


  // Input variables and parsing
  static const char *kwlist[] = {"beam",
                                 "beta",
                                 "alpha",
                                 "gamma",
                                 "lattice_reference",
                                 NULL};

  if (!PyArg_ParseTupleAndKeywords(args, keywds, "|sOOOO",
                                   const_cast<char **>(kwlist),
                                   &Beam,
                                   &List_Beta,
                                   &List_Alpha,
                                   &List_Gamma,
                                   &List_Lattice_Reference
                                  )) {
    return NULL;
  }


  // Check for beta, alpha, gammain the input
  int HasBAG = 0x0;
  if (PyList_Size(List_Beta) != 0) {
    try {
      Beta = OSCARSPY::ListAsTVector2D(List_Beta);
      HasBAG |= 0x4;
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'beta'");
      return NULL;
    }
  }
  if (PyList_Size(List_Alpha) != 0) {
    try {
      Alpha = OSCARSPY::ListAsTVector2D(List_Alpha);
      HasBAG |= 0x2;
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'alpha'");
      return NULL;
    }
  }
  if (PyList_Size(List_Gamma) != 0) {
    try {
      Gamma = OSCARSPY::ListAsTVector2D(List_Gamma);
      HasBAG |= 0x1;
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'gamma'");
      return NULL;
    }
  }

  // Check for Lattice reference in the input
  bool HasReferencePoint = false;
  if (PyList_Size(List_Lattice_Reference) != 0) {
    try {
      Lattice_Reference = OSCARSPY::ListAsTVector3D(List_Lattice_Reference);
      HasReferencePoint = true;
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'lattice_reference'");
      return NULL;
    }
  }

  // Set correct twiss parameters
  switch (HasBAG) {
    case 0x7:
      self->obj->SetTwissParameters(Beam, Beta, Alpha, Gamma, Lattice_Reference, HasReferencePoint);
      break;
    case 0x6:
      self->obj->GetParticleBeam(Beam).SetTwissBetaAlpha(Beta, Alpha, Lattice_Reference, HasReferencePoint);
      break;
    case 0x5:
      self->obj->GetParticleBeam(Beam).SetTwissBetaGamma(Beta, Gamma, Lattice_Reference, HasReferencePoint);
      break;
    case 0x3:
      self->obj->GetParticleBeam(Beam).SetTwissAlphaGamma(Alpha, Gamma, Lattice_Reference, HasReferencePoint);
      break;
    case 0x4:
      self->obj->GetParticleBeam(Beam).SetTwissBetaAlpha(Beta, TVector2D(0, 0), Lattice_Reference, HasReferencePoint);
      break;
    default:
      break;
  }

  // Must return python object None in a special way
  Py_INCREF(Py_None);
  return Py_None;
}






















const char* DOC_OSCARSSR_SetNewParticle = R"docstring(
set_new_particle([, beam, particle])

If no arguments are given sets the current internal particle to a random new particle.  The randomization is based on the weights given for each beam.  This also sets the initial conditions for the particle used in trajectory calculations based on the beam parameters within the randomly sepected beam.  You can specify which beam you want a random particle from using the *beam* parameter.  The *particle* parameter can be 'ideal' if you want the ideal initial conditions for a particle without randomization.

Parameters
----------

beam : str
    The name of the beam from which to get a particle from

particle : str
    'ideal' or 'random', for random, you may omit this.

Returns
-------
None
)docstring";
static PyObject* OSCARSSR_SetNewParticle (OSCARSSRObject* self, PyObject* args, PyObject* keywds)
{
  // Set a new particle within the OSCARSSR object

  char const* Beam_IN = "";
  char const* Particle_IN = "";

  static const char *kwlist[] = {"beam",
                                 "particle",
                                 NULL};

  if (!PyArg_ParseTupleAndKeywords(args, keywds, "|ss",
                                   const_cast<char **>(kwlist),
                                   &Beam_IN,
                                   &Particle_IN)) {
    return NULL;
  }

  // Check if a beam is at least defined
  if (self->obj->GetNParticleBeams() < 1) {
    PyErr_SetString(PyExc_ValueError, "No particle beam defined");
    return NULL;
  }

  std::string Beam = Beam_IN;
  std::string Particle = Particle_IN;

  std::transform(Beam.begin(), Beam.end(), Beam.begin(), ::tolower);
  std::transform(Particle.begin(), Particle.end(), Particle.begin(), ::tolower);

  if (Particle != "" && !(Particle == "ideal" || Particle == "random")) {
    PyErr_SetString(PyExc_ValueError, "'particle' must be 'ideal' or 'random'");
    return NULL;
  }

  try {
    if ((Beam) == "" && (Particle) == "") {
      self->obj->SetNewParticle();
    } else if (Beam != "" && Particle != "") {
      self->obj->SetNewParticle(Beam, Particle);
    } else if (Beam == "" && Particle != "") {
      self->obj->SetNewParticle(Beam, Particle);
    } else if (Beam != "" && Particle == "") {
      self->obj->SetNewParticle(Beam, Particle);
    }
  } catch (std::out_of_range e) {
    PyErr_SetString(PyExc_ValueError, "'beam' name not found");
    return NULL;
  }

  // Must return python object None in a special way
  Py_INCREF(Py_None);
  return Py_None;
}





const char* DOC_OSCARSSR_GetParticleX0 = R"docstring(
get_particle_x0()

Get the initial position for the current particle in [m]

Returns
-------
x0 : [float, float, float]
    Initial particle position
)docstring";
static PyObject* OSCARSSR_GetParticleX0 (OSCARSSRObject* self)
{
  // Get the particle position at particle t0
  return OSCARSPY::TVector3DAsList( self->obj->GetCurrentParticle().GetX0() );
}




const char* DOC_OSCARSSR_GetParticleBeta0 = R"docstring(
get_particle_beta0()

Get the initial :math:`\vec \beta` for the current particle

Returns
-------
b0 : [float, float, float]
    Initial beta of particle
)docstring";
static PyObject* OSCARSSR_GetParticleBeta0 (OSCARSSRObject* self)
{
  // Get the particle beta at particle t0
  return OSCARSPY::TVector3DAsList( self->obj->GetCurrentParticle().GetB0() );
}




const char* DOC_OSCARSSR_GetParticleE0 = R"docstring(
get_particle_e0()

Get the initial energy for the current particle in [GeV]

Returns
-------
e0 : float
    Initial energy of current particle
)docstring";
static PyObject* OSCARSSR_GetParticleE0 (OSCARSSRObject* self)
{
  // Get the particle beta at particle t0
  return Py_BuildValue("f", (self->obj->GetCurrentParticle().GetE0()));
}






const char* DOC_OSCARSSR_AddDriftVolume_Box = R"docstring(
add_drift_box(width, [, rotations, translation, name])

Add a drift volume box.  The particle will drift with constant velocity inside of this volume and the points inside of the volume will not be recorded in the trajectory.  If a dimension of the box is <= 0 that axis is ignored.

Parameters
----------
width: list
    Width of the drift box in x, y, and z

rotations : list, optional
    3-element list representing rotations around x, y, and z axes: [:math:`\theta_x, \theta_y, \theta_z`]

translation : list, optional
    3-element list representing a translation in space [x, y, z]

name : str
    Name of the drift volume

Returns
-------
None
)docstring";
static PyObject* OSCARSSR_AddDriftVolume_Box (OSCARSSRObject* self, PyObject* args, PyObject* keywds)
{
  // Add a magnetic field from a file.
  // UPDATE: add binary file reading


  // Grab the values
  PyObject*   List_Width       = PyList_New(0);
  PyObject*   List_Rotations   = PyList_New(0);
  PyObject*   List_Translation = PyList_New(0);
  char const* Name             = "";

  TVector3D Width(0, 0, 0);
  TVector3D Rotations(0, 0, 0);
  TVector3D Translation(0, 0, 0);


  // Input variables and parsing
  static const char *kwlist[] = {"width",
                                 "rotations",
                                 "translation",
                                 "name",
                                 NULL};

  if (!PyArg_ParseTupleAndKeywords(args, keywds, "O|OOs",
                                   const_cast<char **>(kwlist),
                                   &List_Width,
                                   &List_Rotations,
                                   &List_Translation,
                                   &Name)) {
    return NULL;
  }

  // Get width
  try {
    Width = OSCARSPY::ListAsTVector3D(List_Width);
  } catch (std::length_error e) {
    PyErr_SetString(PyExc_ValueError, "Incorrect format in 'width'");
    return NULL;
  }

  // Check for Rotations in the input
  if (PyList_Size(List_Rotations) != 0) {
    try {
      Rotations = OSCARSPY::ListAsTVector3D(List_Rotations);
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'rotations'");
      return NULL;
    }
  }

  // Check for Translation in the input
  if (PyList_Size(List_Translation) != 0) {
    try {
      Translation = OSCARSPY::ListAsTVector3D(List_Translation);
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'translation'");
      return NULL;
    }
  }

  // Name check
  if (std::string(Name).size() > 0 && Name[0] == '_') {
    PyErr_SetString(PyExc_ValueError, "'name' cannot begin with '_'.  This is reserved for internal use.  Please pick a different name");
    return NULL;
  }

  // Add the drift box to the OSCARSSR object
  try {
    self->obj->AddDriftVolume( new TDriftBox(Width, Translation, Rotations, Name) );
  } catch (...) {
    PyErr_SetString(PyExc_ValueError, "Could not add drift volume.");
    return NULL;
  }

  // Must return python object None in a special way
  Py_INCREF(Py_None);
  return Py_None;
}






const char* DOC_OSCARSSR_RemoveDriftVolume = R"docstring(
remove_drift(name)

Remove all drift volumes with the given name

Parameters
----------
name : str
    Name of the drift volume

Returns
-------
None
)docstring";
static PyObject* OSCARSSR_RemoveDriftVolume (OSCARSSRObject* self, PyObject* args, PyObject* keywds)
{
  // Add a magnetic field from a file.
  // UPDATE: add binary file reading


  // Grab the values
  char const* Name             = "";

  // Input variables and parsing
  static const char *kwlist[] = {"name",
                                 NULL};

  if (!PyArg_ParseTupleAndKeywords(args, keywds, "s",
                                   const_cast<char **>(kwlist),
                                   &Name)) {
    return NULL;
  }

  self->obj->RemoveDriftVolume(Name);

  // Must return python object None in a special way
  Py_INCREF(Py_None);
  return Py_None;
}






const char* DOC_OSCARSSR_ClearDriftVolumes = R"docstring(
clear_drifts()

Remove all drift volumes

Parameters
----------
None

Returns
-------
None
)docstring";
static PyObject* OSCARSSR_ClearDriftVolumes (OSCARSSRObject* self)
{
  // Print all magnetic stored in OSCARSSR

  self->obj->ClearDriftVolumes();

  // Must return python object None in a special way
  Py_INCREF(Py_None);
  return Py_None;
}









const char* DOC_OSCARSSR_PrintDriftVolumes = R"docstring(
print_drifts()

Print information about all drift volumes to the standard out

Parameters
----------
None

Returns
-------
None
)docstring";
static PyObject* OSCARSSR_PrintDriftVolumes (OSCARSSRObject* self)
{
  // Print all magnetic stored in OSCARSSR

  // Out string stream for printing beam information
  std::ostringstream ostream;
  ostream << "*Drift Volumes*\n";
  ostream << self->obj->GetDriftVolumeContainer() << std::endl;

  OSCARSPY::PyPrint_stdout(ostream.str());

  // Must return python object None in a special way
  Py_INCREF(Py_None);
  return Py_None;
}









const char* DOC_OSCARSSR_CorrectTrajectory= R"docstring(
correct_trajectory()

Correct the trajectory to desired position and direction and specified point based on corrector kicks at specified locations.

This function will add magnetic fields to the sr object with names that begin with an underscore.

Parameters
----------
None

Returns
-------
trajectory : list
    A list of points of the form [[[x, y, z], [Beta_x, Beta_y, Beta_z]], ...]
)docstring";
static PyObject* OSCARSSR_CorrectTrajectory (OSCARSSRObject* self)
{
  // Get the CTStop variable from OSCARSSR

  try {
    self->obj->CorrectTrajectory();
  } catch (...) {
    // UPDATE: Do I need t catch something?
    throw;
  }
  // Return the trajectory
  return OSCARSSR_GetTrajectory(self);
}





const char* DOC_OSCARSSR_CalculateTrajectory = R"docstring(
calculate_trajectory()

Calculates the trajectory for the current internal particle.  This calculates the trajectory in 3D from the time set by oscars.sr.set_ctstart() to the time set by oscars.sr.set_ctstop() beginning at the *t0* given by the particle beam from which this particle comes from.  It first does a forward propogation to the stop time, then a backward propogation to the start time.

It is not necessary to call this method before other calculations such as spectrum, power density, or flux calculation methods.

If you have a current particle loaded using *SetNewParticle* this method will calculate the trajectory for that particle.  If no particle is defined one will be randomly selected based on the beam weights and beam parameters.

Parameters
----------
None

Returns
-------
trajectory : list
    A list of points of the form [[[x, y, z], [Beta_x, Beta_y, Beta_z]], ...]
)docstring";
static PyObject* OSCARSSR_CalculateTrajectory (OSCARSSRObject* self)
{
  // Get the CTStop variable from OSCARSSR

  try {
    self->obj->CalculateTrajectory();
  } catch (std::length_error e) {
    PyErr_SetString(PyExc_ValueError, e.what());
    return NULL;
  } catch (std::out_of_range e) {
    PyErr_SetString(PyExc_ValueError, e.what());
    return NULL;
  }

  // Return the trajectory
  return OSCARSSR_GetTrajectory(self);
}





const char* DOC_OSCARSSR_GetTrajectory = R"docstring(
get_trajectory()

Get the current trajectory.  If the trajectory has not been calculated this will return an empty list.  The format of the returned list consists of a list of lists giving you the position and beta (v/c) of the particle at each position.  For a trajectory returnned is of the form: [[[x, y, z], [Beta_x, Beta_y, Beta_z]], ...]

Parameters
----------
None

Returns
-------
trajectory : list
    A list of points of the form [[[x, y, z], [Beta_x, Beta_y, Beta_z]], ...]
)docstring";
static PyObject* OSCARSSR_GetTrajectory (OSCARSSRObject* self)
{
  // Get the Trajectory as 2 3D lists [[x, y, z], [BetaX, BetaY, BetaZ]]

  // Create a python list
  PyObject *PList = PyList_New(0);

  // Grab trajectory
  TParticleTrajectoryPoints const& T = self->obj->GetTrajectory();

  // Number of points in trajectory calculation
  size_t NTPoints = T.GetNPoints();

  PyObject* Value;

  // Loop over all points in trajectory
  for (size_t iT = 0; iT != NTPoints; ++iT) {
    // Create a python list for X and Beta
    PyObject *PList2 = PyList_New(0);

    // Add position and Beta to list
    Value = Py_BuildValue("f", T.GetT(iT));
    PyList_Append(PList2, Value);
    Py_DECREF(Value);

    Value = OSCARSPY::TVector3DAsList(T.GetX(iT));
    PyList_Append(PList2, Value);
    Py_DECREF(Value);

    Value = OSCARSPY::TVector3DAsList(T.GetB(iT));
    PyList_Append(PList2, Value);
    Py_DECREF(Value);

    Value = OSCARSPY::TVector3DAsList(T.GetA(iT));
    PyList_Append(PList2, Value);
    Py_DECREF(Value);

    PyList_Append(PList, PList2);
    Py_DECREF(PList2);
  }

  // Return the python list
  return PList;
}

































const char* DOC_OSCARSSR_CalculateSpectrum = R"docstring(
calculate_spectrum(obs [, npoints, energy_range_eV, energy_points_eV, points_eV, polarization, angle, horizontal_direction=[1, 0, 0], propogation_direction=[0, 0, 1], precision, max_level, nparticles, nthreads, gpu, ngpu, quantity, ofile, bofile])

Calculate the spectrum given a point in space, the range in energy, and the number of points.  The calculation uses the current particle and its initial conditions.  If the trajectory has not been calculated it is calculated first.  The units of this calculation are [:math:`photons / mm^2 / 0.1% bw / s`]

Previously to calling this function you must define a particle beam and define start and stop times at very minimum.

You **must** provide either (*npoints* and *energy_range_eV*) or *points_eV*.

Parameters
----------
obs : list
    Point [x, y, z] where you wish to calculate the spectrum

npoints : int
    Number of points to calculate in the given energy range

energy_range_eV : list
    energy range [min, max] in eV as a list of length 2

points_eV : list
    A list of points to calculate the flux at ie [12.3, 45.6, 78.9, 123.4]

polarization : str
    Which polarization mode to calculate.  Can be 'all', 'linear-horizontal', 'linear-vertical', 'circular-left', 'circular-right', or 'linear' (if linear you must specify the angle parameter).  The short versions are 'lh', 'lv', 'cl', 'cr'.  Theses are not case-sensative.   There is no need to specify 'linear' if you give the angle parameter.

horizontal_direction : list
    The direction you consider to be horizontal.  Should be perpendicular to the photon beam propogation direction.  Default is [1, 0, 0]

propogation_direction : list
    Propogation direction of photon beam. Default is [0, 0, 1]

precision : float
    Calculation precision parameter (typically 0.01 which is 1%)

max_level: int
    Maximum "level" to use for trajectory in the calculation.  Level N corresponds to a total of 2**(N+2) trajectory points.  You cannot go beyond the internal maximum.  You are not guaranteed precision parameter is met if this is used.

max_level_extended: int
    Maximum "level" to use for trajectory in the calculation.  If set to higher than max_level the computation will proceed beyond max_level without creating trajectory arrays in memory (but it will be slower)

angle : float
    The 'angle' is that from the horizontal_direction for the polarization direction you are interested in

nparticles : int
    The number of particles you wish to run for a multi-particle simulation

nthreads : int
    Number of threads to use

gpu : int
    Use the gpu or not (0 or 1).  If 1 will attempt to use ALL gpus available.  This is overridden if you use the input 'ngpu'

ngpu : int or list
    If ngpu is an int, use that number of gpus (if available).
    If ngpu is a list, the list should be a list of gpus you wish to use

quantity: str
    Quantity to return.
    Available are:
        'flux'        (default)
        'precision' - Estimated precision for each point
        'level'     - Trajectory level reached (npoints = 2**(n+1) - 1), if return is -1 the requested precision was not reached

ofile : str
    Output file name

bofile : str
    Binary output file name

Returns
-------
spectrum : list
    A list of 2D lists, each of which is a pair representing the energy [eV] and flux [:math:`photons / mm^2 / 0.1% bw / s`] at that energy.  eg [[energy_0, flux_0], [energy_1, flux_1], ...]

Examples
--------
Calculate the spectrum at a point 30 [m] downstream (assuming beam direction is in the Z direction) in the energy range 100 to 1000 [eV] with 900 points.

    >>> osr.calculate_spectrum(obs=[0, 0, 30], energy_range_eV=[100, 1000], npoints=900)
)docstring";
static PyObject* OSCARSSR_CalculateSpectrum (OSCARSSRObject* self, PyObject* args, PyObject* keywds)
{
  // Calculate the spectrum given an observation point, and energy range


  PyObject*   List_Obs                  = 0x0;
  int         NPoints                   = 0;
  PyObject*   List_EnergyRange_eV       = 0x0;
  PyObject*   List_Points_eV            = 0x0;
  char const* PolarizationIn            = "";
  double      Angle                     = 112188979912321;
  PyObject*   List_HorizontalDirection  = 0x0;
  PyObject*   List_PropogationDirection = 0x0;
  double      Precision                 = 0.01;
  int         MaxLevel                  = -2;
  int         MaxLevelExtended          = 0;
  int         NParticles                = 0;
  int         NThreads                  = 0;
  int         GPU                       = -1;
  PyObject*   NGPU = 0x0;
  char const* ReturnQuantityChars       = "flux";
  const char* OutFileNameText           = "";
  const char* OutFileNameBinary         = "";

  // Input variable list
  static const char *kwlist[] = {"obs",
                                 "npoints",
                                 "energy_range_eV",
                                 "energy_points_eV",
                                 "points_eV", // UPDATE: REMOVE in v2
                                 "polarization",
                                 "angle",
                                 "horizontal_direction",
                                 "propogation_direction",
                                 "precision",
                                 "max_level",
                                 "max_level_extended",
                                 "nparticles",
                                 "nthreads",
                                 "gpu",
                                 "ngpu",
                                 "quantity",
                                 "ofile",
                                 "bofile",
                                 NULL};

  // Parse inputs
  if (!PyArg_ParseTupleAndKeywords(args, keywds, "O|iOOOsdOOdiiiiiOsss",
                                   const_cast<char **>(kwlist),
                                   &List_Obs,
                                   &NPoints,
                                   &List_EnergyRange_eV,
                                   &List_Points_eV,
                                   &List_Points_eV,
                                   &PolarizationIn,
                                   &Angle,
                                   &List_HorizontalDirection,
                                   &List_PropogationDirection,
                                   &Precision,
                                   &MaxLevel,
                                   &MaxLevelExtended,
                                   &NParticles,
                                   &NThreads,
                                   &GPU,
                                   &NGPU,
                                   &ReturnQuantityChars,
                                   &OutFileNameText,
                                   &OutFileNameBinary)) {
    return NULL;
  }

  // Check if a beam is at least defined
  //if (self->obj->GetNParticleBeams() < 1) {
  //  PyErr_SetString(PyExc_ValueError, "No particle beam defined");
  //  return NULL;
  //}


  // Check number of particles
  if (NParticles < 0) {
    PyErr_SetString(PyExc_ValueError, "'nparticles' must be >= 1 (sort of)");
    return NULL;
  }


  // Add all values to a vector
  std::vector<double> VPoints_eV;
  if (List_Points_eV != 0x0) {
    for (int i = 0; i < PyList_Size(List_Points_eV); ++i) {
      VPoints_eV.push_back(PyFloat_AsDouble(PyList_GetItem(List_Points_eV, i)));
    }
  }

  double EStart = 0;
  double EStop = 0;

  if (List_EnergyRange_eV != 0x0) {
    if (PyList_Size(List_EnergyRange_eV) == 2) {
      EStart = PyFloat_AsDouble(PyList_GetItem(List_EnergyRange_eV, 0));
      EStop  = PyFloat_AsDouble(PyList_GetItem(List_EnergyRange_eV, 1));
    } else {
      PyErr_SetString(PyExc_ValueError, "'energy_range_eV' must be a list of length 2");
      return NULL;
    }
  }




  // Observation point
  TVector3D Obs(0, 0, 0);
  try {
    Obs = OSCARSPY::ListAsTVector3D(List_Obs);
  } catch (std::length_error e) {
    PyErr_SetString(PyExc_ValueError, "Incorrect format in 'obs'");
    return NULL;
  }


  // Check for HorizontalDirection in the input
  TVector3D HorizontalDirection(1, 0, 0);
  if (List_HorizontalDirection != 0x0 && PyList_Size(List_HorizontalDirection) != 0) {
    try {
      HorizontalDirection = OSCARSPY::ListAsTVector3D(List_HorizontalDirection);
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'translation'");
      return NULL;
    }
  }


  // Check for PropogationDirection in the input
  TVector3D PropogationDirection(0, 0, 1);
  if (List_PropogationDirection != 0x0 && PyList_Size(List_PropogationDirection) != 0) {
    try {
      PropogationDirection = OSCARSPY::ListAsTVector3D(List_PropogationDirection);
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'translation'");
      return NULL;
    }
  }



  // Check NThreads parameter
  if (NThreads < 0) {
    PyErr_SetString(PyExc_ValueError, "'nthreads' must be > 0");
    return NULL;
  }


  // Check GPU parameter
  if (GPU != 0 && GPU != 1 && GPU != -1) {
    PyErr_SetString(PyExc_ValueError, "'gpu' must be 0 or 1");
    return NULL;
  }


  // Check you are not trying to use threads and GPU
  if (NThreads > 0 && GPU == 1) {
    PyErr_SetString(PyExc_ValueError, "gpu is 1 and nthreads > 0.  Both are not currently allowed.");
    return NULL;
  }

  // Container for spectrum
  TSpectrumContainer SpectrumContainer;

  if (VPoints_eV.size() == 0) {
    // Check NPoints parameter and set minimum if still zero
    if (NPoints < 1) {
      NPoints = fabs(EStop - EStart) + 1 > 100 ? abs((int) (EStop - EStart) + 1) : 100;
    }
    SpectrumContainer.Init(NPoints, EStart, EStop);
  } else {
    SpectrumContainer.Init(VPoints_eV);
  }

  // Check ngpu input
  int NumberOfGPUs = -1;
  std::vector<int> GPUVector;
  if (NGPU != 0x0) {
    if (PyLong_Check(NGPU)) {
      NumberOfGPUs = (int) PyLong_AsLong(NGPU);
    } else if (PyList_Check(NGPU)) {
      OSCARSPY::ListToVectorInt(NGPU, GPUVector);
    }
  }


  int ReturnQuantity = 0;
  std::string ReturnQuantityStr = ReturnQuantityChars;
  std::transform(ReturnQuantityStr.begin(), ReturnQuantityStr.end(), ReturnQuantityStr.begin(), ::toupper);
  if (ReturnQuantityStr == "FLUX") {
    ReturnQuantity = 0;
  } else if (ReturnQuantityStr == "PRECISION") {
    ReturnQuantity = 1;
  } else if (ReturnQuantityStr == "LEVEL") {
    ReturnQuantity = 2;
  } else {
    PyErr_SetString(PyExc_ValueError, "'quantity' must be: 'flux', 'precision', 'level', or blank");
    return NULL;
  }

  // If all is specified make sure no angle is specified
  if (std::string(PolarizationIn) == "all" && Angle != 112188979912321) {
    PyErr_SetString(PyExc_ValueError, "cannot specify 'all' and 'angle'.  This is to save you");
    return NULL;
  }

  // If Angle is specified and polarization is blank set it to linear
  std::string Polarization = "all";
  if (strlen(PolarizationIn) != 0) {
    Polarization = PolarizationIn;
  } else if (Angle != 112188979912321) {
    Polarization = "linear";
  }

  // Actually calculate the spectrum
  try {
    self->obj->CalculateSpectrum(Obs,
                                 SpectrumContainer,
                                 Polarization,
                                 Angle,
                                 HorizontalDirection,
                                 PropogationDirection,
                                 NParticles,
                                 NThreads,
                                 GPU,
                                 NumberOfGPUs,
                                 GPUVector,
                                 Precision,
                                 MaxLevel,
                                 MaxLevelExtended,
                                 ReturnQuantity);

  } catch (std::length_error e) {
    PyErr_SetString(PyExc_ValueError, e.what());
    return NULL;
  } catch (std::out_of_range e) {
    PyErr_SetString(PyExc_ValueError, e.what());
    return NULL;
  } catch (std::invalid_argument e) {
    PyErr_SetString(PyExc_ValueError, e.what());
    return NULL;
  }

  if (!SpectrumContainer.AllConverged()) {
    OSCARSPY::PyPrint_stderr("Not all points converged to desired precision.  Can try increasing 'max_level_extended'\n");
  }


  if (std::string(OutFileNameText) != "") {
    SpectrumContainer.WriteToFileText(OutFileNameText);
  }

  if (std::string(OutFileNameBinary) != "") {
    SpectrumContainer.WriteToFileBinary(OutFileNameBinary);
  }

  // Return the spectrum
  return OSCARSPY::GetSpectrumAsList(SpectrumContainer);
}
















const char* DOC_OSCARSSR_CalculateTotalPower = R"docstring(
calculate_total_power([, precision, max_level, max_level_extended, quantity])

Calculate the total radiated power based on the current particle and that particle's beam current.

See the :doc:`MathematicalNotes` section for the expression used in this calculation.

precision : float
    Calculation precision parameter (typically 0.01 which is 1%)

max_level: int
    Maximum "level" to use for trajectory in the calculation.  Level N corresponds to a total of 2**(N+2) trajectory points.  You cannot go beyond the internal maximum.  You are not guaranteed precision parameter is met if this is used.


max_level_extended: int
    Maximum "level" to use for trajectory in the calculation.  If set to higher than max_level the computation will proceed beyond max_level without creating trajectory arrays in memory (but it will be slower)

quantity: str
    Quantity to return.
    Available are:
        'power'     - (default)
        'precision' - Estimated precision
        'level'     - Trajectory level reached (npoints = 2**(n+1) - 1), if return is -1 the requested precision was not reached

Returns
-------
power : float
    Total power in [W]
)docstring";
static PyObject* OSCARSSR_CalculateTotalPower (OSCARSSRObject* self, PyObject* args, PyObject *keywds)
{
  // Calculate the total power radiated by the current particle

  double      Precision = 0.01;
  int         MaxLevel = -2;
  int         MaxLevelExtended = 0;
  char const* ReturnQuantityChars = "power";


  static const char *kwlist[] = {"precision",
                                 "max_level",
                                 "max_level_extended",
                                 "quantity",
                                 NULL};

  if (!PyArg_ParseTupleAndKeywords(args, keywds, "|diis",
                                   const_cast<char **>(kwlist),
                                   &Precision,
                                   &MaxLevel,
                                   &MaxLevelExtended,
                                   &ReturnQuantityChars)) {
    return NULL;
  }


  int ReturnQuantity = 0;
  std::string ReturnQuantityStr = ReturnQuantityChars;
  std::transform(ReturnQuantityStr.begin(), ReturnQuantityStr.end(), ReturnQuantityStr.begin(), ::toupper);
  if (ReturnQuantityStr == "POWER" || ReturnQuantityStr == "POWER") {
    ReturnQuantity = 0;
  } else if (ReturnQuantityStr == "PRECISION") {
    ReturnQuantity = 1;
  } else if (ReturnQuantityStr == "LEVEL") {
    ReturnQuantity = 2;
  } else {
    PyErr_SetString(PyExc_ValueError, "'quantity' must be: 'power', 'precision', 'level', or blank");
    return NULL;
  }



  double Power = 0;

  // Check if a beam is at least defined
  if (self->obj->GetNParticleBeams() < 1) {
    PyErr_SetString(PyExc_ValueError, "No particle beam defined");
    return NULL;
  }

  // Return the total power
  // UPDATE: This does not fail when no beam defined
  try {
    Power = self->obj->CalculateTotalPower(Precision, MaxLevel, MaxLevelExtended, ReturnQuantity);
  } catch (std::length_error e) {
    PyErr_SetString(PyExc_ValueError, e.what());
    return NULL;
  } catch (std::out_of_range e) {
    PyErr_SetString(PyExc_ValueError, e.what());
    return NULL;
  }

  return Py_BuildValue("f", Power);
}







const char* DOC_OSCARSSR_CalculatePowerDensity = R"docstring(
calculate_power_density(points [, normal, rotations, translation, nparticles, gpu, nthreads, precision, max_level, max_level_extended, quantity, ofile])

Calculate the power density for each point in the list *points*.

See the :doc:`MathematicalNotes` section for the expression used in this calculation.

Parameters
----------

points : list
    A list of points, each point containing a position in 3D (as a list) and a normal vector at that position (also as a 3D list): [[[x, y, z], [nx. ny. nz]], [...], ...]

normal : int
    -1 if you wish to reverse the normal vector, 0 if you wish to ignore the +/- direction in computations, 1 if you with to use the direction of the normal vector as given. 

rotations : list, optional
    3-element list representing rotations around x, y, and z axes: [:math:`\theta_x, \theta_y, \theta_z`]

translation : list, optional
    3-element list representing a translation in space [x, y, z]

nparticles : int
    Number of particles to use for multi-particle calculations

gpu : int
    Use the gpu or not (0 or 1)

nthreads : int
    Number of threads to use

precision : float
    Calculation precision parameter (typically 0.01 which is 1%)

max_level: int
    Maximum "level" to use for trajectory in the calculation.  Level N corresponds to a total of 2**(N+2) trajectory points.  You cannot go beyond the internal maximum.  You are not guaranteed precision parameter is met if this is used.


max_level_extended: int
    Maximum "level" to use for trajectory in the calculation.  If set to higher than max_level the computation will proceed beyond max_level without creating trajectory arrays in memory (but it will be slower)

quantity: str
    Quantity to return.
    Available are:
        'power density' (default)
        'precision' - Estimated precision for each point
        'level'     - Trajectory level reached (npoints = 2**(n+1) - 1), if return is -1 the requested precision was not reached

ofile : str
    Output file name

Returns
-------
power_density : list
    A list, each element of which is a pair representing the position (2D relative (default) or 3D absolute) and power density [:math:`W / mm^2`] at that position.  eg [[[x1_0, x2_0, x3_0], pd_0], [[x1_1, x2_1, x3_1], pd_1]],  ...].  The position is always given as a list of length 3.  For the default (dim=2) the third element is always zero.
)docstring";
static PyObject* OSCARSSR_CalculatePowerDensity (OSCARSSRObject* self, PyObject* args, PyObject *keywds)
{
  // Calculate the power density given a list of points

  PyObject*   List_Translation = PyList_New(0);
  PyObject*   List_Rotations   = PyList_New(0);
  PyObject*   List_Points      = PyList_New(0);
  int         NormalDirection = 0;
  int const   Dim = 3;
  int         NParticles = 0;
  int         GPU = -1;
  PyObject*   NGPU = 0x0;
  int         NThreads = 0;
  double      Precision = 0.01;
  int         MaxLevel = -2;
  int         MaxLevelExtended = 0;
  char const* ReturnQuantityChars = "power density";
  char const* OutFileNameText = "";
  char const* OutFileNameBinary = "";


  static const char *kwlist[] = {"points",
                                 "normal",
                                 "rotations",
                                 "translation",
                                 "nparticles",
                                 "gpu",
                                 "ngpu",
                                 "nthreads",
                                 "precision",
                                 "max_level",
                                 "max_level_extended",
                                 "quantity",
                                 "ofile",
                                 "bofile",
                                 NULL};

  if (!PyArg_ParseTupleAndKeywords(args, keywds, "O|iOOiiOidiisss",
                                   const_cast<char **>(kwlist),
                                   &List_Points,
                                   &NormalDirection,
                                   &List_Rotations,
                                   &List_Translation,
                                   &NParticles,
                                   &GPU,
                                   &NGPU,
                                   &NThreads,
                                   &Precision,
                                   &MaxLevel,
                                   &MaxLevelExtended,
                                   &ReturnQuantityChars,
                                   &OutFileNameText,
                                   &OutFileNameBinary)) {
    return NULL;
  }

  // Check if a beam is at least defined
  if (self->obj->GetNParticleBeams() < 1) {
    PyErr_SetString(PyExc_ValueError, "No particle beam defined");
    return NULL;
  }


  // Check requested dimension
  if (Dim != 2 && Dim != 3) {
    PyErr_SetString(PyExc_ValueError, "'dim' must be 2 or 3");
    return NULL;
  }


  // Vectors for rotations and translations.  Default to 0
  TVector3D Rotations(0, 0, 0);
  TVector3D Translation(0, 0, 0);


  // Check for Rotations in the input
  if (PyList_Size(List_Rotations) != 0) {
    try {
      Rotations = OSCARSPY::ListAsTVector3D(List_Rotations);
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'rotations'");
      return NULL;
    }
  }


  // Check for Translation in the input
  if (PyList_Size(List_Translation) != 0) {
    try {
      Translation = OSCARSPY::ListAsTVector3D(List_Translation);
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'translation'");
      return NULL;
    }
  }

  // Look for arbitrary shape 3D points
  TSurfacePoints_3D Surface;
  for (int i = 0; i < PyList_Size(List_Points); ++i) {
    PyObject* LXN = PyList_GetItem(List_Points, i);
    TVector3D X;
    TVector3D N;
    if (PyList_Size(LXN) == 2) {

      try {
        X = OSCARSPY::ListAsTVector3D(PyList_GetItem(LXN, 0));
        N = OSCARSPY::ListAsTVector3D(PyList_GetItem(LXN, 1));
      } catch (std::length_error e) {
        PyErr_SetString(PyExc_ValueError, "Incorrect format in 'points': Point or Normal does not have 3 elements");
        return NULL;
      }

      // Rotate point and normal
      X.RotateSelfXYZ(Rotations);
      N.RotateSelfXYZ(Rotations);

      // Invert the normal?
      if (NormalDirection == -1) {
        N *= -1;
      }

      // Translate point, normal does not get translated
      X += Translation;

      Surface.AddPoint(X, N);
    } else {
      // input format error
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'points'");
      return NULL;
    }
  }

  // Check number of particles
  if (NParticles < 0) {
    PyErr_SetString(PyExc_ValueError, "'nparticles' must be >= 1 (sort of)");
    return NULL;
  }


  // Check GPU parameter
  if (GPU != 0 && GPU != 1 && GPU != -1) {
    PyErr_SetString(PyExc_ValueError, "'gpu' must be 0 or 1");
    return NULL;
  }

  // Check NThreads parameter
  if (NThreads < 0) {
    PyErr_SetString(PyExc_ValueError, "'nthreads' must be > 0");
    return NULL;
  }

  // Check you are not trying to use threads and GPU
  if (NThreads > 0 && GPU == 1) {
    PyErr_SetString(PyExc_ValueError, "gpu is 1 and nthreads > 0.  Both are not currently allowed.");
    return NULL;
  }

  // Check ngpu input
  int NumberOfGPUs = -1;
  std::vector<int> GPUVector;
  if (NGPU != 0x0) {
    if (PyLong_Check(NGPU)) {
      NumberOfGPUs = (int) PyLong_AsLong(NGPU);
    } else if (PyList_Check(NGPU)) {
      OSCARSPY::ListToVectorInt(NGPU, GPUVector);
    }
  }


  int ReturnQuantity = 0;
  std::string ReturnQuantityStr = ReturnQuantityChars;
  std::transform(ReturnQuantityStr.begin(), ReturnQuantityStr.end(), ReturnQuantityStr.begin(), ::toupper);
  if (ReturnQuantityStr == "POWER DENSITY" || ReturnQuantityStr == "POWERDENSITY") {
    ReturnQuantity = 0;
  } else if (ReturnQuantityStr == "PRECISION") {
    ReturnQuantity = 1;
  } else if (ReturnQuantityStr == "LEVEL") {
    ReturnQuantity = 2;
  } else {
    PyErr_SetString(PyExc_ValueError, "'quantity' must be: 'power density', 'precision', 'level', or blank");
    return NULL;
  }

  // Container for Point plus scalar
  T3DScalarContainer PowerDensityContainer;


  // Actually calculate the spectrum
  bool const Directional = NormalDirection == 0 ? false : true;

  try {
    self->obj->CalculatePowerDensity(Surface,
                                     PowerDensityContainer,
                                     Dim,
                                     Directional,
                                     Precision,
                                     MaxLevel,
                                     MaxLevelExtended,
                                     NParticles,
                                     NThreads,
                                     GPU,
                                     NumberOfGPUs,
                                     GPUVector,
                                     ReturnQuantity);
  } catch (std::length_error e) {
    PyErr_SetString(PyExc_ValueError, e.what());
    return NULL;
  } catch (std::out_of_range e) {
    PyErr_SetString(PyExc_ValueError, e.what());
    return NULL;
  } catch (std::invalid_argument e) {
    PyErr_SetString(PyExc_ValueError, e.what());
    return NULL;
  }

  // If not converged print warning
  if (!PowerDensityContainer.AllConverged()) {
    OSCARSPY::PyPrint_stderr("Not all points converged to desired precision.  Can try increasing 'max_level_extended'\n");
  }

  // Write the output file if requested
  // Text output
  if (std::string(OutFileNameText) != "") {
    PowerDensityContainer.WriteToFileText(OutFileNameText, Dim);
  }

  // Binary output
  if (std::string(OutFileNameBinary) != "") {
    PowerDensityContainer.WriteToFileBinary(OutFileNameBinary, Dim);
  }


  // Build the output list of: [[[x, y, z], PowerDensity], [...]]
  // Create a python list
  PyObject *PList = PyList_New(0);

  // UPDATE: URGENT: PD output ofile, bofile
  size_t const NPoints = PowerDensityContainer.GetNPoints();

  PyObject* Value;

  for (size_t i = 0; i != NPoints; ++i) {
    T3DScalar P = PowerDensityContainer.GetPoint(i);

    // Inner list for each point
    PyObject *PList2 = PyList_New(0);


    // Add position and value to list
    Value = OSCARSPY::TVector3DAsList(P.GetX());
    PyList_Append(PList2, Value);
    Py_DECREF(Value);

    Value = Py_BuildValue("f", P.GetV());
    PyList_Append(PList2, Value);
    Py_DECREF(Value);

    PyList_Append(PList, PList2);
    Py_DECREF(PList2);

  }

  return PList;
}








const char* DOC_OSCARSSR_CalculatePowerDensityRectangle = R"docstring(
calculate_power_density_rectangle(npoints [, plane, width, x0x1x2, rotations, translation, ofile, bofile, normal, nparticles, gpu, ngpu, nthreads, precision, max_level, max_level_extended, dim, quantity])

Calculate the power density in a rectangle either defined by three points, or by defining the plane the rectangle is in and the width, and then rotating and translating it to where it needs be.  The simplest is outlined in the first example below.  By default (dim=2) this returns a list whose position coordinates are in the local coordinate space x1 and x2 (*ie* they do not include the rotations and translation).  if dim=3 the coordinates in the return list are in absolute 3D space.

See the :doc:`MathematicalNotes` section for the expression used in this calculation.

You **must** specify either both (*plane* and *width*) or *x0x1x2*

Parameters
----------
npoints: int
    number of in each dimension for surface

plane : str
    The plane to start in (XY, XZ, YZ, YX, ZX, ZY).  The normal to the surface is defined using the right handed cross product (ie the last three have opposite normal vectors from the first three)

width : list
    Width of rectangle in X1 and X2: [w1, w2]

x0x1x2 : list
    List of three points [[x0, y0, z0], [x1, y1, z1], [x2, y2, z2]] defining a parallelogram (vectors 0->1, and 0->2)

rotations : list, optional
    3-element list representing rotations around x, y, and z axes: [:math:`\theta_x, \theta_y, \theta_z`]

translation : list, optional
    3-element list representing a translation in space [x, y, z]

ofile : str
    Output file name

bofile : str
    Binary output file name

normal : int
    -1 if you wish to reverse the normal vector, 0 if you wish to ignore the +/- direction in computations, 1 if you with to use the direction of the normal vector as given. 

nparticles : int
    Number of particles to use for multi-particle calculations

gpu : int
    Use the gpu or not (0 or 1).  If 1 will attempt to use ALL gpus available.  This is overridden if you use the input 'ngpu'

ngpu : int or list
    If ngpu is an int, use that number of gpus (if available).
    If ngpu is a list, the list should be a list of gpus you wish to use


nthreads : int
    Number of threads to use

precision : float
    Calculation precision parameter (typically 0.01 which is 1%)

max_level: int
    Maximum "level" to use for trajectory in the calculation.  Level N corresponds to a total of 2**(N+2) trajectory points.  You cannot go beyond the internal maximum.  You are not guaranteed precision parameter is met if this is used.

max_level_extended: int
    Maximum "level" to use for trajectory in the calculation.  If set to higher than max_level the computation will proceed beyond max_level without creating trajectory arrays in memory (but it will be slower)

dim : int
    Defaults to 2 where output is in the local plane coordinates X1 and X2.  If you want the return to be given in 3D set dim=3 which will return with X, Y, and Z in absolute coordinates.

quantity: str
    Quantity to return.
    Available are:
        'power density' (default)
        'precision' - Estimated precision for each point
        'level'     - Trajectory level reached (npoints = 2**(n+1) - 1), if return is -1 the requested precision was not reached

Returns
-------
power_density : list
    A list, each element of which is a pair representing the position (2D relative or 3D absolute) and power density [:math:`W / mm^2`] at that position.  eg [[x1_0, x2_0], pd_0, [x1_1, x2_1], pd_1],  ...]

Examples
--------
Calculate the power density within a simple rectangle 1 [cm] x 1 [cm], 30 [m] downstream

    >>> osr.calculate_power_density(plane='XY', width=[0.01, 0.01], npoints=[51, 51], translation=[0, 0, 30])

Calculate the power density within a simple rectangle 1 [cm] x 1 [cm], 30 [m] downstream defined using the three-point method

    >>> osr.calculate_power_density(x0x1x2=[[-0.005, -0.005, 30], [+0.005, -0.005, 30], [-0.005, +0.005, 30]], npoints=[51, 51])

Calculate the power density on a flat surface close to and parallel to the beam direction (beam in z-direction)

    >>> power_density = osr.calculate_power_density(plane='XZ', width=[0.008, 2], npoints=[51, 101], normal=-1)
)docstring";
static PyObject* OSCARSSR_CalculatePowerDensityRectangle (OSCARSSRObject* self, PyObject* args, PyObject *keywds)
{
  // Calculate the spectrum given an observation point, and energy range

  char const* SurfacePlane = "";
  size_t      NX1 = 0;
  size_t      NX2 = 0;
  double      Width_X1 = 0;
  double      Width_X2 = 0;
  PyObject*   List_NPoints     = PyList_New(0);
  PyObject*   List_Width       = PyList_New(0);
  PyObject*   List_Translation = PyList_New(0);
  PyObject*   List_Rotations   = PyList_New(0);
  PyObject*   List_X0X1X2      = PyList_New(0);
  int         NormalDirection = 0;
  int         NParticles = 0;
  int         GPU = -1;
  PyObject*   NGPU = 0x0;
  int         NThreads = 0;
  int         Dim = 2;
  double      Precision = 0.01;
  int         MaxLevel = -2;
  int         MaxLevelExtended = 0;
  char const* ReturnQuantityChars = "power density";
  const char* OutFileNameText = "";
  const char* OutFileNameBinary = "";


  static const char *kwlist[] = {"npoints",
                                 "plane",
                                 "width",
                                 "x0x1x2",
                                 "rotations",
                                 "translation",
                                 "ofile",
                                 "bofile",
                                 "normal",
                                 "nparticles",
                                 "gpu",
                                 "ngpu",
                                 "nthreads",
                                 "precision",
                                 "max_level",
                                 "max_level_extended",
                                 "dim",
                                 "quantity",
                                  NULL};

  if (!PyArg_ParseTupleAndKeywords(args, keywds, "O|sOOOOssiiiOidiiis",
                                   const_cast<char **>(kwlist),
                                   &List_NPoints,
                                   &SurfacePlane,
                                   &List_Width,
                                   &List_X0X1X2,
                                   &List_Rotations,
                                   &List_Translation,
                                   &OutFileNameText,
                                   &OutFileNameBinary,
                                   &NormalDirection,
                                   &NParticles,
                                   &GPU,
                                   &NGPU,
                                   &NThreads,
                                   &Precision,
                                   &MaxLevel,
                                   &MaxLevelExtended,
                                   &Dim,
                                   &ReturnQuantityChars)) {
    return NULL;
  }

  // Check if a beam is at least defined
  if (self->obj->GetNParticleBeams() < 1) {
    PyErr_SetString(PyExc_ValueError, "No particle beam defined");
    return NULL;
  }


  // Check requested dimension
  if (Dim != 2 && Dim != 3) {
    PyErr_SetString(PyExc_ValueError, "'dim' must be 2 or 3");
    return NULL;
  }


  // The rectangular surface object we'll use
  TSurfacePoints_Rectangle Surface;

  if (PyList_Size(List_NPoints) == 2) {
    // NPoints in [m]
    NX1 = PyLong_AsSsize_t(PyList_GetItem(List_NPoints, 0));
    NX2 = PyLong_AsSsize_t(PyList_GetItem(List_NPoints, 1));
  } else {
    PyErr_SetString(PyExc_ValueError, "'npoints' must be [int, int]");
    return NULL;
  }



  if (NX1 <= 0 || NX2 <= 0) {
    PyErr_SetString(PyExc_ValueError, "an entry in 'npoints' is <= 0");
    return NULL;
  }

  // Vectors for rotations and translations.  Default to 0
  TVector3D Rotations(0, 0, 0);
  TVector3D Translation(0, 0, 0);

  // Check for Rotations in the input
  if (PyList_Size(List_Rotations) != 0) {
    try {
      Rotations = OSCARSPY::ListAsTVector3D(List_Rotations);
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'rotations'");
      return NULL;
    }
  }


  // Check for Translation in the input
  if (PyList_Size(List_Translation) != 0) {
    try {
      Translation = OSCARSPY::ListAsTVector3D(List_Translation);
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'translation'");
      return NULL;
    }
  }

  if (PyList_Size(List_Width) == 2) {
    // Width in [m]
    Width_X1 = PyFloat_AsDouble(PyList_GetItem(List_Width, 0));
    Width_X2 = PyFloat_AsDouble(PyList_GetItem(List_Width, 1));
  }



  // If you are requesting a simple surface plane, check that you have widths
  if (std::strlen(SurfacePlane) != 0 && Width_X1 > 0 && Width_X2 > 0) {
    try {
      Surface.Init(SurfacePlane, (int) NX1, (int) NX2, Width_X1, Width_X2, Rotations, Translation, NormalDirection);
    } catch (std::invalid_argument e) {
      PyErr_SetString(PyExc_ValueError, e.what());
      return NULL;
    }
  }



  // If X0X1X2 defined
  std::vector<TVector3D> X0X1X2;

  if (PyList_Size(List_X0X1X2) != 0) {
    if (PyList_Size(List_X0X1X2) == 3) {
      for (int i = 0; i != 3; ++i) {
        PyObject* List_X = PyList_GetItem(List_X0X1X2, i);

        try {
          X0X1X2.push_back(OSCARSPY::ListAsTVector3D(List_X));
        } catch (std::length_error e) {
          PyErr_SetString(PyExc_ValueError, "Incorrect format in 'x0x1x2'");
          return NULL;
        }
      }
    } else {
      PyErr_SetString(PyExc_ValueError, "'x0x1x2' must have 3 XYZ points defined correctly");
      return NULL;
    }

    for (std::vector<TVector3D>::iterator it = X0X1X2.begin(); it != X0X1X2.end(); ++it) {
      it->RotateSelfXYZ(Rotations);
      *it += Translation;
    }

    // UPDATE: Check for orthogonality
    Surface.Init((int) NX1, (int) NX2, X0X1X2[0], X0X1X2[1], X0X1X2[2], NormalDirection);
  }


  // Check number of particles
  if (NParticles < 0) {
    PyErr_SetString(PyExc_ValueError, "'nparticles' must be >= 1 (sort of)");
    return NULL;
  }


  // Check GPU parameter
  if (GPU != 0 && GPU != 1 && GPU != -1) {
    PyErr_SetString(PyExc_ValueError, "'gpu' must be 0 or 1");
    return NULL;
  }

  // Check NThreads parameter
  if (NThreads < 0) {
    PyErr_SetString(PyExc_ValueError, "'nthreads' must be > 0");
    return NULL;
  }

  // Check you are not trying to use threads and GPU
  if (NThreads > 0 && GPU == 1) {
    PyErr_SetString(PyExc_ValueError, "gpu is 1 and nthreads > 0.  Both are not currently allowed.");
    return NULL;
  }

  // Check ngpu input
  int NumberOfGPUs = -1;
  std::vector<int> GPUVector;
  if (NGPU != 0x0) {
    if (PyLong_Check(NGPU)) {
      NumberOfGPUs = (int) PyLong_AsLong(NGPU);
    } else if (PyList_Check(NGPU)) {
      OSCARSPY::ListToVectorInt(NGPU, GPUVector);
    }
  }

  int ReturnQuantity = 0;
  std::string ReturnQuantityStr = ReturnQuantityChars;
  std::transform(ReturnQuantityStr.begin(), ReturnQuantityStr.end(), ReturnQuantityStr.begin(), ::toupper);
  if (ReturnQuantityStr == "POWER DENSITY" || ReturnQuantityStr == "POWERDENSITY") {
    ReturnQuantity = 0;
  } else if (ReturnQuantityStr == "PRECISION") {
    ReturnQuantity = 1;
  } else if (ReturnQuantityStr == "LEVEL") {
    ReturnQuantity = 2;
  } else {
    PyErr_SetString(PyExc_ValueError, "'quantity' must be: 'power density', 'precision', 'level', or blank");
    return NULL;
  }


  // Container for Point plus scalar
  T3DScalarContainer PowerDensityContainer;

  // Actually calculate the spectrum
  bool const Directional = NormalDirection == 0 ? false : true;
  try {
    self->obj->CalculatePowerDensity(Surface,
                                     PowerDensityContainer,
                                     Dim,
                                     Directional,
                                     Precision,
                                     MaxLevel,
                                     MaxLevelExtended,
                                     NParticles,
                                     NThreads,
                                     GPU,
                                     NumberOfGPUs,
                                     GPUVector,
                                     ReturnQuantity);

  } catch (std::length_error e) {
    PyErr_SetString(PyExc_ValueError, e.what());
    return NULL;
  } catch (std::out_of_range e) {
    PyErr_SetString(PyExc_ValueError, e.what());
    return NULL;
  } catch (std::invalid_argument e) {
    PyErr_SetString(PyExc_ValueError, e.what());
    return NULL;
  }

  // If not converged print warning
  if (!PowerDensityContainer.AllConverged()) {
    OSCARSPY::PyPrint_stderr("Not all points converged to desired precision.  Can try increasing 'max_level_extended'\n");
  }

  // Write the output file if requested
  // Text output
  if (std::string(OutFileNameText) != "") {
    PowerDensityContainer.WriteToFileText(OutFileNameText, Dim);
  }

  // Binary output
  if (std::string(OutFileNameBinary) != "") {
    PowerDensityContainer.WriteToFileBinary(OutFileNameBinary, Dim);
  }


  // Build the output list of: [[[x, y, z], PowerDensity], [...]]
  // Create a python list
  PyObject *PList = PyList_New(0);

  // UPDATE: URGENT: PD output ofile, bofile
  size_t const NPoints = PowerDensityContainer.GetNPoints();

  PyObject* Value;

  for (size_t i = 0; i != NPoints; ++i) {
    T3DScalar P = PowerDensityContainer.GetPoint(i);

    // Inner list for each point
    PyObject *PList2 = PyList_New(0);


    // Add position and value to list
    Value = OSCARSPY::TVector3DAsList(P.GetX());
    PyList_Append(PList2, Value);
    Py_DECREF(Value);

    Value = Py_BuildValue("f", P.GetV());
    PyList_Append(PList2, Value);
    Py_DECREF(Value);

    PyList_Append(PList, PList2);
    Py_DECREF(PList2);

  }

  return PList;
}




















const char* DOC_OSCARSSR_CalculatePowerDensitySTL = R"docstring(
calculate_power_density_stl(ifiles [, rotations, translation, ofile, bofile, stlofile, normal, scale, nparticles, gpu, ngpu, nthreads, precision, max_level, max_level_extended, quantity])

Calculate the power density on surfaces described in STL format input files

See the :doc:`MathematicalNotes` section for the expression used in this calculation.

You **must** specify either both (*plane* and *width*) or *x0x1x2*

Parameters
----------
ifiles : list[str]
    A list of STL files to import.

rotations : list, optional
    3-element list representing rotations around x, y, and z axes: [:math:`\theta_x, \theta_y, \theta_z`]

translation : list, optional
    3-element list representing a translation in space [x, y, z]

ofile : str
    Output file name

bofile : str
    Binary output file name

stlofile : str
    STL output file name.  Will include RGB color information in the attribute byte count, maybe

normal : int
    -1 if you wish to reverse the normal vector, 0 if you wish to ignore the +/- direction in computations, 1 if you with to use the direction of the normal vector as given. 

scale : float
    What to scale the dimensions by.  Default input is in meters.

nparticles : int
    Number of particles to use for multi-particle calculations

gpu : int
    Use the gpu or not (0 or 1).  If 1 will attempt to use ALL gpus available.  This is overridden if you use the input 'ngpu'

ngpu : int or list
    If ngpu is an int, use that number of gpus (if available).
    If ngpu is a list, the list should be a list of gpus you wish to use


nthreads : int
    Number of threads to use

precision : float
    Calculation precision parameter (typically 0.01 which is 1%)

max_level: int
    Maximum "level" to use for trajectory in the calculation.  Level N corresponds to a total of 2**(N+2) trajectory points.  You cannot go beyond the internal maximum.  You are not guaranteed precision parameter is met if this is used.

max_level_extended: int
    Maximum "level" to use for trajectory in the calculation.  If set to higher than max_level the computation will proceed beyond max_level without creating trajectory arrays in memory (but it will be slower)

quantity: str
    Quantity to return.
    Available are:
        'power density' (default)
        'precision' - Estimated precision for each point
        'level'     - Trajectory level reached (npoints = 2**(n+1) - 1), if return is -1 the requested precision was not reached

Returns
-------
power_density : list
    This return list is NOT the same as other power density calculations.  The return format is as follows:
    [ [[[T0x, T0y, T0z], [T1x, T1y, T1z], [T2x, T2y, T2z]], PD], [...], ...]
    where the Ts are the triangle vertices, and PD is the power density in the center of that triangle.

Examples
--------
coming
)docstring";
static PyObject* OSCARSSR_CalculatePowerDensitySTL (OSCARSSRObject* self, PyObject* args, PyObject *keywds)
{
  // Calculate the spectrum given an observation point, and energy range

  PyObject*   List_Files = PyList_New(0);
  const char* InFileName = "";
  PyObject*   List_Translation = PyList_New(0);
  PyObject*   List_Rotations   = PyList_New(0);
  int         NormalDirection = -1;    // Default to -1 for STL
  double      Scale = 1;
  int         NParticles = 0;
  int         GPU = -1;
  PyObject*   NGPU = 0x0;
  int         NThreads = 0;
  double      Precision = 0.01;
  int         MaxLevel = -2;
  int         MaxLevelExtended = 0;
  char const* ReturnQuantityChars = "power density";
  const char* OutFileNameText = "";
  const char* OutFileNameBinary = "";
  const char* OutFileNameSTL = "";

  int const Dim = 3;


  static const char *kwlist[] = {"ifiles",
                                 "ifile",
                                 "rotations",
                                 "translation",
                                 "ofile",
                                 "bofile",
                                 "stlofile",
                                 "normal",
                                 "scale",
                                 "nparticles",
                                 "gpu",
                                 "ngpu",
                                 "nthreads",
                                 "precision",
                                 "max_level",
                                 "max_level_extended",
                                 "quantity",
                                  NULL};

  if (!PyArg_ParseTupleAndKeywords(args, keywds, "|OsOOsssidiiOidiii",
                                   const_cast<char **>(kwlist),
                                   &List_Files,
                                   &InFileName,
                                   &List_Rotations,
                                   &List_Translation,
                                   &OutFileNameText,
                                   &OutFileNameBinary,
                                   &OutFileNameSTL,
                                   &NormalDirection,
                                   &Scale,
                                   &NParticles,
                                   &GPU,
                                   &NGPU,
                                   &NThreads,
                                   &Precision,
                                   &MaxLevel,
                                   &MaxLevelExtended,
                                   &ReturnQuantityChars)) {
    return NULL;
  }

  // Check if a beam is at least defined
  if (self->obj->GetNParticleBeams() < 1) {
    PyErr_SetString(PyExc_ValueError, "No particle beam defined");
    return NULL;
  }



  // Vectors for rotations and translations.  Default to 0
  TVector3D Rotations(0, 0, 0);
  TVector3D Translation(0, 0, 0);

  // Check for Rotations in the input
  if (PyList_Size(List_Rotations) != 0) {
    try {
      Rotations = OSCARSPY::ListAsTVector3D(List_Rotations);
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'rotations'");
      return NULL;
    }
  }


  // Check for Translation in the input
  if (PyList_Size(List_Translation) != 0) {
    try {
      Translation = OSCARSPY::ListAsTVector3D(List_Translation);
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'translation'");
      return NULL;
    }
  }



  // Check number of particles
  if (NParticles < 0) {
    PyErr_SetString(PyExc_ValueError, "'nparticles' must be >= 1 (sort of)");
    return NULL;
  }


  // Check GPU parameter
  if (GPU != 0 && GPU != 1 && GPU != -1) {
    PyErr_SetString(PyExc_ValueError, "'gpu' must be 0 or 1");
    return NULL;
  }

  // Check NThreads parameter
  if (NThreads < 0) {
    PyErr_SetString(PyExc_ValueError, "'nthreads' must be > 0");
    return NULL;
  }

  // Check you are not trying to use threads and GPU
  if (NThreads > 0 && GPU == 1) {
    PyErr_SetString(PyExc_ValueError, "gpu is 1 and nthreads > 0.  Both are not currently allowed.");
    return NULL;
  }

  // Check ngpu input
  int NumberOfGPUs = -1;
  std::vector<int> GPUVector;
  if (NGPU != 0x0) {
    if (PyLong_Check(NGPU)) {
      NumberOfGPUs = (int) PyLong_AsLong(NGPU);
    } else if (PyList_Check(NGPU)) {
      OSCARSPY::ListToVectorInt(NGPU, GPUVector);
    }
  }

  int ReturnQuantity = 0;
  std::string ReturnQuantityStr = ReturnQuantityChars;
  std::transform(ReturnQuantityStr.begin(), ReturnQuantityStr.end(), ReturnQuantityStr.begin(), ::toupper);
  if (ReturnQuantityStr == "POWER DENSITY" || ReturnQuantityStr == "POWERDENSITY") {
    ReturnQuantity = 0;
  } else if (ReturnQuantityStr == "PRECISION") {
    ReturnQuantity = 1;
  } else if (ReturnQuantityStr == "LEVEL") {
    ReturnQuantity = 2;
  } else {
    PyErr_SetString(PyExc_ValueError, "'quantity' must be: 'power density', 'precision', 'level', or blank");
    return NULL;
  }

  TTriangle3DContainer STLContainer;
  try {
    STLContainer.ReadSTLFile(InFileName, Scale);
  } catch (...) {
    PyErr_SetString(PyExc_ValueError, "Cannot read STL file");
    return NULL;
  }
  STLContainer.RotateSelfXYZ(Rotations);
  STLContainer.TranslateSelf(Translation);


  TSurfacePoints_3D Surface;

  for (size_t istl = 0; istl != STLContainer.GetNPoints(); ++istl) {
    TVector3D Center = STLContainer.GetPoint(istl).GetCenter();
    if (NormalDirection < 0) {
      Surface.AddPoint(Center, -STLContainer.GetPoint(istl).GetNormal());
    } else {
      Surface.AddPoint(Center, STLContainer.GetPoint(istl).GetNormal());
    }
  }

  // Container for Point plus scalar
  T3DScalarContainer PowerDensityContainer;

  // Actually calculate the spectrum
  bool const Directional = NormalDirection == 0 ? false : true;
  try {
    self->obj->CalculatePowerDensity(Surface,
                                     PowerDensityContainer,
                                     Dim,
                                     Directional,
                                     Precision,
                                     MaxLevel,
                                     MaxLevelExtended,
                                     NParticles,
                                     NThreads,
                                     GPU,
                                     NumberOfGPUs,
                                     GPUVector,
                                     ReturnQuantity);

  } catch (std::length_error e) {
    PyErr_SetString(PyExc_ValueError, e.what());
    return NULL;
  } catch (std::out_of_range e) {
    PyErr_SetString(PyExc_ValueError, e.what());
    return NULL;
  } catch (std::invalid_argument e) {
    PyErr_SetString(PyExc_ValueError, e.what());
    return NULL;
  }

  // Write the output file if requested
  // Text output
  if (std::string(OutFileNameText) != "") {
    PowerDensityContainer.WriteToFileText(OutFileNameText, Dim);
  }

  // Binary output
  if (std::string(OutFileNameBinary) != "") {
    PowerDensityContainer.WriteToFileBinary(OutFileNameBinary, Dim);
  }

  if (std::string(OutFileNameSTL) != "") {
    STLContainer.WriteSTLFile(OutFileNameSTL);
  }
      


  // Build the output list of: [[[x, y, z], PowerDensity], [...]]
  // Create a python list
  PyObject *PList = PyList_New(0);

  // UPDATE: URGENT: PD output ofile, bofile
  size_t const NPoints = PowerDensityContainer.GetNPoints();

  PyObject* Value;

  for (size_t i = 0; i != NPoints; ++i) {
    T3DScalar P = PowerDensityContainer.GetPoint(i);
    TTriangle3D T = STLContainer.GetPoint(i);

    // Inner list for each point
    PyObject *PList2 = PyList_New(0);

    PyObject *PListT = PyList_New(0);

    Value = OSCARSPY::TVector3DAsList(T[0]);
    PyList_Append(PListT, Value);
    Py_DECREF(Value);

    Value = OSCARSPY::TVector3DAsList(T[1]);
    PyList_Append(PListT, Value);
    Py_DECREF(Value);

    Value = OSCARSPY::TVector3DAsList(T[2]);
    PyList_Append(PListT, Value);
    Py_DECREF(Value);

    // Add position and value to list
    PyList_Append(PList2, PListT);
    Py_DECREF(PListT);

    Value = Py_BuildValue("f", P.GetV());
    PyList_Append(PList2, Value);
    Py_DECREF(Value);

    PyList_Append(PList, PList2);
    Py_DECREF(PList2);

  }

  return PList;
}




















const char* DOC_OSCARSSR_CalculatePowerDensityLine = R"docstring(
calculate_power_density_line(x1, x2, normal_direction, npoints [, ofile, bofile, normal, nparticles, gpu, nthreads, precision, max_level])

See the :doc:`MathematicalNotes` section for the expression used in this calculation.

Parameters
----------
x1 : list[3]
    Point in space to start [x0, y0, z0]

x2 : list[3]
    Point in space to end [x1, y1, z1]

normal_direction : list[3]
    Direction of the normal vector for power density calculation.  Typically this is 

npoints: int
    number of in each dimension for surface

ofile : str
    Output file name

bofile : str
    Binary output file name

normal : int
    -1 if you wish to reverse the normal vector, 0 if you wish to ignore the +/- direction in computations, 1 if you with to use the direction of the normal vector as given. 

nparticles : int
    Number of particles to use for multi-particle calculations

gpu : int
    Use the gpu or not (0 or 1)

nthreads : int
    Number of threads to use

precision : float
    Calculation precision parameter (typically 0.01 which is 1%)

max_level: int
    Maximum "level" to us for trajectory in the calculation.  Level N corresponds to a total of 2**(N+2) trajectory points.  You cannot go beyond the internal maximum.  You are not guaranteed precision parameter is met if this is used.

max_level_extended: int
    Maximum "level" to use for trajectory in the calculation.  If set to higher than max_level the computation will proceed beyond max_level without creating trajectory arrays in memory (but it will be slower)

Returns
-------
power_density_1d : list
    A list, each element of which is a pair representing the length along the line specified and power density [:math:`W / mm^2`] at that position.  eg [[x0, p0], [x1, p1], ...]

)docstring";
static PyObject* OSCARSSR_CalculatePowerDensityLine (OSCARSSRObject* self, PyObject* args, PyObject *keywds)
{
  // Calculate the spectrum given an observation point, and energy range

  int         NPoints         = 0;
  PyObject*   List_x1       = PyList_New(0);
  PyObject*   List_x2         = PyList_New(0);
  int         NormalDirection = 0;
  int         NParticles = 0;
  int         GPU = -1;
  int         NThreads = 0;
  const char* OutFileNameText = "";
  const char* OutFileNameBinary = "";
  double      Precision = 0.01;
  int         MaxLevel = -2;
  int         MaxLevelExtended = 0;
  int         Dim = 1;


  static const char *kwlist[] = {"x1",
                                 "x2",
                                 "npoints",
                                 "ofile",
                                 "bofile",
                                 "normal",
                                 "nparticles",
                                 "gpu",
                                 "nthreads",
                                 "precision",
                                 "max_level",
                                 "max_level_extended",
                                 "dim",
                                  NULL};

  if (!PyArg_ParseTupleAndKeywords(args, keywds, "OO|issiiiidiii",
                                   const_cast<char **>(kwlist),
                                   &List_x1,
                                   &List_x2,
                                   &NPoints,
                                   &OutFileNameText,
                                   &OutFileNameBinary,
                                   &NormalDirection,
                                   &NParticles,
                                   &GPU,
                                   &NThreads,
                                   &Precision,
                                   &MaxLevel,
                                   &MaxLevelExtended,
                                   &Dim
                                   )) {
    return NULL;
  }

  // Check if a beam is at least defined
  if (self->obj->GetNParticleBeams() < 1) {
    PyErr_SetString(PyExc_ValueError, "No particle beam defined");
    return NULL;
  }


  // Check requested dimension
  if (Dim != 1 && Dim != 3) {
    PyErr_SetString(PyExc_ValueError, "'dim' must be 1 or 3");
    return NULL;
  }


  if (NPoints < 2) {
    PyErr_SetString(PyExc_ValueError, "'npoints' must be >= 2");
    return NULL;
  }



  // Vectors for from and to
  TVector3D x1(0, 0, 0);
  TVector3D x2(0, 0, 0);

  // Check for from in the input
  if (PyList_Size(List_x1) != 0) {
    try {
      x1 = OSCARSPY::ListAsTVector3D(List_x1);
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'x1'");
      return NULL;
    }
  }

  // Check for to in the input
  if (PyList_Size(List_x2) != 0) {
    try {
      x2 = OSCARSPY::ListAsTVector3D(List_x2);
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'x2'");
      return NULL;
    }
  }

  TVector3D const Step = (x2 - x1) / (NPoints - 1);

  // Check number of particles
  if (NParticles < 0) {
    PyErr_SetString(PyExc_ValueError, "'nparticles' must be >= 1 (sort of)");
    return NULL;
  }


  // Check GPU parameter
  if (GPU != 0 && GPU != 1 && GPU != -1) {
    PyErr_SetString(PyExc_ValueError, "'gpu' must be 0 or 1");
    return NULL;
  }

  // Check NThreads parameter
  if (NThreads < 0) {
    PyErr_SetString(PyExc_ValueError, "'nthreads' must be > 0");
    return NULL;
  }

  // Check you are not trying to use threads and GPU
  if (NThreads > 0 && GPU == 1) {
    PyErr_SetString(PyExc_ValueError, "gpu is 1 and nthreads > 0.  Both are not currently allowed.");
    return NULL;
  }

  TSurfacePoints_3D Surface;
  for (int i = 0; i != NPoints; ++i) {
    Surface.AddPoint(x1 + i * Step);
  }

  // Container for Point plus scalar
  T3DScalarContainer PowerDensityContainer;

  // Actually calculate the spectrum
  bool const Directional = NormalDirection == 0 ? false : true;
  try {
    self->obj->CalculatePowerDensity(Surface, PowerDensityContainer, Dim, Directional, Precision, MaxLevel, MaxLevelExtended, NParticles, NThreads, GPU);
  } catch (std::length_error e) {
    PyErr_SetString(PyExc_ValueError, e.what());
    return NULL;
  } catch (std::out_of_range e) {
    PyErr_SetString(PyExc_ValueError, e.what());
    return NULL;
  } catch (std::invalid_argument e) {
    PyErr_SetString(PyExc_ValueError, e.what());
    return NULL;
  }

  // Write the output file if requested
  // Text output
  if (std::string(OutFileNameText) != "") {
    PowerDensityContainer.WriteToFileText(OutFileNameText, Dim);
  }

  // Binary output
  if (std::string(OutFileNameBinary) != "") {
    PowerDensityContainer.WriteToFileBinary(OutFileNameBinary, Dim);
  }


  // Build the output list of: [[[x, y, z], PowerDensity], [...]]
  // Create a python list
  PyObject *PList = PyList_New(0);

  PyObject* Value;

  for (size_t i = 0; i != NPoints; ++i) {
    T3DScalar P = PowerDensityContainer.GetPoint(i);

    // Inner list for each point
    PyObject *PList2 = PyList_New(0);


    // Add position and value to list
    Value = OSCARSPY::TVector3DAsList(P.GetX());
    PyList_Append(PList2, Value);
    Py_DECREF(Value);

    Value = Py_BuildValue("f", P.GetV());
    PyList_Append(PList2, Value);
    Py_DECREF(Value);

    PyList_Append(PList, PList2);
    Py_DECREF(PList2);

  }

  return PList;
}




















const char* DOC_OSCARSSR_CalculateFlux = R"docstring(
calculate_flux(energy_eV, points [, normal, rotations, translation, nparticles, nthreads, gpu, ngpu, precision, max_level, max_level_extended, ofile, bofile, quantity])

Calculates the flux at a given set of points

See the :doc:`MathematicalNotes` section for the expression used in this calculation.

Parameters
----------

energy_eV : float
    Photon energy of interest

points : list
    A list of points, each point containing a position in 3D (as a list) and a normal vector at that position (also as a 3D list): [[[x, y, z], [nx. ny. nz]], [...], ...]

normal : int
    -1 if you wish to reverse the normal vector, 0 if you wish to ignore the +/- direction in computations, 1 if you with to use the direction of the normal vector as given. 

rotations : list, optional
    3-element list representing rotations around x, y, and z axes: [:math:`\theta_x, \theta_y, \theta_z`]

translation : list, optional
    3-element list representing a translation in space [x, y, z]

nparticles : int
    Number of particles to use for multi-particle calculations

nthreads : int
    Number of threads to use

gpu : int
    Use the gpu or not (0 or 1).  If 1 will attempt to use ALL gpus available.  This is overridden if you use the input 'ngpu'

ngpu : int or list
    If ngpu is an int, use that number of gpus (if available).
    If ngpu is a list, the list should be a list of gpus you wish to use

precision : float
    Calculation precision parameter (typically 0.01 which is 1%)

max_level: int
    Maximum "level" to us for trajectory in the calculation.  Level N corresponds to a total of 2**(N+2) trajectory points.  You cannot go beyond the internal maximum.  You are not guaranteed precision parameter is met if this is used.

max_level_extended: int
    Maximum "level" to us for trajectory in the calculation.  If set to higher than max_level the computation will proceed beyond max_level without creating trajectory arrays in memory (but it will be slower)

ofile : str
    Output file name

bofile : str
    Binary output file name

quantity: str
    Quantity to return.
    Available are:
        'power density' (default)
        'precision' - Estimated precision for each point
        'level'     - Trajectory level reached (npoints = 2**(n+1) - 1), if return is -1 the requested precision was not reached

Returns
-------
power_density : list
    A list, each element of which is a pair representing the position (2D relative (default) or 3D absolute) and power density [:math:`W / mm^2`] at that position.  eg [[[x1_0, x2_0, x3_0], pd_0], [[x1_1, x2_1, x3_1], pd_1]],  ...].  The position is always given as a list of length 3.  For the default (dim=2) the third element is always zero.
)docstring";
static PyObject* OSCARSSR_CalculateFlux (OSCARSSRObject* self, PyObject* args, PyObject *keywds)
{
  // Calculate the flux on a surface given an energy and list of points in 3D

  double      Energy_eV = 0;
  PyObject*   List_Translation = PyList_New(0);
  PyObject*   List_Rotations   = PyList_New(0);
  PyObject*   List_Points      = PyList_New(0);
  int         NormalDirection = 0;
  int         Dim = 3;
  int         NParticles = 0;
  int         NThreads = 0;
  int         GPU = -1;
  PyObject*   NGPU = 0x0;
  double      Precision = 0.01;
  int         MaxLevel = -2;
  int         MaxLevelExtended = 0;
  char const* OutFileNameText = "";
  char const* OutFileNameBinary = "";
  char const* ReturnQuantityChars = "flux";


  static const char *kwlist[] = {"energy_eV",
                                 "points",
                                 "normal",
                                 "rotations",
                                 "translation",
                                 "nparticles",
                                 "nthreads",
                                 "gpu",
                                 "ngpu",
                                 "precision",
                                 "max_level",
                                 "max_level_extended",
                                 "ofile",
                                 "bofile",
                                 "quantity",
                                 NULL};

  if (!PyArg_ParseTupleAndKeywords(args, keywds, "|dOiOOiiiOdiisss",
                                   const_cast<char **>(kwlist),
                                   &Energy_eV,
                                   &List_Points,
                                   &NormalDirection,
                                   &List_Rotations,
                                   &List_Translation,
                                   &NParticles,
                                   &NThreads,
                                   &GPU,
                                   &NGPU,
                                   &Precision,
                                   &MaxLevel,
                                   &MaxLevelExtended,
                                   &OutFileNameText,
                                   &OutFileNameBinary,
                                   &ReturnQuantityChars)) {
    return NULL;
  }

  // Check if a beam is at least defined
  if (self->obj->GetNParticleBeams() < 1) {
    PyErr_SetString(PyExc_ValueError, "No particle beam defined");
    return NULL;
  }


  // Check requested dimension
  if (Dim != 2 && Dim != 3) {
    PyErr_SetString(PyExc_ValueError, "'dim' must be 2 or 3");
    return NULL;
  }


  // Vectors for rotations and translations.  Default to 0
  TVector3D Rotations(0, 0, 0);
  TVector3D Translation(0, 0, 0);


  // Check for Rotations in the input
  if (PyList_Size(List_Rotations) != 0) {
    try {
      Rotations = OSCARSPY::ListAsTVector3D(List_Rotations);
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'rotations'");
      return NULL;
    }
  }


  // Check for Translation in the input
  if (PyList_Size(List_Translation) != 0) {
    try {
      Translation = OSCARSPY::ListAsTVector3D(List_Translation);
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'translation'");
      return NULL;
    }
  }

  // Look for arbitrary shape 3D points
  TSurfacePoints_3D Surface;
  for (int i = 0; i < PyList_Size(List_Points); ++i) {
    PyObject* LXN = PyList_GetItem(List_Points, i);
    TVector3D X;
    TVector3D N;
    if (PyList_Size(LXN) == 2) {

      try {
        X = OSCARSPY::ListAsTVector3D(PyList_GetItem(LXN, 0));
        N = OSCARSPY::ListAsTVector3D(PyList_GetItem(LXN, 1));
      } catch (std::length_error e) {
        PyErr_SetString(PyExc_ValueError, "Incorrect format in 'points': Point or Normal does not have 3 elements");
        return NULL;
      }

      // Rotate point and normal
      X.RotateSelfXYZ(Rotations);
      N.RotateSelfXYZ(Rotations);

      // Translate point, normal does not get translated
      X += Translation;

      Surface.AddPoint(X, N);
    } else {
      // input format error
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'points'");
      return NULL;
    }
  }

  // Check number of particles
  if (NParticles < 0) {
    PyErr_SetString(PyExc_ValueError, "'nparticles' must be >= 1 (sort of)");
    return NULL;
  }

  // Check NThreads parameter
  if (NThreads < 0) {
    PyErr_SetString(PyExc_ValueError, "'nthreads' must be > 0");
    return NULL;
  }

  // Check GPU parameter
  if (GPU != 0 && GPU != 1 && GPU != -1) {
    PyErr_SetString(PyExc_ValueError, "'gpu' must be 0 or 1");
    return NULL;
  }

  // Check you are not trying to use threads and GPU
  if (NThreads > 0 && GPU == 1) {
    PyErr_SetString(PyExc_ValueError, "gpu is 1 and nthreads > 0.  Both are not currently allowed.");
    return NULL;
  }

  // Check ngpu input
  int NumberOfGPUs = -1;
  std::vector<int> GPUVector;
  if (NGPU != 0x0) {
    if (PyLong_Check(NGPU)) {
      NumberOfGPUs = (int) PyLong_AsLong(NGPU);
    } else if (PyList_Check(NGPU)) {
      OSCARSPY::ListToVectorInt(NGPU, GPUVector);
    }
  }

  int ReturnQuantity = 0;
  std::string ReturnQuantityStr = ReturnQuantityChars;
  std::transform(ReturnQuantityStr.begin(), ReturnQuantityStr.end(), ReturnQuantityStr.begin(), ::toupper);
  if (ReturnQuantityStr == "FLUX") {
    ReturnQuantity = 0;
  } else if (ReturnQuantityStr == "PRECISION") {
    ReturnQuantity = 1;
  } else if (ReturnQuantityStr == "LEVEL") {
    ReturnQuantity = 2;
  } else {
    PyErr_SetString(PyExc_ValueError, "'quantity' must be: 'flux', 'precision', 'level', or blank");
    return NULL;
  }


  // Container for Point plus scalar
  T3DScalarContainer FluxContainer;

  try {
    throw;
    // UPDATE: Must fix single flux to accept polarizaton and angle
    self->obj->CalculateFlux(Surface, Energy_eV, FluxContainer, "all", 0, TVector3D(1, 0, 0), TVector3D(0, 1, 0), NParticles, NThreads, GPU, NumberOfGPUs, GPUVector, Precision, MaxLevel, MaxLevelExtended, Dim, ReturnQuantity);
  } catch (std::length_error e) {
    PyErr_SetString(PyExc_ValueError, e.what());
    return NULL;
  } catch (std::out_of_range e) {
    PyErr_SetString(PyExc_ValueError, e.what());
    return NULL;
  } catch (std::invalid_argument e) {
    PyErr_SetString(PyExc_ValueError, e.what());
    return NULL;
  }

  if (!FluxContainer.AllConverged()) {
    OSCARSPY::PyPrint_stderr("Not all points converged to desired precision.  Can try increasing 'max_level_extended'\n");
  }

  // Write the output file if requested
  // Text output
  if (std::string(OutFileNameText) != "") {
    FluxContainer.WriteToFileText(OutFileNameText, Dim);
  }

  // Binary output
  if (std::string(OutFileNameBinary) != "") {
    FluxContainer.WriteToFileBinary(OutFileNameBinary, Dim);
  }

  // Build the output list of: [[[x, y, z], Flux], [...]]
  // Create a python list
  PyObject *PList = PyList_New(0);

  size_t const NPoints = FluxContainer.GetNPoints();

  PyObject* Value;

  for (size_t i = 0; i != NPoints; ++i) {
    T3DScalar P = FluxContainer.GetPoint(i);

    // Inner list for each point
    PyObject *PList2 = PyList_New(0);


    // Add position and value to list
    Value = OSCARSPY::TVector3DAsList(P.GetX());
    PyList_Append(PList2, Value);
    Py_DECREF(Value);

    Value = Py_BuildValue("f", P.GetV());
    PyList_Append(PList2, Value);
    Py_DECREF(Value);

    PyList_Append(PList, PList2);
    Py_DECREF(PList2);

  }

  return PList;
}











const char* DOC_OSCARSSR_CalculateFluxRectangle = R"docstring(
calculate_flux_rectangle(energy_eV, npoints [, plane, normal, dim, width, rotations, translation, x0x1x2, polarization, angle, horizontal_direction=[1, 0, 0], propogation_direction=[0, 0, 1], nparticles, nthreads, gpu, ngpu, precision, max_level, max_level_extended, quantity, ofile, bofile])

Calculate the flux density in a rectangle either defined by three points, or by defining the plane the rectangle is in and the width, and then rotating and translating it to where it needs be.  The simplest is outlined in the first example below.  By default (dim=2) this returns a list whose position coordinates are in the local coordinate space x1 and x2 (*ie* they do not include the rotations and translation).  if dim=3 the coordinates in the return list are in absolute 3D space.

You **must** specify either both (*plane* and *width*) or *x0x1x2*

See the :doc:`MathematicalNotes` section for the expression used in this calculation.

Parameters
----------
energy_eV : float
    Photon energy of interest

npoints : list [int, int]
    Number of points in X1 and X2 dimension [n1, n2]

plane : str
    The plane to start in (XY, XZ, YZ, YX, ZX, ZY).  The normal to the surface is defined using the right handed cross product (ie the last three have opposite normal vectors from the first three)

normal : int
    -1 if you wish to reverse the normal vector, 0 if you wish to ignore the +/- direction in computations, 1 if you with to use the direction of the normal vector as given. 

dim : int
    Defaults to 2 where output is in the local plane coordinates X1 and X2.  If you want the return to be given in 3D set dim=3 which will return with X, Y, and Z in absolute coordinates.

width : list
    Width of rectangle in X1 and X2: [w1, w2]

rotations : list, optional
    3-element list representing rotations around x, y, and z axes: [:math:`\theta_x, \theta_y, \theta_z`]

translation : list, optional
    3-element list representing a translation in space [x, y, z]

x0x1x2 : list
    List of three points [[x0, y0, z0], [x1, y1, z1], [x2, y2, z2]] defining a parallelogram (vectors 0->1, and 0->2)

polarization : str
    Which polarization mode to calculate.  Can be 'all', 'linear-horizontal', 'linear-vertical', 'circular-left', 'circular-right', or 'linear' (if linear you must specify the angle parameter).  The short versions are 'lh', 'lv', 'cl', 'cr'.  Theses are not case-sensative.   There is no need to specify 'linear' if you give the angle parameter.

angle : float
    Only used if polarization='linear' is specified.  The 'angle' is that from the horizontal_direction for the polarization directino you are interested in

horizontal_direction : list
    The direction you consider to be horizontal.  Should be perpendicular to the photon beam propogation direction.  Default is [1, 0, 0]

propogation_direction : list
    Propogation direction of photon beam. Default is [0, 0, 1]

nparticles : int
    The number of particles you wish to run for a multi-particle simulation

nthreads : int
    Number of threads to use

gpu : int
    Use the gpu or not (0 or 1).  If 1 will attempt to use ALL gpus available.  This is overridden if you use the input 'ngpu'

ngpu : int or list
    If ngpu is an int, use that number of gpus (if available).
    If ngpu is a list, the list should be a list of gpus you wish to use

precision : float
    Calculation precision parameter (typically 0.01 which is 1%)

max_level: int
    Maximum "level" to use for trajectory in the calculation.  Level N corresponds to a total of 2**(N+2) trajectory points.  You cannot go beyond the internal maximum.  You are not guaranteed precision parameter is met if this is used.

max_level_extended: int
    Maximum "level" to use for trajectory in the calculation.  If set to higher than max_level the computation will proceed beyond max_level without creating trajectory arrays in memory (but it will be slower)

quantity: str
    Quantity to return.
    Available are:
        'flux'        (default)
        'precision' - Estimated precision for each point
        'level'     - Trajectory level reached (npoints = 2**(n+1) - 1), if return is -1 the requested precision was not reached

ofile : str
    Output file name

bofile : str
    Binary output file name

Returns
-------
flux : list
    A list, each element of which is a pair representing the position (2D relative (default) or 3D absolute) and flux [:math:`W / mm^2`] at that position.  eg [[[x1_0, x2_0, x3_0], f_0], [[x1_1, x2_1, x3_1], f_1]],  ...].  The position is always given as a list of length 3.  For the default (dim=2) the third element is always zero.
)docstring";
static PyObject* OSCARSSR_CalculateFluxRectangle (OSCARSSRObject* self, PyObject* args, PyObject *keywds)
{
  // Calculate the spectrum given an observation point, and energy range

  char const* SurfacePlane = "";
  size_t      NX1 = 0;
  size_t      NX2 = 0;
  double      Width_X1 = 0;
  double      Width_X2 = 0;
  PyObject*   List_NPoints= PyList_New(0);
  PyObject*   List_Width= PyList_New(0);
  PyObject*   List_Translation = PyList_New(0);
  PyObject*   List_Rotations = PyList_New(0);
  PyObject*   List_X0X1X2 = PyList_New(0);
  int         NormalDirection = 0;
  int         Dim = 2;
  double      Energy_eV = 0;
  char const* PolarizationIn = "";
  double      Angle                     = 112188979912321;
  PyObject*   List_HorizontalDirection  = 0x0;
  PyObject*   List_PropogationDirection = 0x0;
  int         NParticles = 0;
  int         NThreads = 0;
  int         GPU = -1;
  PyObject*   NGPU = 0x0;
  double      Precision = 0.01;
  int         MaxLevel = -2;
  int         MaxLevelExtended = 0;
  char const* ReturnQuantityChars = "flux";
  char const* OutFileNameText = "";
  char const* OutFileNameBinary = "";


  static const char *kwlist[] = {"energy_eV",
                                 "npoints",
                                 "plane",
                                 "normal",
                                 "dim",
                                 "width",
                                 "rotations",
                                 "translation",
                                 "x0x1x2",
                                 "polarization",
                                 "angle",
                                 "horizontal_direction",
                                 "propogation_direction",
                                 "nparticles",
                                 "nthreads",
                                 "gpu",
                                 "ngpu",
                                 "precision",
                                 "max_level",
                                 "max_level_extended",
                                 "quantity",
                                 "ofile",
                                 "bofile",
                                 NULL};

  if (!PyArg_ParseTupleAndKeywords(args, keywds, "dO|siiOOOOsdOOiiiOdiisss",
                                   const_cast<char **>(kwlist),
                                   &Energy_eV,
                                   &List_NPoints,
                                   &SurfacePlane,
                                   &NormalDirection,
                                   &Dim,
                                   &List_Width,
                                   &List_Rotations,
                                   &List_Translation,
                                   &List_X0X1X2,
                                   &PolarizationIn,
                                   &Angle,
                                   &List_HorizontalDirection,
                                   &List_PropogationDirection,
                                   &NParticles,
                                   &NThreads,
                                   &GPU,
                                   &NGPU,
                                   &Precision,
                                   &MaxLevel,
                                   &MaxLevelExtended,
                                   &ReturnQuantityChars,
                                   &OutFileNameText,
                                   &OutFileNameBinary)) {
    return NULL;
  }

  // Check if a beam is at least defined
  if (self->obj->GetNParticleBeams() < 1) {
    PyErr_SetString(PyExc_ValueError, "No particle beam defined");
    return NULL;
  }


  // Check requested dimension
  if (Dim != 2 && Dim != 3) {
    PyErr_SetString(PyExc_ValueError, "'dim' must be 2 or 3");
    return NULL;
  }


  // The rectangular surface object we'll use
  TSurfacePoints_Rectangle Surface;

  if (PyList_Size(List_NPoints) == 2) {
    // NPoints in [m]
    NX1 = PyLong_AsSsize_t(PyList_GetItem(List_NPoints, 0));
    NX2 = PyLong_AsSsize_t(PyList_GetItem(List_NPoints, 1));
  } else {
    PyErr_SetString(PyExc_ValueError, "'npoints' must be [int, int]");
    return NULL;
  }



  if (NX1 <= 0 || NX2 <= 0) {
    PyErr_SetString(PyExc_ValueError, "an entry in 'npoints' is <= 0");
    return NULL;
  }

  // Vectors for rotations and translations.  Default to 0
  TVector3D Rotations(0, 0, 0);
  TVector3D Translation(0, 0, 0);

  // Check for Rotations in the input
  if (PyList_Size(List_Rotations) != 0) {
    try {
      Rotations = OSCARSPY::ListAsTVector3D(List_Rotations);
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'rotations'");
      return NULL;
    }
  }


  // Check for Translation in the input
  if (PyList_Size(List_Translation) != 0) {
    try {
      Translation = OSCARSPY::ListAsTVector3D(List_Translation);
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'translation'");
      return NULL;
    }
  }

  if (PyList_Size(List_Width) == 2) {
    // Width in [m]
    Width_X1 = PyFloat_AsDouble(PyList_GetItem(List_Width, 0));
    Width_X2 = PyFloat_AsDouble(PyList_GetItem(List_Width, 1));
  }

  // Check normal
  if (abs(NormalDirection) > 1) {
    PyErr_SetString(PyExc_ValueError, "'normal' must be -1, 0, or 1");
    return NULL;
  }


  // If you are requesting a simple surface plane, check that you have widths
  if (std::strlen(SurfacePlane) != 0 && Width_X1 > 0 && Width_X2 > 0) {
    try {
      Surface.Init(SurfacePlane, (int) NX1, (int) NX2, Width_X1, Width_X2, Rotations, Translation, NormalDirection);
    } catch (std::invalid_argument e) {
      PyErr_SetString(PyExc_ValueError, e.what());
      return NULL;
    }
  }



  // If X0X1X2 defined
  std::vector<TVector3D> X0X1X2;

  if (PyList_Size(List_X0X1X2) != 0) {
    if (PyList_Size(List_X0X1X2) == 3) {
      for (int i = 0; i != 3; ++i) {
        PyObject* List_X = PyList_GetItem(List_X0X1X2, i);

        try {
          X0X1X2.push_back(OSCARSPY::ListAsTVector3D(List_X));
        } catch (std::length_error e) {
          PyErr_SetString(PyExc_ValueError, "Incorrect format in 'x0x1x2'");
          return NULL;
        }
      }
    } else {
      PyErr_SetString(PyExc_ValueError, "'x0x1x2' must have 3 XYZ points defined correctly");
      return NULL;
    }

    for (std::vector<TVector3D>::iterator it = X0X1X2.begin(); it != X0X1X2.end(); ++it) {
      it->RotateSelfXYZ(Rotations);
      *it += Translation;
    }

    // UPDATE: Check for orthogonality
    Surface.Init((int) NX1, (int) NX2, X0X1X2[0], X0X1X2[1], X0X1X2[2], NormalDirection);
  }



  // Check for HorizontalDirection in the input
  TVector3D HorizontalDirection(1, 0, 0);
  if (List_HorizontalDirection != 0x0 && PyList_Size(List_HorizontalDirection) != 0) {
    try {
      HorizontalDirection = OSCARSPY::ListAsTVector3D(List_HorizontalDirection);
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'translation'");
      return NULL;
    }
  }


  // Check for PropogationDirection in the input
  TVector3D PropogationDirection(0, 0, 1);
  if (List_PropogationDirection != 0x0 && PyList_Size(List_PropogationDirection) != 0) {
    try {
      PropogationDirection = OSCARSPY::ListAsTVector3D(List_PropogationDirection);
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'translation'");
      return NULL;
    }
  }


  // Check number of particles
  if (NParticles < 0) {
    PyErr_SetString(PyExc_ValueError, "'nparticles' must be >= 1 (sort of)");
    return NULL;
  }

  // Check NThreads parameter
  if (NThreads < 0) {
    PyErr_SetString(PyExc_ValueError, "'nthreads' must be > 0");
    return NULL;
  }

  // Check GPU parameter
  if (GPU != 0 && GPU != 1 && GPU != -1) {
    PyErr_SetString(PyExc_ValueError, "'gpu' must be 0 or 1");
    return NULL;
  }

  // Check you are not trying to use threads and GPU
  if (NThreads > 0 && GPU == 1) {
    PyErr_SetString(PyExc_ValueError, "gpu is 1 and nthreads > 0.  Both are not currently allowed.");
    return NULL;
  }

  // Check ngpu input
  int NumberOfGPUs = -1;
  std::vector<int> GPUVector;
  if (NGPU != 0x0) {
    if (PyLong_Check(NGPU)) {
      NumberOfGPUs = (int) PyLong_AsLong(NGPU);
    } else if (PyList_Check(NGPU)) {
      OSCARSPY::ListToVectorInt(NGPU, GPUVector);
    }
  }

  int ReturnQuantity = 0;
  std::string ReturnQuantityStr = ReturnQuantityChars;
  std::transform(ReturnQuantityStr.begin(), ReturnQuantityStr.end(), ReturnQuantityStr.begin(), ::toupper);
  if (ReturnQuantityStr == "FLUX") {
    ReturnQuantity = 0;
  } else if (ReturnQuantityStr == "PRECISION") {
    ReturnQuantity = 1;
  } else if (ReturnQuantityStr == "LEVEL") {
    ReturnQuantity = 2;
  } else {
    PyErr_SetString(PyExc_ValueError, "'quantity' must be: 'flux', 'precision', 'level', or blank");
    return NULL;
  }


  // Container for Point plus scalar
  T3DScalarContainer FluxContainer;

  // UPDATE: Needed, directional?
  //bool const Directional = NormalDirection == 0 ? false : true;

  // If all is specified make sure no angle is specified
  if (std::string(PolarizationIn) == "all" && Angle != 112188979912321) {
    PyErr_SetString(PyExc_ValueError, "cannot specify 'all' and 'angle'.  This is to save you");
    return NULL;
  }

  // If Angle is specified and polarization is blank set it to linear
  std::string Polarization = "all";
  if (strlen(PolarizationIn) != 0) {
    Polarization = PolarizationIn;
  } else if (Angle != 112188979912321) {
    Polarization = "linear";
  }

  try {
    self->obj->CalculateFlux(Surface,
                             Energy_eV,
                             FluxContainer,
                             Polarization,
                             Angle,
                             HorizontalDirection,
                             PropogationDirection,
                             NParticles,
                             NThreads,
                             GPU,
                             NumberOfGPUs,
                             GPUVector,
                             Precision,
                             MaxLevel,
                             MaxLevelExtended,
                             Dim,
                             ReturnQuantity);

  } catch (std::length_error e) {
    PyErr_SetString(PyExc_ValueError, e.what());
    return NULL;
  } catch (std::out_of_range e) {
    PyErr_SetString(PyExc_ValueError, e.what());
    return NULL;
  } catch (std::invalid_argument e) {
    PyErr_SetString(PyExc_ValueError, e.what());
    return NULL;
  }

  if (!FluxContainer.AllConverged()) {
    OSCARSPY::PyPrint_stderr("Not all points converged to desired precision.  Can try increasing 'max_level_extended'\n");
  }

  // Write the output file if requested
  // Text output
  if (std::string(OutFileNameText) != "") {
    FluxContainer.WriteToFileText(OutFileNameText, Dim);
  }

  // Binary output
  if (std::string(OutFileNameBinary) != "") {
    FluxContainer.WriteToFileBinary(OutFileNameBinary, Dim);
  }


  // Build the output list of: [[[x, y, z], Flux], [...]]
  // Create a python list
  PyObject *PList = PyList_New(0);

  size_t const NPoints = FluxContainer.GetNPoints();

  PyObject* Value;
  for (size_t i = 0; i != NPoints; ++i) {
    T3DScalar P = FluxContainer.GetPoint(i);

    // Inner list for each point
    PyObject *PList2 = PyList_New(0);


    // Add position and value to list
    Value = OSCARSPY::TVector3DAsList(P.GetX());
    PyList_Append(PList2, Value);
    Py_DECREF(Value);

    Value = Py_BuildValue("f", P.GetV());
    PyList_Append(PList2, Value);
    Py_DECREF(Value);

    PyList_Append(PList, PList2);
    Py_DECREF(PList2);

  }

  return PList;
}













const char* DOC_OSCARSSR_AverageSpectra = R"docstring(
average_spectra([, ifiles, bifiles, ofile, bofile])

Average spectra from different files.  The input files must have the same format.

Parameters
----------
ifiles : list
    The input file names as strings: ['f0.txt', 'f1.txt', ...]

bifiles : list
    The binary input file names as strings: ['f0.dat', 'f1.dat', ...]

ofile : str
    The output file name

bofile : str
    The binary output file name

Returns
-------
spectrum : list
    A list of 2D lists, each of which is a pair representing the energy and flux at that energy.  eg [[energy_0, flux_0], [energy_1, flux_1], ...]
)docstring";
static PyObject* OSCARSSR_AverageSpectra (OSCARSSRObject* self, PyObject* args, PyObject *keywds)
{
  // Calculate the flux on a surface given an energy and list of points in 3D

  PyObject*   List_InFileNamesText = PyList_New(0);
  PyObject*   List_InFileNamesBinary = PyList_New(0);
  char const* OutFileNameText = "";
  char const* OutFileNameBinary = "";


  static const char *kwlist[] = {"ifiles",
                                 "bifiles",
                                 "ofile",
                                 "bofile",
                                 NULL};

  if (!PyArg_ParseTupleAndKeywords(args, keywds, "|OOss",
                                   const_cast<char **>(kwlist),
                                   &List_InFileNamesText,
                                   &List_InFileNamesBinary,
                                   &OutFileNameText,
                                   &OutFileNameBinary)) {
    return NULL;
  }

  // Grab the number of input files for both text and binary lists
  size_t const NFilesText = PyList_Size(List_InFileNamesText);
  size_t const NFilesBinary = PyList_Size(List_InFileNamesBinary);

  // Doesn't allow for both binary and text input at the same time
  if (NFilesText != 0 && NFilesBinary !=0) {
    PyErr_SetString(PyExc_ValueError, "either text or binary files may be added, but not both.");
    return NULL;
  }

  // Check that there is at least one file
  if (NFilesText + NFilesBinary < 1) {
    PyErr_SetString(PyExc_ValueError, "No files given.  You need at least one file as input in a list.");
    return NULL;
  }

  // Add file names to vector
  std::vector<std::string> FileNames;
  std::vector<std::string> FileNamesBinary;
  for (size_t i = 0; i != NFilesText; ++i) {
    FileNames.push_back( OSCARSPY::GetAsString(PyList_GetItem(List_InFileNamesText, i)) );
  }
  for (size_t i = 0; i != NFilesBinary; ++i) {
    FileNamesBinary.push_back( OSCARSPY::GetAsString(PyList_GetItem(List_InFileNamesBinary, i)) );
  }

  // Container for flux average
  TSpectrumContainer Container;

  // Check that only one type is added.  Can update to both if useful later
  if (NFilesText > 0 && NFilesBinary > 0) {
    PyErr_SetString(PyExc_ValueError, "Currently adding mixed types of binary and text files is not supported.");
    return NULL;
  }

  // Either they are text files or binary files
  if (NFilesText > 0) {
    try {
      Container.AverageFromFilesText(FileNames);
    } catch (std::invalid_argument e) {
      PyErr_SetString(PyExc_ValueError, e.what());
      return NULL;
    }
  } else {
    try {
      Container.AverageFromFilesBinary(FileNamesBinary);
    } catch (std::invalid_argument e) {
      PyErr_SetString(PyExc_ValueError, e.what());
      return NULL;
    }
  }

  // Text output
  if (std::string(OutFileNameText) != "") {
    Container.WriteToFileText(OutFileNameText);
  }

  // Binary output
  if (std::string(OutFileNameBinary) != "") {
    Container.WriteToFileBinary(OutFileNameBinary);
  }


  return OSCARSPY::GetSpectrumAsList(Container);
}








const char* DOC_OSCARSSR_AddToSpectrum = R"docstring(
add_to_spectrum(spectrum [, weight])

Add a spectrum to the current spectrum with a given weight.

Parameters
----------
spectrum : list
    A list of pairs of numbers (spectrum format)

weight : float
    Weight for *this* spectrum

Returns
-------
None
)docstring";
static PyObject* OSCARSSR_AddToSpectrum (OSCARSSRObject* self, PyObject* args, PyObject *keywds)
{
  // Calculate the flux on a surface given an energy and list of points in 3D

  PyObject*   List_Spectrum = PyList_New(0);
  double Weight = 1;


  static const char *kwlist[] = {"spectrum",
                                 "weight",
                                 NULL};

  if (!PyArg_ParseTupleAndKeywords(args, keywds, "O|d",
                                   const_cast<char **>(kwlist),
                                   &List_Spectrum,
                                   &Weight)) {
    return NULL;
  }


  // Check if there is an input spectrum

  if (PyList_Size(List_Spectrum) < 1) {
    PyErr_SetString(PyExc_ValueError, "No points in spectrum.");
    return NULL;
  }
  TSpectrumContainer S = OSCARSPY::GetSpectrumFromList(List_Spectrum);

  self->obj->AddToSpectrum(S, Weight);

  // Must return python object None in a special way
  Py_INCREF(Py_None);
  return Py_None;
}




const char* DOC_OSCARSSR_GetSpectrum = R"docstring(
get_spectrum()

Get the current spectrum

Parameters
----------
None

Returns
-------
spectrum : list
    A list of 2D lists, each of which is a pair representing the energy and flux at that energy.  eg [[energy_0, flux_0], [energy_1, flux_1], ...]
)docstring";
static PyObject* OSCARSSR_GetSpectrum (OSCARSSRObject* self)
{
  // Calculate the flux on a surface given an energy and list of points in 3D

  return OSCARSPY::GetSpectrumAsList(self->obj->GetSpectrum());
}





const char* DOC_OSCARSSR_AverageT3DScalars_Flux = R"docstring(
average_flux([, ifiles, bifiles, ofile, bofile, dim])

Average from different files and output to specified file.  The input files must have the same format.

Parameters
----------
ifiles : list
    The input file names as strings: ['f0.txt', 'f1.txt', ...]

bifiles : list
    The binary input file names as strings: ['f0.dat', 'f1.dat', ...]

ofile : str
    The output file name

bofile : str
    The binary output file name

dim : int
    in 2 or 3 dimensions (default is 2)

Returns
-------
flux : list
    A list, each element of which is a pair representing the position (x) and value (v) at that position.  eg [[[x1_0, x2_0, x3_0], v_0], [[x1_1, x2_1, x3_1], v_1]],  ...].  The position is always given as a list of length 3.
)docstring";
const char* DOC_OSCARSSR_AverageT3DScalars_PowerDensity = R"docstring(
average_power_density([, ifiles, bifiles, ofile, bofile, dim])

Average from different files and output to specified file.  The input files must have the same format.

Parameters
----------
ifiles : list
    The input file names as strings: ['f0.txt', 'f1.txt', ...]

bifiles : list
    The binary input file names as strings: ['f0.dat', 'f1.dat', ...]

ofile : str
    The output file name

bofile : str
    The binary output file name

dim : int
    in 2 or 3 dimensions (default is 2)

Returns
-------
power_density : list
    A list, each element of which is a pair representing the position (x) and value (v) at that position.  eg [[[x1_0, x2_0, x3_0], v_0], [[x1_1, x2_1, x3_1], v_1]],  ...].  The position is always given as a list of length 3.
)docstring";
static PyObject* OSCARSSR_AverageT3DScalars (OSCARSSRObject* self, PyObject* args, PyObject *keywds)
{
  // Calculate the flux on a surface given an energy and list of points in 3D

  PyObject*   List_InFileNamesText = PyList_New(0);
  PyObject*   List_InFileNamesBinary = PyList_New(0);
  int         Dim = 2;
  char const* OutFileNameText = "";
  char const* OutFileNameBinary = "";


  static const char *kwlist[] = {"ifiles",
                                 "bifiles",
                                 "ofile",
                                 "bofile",
                                 "dim",
                                 NULL};

  if (!PyArg_ParseTupleAndKeywords(args, keywds, "|OOssi",
                                   const_cast<char **>(kwlist),
                                   &List_InFileNamesText,
                                   &List_InFileNamesBinary,
                                   &OutFileNameText,
                                   &OutFileNameBinary,
                                   &Dim)) {
    return NULL;
  }

  // Grab the number of input files for both text and binary lists
  size_t const NFilesText = PyList_Size(List_InFileNamesText);
  size_t const NFilesBinary = PyList_Size(List_InFileNamesBinary);

  // Doesn't allow for both binary and text input at the same time
  if (NFilesText != 0 && NFilesBinary !=0) {
    PyErr_SetString(PyExc_ValueError, "either text or binary files may be added, but not both.");
    return NULL;
  }

  // Check that there is at least one file
  if (NFilesText + NFilesBinary < 1) {
    PyErr_SetString(PyExc_ValueError, "No files given.  You need at least one file as input in a list.");
    return NULL;
  }

  // Add file names to vector
  std::vector<std::string> FileNames;
  for (size_t i = 0; i != NFilesText; ++i) {
    FileNames.push_back( OSCARSPY::GetAsString(PyList_GetItem(List_InFileNamesText, i)) );
  }
  for (size_t i = 0; i != NFilesBinary; ++i) {
    FileNames.push_back( OSCARSPY::GetAsString(PyList_GetItem(List_InFileNamesBinary, i)) );
  }


  // Container for flux average
  T3DScalarContainer Container;

  // Either they are text files or binary files
  try {
    if (NFilesText > 0) {
      Container.AverageFromFilesText(FileNames, Dim);
    } else {
      Container.AverageFromFilesBinary(FileNames, Dim);
    }
  } catch (std::length_error e) {
    PyErr_SetString(PyExc_ValueError, e.what());
    return NULL;
  } catch (std::invalid_argument e) {
    PyErr_SetString(PyExc_ValueError, e.what());
    return NULL;
  }

  // Build the output list of: [[[x, y, z], Value], [...]]
  // Create a python list
  PyObject *PList = PyList_New(0);

  // Number of points in container
  size_t const NPoints = Container.GetNPoints();

  PyObject* Value;

  for (size_t i = 0; i != NPoints; ++i) {

    // This point in container
    T3DScalar P = Container.GetPoint(i);

    // Inner list for each point
    PyObject *PList2 = PyList_New(0);


    // Add position and value to list
    Value = OSCARSPY::TVector3DAsList(P.GetX());
    PyList_Append(PList2, Value);
    Py_DECREF(Value);

    Value = Py_BuildValue("f", P.GetV());
    PyList_Append(PList2, Value);
    Py_DECREF(Value);

    PyList_Append(PList, PList2);
    Py_DECREF(PList2);

  }

  // Text output
  if (std::string(OutFileNameText) != "") {
    Container.WriteToFileText(OutFileNameText, Dim);
  }

  // Binary output
  if (std::string(OutFileNameBinary) != "") {
    Container.WriteToFileBinary(OutFileNameBinary, Dim);
  }

  return PList;
}




const char* DOC_OSCARSSR_AddToFlux = R"docstring(
add_to_flux(flux [, weight])

Add flux map to the current flux map with a given weight.

Parameters
----------
flux : list
    A list of points and fluxes (in flux format)

weight : float
    Weight for *this* flux map

Returns
-------
None
)docstring";
static PyObject* OSCARSSR_AddToFlux (OSCARSSRObject* self, PyObject* args, PyObject *keywds)
{
  // Calculate the flux on a surface given an energy and list of points in 3D

  PyObject*   List_Flux = PyList_New(0);
  double Weight = 1;


  static const char *kwlist[] = {"flux",
                                 "weight",
                                 NULL};

  if (!PyArg_ParseTupleAndKeywords(args, keywds, "O|d",
                                   const_cast<char **>(kwlist),
                                   &List_Flux,
                                   &Weight)) {
    return NULL;
  }


  // Check if there is an input spectrum

  if (PyList_Size(List_Flux) < 1) {
    PyErr_SetString(PyExc_ValueError, "No points in flux.");
    return NULL;
  }
  T3DScalarContainer F = OSCARSPY::GetT3DScalarContainerFromList(List_Flux);

  self->obj->AddToFlux(F, Weight);

  // Must return python object None in a special way
  Py_INCREF(Py_None);
  return Py_None;
}




const char* DOC_OSCARSSR_GetFlux = R"docstring(
get_flux()

Get the current flux map stored in memory

Parameters
----------
None

Returns
-------
flux : list
    list of [[[x, y, z], flux], ...]
)docstring";
static PyObject* OSCARSSR_GetFlux (OSCARSSRObject* self)
{
  // Return flux list

  return OSCARSSR_GetT3DScalarAsList(self->obj->GetFlux());
}





const char* DOC_OSCARSSR_AddToPowerDensity = R"docstring(
add_to_power_density(power_density [, weight])

Add power_density map to the current power_density map with a given weight.

Parameters
----------
power_density : list
    A list of points and power_densities (in power_density format)

weight : float
    Weight for *this* power_density map

Returns
-------
None
)docstring";
static PyObject* OSCARSSR_AddToPowerDensity (OSCARSSRObject* self, PyObject* args, PyObject *keywds)
{
  // Calculate the flux on a surface given an energy and list of points in 3D

  PyObject*   List_PowerDensity = PyList_New(0);
  double Weight = 1;


  static const char *kwlist[] = {"power_density",
                                 "weight",
                                 NULL};

  if (!PyArg_ParseTupleAndKeywords(args, keywds, "O|d",
                                   const_cast<char **>(kwlist),
                                   &List_PowerDensity,
                                   &Weight)) {
    return NULL;
  }


  // Check if there is an input spectrum

  if (PyList_Size(List_PowerDensity) < 1) {
    PyErr_SetString(PyExc_ValueError, "No points in flux.");
    return NULL;
  }
  T3DScalarContainer F = OSCARSPY::GetT3DScalarContainerFromList(List_PowerDensity);

  self->obj->AddToPowerDensity(F, Weight);

  // Must return python object None in a special way
  Py_INCREF(Py_None);
  return Py_None;
}




const char* DOC_OSCARSSR_GetPowerDensity = R"docstring(
get_power_density()

Get the current power_density map stored in memory

Parameters
----------
None

Returns
-------
power_density : list
    list of [[[x, y, z], power_density], ...]
)docstring";
static PyObject* OSCARSSR_GetPowerDensity (OSCARSSRObject* self)
{
  // Return flux list

  return OSCARSSR_GetT3DScalarAsList(self->obj->GetPowerDensity());
}













const char* DOC_OSCARSSR_CalculateElectricFieldTimeDomain = R"docstring(
calculate_efield_vs_time(obs [, ofile])

Calculate the electric field in the time domain for a single particle

See the :doc:`MathematicalNotes` section for the expression used in this calculation.

Parameters
----------
obs : lsit
    Point where you wish to calculate the electric field [x, y, z]

ofile : str
    Output file name

Returns
-------
efield : list
    A list, each element of which has a time (in [s]) and a 3-dimensional list representing the x, y, and z componemts of the electric field at that time: [[t, [Ex, Ey, Ez]], ...]
)docstring";
static PyObject* OSCARSSR_CalculateElectricFieldTimeDomain (OSCARSSRObject* self, PyObject* args, PyObject *keywds)
{
  // Calculate the electric field in the proper time domain.
  // The warning for using this function is that it returns unevenly spaced
  // time steps in the lab frame.  The calculation is based on the time stamps/steps
  // in the Trajectory object.  ie don't blindly employ a DFT or FFT.

  PyObject*   List_Obs = PyList_New(0);
  char const* OutFileName = "";


  static const char *kwlist[] = {"obs",
                                 "ofile",
                                 NULL};

  if (!PyArg_ParseTupleAndKeywords(args, keywds, "O|O",
                                   const_cast<char **>(kwlist),
                                   &List_Obs,
                                   &OutFileName)) {
    return NULL;
  }

  // Check if a beam is at least defined
  if (self->obj->GetNParticleBeams() < 1) {
    PyErr_SetString(PyExc_ValueError, "No particle beam defined");
    return NULL;
  }


  // Observation point
  TVector3D Obs(0, 0, 0);
  try {
      Obs = OSCARSPY::ListAsTVector3D(List_Obs);
  } catch (std::length_error e) {
    PyErr_SetString(PyExc_ValueError, "Incorrect format in 'obs'");
    return NULL;
  }

  T3DScalarContainer XYZT;
  self->obj->CalculateElectricFieldTimeDomain(Obs, XYZT);

  // UPDATE: Format is not great for XYZT output
  if (std::string(OutFileName) != "") {
    XYZT.WriteToFileText(OutFileName, 3);
  }

  // Build the output list of: [[[x, y, z], Flux], [...]]
  // Create a python list
  PyObject *PList = PyList_New(0);

  size_t const NPoints = XYZT.GetNPoints();

  PyObject* Value;

  for (size_t i = 0; i != NPoints; ++i) {
    T3DScalar P = XYZT.GetPoint(i);

    // Inner list for each point
    PyObject *PList2 = PyList_New(0);


    // Add position and value to list
    Value = Py_BuildValue("f", P.GetV());
    PyList_Append(PList2, Value);
    Py_DECREF(Value);

    Value = OSCARSPY::TVector3DAsList(P.GetX());
    PyList_Append(PList2, Value);
    Py_DECREF(Value);

    PyList_Append(PList, PList2);
    Py_DECREF(PList2);

  }

  return PList;
}





const char* DOC_OSCARSSR_PrintGPU = R"docstring(
print_gpu()

Print information about all gpus to standard out

Parameters
----------
None

Returns
-------
None
)docstring";
static PyObject* OSCARSSR_PrintGPU (OSCARSSRObject* self)
{
  // Print all magnetic stored in OSCARSSR

  int const NGPU = self->obj->CheckGPU();
  // Out string stream for printing beam information
  std::ostringstream ostream;
  ostream << "*GPUs*\n";
  ostream << "Use GPU Globally: " << self->obj->GetUseGPUGlobal() << "\n";
  ostream << "Number of GPUs: " << NGPU << "\n" << std::endl;

  if (NGPU == -1) {
    ostream << " GPU not enabled in this compiled binary\n";
  }

  for (int i = 0; i < NGPU; ++i) {
    ostream << "GPU " << i << "\n";
    ostream << self->obj->GetGPUInfo(i) << "\n";
  }
  ostream << std::endl;

  OSCARSPY::PyPrint_stdout(ostream.str());

  // Must return python object None in a special way
  Py_INCREF(Py_None);
  return Py_None;
}







const char* DOC_OSCARSSR_PrintNThreads = R"docstring(
print_gpu()

Print information about all gpus to standard out

Parameters
----------
None

Returns
-------
None
)docstring";
static PyObject* OSCARSSR_PrintNThreads (OSCARSSRObject* self)
{
  // Print all magnetic stored in OSCARSSR

  int const NThreads = self->obj->GetNThreadsGlobal();
  // Out string stream for printing beam information
  std::ostringstream ostream;
  ostream << "*NThreads Globals*\n";
  ostream << "Number of Threads to use: " << NThreads << "\n" << std::endl;

  OSCARSPY::PyPrint_stdout(ostream.str());

  // Must return python object None in a special way
  Py_INCREF(Py_None);
  return Py_None;
}







const char* DOC_OSCARSSR_PrintAll = R"docstring(
print_all()

Print all internal information related to beams and fields

Parameters
----------
None

Returns
-------
None
)docstring";
static PyObject* OSCARSSR_PrintAll (OSCARSSRObject* self)
{
  // Print beams and fields

  OSCARSSR_PrintParticleBeams(self);
  OSCARSSR_PrintMagneticFields(self);
  OSCARSSR_PrintElectricFields(self);
  OSCARSSR_PrintDriftVolumes(self);
  OSCARSSR_PrintGPU(self);
  OSCARSSR_PrintNThreads(self);

  // Must return python object None in a special way
  Py_INCREF(Py_None);
  return Py_None;
}















static PyObject* OSCARSSR_Fake (OSCARSSRObject* self, PyObject* args, PyObject *keywds)
{
    PyErr_SetString(PyExc_RuntimeError, "You must create an object to use this function: osr = oscars.sr.sr()");
    return NULL;
}


static PyMethodDef OSCARSSR_methods_fake[] = {
  {"version",                           (PyCFunction) OSCARSSR_Version, METH_NOARGS,               DOC_OSCARSSR_Version},
  {"pi",                                (PyCFunction) OSCARSSR_Fake, METH_NOARGS,                  DOC_OSCARSSR_Pi},
  {"qe",                                (PyCFunction) OSCARSSR_Fake, METH_NOARGS,                  DOC_OSCARSSR_Qe},
  {"me",                                (PyCFunction) OSCARSSR_Fake, METH_NOARGS,                  DOC_OSCARSSR_Me},
  {"rand",                              (PyCFunction) OSCARSSR_Fake, METH_NOARGS,                  DOC_OSCARSSR_Random},
  {"norm",                              (PyCFunction) OSCARSSR_Fake, METH_NOARGS,                  DOC_OSCARSSR_RandomNormal},
  {"set_seed",                          (PyCFunction) OSCARSSR_Fake, METH_O,                       DOC_OSCARSSR_SetSeed},
  {"set_gpu_global",                    (PyCFunction) OSCARSSR_Fake, METH_O,                       DOC_OSCARSSR_SetGPUGlobal},
  {"check_gpu",                         (PyCFunction) OSCARSSR_Fake, METH_NOARGS,                  DOC_OSCARSSR_CheckGPU},
  {"set_nthreads_global",               (PyCFunction) OSCARSSR_Fake, METH_O,                       DOC_OSCARSSR_SetNThreadsGlobal},
                                                                                                                            
  {"get_ctstart",                       (PyCFunction) OSCARSSR_Fake, METH_NOARGS,                  DOC_OSCARSSR_GetCTStart},
  {"get_ctstop",                        (PyCFunction) OSCARSSR_Fake, METH_NOARGS,                  DOC_OSCARSSR_GetCTStop},
  {"set_ctstartstop",                   (PyCFunction) OSCARSSR_Fake, METH_VARARGS,                 DOC_OSCARSSR_SetCTStartStop},
  {"get_npoints_trajectory",            (PyCFunction) OSCARSSR_Fake, METH_NOARGS,                  DOC_OSCARSSR_GetNPointsTrajectory},
  {"set_npoints_trajectory",            (PyCFunction) OSCARSSR_Fake, METH_O,                       DOC_OSCARSSR_SetNPointsTrajectory},
  {"set_npoints_per_meter_trajectory",  (PyCFunction) OSCARSSR_Fake, METH_O,                       DOC_OSCARSSR_SetNPointsPerMeterTrajectory},
                                                                                          
  {"add_bfield_file",                   (PyCFunction) OSCARSSR_Fake, METH_VARARGS | METH_KEYWORDS, DOC_OSCARSSR_AddMagneticField},
  {"add_bfield_interpolated",           (PyCFunction) OSCARSSR_Fake, METH_VARARGS | METH_KEYWORDS, DOC_OSCARSSR_AddMagneticFieldInterpolated},
  {"add_bfield_function",               (PyCFunction) OSCARSSR_Fake, METH_VARARGS | METH_KEYWORDS, DOC_OSCARSSR_AddMagneticFieldFunction},
  {"add_bfield_gaussian",               (PyCFunction) OSCARSSR_Fake, METH_VARARGS | METH_KEYWORDS, DOC_OSCARSSR_AddMagneticFieldGaussian},
  {"add_bfield_uniform",                (PyCFunction) OSCARSSR_Fake, METH_VARARGS | METH_KEYWORDS, DOC_OSCARSSR_AddMagneticFieldUniform},
  {"add_bfield_undulator",              (PyCFunction) OSCARSSR_Fake, METH_VARARGS | METH_KEYWORDS, DOC_OSCARSSR_AddMagneticFieldIdealUndulator},
  {"add_bfield_halbach",                (PyCFunction) OSCARSSR_Fake, METH_VARARGS | METH_KEYWORDS, DOC_OSCARSSR_AddMagneticFieldHalbach},
  {"add_bfield_quadrupole",             (PyCFunction) OSCARSSR_Fake, METH_VARARGS | METH_KEYWORDS, DOC_OSCARSSR_AddMagneticFieldQuadrupole},
  {"remove_bfield",                     (PyCFunction) OSCARSSR_Fake, METH_VARARGS | METH_KEYWORDS, DOC_OSCARSSR_RemoveMagneticField},
  {"get_bfield",                        (PyCFunction) OSCARSSR_Fake, METH_VARARGS | METH_KEYWORDS, DOC_OSCARSSR_GetBField},
  {"clear_bfields",                     (PyCFunction) OSCARSSR_Fake, METH_NOARGS,                  DOC_OSCARSSR_ClearMagneticFields},
  {"print_bfields",                     (PyCFunction) OSCARSSR_Fake, METH_NOARGS,                  DOC_OSCARSSR_PrintMagneticFields},

  {"add_efield_file",                   (PyCFunction) OSCARSSR_Fake, METH_VARARGS | METH_KEYWORDS, DOC_OSCARSSR_AddElectricField},
  {"add_efield_interpolated",           (PyCFunction) OSCARSSR_Fake, METH_VARARGS | METH_KEYWORDS, DOC_OSCARSSR_AddElectricFieldInterpolated},
  {"add_efield_function",               (PyCFunction) OSCARSSR_Fake, METH_VARARGS | METH_KEYWORDS, DOC_OSCARSSR_AddElectricFieldFunction},
  {"add_efield_gaussian",               (PyCFunction) OSCARSSR_Fake, METH_VARARGS | METH_KEYWORDS, DOC_OSCARSSR_AddElectricFieldGaussian},
  {"add_efield_uniform",                (PyCFunction) OSCARSSR_Fake, METH_VARARGS | METH_KEYWORDS, DOC_OSCARSSR_AddElectricFieldUniform},
  {"add_efield_undulator",              (PyCFunction) OSCARSSR_Fake, METH_VARARGS | METH_KEYWORDS, DOC_OSCARSSR_AddElectricFieldIdealUndulator},
  {"remove_efield",                     (PyCFunction) OSCARSSR_Fake, METH_VARARGS | METH_KEYWORDS, DOC_OSCARSSR_RemoveElectricField},
  {"get_efield",                        (PyCFunction) OSCARSSR_Fake, METH_VARARGS | METH_KEYWORDS, DOC_OSCARSSR_GetEField},
  {"clear_efields",                     (PyCFunction) OSCARSSR_Fake, METH_NOARGS,                  DOC_OSCARSSR_ClearElectricFields},
  {"print_efields",                     (PyCFunction) OSCARSSR_Fake, METH_NOARGS,                  DOC_OSCARSSR_PrintElectricFields},
 

  {"write_bfield",                      (PyCFunction) OSCARSSR_Fake, METH_VARARGS | METH_KEYWORDS, DOC_OSCARSSR_WriteMagneticField},
  {"write_efield",                      (PyCFunction) OSCARSSR_Fake, METH_VARARGS | METH_KEYWORDS, DOC_OSCARSSR_WriteElectricField},

                                                                                          
  {"set_particle_beam",                 (PyCFunction) OSCARSSR_Fake, METH_VARARGS | METH_KEYWORDS, DOC_OSCARSSR_SetParticleBeam},
  {"add_particle_beam",                 (PyCFunction) OSCARSSR_Fake, METH_VARARGS | METH_KEYWORDS, DOC_OSCARSSR_AddParticleBeam},
  {"set_particle_beam_size",            (PyCFunction) OSCARSSR_Fake, METH_VARARGS | METH_KEYWORDS, DOC_OSCARSSR_SetParticleBeamSize},
  {"clear_particle_beams",              (PyCFunction) OSCARSSR_Fake, METH_NOARGS,                  DOC_OSCARSSR_ClearParticleBeams},
  {"print_particle_beams",              (PyCFunction) OSCARSSR_Fake, METH_NOARGS,                  DOC_OSCARSSR_PrintParticleBeams},
                                                                                          
  {"set_twiss_parameters",              (PyCFunction) OSCARSSR_Fake, METH_NOARGS,                  DOC_OSCARSSR_SetTwissParameters},
                                                                                          
  {"set_new_particle",                  (PyCFunction) OSCARSSR_Fake, METH_VARARGS | METH_KEYWORDS, DOC_OSCARSSR_SetNewParticle},
  {"get_particle_x0",                   (PyCFunction) OSCARSSR_Fake, METH_NOARGS,                  DOC_OSCARSSR_GetParticleX0},
  {"get_particle_beta0",                (PyCFunction) OSCARSSR_Fake, METH_NOARGS,                  DOC_OSCARSSR_GetParticleBeta0},
  {"get_particle_e0",                   (PyCFunction) OSCARSSR_Fake, METH_NOARGS,                  DOC_OSCARSSR_GetParticleE0},
 
  {"add_drift_box",                     (PyCFunction) OSCARSSR_Fake, METH_VARARGS | METH_KEYWORDS, DOC_OSCARSSR_AddDriftVolume_Box},
  {"remove_drift",                      (PyCFunction) OSCARSSR_Fake, METH_VARARGS | METH_KEYWORDS, DOC_OSCARSSR_RemoveDriftVolume},
  {"clear_drifts",                      (PyCFunction) OSCARSSR_Fake, METH_NOARGS,                  DOC_OSCARSSR_ClearDriftVolumes},
  {"print_drifts",                      (PyCFunction) OSCARSSR_Fake, METH_NOARGS,                  DOC_OSCARSSR_PrintDriftVolumes},

  {"correct_trajectory",                (PyCFunction) OSCARSSR_Fake, METH_NOARGS,                  DOC_OSCARSSR_CorrectTrajectory},
  {"calculate_trajectory",              (PyCFunction) OSCARSSR_Fake, METH_NOARGS,                  DOC_OSCARSSR_CalculateTrajectory},
  {"get_trajectory",                    (PyCFunction) OSCARSSR_Fake, METH_NOARGS,                  DOC_OSCARSSR_GetTrajectory},

  {"calculate_spectrum",                (PyCFunction) OSCARSSR_Fake, METH_VARARGS | METH_KEYWORDS, DOC_OSCARSSR_CalculateSpectrum},

  {"calculate_total_power",             (PyCFunction) OSCARSSR_Fake, METH_VARARGS | METH_KEYWORDS, DOC_OSCARSSR_CalculateTotalPower},
  {"calculate_power_density",           (PyCFunction) OSCARSSR_Fake, METH_VARARGS | METH_KEYWORDS, DOC_OSCARSSR_CalculatePowerDensity},
  {"calculate_power_density_rectangle", (PyCFunction) OSCARSSR_Fake, METH_VARARGS | METH_KEYWORDS, DOC_OSCARSSR_CalculatePowerDensityRectangle},
  {"calculate_power_density_stl",       (PyCFunction) OSCARSSR_Fake, METH_VARARGS | METH_KEYWORDS, DOC_OSCARSSR_CalculatePowerDensitySTL},
  {"calculate_power_density_line",      (PyCFunction) OSCARSSR_Fake, METH_VARARGS | METH_KEYWORDS, DOC_OSCARSSR_CalculatePowerDensityLine},

  //{"calculate_flux",                    (PyCFunction) OSCARSSR_Fake, METH_VARARGS | METH_KEYWORDS, DOC_OSCARSSR_CalculateFlux},
  {"calculate_flux_rectangle",          (PyCFunction) OSCARSSR_Fake, METH_VARARGS | METH_KEYWORDS, DOC_OSCARSSR_CalculateFluxRectangle},

  {"average_spectra",                   (PyCFunction) OSCARSSR_Fake, METH_VARARGS | METH_KEYWORDS, DOC_OSCARSSR_AverageSpectra},
  {"add_to_spectrum",                   (PyCFunction) OSCARSSR_Fake, METH_VARARGS | METH_KEYWORDS, DOC_OSCARSSR_AddToSpectrum},
  {"get_spectrum",                      (PyCFunction) OSCARSSR_Fake, METH_VARARGS | METH_KEYWORDS, DOC_OSCARSSR_GetSpectrum},

  {"average_flux",                      (PyCFunction) OSCARSSR_Fake, METH_VARARGS | METH_KEYWORDS, DOC_OSCARSSR_AverageT3DScalars_Flux},
  {"average_power_density",             (PyCFunction) OSCARSSR_Fake, METH_VARARGS | METH_KEYWORDS, DOC_OSCARSSR_AverageT3DScalars_PowerDensity},

  {"add_to_flux",                       (PyCFunction) OSCARSSR_Fake, METH_VARARGS | METH_KEYWORDS, DOC_OSCARSSR_AddToFlux},
  {"get_flux",                          (PyCFunction) OSCARSSR_Fake, METH_VARARGS | METH_KEYWORDS, DOC_OSCARSSR_GetFlux},

  {"add_to_power_density",              (PyCFunction) OSCARSSR_Fake, METH_VARARGS | METH_KEYWORDS, DOC_OSCARSSR_AddToPowerDensity},
  {"get_power_density",                 (PyCFunction) OSCARSSR_Fake, METH_VARARGS | METH_KEYWORDS, DOC_OSCARSSR_GetPowerDensity},

  {"calculate_efield_vs_time",          (PyCFunction) OSCARSSR_Fake, METH_VARARGS | METH_KEYWORDS, DOC_OSCARSSR_CalculateElectricFieldTimeDomain},

  {"print_gpu",                         (PyCFunction) OSCARSSR_Fake, METH_NOARGS,                  DOC_OSCARSSR_PrintGPU},

  {"print_all",                         (PyCFunction) OSCARSSR_Fake, METH_NOARGS,                  DOC_OSCARSSR_PrintAll},

  {NULL, NULL, 0, NULL}  /* Sentinel */
};


static PyMethodDef OSCARSSR_methods[] = {
  // We must tell python about the function we allow access as well as give them nice
  // python names, and tell python the method of input parameters.

  {"version",                           (PyCFunction) OSCARSSR_Version,                         METH_NOARGS,                  DOC_OSCARSSR_Version},
  {"pi",                                (PyCFunction) OSCARSSR_Pi,                              METH_NOARGS,                  DOC_OSCARSSR_Pi},
  {"qe",                                (PyCFunction) OSCARSSR_Qe,                              METH_NOARGS,                  DOC_OSCARSSR_Qe},
  {"me",                                (PyCFunction) OSCARSSR_Me,                              METH_NOARGS,                  DOC_OSCARSSR_Me},
  {"rand",                              (PyCFunction) OSCARSSR_Random,                          METH_NOARGS,                  DOC_OSCARSSR_Random},
  {"norm",                              (PyCFunction) OSCARSSR_RandomNormal,                    METH_NOARGS,                  DOC_OSCARSSR_RandomNormal},
  {"set_seed",                          (PyCFunction) OSCARSSR_SetSeed,                         METH_O,                       DOC_OSCARSSR_SetSeed},
  {"set_gpu_global",                    (PyCFunction) OSCARSSR_SetGPUGlobal,                    METH_O,                       DOC_OSCARSSR_SetGPUGlobal},
  {"check_gpu",                         (PyCFunction) OSCARSSR_CheckGPU,                        METH_NOARGS,                  DOC_OSCARSSR_CheckGPU},
  {"set_nthreads_global",               (PyCFunction) OSCARSSR_SetNThreadsGlobal,               METH_O,                       DOC_OSCARSSR_SetNThreadsGlobal},
                                                                                                                            
  {"get_ctstart",                       (PyCFunction) OSCARSSR_GetCTStart,                      METH_NOARGS,                  DOC_OSCARSSR_GetCTStart},
  {"get_ctstop",                        (PyCFunction) OSCARSSR_GetCTStop,                       METH_NOARGS,                  DOC_OSCARSSR_GetCTStop},
  {"set_ctstartstop",                   (PyCFunction) OSCARSSR_SetCTStartStop,                  METH_VARARGS,                 DOC_OSCARSSR_SetCTStartStop},
  {"get_npoints_trajectory",            (PyCFunction) OSCARSSR_GetNPointsTrajectory,            METH_NOARGS,                  DOC_OSCARSSR_GetNPointsTrajectory},
  {"set_npoints_trajectory",            (PyCFunction) OSCARSSR_SetNPointsTrajectory,            METH_O,                       DOC_OSCARSSR_SetNPointsTrajectory},
  {"set_npoints_per_meter_trajectory",  (PyCFunction) OSCARSSR_SetNPointsPerMeterTrajectory,    METH_O,                       DOC_OSCARSSR_SetNPointsPerMeterTrajectory},
                                                                                          
  {"add_bfield_file",                   (PyCFunction) OSCARSSR_AddMagneticField,                METH_VARARGS | METH_KEYWORDS, DOC_OSCARSSR_AddMagneticField},
  {"add_bfield_interpolated",           (PyCFunction) OSCARSSR_AddMagneticFieldInterpolated,    METH_VARARGS | METH_KEYWORDS, DOC_OSCARSSR_AddMagneticFieldInterpolated},
  {"add_bfield_function",               (PyCFunction) OSCARSSR_AddMagneticFieldFunction,        METH_VARARGS | METH_KEYWORDS, DOC_OSCARSSR_AddMagneticFieldFunction},
  {"add_bfield_gaussian",               (PyCFunction) OSCARSSR_AddMagneticFieldGaussian,        METH_VARARGS | METH_KEYWORDS, DOC_OSCARSSR_AddMagneticFieldGaussian},
  {"add_bfield_uniform",                (PyCFunction) OSCARSSR_AddMagneticFieldUniform,         METH_VARARGS | METH_KEYWORDS, DOC_OSCARSSR_AddMagneticFieldUniform},
  {"add_bfield_undulator",              (PyCFunction) OSCARSSR_AddMagneticFieldIdealUndulator,  METH_VARARGS | METH_KEYWORDS, DOC_OSCARSSR_AddMagneticFieldIdealUndulator},
  {"add_bfield_halbach",                (PyCFunction) OSCARSSR_AddMagneticFieldHalbach,         METH_VARARGS | METH_KEYWORDS, DOC_OSCARSSR_AddMagneticFieldHalbach},
  {"add_bfield_quadrupole",             (PyCFunction) OSCARSSR_AddMagneticFieldQuadrupole,      METH_VARARGS | METH_KEYWORDS, DOC_OSCARSSR_AddMagneticFieldQuadrupole},
  {"remove_bfield",                     (PyCFunction) OSCARSSR_RemoveMagneticField,             METH_VARARGS | METH_KEYWORDS, DOC_OSCARSSR_RemoveMagneticField},
  {"get_bfield",                        (PyCFunction) OSCARSSR_GetBField,                       METH_VARARGS | METH_KEYWORDS, DOC_OSCARSSR_GetBField},
  {"clear_bfields",                     (PyCFunction) OSCARSSR_ClearMagneticFields,             METH_NOARGS,                  DOC_OSCARSSR_ClearMagneticFields},
  {"print_bfields",                     (PyCFunction) OSCARSSR_PrintMagneticFields,             METH_NOARGS,                  DOC_OSCARSSR_PrintMagneticFields},

  {"add_efield_file",                   (PyCFunction) OSCARSSR_AddElectricField,                METH_VARARGS | METH_KEYWORDS, DOC_OSCARSSR_AddElectricField},
  {"add_efield_interpolated",           (PyCFunction) OSCARSSR_AddElectricFieldInterpolated,    METH_VARARGS | METH_KEYWORDS, DOC_OSCARSSR_AddElectricFieldInterpolated},
  {"add_efield_function",               (PyCFunction) OSCARSSR_AddElectricFieldFunction,        METH_VARARGS | METH_KEYWORDS, DOC_OSCARSSR_AddElectricFieldFunction},
  {"add_efield_gaussian",               (PyCFunction) OSCARSSR_AddElectricFieldGaussian,        METH_VARARGS | METH_KEYWORDS, DOC_OSCARSSR_AddElectricFieldGaussian},
  {"add_efield_uniform",                (PyCFunction) OSCARSSR_AddElectricFieldUniform,         METH_VARARGS | METH_KEYWORDS, DOC_OSCARSSR_AddElectricFieldUniform},
  {"add_efield_undulator",              (PyCFunction) OSCARSSR_AddElectricFieldIdealUndulator,  METH_VARARGS | METH_KEYWORDS, DOC_OSCARSSR_AddElectricFieldIdealUndulator},
  {"remove_efield",                     (PyCFunction) OSCARSSR_RemoveElectricField,             METH_VARARGS | METH_KEYWORDS, DOC_OSCARSSR_RemoveElectricField},
  {"get_efield",                        (PyCFunction) OSCARSSR_GetEField,                       METH_VARARGS | METH_KEYWORDS, DOC_OSCARSSR_GetEField},
  {"clear_efields",                     (PyCFunction) OSCARSSR_ClearElectricFields,             METH_NOARGS,                  DOC_OSCARSSR_ClearElectricFields},
  {"print_efields",                     (PyCFunction) OSCARSSR_PrintElectricFields,             METH_NOARGS,                  DOC_OSCARSSR_PrintElectricFields},
 

  {"write_bfield",                      (PyCFunction) OSCARSSR_WriteMagneticField,              METH_VARARGS | METH_KEYWORDS, DOC_OSCARSSR_WriteMagneticField},
  {"write_efield",                      (PyCFunction) OSCARSSR_WriteElectricField,              METH_VARARGS | METH_KEYWORDS, DOC_OSCARSSR_WriteElectricField},

                                                                                          
  {"set_particle_beam",                 (PyCFunction) OSCARSSR_SetParticleBeam,                 METH_VARARGS | METH_KEYWORDS, DOC_OSCARSSR_SetParticleBeam},
  {"add_particle_beam",                 (PyCFunction) OSCARSSR_AddParticleBeam,                 METH_VARARGS | METH_KEYWORDS, DOC_OSCARSSR_AddParticleBeam},
  {"set_particle_beam_size",            (PyCFunction) OSCARSSR_SetParticleBeamSize,             METH_VARARGS | METH_KEYWORDS, DOC_OSCARSSR_SetParticleBeamSize},
  {"clear_particle_beams",              (PyCFunction) OSCARSSR_ClearParticleBeams,              METH_NOARGS,                  DOC_OSCARSSR_ClearParticleBeams},
  {"print_particle_beams",              (PyCFunction) OSCARSSR_PrintParticleBeams,              METH_NOARGS,                  DOC_OSCARSSR_PrintParticleBeams},
     
  {"set_twiss_parameters",              (PyCFunction) OSCARSSR_SetTwissParameters,              METH_VARARGS | METH_KEYWORDS, DOC_OSCARSSR_SetTwissParameters},

  {"set_new_particle",                  (PyCFunction) OSCARSSR_SetNewParticle,                  METH_VARARGS | METH_KEYWORDS, DOC_OSCARSSR_SetNewParticle},
  {"get_particle_x0",                   (PyCFunction) OSCARSSR_GetParticleX0,                   METH_NOARGS,                  DOC_OSCARSSR_GetParticleX0},
  {"get_particle_beta0",                (PyCFunction) OSCARSSR_GetParticleBeta0,                METH_NOARGS,                  DOC_OSCARSSR_GetParticleBeta0},
  {"get_particle_e0",                   (PyCFunction) OSCARSSR_GetParticleE0,                   METH_NOARGS,                  DOC_OSCARSSR_GetParticleE0},

  {"add_drift_box",                     (PyCFunction) OSCARSSR_AddDriftVolume_Box,              METH_VARARGS | METH_KEYWORDS, DOC_OSCARSSR_AddDriftVolume_Box},
  {"remove_drift",                      (PyCFunction) OSCARSSR_RemoveDriftVolume,               METH_VARARGS | METH_KEYWORDS, DOC_OSCARSSR_RemoveDriftVolume},
  {"clear_drifts",                      (PyCFunction) OSCARSSR_ClearDriftVolumes,               METH_NOARGS,                  DOC_OSCARSSR_ClearDriftVolumes},
  {"print_drifts",                      (PyCFunction) OSCARSSR_PrintDriftVolumes,               METH_NOARGS,                  DOC_OSCARSSR_PrintDriftVolumes},
                                                                                          
  {"correct_trajectory",                (PyCFunction) OSCARSSR_CorrectTrajectory,               METH_NOARGS,                  DOC_OSCARSSR_CorrectTrajectory},
  {"calculate_trajectory",              (PyCFunction) OSCARSSR_CalculateTrajectory,             METH_NOARGS,                  DOC_OSCARSSR_CalculateTrajectory},
  {"get_trajectory",                    (PyCFunction) OSCARSSR_GetTrajectory,                   METH_NOARGS,                  DOC_OSCARSSR_GetTrajectory},

  {"calculate_spectrum",                (PyCFunction) OSCARSSR_CalculateSpectrum,               METH_VARARGS | METH_KEYWORDS, DOC_OSCARSSR_CalculateSpectrum},

  {"calculate_total_power",             (PyCFunction) OSCARSSR_CalculateTotalPower,             METH_VARARGS | METH_KEYWORDS, DOC_OSCARSSR_CalculateTotalPower},
  {"calculate_power_density",           (PyCFunction) OSCARSSR_CalculatePowerDensity,           METH_VARARGS | METH_KEYWORDS, DOC_OSCARSSR_CalculatePowerDensity},
  {"calculate_power_density_rectangle", (PyCFunction) OSCARSSR_CalculatePowerDensityRectangle,  METH_VARARGS | METH_KEYWORDS, DOC_OSCARSSR_CalculatePowerDensityRectangle},
  {"calculate_power_density_stl",       (PyCFunction) OSCARSSR_CalculatePowerDensitySTL,        METH_VARARGS | METH_KEYWORDS, DOC_OSCARSSR_CalculatePowerDensitySTL},
  {"calculate_power_density_line",      (PyCFunction) OSCARSSR_CalculatePowerDensityLine,       METH_VARARGS | METH_KEYWORDS, DOC_OSCARSSR_CalculatePowerDensityLine},

  //{"calculate_flux",                    (PyCFunction) OSCARSSR_CalculateFlux,                   METH_VARARGS | METH_KEYWORDS, DOC_OSCARSSR_CalculateFlux},
  {"calculate_flux_rectangle",          (PyCFunction) OSCARSSR_CalculateFluxRectangle,          METH_VARARGS | METH_KEYWORDS, DOC_OSCARSSR_CalculateFluxRectangle},

  {"average_spectra",                   (PyCFunction) OSCARSSR_AverageSpectra,                  METH_VARARGS | METH_KEYWORDS, DOC_OSCARSSR_AverageSpectra},
  {"add_to_spectrum",                   (PyCFunction) OSCARSSR_AddToSpectrum,                   METH_VARARGS | METH_KEYWORDS, DOC_OSCARSSR_AddToSpectrum},
  {"get_spectrum",                      (PyCFunction) OSCARSSR_GetSpectrum,                     METH_VARARGS | METH_KEYWORDS, DOC_OSCARSSR_GetSpectrum},

  {"average_flux",                      (PyCFunction) OSCARSSR_AverageT3DScalars,               METH_VARARGS | METH_KEYWORDS, DOC_OSCARSSR_AverageT3DScalars_Flux},
  {"average_power_density",             (PyCFunction) OSCARSSR_AverageT3DScalars,               METH_VARARGS | METH_KEYWORDS, DOC_OSCARSSR_AverageT3DScalars_PowerDensity},

  {"add_to_flux",                       (PyCFunction) OSCARSSR_AddToFlux,                       METH_VARARGS | METH_KEYWORDS, DOC_OSCARSSR_AddToFlux},
  {"get_flux",                          (PyCFunction) OSCARSSR_GetFlux,                         METH_VARARGS | METH_KEYWORDS, DOC_OSCARSSR_GetFlux},

  {"add_to_power_density",              (PyCFunction) OSCARSSR_AddToPowerDensity,               METH_VARARGS | METH_KEYWORDS, DOC_OSCARSSR_AddToPowerDensity},
  {"get_power_density",                 (PyCFunction) OSCARSSR_GetPowerDensity,                 METH_VARARGS | METH_KEYWORDS, DOC_OSCARSSR_GetPowerDensity},

  {"calculate_efield_vs_time",          (PyCFunction) OSCARSSR_CalculateElectricFieldTimeDomain,METH_VARARGS | METH_KEYWORDS, DOC_OSCARSSR_CalculateElectricFieldTimeDomain},

  {"print_gpu",                         (PyCFunction) OSCARSSR_PrintGPU,                           METH_NOARGS,                  DOC_OSCARSSR_PrintGPU},

  {"print_all",                         (PyCFunction) OSCARSSR_PrintAll,                           METH_NOARGS,                  DOC_OSCARSSR_PrintAll},

  {NULL, NULL, 0, NULL}  /* Sentinel */
};


#if PY_MAJOR_VERSION >= 3
PyMODINIT_FUNC PyInit_sr(void);
#else
PyMODINIT_FUNC initsr(OSCARSSRObject* self, PyObject* args, PyObject* kwds);
#endif





#if PY_MAJOR_VERSION >= 3
static PyTypeObject OSCARSSRType = {
  PyVarObject_HEAD_INIT(NULL, 0)
  "sr",            /* tp_name */
  sizeof(OSCARSSRObject),       /* tp_basicsize */
  0,                          /* tp_itemsize */
  (destructor)OSCARSSR_dealloc, /* tp_dealloc */
  0,                          /* tp_print */
  0,                          /* tp_getattr */
  0,                          /* tp_setattr */
  0,                          /* tp_reserved */
  0,                          /* tp_repr */
  0,                          /* tp_as_number */
  0,                          /* tp_as_sequence */
  0,                          /* tp_as_mapping */
  0,                          /* tp_hash  */
  0,                          /* tp_call */
  0,                          /* tp_str */
  0,                          /* tp_getattro */
  0,                          /* tp_setattro */
  0,                          /* tp_as_buffer */
  Py_TPFLAGS_DEFAULT |
  Py_TPFLAGS_BASETYPE,        /* tp_flags */
  "oscars sr class",           /* tp_doc */
  0,                          /* tp_traverse */
  0,                          /* tp_clear */
  0,                          /* tp_richcompare */
  0,                          /* tp_weaklistoffset */
  0,                          /* tp_iter */
  0,                          /* tp_iternext */
  OSCARSSR_methods,             /* tp_methods */
  0,                          /* tp_members */
  0,                          /* tp_getset */
  0,                          /* tp_base */
  0,                          /* tp_dict */
  0,                          /* tp_descr_get */
  0,                          /* tp_descr_set */
  0,                          /* tp_dictoffset */
  0,      /* tp_init */
  0,                          /* tp_alloc */
  OSCARSSR_new,                 /* tp_new */
};
#else
static PyTypeObject OSCARSSRType = {
  // The python object.  Fully defined elsewhere.  only put here what you need,
  // otherwise default values

  PyObject_HEAD_INIT(NULL)
  0,                                        /* ob_size */
  "sr",                                 /* tp_name */
  sizeof(OSCARSSRObject),                         /* tp_basicsize */
  0,                                        /* tp_itemsize */
  (destructor) OSCARSSR_dealloc,                 /* tp_dealloc */
  0,                                        /* tp_print */
  0,                                        /* tp_getattr */
  0,                                        /* tp_setattr */
  0,                                        /* tp_compare */
  0,                                        /* tp_repr */
  0,                                        /* tp_as_number */
  0,                                        /* tp_as_sequence */
  0,                                        /* tp_as_mapping */
  0,                                        /* tp_hash */
  0,                                        /* tp_call */
  0,                                        /* tp_str */
  0,                                        /* tp_getattro */
  0,                                        /* tp_setattro */
  0,                                        /* tp_as_buffer */
  Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, /* tp_flags */
  "oscars sr class",                              /* tp_doc */
  0,                                        /* tp_traverse */
  0,                                        /* tp_clear */
  0,                                        /* tp_richcompare */
  0,                                        /* tp_weaklistoffset */
  0,                                        /* tp_iter */
  0,                                        /* tp_iternext */
  OSCARSSR_methods,                             /* tp_methods */
  0,                                        /* tp_members */
  0,                                        /* tp_getset */
  0,                                        /* tp_base */
  0,                                        /* tp_dict */
  0,                                        /* tp_descr_get */
  0,                                        /* tp_descr_set */
  0,                                        /* tp_dictoffset */
  0,                                        /* tp_init */
  0,                                        /* tp_alloc */
  OSCARSSR_new,                                  /* tp_new */
};
#endif




//static PyMethodDef module_methods[] = {
//  // I do not need
//  {NULL}  /* Sentinel */
//};


#if PY_MAJOR_VERSION >= 3
static PyModuleDef OSCARSSRmodule = {
  PyModuleDef_HEAD_INIT,
  "sr",
  "OSCARSSR module extension.",
  -1,
  OSCARSSR_methods_fake,
  NULL, NULL, NULL, NULL
};
#endif



#if PY_MAJOR_VERSION >= 3
PyMODINIT_FUNC PyInit_sr(void)
{
  if (PyType_Ready(&OSCARSSRType) < 0) {
    return NULL;
  }
  PyObject* m = PyModule_Create(&OSCARSSRmodule);
  if (m == NULL) {
    return NULL;
  }
  Py_INCREF(&OSCARSSRType);
  PyModule_AddObject(m, "sr", (PyObject *)&OSCARSSRType);


  std::string Message = "OSCARS v" + OSCARSPY::GetVersionString() + " - Open Source Code for Advanced Radiation Simulation\nBrookhaven National Laboratory, Upton NY, USA\nhttp://oscars.bnl.gov\noscars@bnl.gov\n";
  OSCARSPY::PyPrint_stdout(Message);

  return m;
}
#else
PyMODINIT_FUNC initsr(OSCARSSRObject* self, PyObject* args, PyObject* kwds)
{
  if (PyType_Ready(&OSCARSSRType) < 0) {
    return;
  }
  PyObject *m = Py_InitModule("oscars.sr", OSCARSSR_methods);
  if (m == NULL) {
    return;
  }
  Py_INCREF(&OSCARSSRType);
  PyModule_AddObject(m, "sr", (PyObject *)&OSCARSSRType);

  std::string Message = "OSCARS v" + OSCARSPY::GetVersionString() + " - Open Source Code for Advanced Radiation Simulation\nBrookhaven National Laboratory, Upton NY, USA\nhttp://oscars.bnl.gov\noscars@bnl.gov\n";
  OSCARSPY::PyPrint_stdout(Message);

  return;
}
#endif












static PyObject* OSCARSSR_GetT3DScalarAsList (T3DScalarContainer const& C)
{
  // Get the spectrum as a list format for python output

  // Create a python list
  PyObject *PList = PyList_New(0);

  // Number of points in trajectory calculation
  size_t NPoints = C.GetNPoints();

  PyObject* Value;

  // Loop over all points
  for (size_t i = 0; i != NPoints; ++i) {
    // Create a python list
    PyObject *PList2 = PyList_New(0);

    double const V = C.GetPoint(i).GetV();

    // Add position and Beta to list
    Value = OSCARSPY::TVector3DAsList(C.GetPoint(i).GetX());
    PyList_Append(PList2, Value);
    Py_DECREF(Value);

    Value = Py_BuildValue("f", V);
    PyList_Append(PList2, Value);
    Py_DECREF(Value);

    PyList_Append(PList, PList2);
    Py_DECREF(PList2);
  }

  // Return the python list
  return PList;
}






