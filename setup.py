import os
import sys
try:
    from setuptools import setup, Extension
except ImportError:
    from distutils.core import setup, Extension

# Get Version numbering from Version.h
v_major = ''
v_minor = ''
v_rev = ''
v_rel = ''
with open('include/Version.h', 'r') as f:
    for l in f:
        if "OSCARS_VMAJOR" in l:
            s = l.split()
            v_major = s[-1]
        elif "OSCARS_VMINOR" in l:
            s = l.split()
            v_minor = s[-1]
        elif "OSCARS_REVISION" in l:
            s = l.split()
            v_rev = s[-1]
        elif "OSCARS_RELEASE" in l:
            s = l.split()
            if s[-1] != 'NULL':
                v_rev = v_rev + '.' + s[-1].strip('"')

VERSION=v_major+'.'+v_minor+'.'+v_rev

#os.environ["CC"] = "nvcc"
#os.environ["CXX"] = "nvcc"


extra_compile_args=[]
extra_objects_sr=[]
extra_objects_th=[]
library_dirs=[]
libraries=[]

# Check distribution for flags and libs
if sys.platform == "linux" or sys.platform == "linux2":
    extra_compile_args.append('-std=c++11')
    extra_compile_args.append('-O3')
    extra_compile_args.append('-fPIC')
    extra_compile_args.append('-pthread')

    if os.path.exists('lib/OSCARSSR_Cuda.o') and os.path.exists('lib/OSCARSTH_Cuda.o'):
        library_dirs.append('/usr/local/cuda/lib64')
        library_dirs.append('/lib64')
        library_dirs.append('/usr/lib64')

        extra_compile_args.append('-DCUDA')
        libraries.append('cudart_static')
        libraries.append('rt')
        extra_objects_sr.append('lib/OSCARSSR_Cuda.o')
        extra_objects_th.append('lib/OSCARSTH_Cuda.o')

elif sys.platform == 'darwin':
    extra_compile_args.append('-std=c++11')
    extra_compile_args.append('-O3')
    extra_compile_args.append('-fPIC')
    extra_compile_args.append('-pthread')

    if 'conda' not in sys.version:
        extra_compile_args.append('-mmacosx-version-min=10.9')

    if os.path.exists('lib/OSCARSSR_Cuda.o') and os.path.exists('lib/OSCARSTH_Cuda.o'):
        library_dirs.append('/usr/local/cuda/lib')

        extra_compile_args.append('-DCUDA')
        libraries.append('cudart_static')
        extra_objects_sr.append('lib/OSCARSSR_Cuda.o')
        extra_objects_th.append('lib/OSCARSTH_Cuda.o')

elif sys.platform == 'win32':
    library_dirs.append('C:\\Program Files\\NVIDIA GPU Computing Toolkit\\CUDA\\v8.0\\lib\\x64')
    if sys.version_info >= (3, 6):
        ccbin=''
    elif sys.version_info >= (3, 5):
        pass
    elif sys.version_info >= (3, 4):
        pass
    elif sys.version_info >= (2, 8):
        pass
    elif sys.version_info >= (2, 7):
        pass
    elif sys.version_info < (2, 7):
        pass

    if os.path.exists('lib/OSCARSSR_Cuda.o'):
        extra_compile_args.append('/DCUDA')
        libraries.append('cudart_static')
        extra_objects_sr.append('lib/OSCARSSR_Cuda.o')
        extra_objects_th.append('lib/OSCARSTH_Cuda.o')
    




moduleOSCARSSR = Extension('oscars.sr',
                      include_dirs = ['include'],
                      sources = ['src/OSCARSSR.cc',
                                 'src/OSCARSSR_Python.cc',
                                 'src/T3DScalarContainer.cc',
                                 'src/TField3D_Grid.cc',
                                 'src/TField3D_Gaussian.cc',
                                 'src/TFieldContainer.cc',
                                 'src/TField3D_IdealUndulator.cc',
                                 'src/TField3D_Halbach.cc',
                                 'src/TField3D_UniformBox.cc',
                                 'src/TFieldPythonFunction.cc',
                                 'src/TParticleA.cc',
                                 'src/TParticleBeam.cc',
                                 'src/TParticleBeamContainer.cc',
                                 'src/TParticleTrajectoryPoint.cc',
                                 'src/TParticleTrajectoryPoints.cc',
                                 'src/TParticleTrajectoryInterpolated.cc',
                                 'src/TParticleTrajectoryInterpolatedPoints.cc',
                                 'src/TRandomA.cc',
                                 'src/TSpectrumContainer.cc',
                                 'src/TSurfaceOfPoints.cc',
                                 'src/TSurfacePoint.cc',
                                 'src/TSurfacePoints_3D.cc',
                                 'src/TSurfacePoints_Rectangle.cc',
                                 'src/TVector2D.cc',
                                 'src/TVector3D.cc',
                                 'src/TVector3DC.cc',
                                 'src/TVector4D.cc',
                                 'src/TField3D_Quadrupole.cc',
                                 'src/TOMATH.cc',
                                 'src/TDriftBox.cc',
                                 'src/TTriangle3D.cc',
                                 'src/TTriangle3DContainer.cc',
                                 'src/TDriftVolumeContainer.cc',
                                 'src/OSCARSPY.cc'],
                      extra_compile_args=extra_compile_args,
                      libraries=libraries,
                      library_dirs=library_dirs,
                      extra_objects=extra_objects_sr
                     )


moduleOSCARSTH = Extension('oscars.th',
                      include_dirs = ['include'],
                      sources = ['src/OSCARSTH.cc',
                                 'src/OSCARSTH_Python.cc',
                                 'src/T3DScalarContainer.cc',
                                 'src/TField3D_Grid.cc',
                                 'src/TField3D_Gaussian.cc',
                                 'src/TFieldContainer.cc',
                                 'src/TField3D_IdealUndulator.cc',
                                 'src/TField3D_Halbach.cc',
                                 'src/TField3D_UniformBox.cc',
                                 'src/TFieldPythonFunction.cc',
                                 'src/TParticleA.cc',
                                 'src/TParticleBeam.cc',
                                 'src/TParticleBeamContainer.cc',
                                 'src/TParticleTrajectoryPoint.cc',
                                 'src/TParticleTrajectoryPoints.cc',
                                 'src/TParticleTrajectoryInterpolated.cc',
                                 'src/TParticleTrajectoryInterpolatedPoints.cc',
                                 'src/TRandomA.cc',
                                 'src/TSpectrumContainer.cc',
                                 'src/TSurfaceOfPoints.cc',
                                 'src/TSurfacePoint.cc',
                                 'src/TSurfacePoints_3D.cc',
                                 'src/TSurfacePoints_Rectangle.cc',
                                 'src/TVector2D.cc',
                                 'src/TVector3D.cc',
                                 'src/TVector3DC.cc',
                                 'src/TVector4D.cc',
                                 'src/TField3D_Quadrupole.cc',
                                 'src/TOMATH.cc',
                                 'src/TDriftBox.cc',
                                 'src/TTriangle3D.cc',
                                 'src/TTriangle3DContainer.cc',
                                 'src/TDriftVolumeContainer.cc',
                                 'src/OSCARSPY.cc'],
                      extra_compile_args=extra_compile_args,
                      libraries=libraries,
                      library_dirs=library_dirs,
                      extra_objects=extra_objects_th
                     )




setup(
  name="oscars",
  version=VERSION,
  description = 'Open Source Code for Advanced Radiation Simulation',
  author = 'Dean Andrew Hidas',
  author_email = 'dhidas@bnl.gov',
  url = 'http://oscars.bnl.gov/',
  license = 'LICENSE.txt',
  long_description = '''The OSCARS Package.''',
  ext_modules = [moduleOSCARSSR, moduleOSCARSTH],
  data_files=[('oscars', ['LICENSE.txt', 'COPYRIGHT.txt'])],
  package_data = {'' : ['LICENSE.txt']},
  #package_dir = {'oscars': 'python'},
  #include_package_data=True,
  py_modules = ['oscars.plots_mpl', 'oscars.plots3d_mpl', 'oscars.parametric_surfaces', 'oscars.util', 'oscars.fit', 'oscars.bl', 'oscars.lut', 'oscars.brightness', 'oscars.twiss', 'oscars.me']
)
