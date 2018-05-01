# OSCARS - Open Source Code for Advanced Radiation Simulation

[![Build Status](https://travis-ci.org/dhidas/OSCARS.svg?branch=master)](https://travis-ci.org/dhidas/OSCARS)
[![Build status](https://ci.appveyor.com/api/projects/status/485457m1fdd38dar/branch/master?svg=true)](https://ci.appveyor.com/project/dhidas/oscars/branch/master)

This is the OSCARS developer repository.  The main OSCARS website is <http://oscars.bnl.gov>

## Installation - Basic

The easiest way to install OSCARS is:
```
pip install oscars
```
That is all you need to do.  It is recommended that you also install jupyter, but it is not necessary.




## Installation - Intermediate : Python virtualenv

One way to install OSCARS is using a python virtual environment:

```
# Create virtual env
python3 -m venv oscarsenv

# Activate the oscars conda environment
source oscarsenv/bin/activate

# Install OSCARS
pip install oscars

# Install other useful utilities
pip install numpy
pip install scipy
pip install matplotlib
pip install jupyter
```

At this point you should be able to run jupyter notebook (not necessairy, but nice)
```
jupyter notebook
```
and point your browser to <http://localhost:8888/> if it hasn't taken you there already.  Next, you are likely ready to try some of the examples on <http://oscars.bnl.gov/examples.php>





## Installation - Intermediate : Conda

One way to install OSCARS is using a conda environment.  First:

* [Install conda](http://conda.pydata.org/miniconda.html)

```
# Fix for temporary conda bug for some distributions
conda install -n root pyyaml

# Create oscars environment (python=3.6 is optional, but recommended)
conda create --name oscars python=3.6

# Activate the oscars conda environment
source activate oscars

# Install OSCARS
pip install oscars

# Install other useful utilities
pip install numpy
pip install scipy
pip install matplotlib
pip install jupyter
```

At this point you should be able to run jupyter notebook (not necessairy, but nice)
```
jupyter notebook
```
and point your browser to <http://localhost:8888/> if it hasn't taken you there already.  Next, you are likely ready to try some of the examples on <http://oscars.bnl.gov/examples.php>



## Installation Advanced - GPU Utilities

Compiling for the GPU requires the nvidia compiler nvcc:
* https://developer.nvidia.com/cuda-downloads

Once this is installed you need to download OSCARS and run the following:
```
make
python setup.py install
```
The first command will build the gpu library lib/OSCARSSR_Cuda.o and the second will compile the rest of OSCARS and install it.  The complete picture for this using conda is given below:
* [Install conda](http://conda.pydata.org/miniconda.html)

```
# Install git in your root conda environment (if you don't already have it).
conda install -n root git

# Fix for temporary conda bug for some distributions
conda install -n root pyyaml

# Download OSCARS
git clone https://github.com/dhidas/OSCARS -b 2.0.20

# Create a new "conda environment" and install the required Python packages.
cd OSCARS
conda env create -f environment.yml

# Activate the oscars conda environment
source activate oscars

# Build and install OSCARS with GPU support
make
python setup.py install
```
