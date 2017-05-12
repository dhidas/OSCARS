# OSCARS - Open Source Code for Advanced Radiation Simulation

This is the OSCARS developer repository.  The main OSCARS website is <http://oscars.bnl.gov>

## Installation - Basic

The easiest way to install OSCARS is:
```
pip install oscars -U
```
That is all you need to do.  It is recommended that you also install jupyter, but it is not necessary.



## Installation - Intermediate (if you want to also create a conda environment)

The best way to install OSCARS is using a conda environment, so first

* [Install conda](http://conda.pydata.org/miniconda.html)

```
# Install git in your root conda environment (if you don't already have it).
conda install -n root git

# Fix for temporary conda bug for some distributions
conda install -n root pyyaml

# Download OSCARS
git clone https://github.com/dhidas/OSCARS -b 1.36.03

# Create a new "conda environment" and install the required Python packages.
cd OSCARS
conda env create -f environment.yml

# Activate the oscars conda environment
source activate oscars

# Build and install OSCARS
python setup.py install
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
python setup_gpu.py install
```
The first command will build the gpu library lib/OSCARSSR_Cuda.o and the second will compile the rest of OSCARS and install it.  The complete picture for this using conda is given below:
* [Install conda](http://conda.pydata.org/miniconda.html)

```
# Install git in your root conda environment (if you don't already have it).
conda install -n root git

# Fix for temporary conda bug for some distributions
conda install -n root pyyaml

# Download OSCARS
git clone https://github.com/dhidas/OSCARS -b 1.36.03

# Create a new "conda environment" and install the required Python packages.
cd OSCARS
conda env create -f environment.yml

# Activate the oscars conda environment
source activate oscars

# Build and install OSCARS
python setup_gpu.py install
```
