# OSCARS - Open Source Code for Advanced Radiation Simulation


## Installation

The best way to install OSCARS is using a conda environment, so first

* [Install conda](http://conda.pydata.org/miniconda.html)


## Getting Started

```
# Install git in your root conda environment (if you don't already have it).
conda install -n root git

# Download the materials.
git clone https://github.com/dhidas/OSCARS

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
and point your browser to <http://localhost:8888/> if it hasn't taken you there already.  Next, you are likely ready to try some of the examples on <http://oscars.bnl.gov/examples>
