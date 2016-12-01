# OSCARS

## Open Source Code for Advanced Radiation Simulation


## Installation

## Before you begin

* [Install conda](http://conda.pydata.org/miniconda.html)

## Getting Started

### One-time installation of tutorial software and scripts

On Windows, open "Anaconda Prompt". (You can search for it from the Start Menu.)
On Mac or Linux, open any terminal.

```
# Update conda (not necessary if you *just* installed it)
conda update -n root conda

# Install git in your root conda environment (if you don't already have it).
conda install -n root git

# Download the materials.
git clone https://github.com/dhidas/OSCARS

# Create a new "conda environment" and install the required Python packages.
cd OSCARS
conda env create -f environment.yml
```


Once you have a proper environment you can build and intall the basic OSCARS package (without GPU utilities) by simply:
```
# Build and install OSCARS
python setup.py install
```
