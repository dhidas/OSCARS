# Instructions for running on OSG (Open Science Grid)


## Compilation for OSG
In order to compile the oscars package for osg we must setup the correct environment as well as compile and install oscars locally.

```
# Load the desired modules
module load gcc/4.9.2
module load python/3.5.2

# Compile OSCARS locally
python setup.py build_ext --inplace
```
