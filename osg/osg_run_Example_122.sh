#!/bin/bash

Process=$1
NAME=$2

echo $PWD


tar zxf osg.tgz

find .



# Load same modules you used to compile oscars
module load gcc/4.9.2
module load python/3.5.2

python3 examples/Example_122_UndulatorSpectrum_MultiParticle_OSG.py $Process $NAME
