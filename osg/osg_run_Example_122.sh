#!/bin/bash

Process=$1
NAME=$2


# Load same modules you used to compile oscars
module load gcc/4.9.2
module load python/3.5.2

python examples/Example_122_OSG.py $Process $NAME
