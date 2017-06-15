#!/bin/bash

# This script should be run using the standard centos5 manylinux docker image

PYALL="cp27-cp27m cp27-cp27mu cp33-cp33m cp34-cp34m cp35-cp35m cp36-cp36m"

rm -f /wheelhouse/*whl

#for PYBIN in /opt/python/*/bin; do
for PYBIN in $PYALL; do
    echo "/opt/python/${PYBIN}/python"
    "/opt/python/${PYBIN}/bin/python" setup.py bdist_wheel --dist-dir wheelhouse
done


# Bundle external shared libraries into the wheels
for whl in wheelhouse/*.whl; do
    auditwheel repair "$whl" -w /wheelhouse/
    rm "$whl"
done

/opt/python/cp35-cp35m/bin/twine upload -r pypi /wheelhouse/*manylinux*whl
