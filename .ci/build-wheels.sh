#!/bin/bash
set -e -x


PYALL="cp27-cp27m cp27-cp27mu cp33-cp33m cp34-cp34m cp35-cp35m cp36-cp36m"

echo In build-wheels.sh
pwd
ls .

mkdir wheelhouse
ls /io

cd /io

mkdir wheelhouse

ls -la /opt/python/
# Compile wheels
#for PYBIN in $PYALL; do
for PYBIN in /opt/python/*/bin; do
    #echo "/opt/python/${PYBIN}/bin/pip"
    #"/opt/python/${PYBIN}/bin/pip" install -r /io/requirements.txt
    "${PYBIN}/pip" install -r /io/requirements.txt
    #"/opt/python/${PYBIN}/bin/pip" wheel /io/ -w wheelhouse/
    #"/opt/python/${PYBIN}/bin/python" setup.py bdist_wheel --dist-dir wheelhouse
    "${PYBIN}/python" setup.py bdist_wheel --dist-dir wheelhouse
done

# Bundle external shared libraries into the wheels
for whl in wheelhouse/*.whl; do
    auditwheel repair "$whl" -w /io/wheelhouse/
done

# Install packages and test
#for PYBIN in /opt/python/*/bin/; do
#    "${PYBIN}/pip" install python-manylinux-demo --no-index -f /io/wheelhouse
#    (cd "$HOME"; "${PYBIN}/nosetests" pymanylinuxdemo)
#done
