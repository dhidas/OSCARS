#!/bin/bash
set -e -x

# Install a system package required by our library
yum install -y atlas-devel

PYALL="cp27-cp27m cp27-cp27mu cp33-cp33m cp34-cp34m cp35-cp35m cp36-cp36m"

# Compile wheels
for PYBIN in $PYALL; do
    echo "/opt/python/${PYBIN}/bin/pip"
    "/opt/python/${PYBIN}/bin/pip" install -r /io/dev-requirements.txt
    "/opt/python/${PYBIN}/bin/pip" wheel /io/ -w wheelhouse/
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
