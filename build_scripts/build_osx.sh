#!/bin/bash

# Make the gpu cuda library
make

source ~/venv/py2.7/bin/activate
python2.7 setup.py bdist_wheel --dist-dir wheelhouse
deactivate

source ~/venv/py3.4/bin/activate
python setup.py bdist_wheel --dist-dir wheelhouse
deactivate

source ~/venv/py3.5/bin/activate
python setup.py bdist_wheel --dist-dir wheelhouse
deactivate

source ~/venv/py3.6/bin/activate
python setup.py bdist_wheel --dist-dir wheelhouse
deactivate


source ~/venv/py3.6/bin/activate
twine upload -r pypi wheelhouse/*macosx*whl
deactivate
