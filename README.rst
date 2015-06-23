simphony-jyulb
==============

Internal and file-IO based JYU-LB engine-wrappers for the SimPhoNy framework.

.. image:: https://travis-ci.org/simphony/simphony-jyulb.svg?branch=master
    :target: https://travis-ci.org/simphony/simphony-jyulb

.. image:: https://coveralls.io/repos/simphony/simphony-jyulb/badge.svg
   :target: https://coveralls.io/r/simphony/simphony-jyulb

Repository
----------

Simphony-jyulb is hosted on github: https://github.com/simphony/simphony-jyulb

Requirements
------------
- `simphony-common`_ >= 0.1.1

.. _simphony-common: https://github.com/simphony/simphony-common

- `cython`_

.. _cython: https://pypi.python.org/pypi/Cython/

- `JYU-LB`_ >= 0.1.1

.. _JYU-LB: https://github.com/simphony/JYU-LB

.. note::
  simphony-jyulb uses keywords (COLLISION_OPERATOR, REFERENCE_DENSITY, GRAVITY,
  FLOW_TYPE, EXTERNAL_FORCING) which currently do not belong to the CUBA-keywords.

  simphony-jyulb treats the BC data component as a dummy variable (currently there
  is no SimPhoNy definition for BC).

Installation
------------

The package requires python 2.7.x, cython and JYU-LB header and source files along with
`simphony-common` package to build and install cython extension. To get the
required files clone the repository and install dependencies::

    # Clone recursively
    git clone --recursive https://github.com/simphony/simphony-jyulb.git

    # Install required dependencies
    pip install -r requirements.txt

The next step is straight forward and is based on setuptools::

    # build and install
    python setup.py install

or::

    # build for in-place development
    python setup.py develop

.. note::
  In case of changes to JYU-LB repository run the below command to fetch the changes::

    # Get submodule's changes
    git submodule update --remote

  For more information on git submodules see: https://git-scm.com/book/en/v2/Git-Tools-Submodules

JYU-LB installation
~~~~~~~~~~~~~~~~~~~

These engine-wrappers use JYU-LB fluid flow simulators and the related source code.

The JYU-LB source code is available on github: https://github.com/simphony/JYU-LB

The file-IO based engine wrapper assumes that there is an executable called
"jyu_lb_isothermal.exe" that can be found in the PATH and an exception is thrown
if this is not the case.

For convenience there is a shell script called `install_external.sh` which
is included in this repository. Running this script will install the executable.
Default installation path points to `/usr/local/bin/` but it is possible to pass it via
`--prefix` parameter (should be absolute path)::

  # Install jyu_lb_isothermal.exe to `/usr/local/bin/` with super user rights
  sudo ./install_external.sh

  # Install to a custom path
  ./install_external.sh --prefix=/home/username/.local/bin

One can add the custom path to the environment (`~/.bashrc`)::

  export PATH=$PATH:/home/username/.local/bin



Testing
-------

To run the full test-suite run::

    python -m unittest discover

Directory structure
-------------------

- jyulb -- to hold the wrapper implementations
- jyulb/fileio -- to hold the file-IO based wrapper implementations
- jyulb/internal -- to hold the internal wrapper implementations
