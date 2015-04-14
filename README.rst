simphony-jyulb
==============

Internal and file-IO based JYU-LB engine-wrappers for the SimPhoNy framework.

.. image:: https://travis-ci.org/simphony/simphony-jyulb.svg?branch=master
    :target: https://travis-ci.org/simphony/simphony-jyulb
      :alt: Build status

.. image:: https://coveralls.io/repos/simphony/simphony-jyulb/badge.svg
   :target: https://coveralls.io/r/simphony/simphony-jyulb
      :alt: Test coverage

Repository
----------

simphony-jyulb is hosted on github: https://github.com/simphony/simphony-jyulb

Requirements
------------
 - `simphony-common`_ >= 0.1.1 

.. _simphony-common: https://github.com/simphony/simphony-common

.. note::
  simphony-jyulb uses keywords (COLLISION_OPERATOR, REFERENCE_DENSITY, GRAVITY,
  FLOW_TYPE, EXTERNAL_FORCING) which currently do not belong to the CUBA-keywords.

  simphony-jyulb treats the BC data component as a dummy variable (currently there
  is no SimPhoNy definition for BC).  

Installation
------------

The package requires python 2.7.x. Installation is based on setuptools::

    # build and install
    python setup.py install

or::

    # build for in-place development
    python setup.py develop

JYU-LB installation
~~~~~~~~~~~~~~~~~~~

These engine-wrappers use JYU-LB fluid flow simulators and the related source code.

The JYU-LB source code is available on github: https://github.com/simphony/JYU-LB

The file-IO based engine wrapper assumes that there is an executable called
"jyu_lb_isothermal.exe" that can be found in the PATH and an exception is thrown
if this is not the case.  

Testing
-------

To run the full test-suite run::

    python -m unittest discover

Directory structure
-------------------

- jyulb -- to hold the wrapper implementations
- jyulb/fileio -- to hold the file-IO based wrapper implementations
- jyulb/internal -- to hold the internal wrapper implementations
