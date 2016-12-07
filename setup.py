import os
from setuptools import setup, find_packages, Extension
from Cython.Build import cythonize
import numpy

VERSION = '0.2.1'


def write_version_py(filename=None):
    if filename is None:
        filename = os.path.join(
            os.path.dirname(__file__), 'jyulb', 'version.py')
    ver = """\
version = '%s'
"""
    fh = open(filename, 'wb')
    try:
        fh.write(ver % VERSION)
    finally:
        fh.close()


write_version_py()

extensions = [
    Extension(
        name="jyulb.internal.isothermal.solver",
        sources=["jyulb/internal/isothermal/solver.pyx",
                 "JYU-LB/include/common/node.cpp",
                 "JYU-LB/include/common/filter.cpp",
                 "JYU-LB/include/collision/collision.cpp",
                 "JYU-LB/include/kernel/kernel.cpp",
                 "JYU-LB/include/solver/solver.cpp"],
        include_dirs=['JYU-LB/include/common/',
                      'JYU-LB/include/dvs/',
                      'JYU-LB/include/collision/',
                      'JYU-LB/include/kernel/',
                      'JYU-LB/include/solver/',
                      numpy.get_include()],
        extra_compile_args=['-fopenmp', '-O3'],
        extra_link_args=['-fopenmp', '-O3'],
        language="c++",
    ),
    Extension(
        name="jyulb.internal.common.domain",
        sources=["jyulb/internal/common/domain.pyx",
                 "JYU-LB/include/common/node.cpp"],
        include_dirs=['JYU-LB/include/common/',
                      numpy.get_include()],
        extra_compile_args=['-fopenmp', '-O3'],
        extra_link_args=['-fopenmp', '-O3'],
        language="c++",
    )
]

setup(
    name='jyulb_engine',
    version=VERSION,
    author='SimPhoNy FP7 European Project',
    description='Implementation of JYU-LB wrappers',
    packages=find_packages(),
    install_requires=['simphony~=0.4'],
    ext_modules=cythonize(extensions),
    entry_points={
        'simphony.engine': [
            'jyulb_fileio_isothermal =' +
            'jyulb.fileio.isothermal.jyulb_engine',
            'jyulb_internal_isothermal =' +
            'jyulb.internal.isothermal.jyulb_engine'
        ]
    },
)
