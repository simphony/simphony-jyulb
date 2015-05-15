from setuptools import find_packages
from distutils.core import setup, Extension
from Cython.Build import cythonize
import numpy

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
    name='jyu_engine',
    version='0.1.3.dev0',
    author='SimPhoNy FP7 European Project',
    description='Implementation of JYU-LB wrappers',
    packages=find_packages(),
    install_requires=['simphony'],
    ext_modules=cythonize(extensions),
    entry_points={
        'simphony.engine': [
            'jyulb_fileio_isothermal = jyulb.fileio.isothermal.jyu_engine',
            'jyulb_internal_isothermal = jyulb.internal.isothermal.jyu_engine',
        ]
    }
)
