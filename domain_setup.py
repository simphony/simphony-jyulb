from distutils.core import setup, Extension
from Cython.Build import cythonize
import numpy

setup(ext_modules=cythonize(Extension(
    "domain",                                     # extesion name
    sources=["jyulb/internal/common/domain.pyx",  # Cython/C++ source files
             "JYU-LB/include/common/node.cpp"],
    include_dirs=['JYU-LB/include/common/',       # include paths
                  numpy.get_include()],
    extra_compile_args=['-fopenmp', '-O3'],       # compiler flags
    extra_link_args=['-fopenmp', '-O3'],          # linker flags
    language="c++",                               # generate/compile C++ code
)))
