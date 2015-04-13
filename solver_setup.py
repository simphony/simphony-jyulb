from distutils.core import setup, Extension
from Cython.Build import cythonize
import numpy

setup(ext_modules = cythonize(Extension(
           "solver",                                           # the extesion name
           sources=["jyulb/internal/isothermal/solver.pyx",    # the Cython and C++ source files
                    "JYU-LB/include/common/node.cpp",
                    "JYU-LB/include/common/filter.cpp",
                    "JYU-LB/include/collision/collision.cpp",
                    "JYU-LB/include/kernel/kernel.cpp",
                    "JYU-LB/include/solver/solver.cpp"],
           include_dirs=['JYU-LB/include/common/',             # include paths
                         'JYU-LB/include/dvs/',     
                         'JYU-LB/include/collision/',     
                         'JYU-LB/include/kernel/',     
                         'JYU-LB/include/solver/',     
                         numpy.get_include()],                 
           extra_compile_args=['-fopenmp','-O3'],              # compiler flags
           extra_link_args=['-fopenmp','-O3'],                 # linker flags
           language="c++",                                     # generate and compile C++ code
)))
