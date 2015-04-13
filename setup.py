from setuptools import setup, find_packages
import os

cmmnd1 = "python domain_setup.py build_ext " \
         "--build-lib=jyulb/internal/common/"

cmmnd2 = "python solver_setup.py build_ext " \
         "--build-lib=jyulb/internal/isothermal/"

os.system(cmmnd1)
os.system(cmmnd2)

setup(
    name='jyu_engine',
    version='0.1.2.dev0',
    author='SimPhoNy FP7 European Project',
    description='Implementation of JYU-LB wrappers',
    packages=find_packages(),
    install_requires=['simphony'],
    entry_points={
        'simphony.engine': [
            'jyulb_fileio_isothermal = jyulb.fileio.isothermal',
            'jyulb_internal_isothermal = jyulb.internal.isothermal',
        ]
    })
