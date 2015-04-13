from setuptools import setup, find_packages
import os

os.system("python domain_setup.py build_ext --build-lib=jyulb/internal/common/")
os.system("python solver_setup.py build_ext --build-lib=jyulb/internal/isothermal/")

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
