from setuptools import setup, find_packages

setup(
    name='jyu_engine',
    version='0.1.1.dev0',
    author='SimPhoNy FP7 European Project',
    description='Implementation of JYU-LB wrappers',
    packages=find_packages(),
    install_requires=['simphony'],
    entry_points={
        'simphony.engine': ['jyulb = jyulb.fileio.isothermal3D']
    })
