import os
from setuptools import setup, find_packages

VERSION = '0.3.0'

def write_version_py(filename=None):
    if filename is None:
        filename = os.path.join(
            os.path.dirname(__file__), 'wrapper', 'version.py')

    ver = """version = '%s'"""
    fh = open(filename, 'wb')
    try:
        fh.write(ver % VERSION)
    finally:
        fh.close()

write_version_py()

setup(
    name='simphony_jyulb',
    version=VERSION,
    author='SimPhoNy, EU FP7 Project (Nr. 604005) www.simphony-project.eu',
    description='JYU-LB wrapper for the SimPhoNy framework',
    packages=find_packages(),
    install_requires=['numpy>=1.9.1','simphony>=0.6','jyulb>=0.2.0'],
    entry_points={'simphony.engine': ['jyulb = wrapper']}
)
