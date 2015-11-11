"""Testing module for jyu-lb proxy-lattice data class."""
import unittest
import numpy as np

from simphony.core.cuba import CUBA
from simphony.core.data_container import DataContainer
from jyulb.fileio.common.proxy_lattice import ProxyLattice
from simphony.testing.abc_check_lattice import (
    CheckLatticeContainer, CheckLatticeNodeCoordinates)
from jyulb.testing.jyulb_check_proxy_lattice import (
    ProxyLatticeNodeOperations)


def _create_zeroed_lattice(name, primitive_cell, size, origin):
    """ Returns a lattice where the node-data only contains null values

    """
    ext_ndata = DataContainer()
    ext_ndata[CUBA.MATERIAL_ID] = np.zeros(size[::-1], dtype=np.uint8)
    ext_ndata[CUBA.DENSITY] = np.zeros(size[::-1], dtype=np.float64)
    ext_ndata[CUBA.VELOCITY] = np.zeros(size[::-1] + (3,), dtype=np.float64)
    ext_ndata[CUBA.FORCE] = np.zeros(size[::-1] + (3,), dtype=np.float64)
    return ProxyLattice(name, primitive_cell, size, origin, ext_ndata)


class TestProxyLatticeContainer(CheckLatticeContainer, unittest.TestCase):

    """Test case for ProxyLattice class."""
    def container_factory(self, name, primitive_cell, size, origin):
        return _create_zeroed_lattice(name, primitive_cell, size, origin)

    def supported_cuba(self):
        return set(CUBA)


class TestProxyLatticeNodeOperations(ProxyLatticeNodeOperations,
                                     unittest.TestCase):

    def container_factory(self, name, primitive_cell, size, origin):
        return _create_zeroed_lattice(name, primitive_cell, size, origin)

    def supported_cuba(self):
        return [CUBA.MATERIAL_ID, CUBA.DENSITY, CUBA.VELOCITY, CUBA.FORCE]

    def _create_data_with_zero_values(self):
        """ Return a DataContainer containing ZERO values for the
        supported CUBA keys

        """
        return DataContainer(MATERIAL_ID=0, DENSITY=0.0,
                             VELOCITY=[0, 0, 0], FORCE=[0, 0, 0])


class TestProxyLatticeNodeCoordinates(CheckLatticeNodeCoordinates,
                                      unittest.TestCase):

    def container_factory(self, name, primitive_cell, size, origin):
        return _create_zeroed_lattice(name, primitive_cell, size, origin)

    def supported_cuba(self):
        return [CUBA.MATERIAL_ID, CUBA.DENSITY, CUBA.VELOCITY, CUBA.FORCE]

if __name__ == '__main__':
    unittest.main()
