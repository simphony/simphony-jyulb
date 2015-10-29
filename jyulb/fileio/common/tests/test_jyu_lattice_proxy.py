"""Testing module for jyu-lb proxy-lattice data class."""
import unittest
import numpy as np
import numpy.testing as np_test
from simphony.core.cuba import CUBA
from simphony.core.data_container import DataContainer
from simphony.cuds.lattice import LatticeNode
from simphony.cuds.abc_lattice import ABCLattice
from jyulb.fileio.common.jyu_lattice_proxy import JYULatticeProxy

from simphony.testing.abc_check_lattice import (
    CheckLatticeContainer, CheckLatticeNodeOperations,
    CheckLatticeNodeCoordinates)


class TestJYULatticeProxyContainer(CheckLatticeContainer,
                                   unittest.TestCase):

    """Test case for JYULatticeProxy class."""
    def container_factory(self, name, primitive_cell, size, origin):
        self.ext_ndata = DataContainer()
        self.ext_ndata[CUBA.MATERIAL_ID] = np.zeros(size[::-1], dtype=np.uint8)
        self.ext_ndata[CUBA.DENSITY] = np.zeros(size[::-1], dtype=np.float64)
        self.ext_ndata[CUBA.VELOCITY] = np.zeros(size[::-1] + (3,), dtype=np.float64)
        self.ext_ndata[CUBA.FORCE] = np.zeros(size[::-1] + (3,), dtype=np.float64)
        return JYULatticeProxy(name, primitive_cell, size, origin,
                               self.ext_ndata)

    def supported_cuba(self):
        return set(CUBA)


class TestJYULatticeProxyNodeOperations(CheckLatticeNodeOperations,
                                        unittest.TestCase):

    def container_factory(self, name, primitive_cell, size, origin):
        self.ext_ndata = DataContainer()
        self.ext_ndata[CUBA.MATERIAL_ID] = np.zeros(size[::-1], dtype=np.uint8)
        self.ext_ndata[CUBA.DENSITY] = np.zeros(size[::-1], dtype=np.float64)
        self.ext_ndata[CUBA.VELOCITY] = np.zeros(size[::-1] + (3,), dtype=np.float64)
        self.ext_ndata[CUBA.FORCE] = np.zeros(size[::-1] + (3,), dtype=np.float64)
        return JYULatticeProxy(name, primitive_cell, size, origin,
                               self.ext_ndata)

    def supported_cuba(self):
        return set([CUBA.MATERIAL_ID, CUBA.DENSITY, CUBA.VELOCITY,
                    CUBA.FORCE])



if __name__ == '__main__':
    unittest.main()
