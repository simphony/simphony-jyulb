"""Testing module for jyu-lb proxy-lattice data class."""
import unittest
import numpy as np
import numpy.testing as np_test
from simphony.core.cuba import CUBA
from simphony.core.data_container import DataContainer
from simphony.cuds.lattice import LatticeNode
from simphony.cuds.abstractlattice import ABCLattice
from jyulb.fileio.common.jyu_lattice_proxy import JYULatticeProxy


class JYULatticeProxyTestCase(unittest.TestCase):

    """Test case for JYULatticeProxy class."""

    def setUp(self):
        h = 0.1
        nx = 10
        ny = 20
        nz = 30
        ox = 0
        oy = 0
        oz = 0

        self.h = h
        self.nx = nx
        self.ny = ny
        self.nz = nz
        self.ox = ox
        self.oy = oy
        self.oz = oz
        self.name = 'Test'
        self.type = 'Cubic'
        self.bvec = (h, h, h)
        self.size = (nx, ny, nz)
        self.orig = (ox, oy, oz)

        self.data = DataContainer()
        self.data[CUBA.STATUS] = 0

        geom_arr = np.zeros((nz, ny, nx), dtype=np.uint8)
        den_arr = np.zeros((nz, ny, nx), dtype=np.float64)
        vel_arr = np.zeros((nz, ny, nx, 3), dtype=np.float64)
        frc_arr = np.zeros((nz, ny, nx, 3), dtype=np.float64)

        self.ext_ndata = {}
        self.ext_ndata[CUBA.MATERIAL_ID] = geom_arr
        self.ext_ndata[CUBA.DENSITY] = den_arr
        self.ext_ndata[CUBA.VELOCITY] = vel_arr
        self.ext_ndata[CUBA.FORCE] = frc_arr

    def test_construct_lattice(self):
        """Construction of a lattice."""
        lat = JYULatticeProxy(self.name, self.type, self.bvec, self.size,
                              self.orig, self.ext_ndata)
        lat.data = self.data

        self.assertIsInstance(lat, ABCLattice, "Error: not a ABCLattice!")

        self.assertEqual(lat.name, self.name)
        self.assertEqual(lat.type, self.type)
        self.assertEqual(lat._external_node_data, self.ext_ndata)
        self.assertEqual(lat.data, self.data)
        np_test.assert_array_equal(lat.size, self.size)
        np_test.assert_array_equal(lat.origin, self.orig)
        np_test.assert_array_equal(lat.base_vect, self.bvec)

        try:
            lat.get_node((-1, 0, 0))
        except IndexError:
            pass
        else:
            raise AssertionError('Negative indices incorrectly accepted.')

        test_node = LatticeNode((0, -1, 0))
        try:
            lat.update_node(test_node)
        except IndexError:
            pass
        else:
            raise AssertionError('Negative indices incorrectly accepted.')

    def test_set_get_iter_lattice_nodes(self):
        """Creation of lattices using the factory functions."""
        lat = JYULatticeProxy(self.name, self.type, self.bvec, self.size,
                              self.orig, self.ext_ndata)

        SOLID = 1
        FLUID = 255

        geom2 = np.zeros((self.nz, self.ny, self.nx), dtype=np.uint8)
        geom2[:, :, 0] = SOLID
        geom2[:, :, self.nx-1] = SOLID
        geom2[:, :, 1:self.nx-1] = FLUID

        for node in lat.iter_nodes():
            if node.index[0] == 0 or node.index[0] == self.nx-1:
                node.data[CUBA.MATERIAL_ID] = SOLID
            else:
                node.data[CUBA.MATERIAL_ID] = FLUID
            lat.update_node(node)

        np_test.assert_array_equal(self.ext_ndata[CUBA.MATERIAL_ID], geom2)

        for node in lat.iter_nodes():
            if node.data[CUBA.MATERIAL_ID] == FLUID:
                node.data[CUBA.VELOCITY] = lat.get_coordinate(node.index)
                node.data[CUBA.DENSITY] = 1.0
            lat.update_node(node)

        self.assertEqual(self.ext_ndata[CUBA.DENSITY].sum(),
                         self.nz*self.ny*(self.nx-2))

        for node in lat.iter_nodes(np.ndindex(7, 6, 5)):
            if node.data[CUBA.MATERIAL_ID] == FLUID:
                np_test.assert_array_equal(node.data[CUBA.VELOCITY],
                                           lat.get_coordinate(node.index))

        node = lat.get_node((self.nx/2, self.ny/2, self.nz/2))
        node.data[CUBA.FORCE] = (-3, -2, -1)
        lat.update_node(node)

        self.assertEqual(self.ext_ndata[CUBA.FORCE][:, :, :, 0].sum(), -3)
        self.assertEqual(self.ext_ndata[CUBA.FORCE][:, :, :, 1].sum(), -2)
        self.assertEqual(self.ext_ndata[CUBA.FORCE][:, :, :, 2].sum(), -1)

if __name__ == '__main__':
    unittest.main()
