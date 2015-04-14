"""Testing module for jyu-lb proxy-lattice data class."""
import unittest
import numpy as np
import numpy.testing as np_test
from simphony.core.cuba import CUBA
from simphony.cuds.lattice import LatticeNode
from simphony.cuds.abstractlattice import ABCLattice
from simphony.core.data_container import DataContainer
from jyulb.internal.common.proxy_lattice import ProxyLattice
from jyulb.internal.common.domain import PyLattice, PyGeometry
from jyulb.internal.common.domain import PyIsothermalData


class ProxyLatticeTestCase(unittest.TestCase):

    """Test case for ProxyLattice class."""

    def setUp(self):
        self.h = 0.1
        self.nx = 10
        self.ny = 20
        self.nz = 30
        self.name = 'Test'
        self.type = 'Cubic'
        self.bvec = (self.h, self.h, self.h)
        self.size = np.array((self.nx, self.ny, self.nz), dtype=np.uint32)
        self.orig = np.array((0, 0, 0), dtype=np.float64)

        self.hl_data = DataContainer()
        self.hl_data[CUBA.STATUS] = 0

        self.lat = PyLattice(self.size, self.orig)
        self.geom = PyGeometry(self.lat)

        for index in np.ndindex(tuple(self.size)):
            if index[0] > 0 and index[0] < self.nx-1:
                ijk = np.array(index, dtype=np.uint32)
                self.geom.set_material_ijk(ijk, ProxyLattice.FLUID_ENUM)

        ivel = np.array((0.0, 0.0, 0.0), dtype=np.float64)
        ifrc = np.array((0.0, 0.0, 0.0), dtype=np.float64)
        self.fdata = PyIsothermalData(self.lat, 0.0, ivel, ifrc)

    def test_construct_lattice(self):
        """Construction of a lattice."""
        proxy = ProxyLattice(self.name, self.type,
                             self.bvec, self.geom, self.fdata)

        proxy.data = self.hl_data

        self.assertIsInstance(proxy, ABCLattice, "Error: not a ABCLattice!")

        self.assertEqual(proxy.name, self.name)
        self.assertEqual(proxy.type, self.type)
        self.assertEqual(proxy._geom, self.geom)
        self.assertEqual(proxy._fdata, self.fdata)
        self.assertEqual(proxy.data, self.hl_data)
        np_test.assert_array_equal(proxy.size, self.size)
        np_test.assert_array_equal(proxy.origin, self.orig)
        np_test.assert_array_equal(proxy.base_vect, self.bvec)

        with self.assertRaises(IndexError):
            proxy.get_node((-1, 0, 0))

        test_node = LatticeNode((0, -1, 0))
        with self.assertRaises(IndexError):
            proxy.update_node(test_node)

    def test_set_get_iter_lattice_nodes(self):
        """Creation of lattices using the factory functions."""
        proxy = ProxyLattice(self.name, self.type,
                             self.bvec, self.geom, self.fdata)

        iden = 1.0
        ivel = np.array((0.0, 0.0, 1.0), dtype=np.float64)
        ifrc = np.array((0.0, 1.0, 0.0), dtype=np.float64)

        for node in proxy.iter_nodes():
            node.data[CUBA.DENSITY] = iden
            node.data[CUBA.VELOCITY] = ivel
            node.data[CUBA.FORCE] = ifrc
            proxy.update_node(node)

        check_iden = 0.0
        check_ivelx = 0.0
        check_ively = 0.0
        check_ivelz = 0.0
        check_ifrcx = 0.0
        check_ifrcy = 0.0
        check_ifrcz = 0.0
        check_sum = (self.nx-2)*self.ny*self.nz

        for node in proxy.iter_nodes():
            if node.data[CUBA.MATERIAL_ID] == ProxyLattice.FLUID_ENUM:
                check_iden = check_iden + node.data[CUBA.DENSITY]
                check_ivelx = check_ivelx + node.data[CUBA.VELOCITY][0]
                check_ively = check_ively + node.data[CUBA.VELOCITY][1]
                check_ivelz = check_ivelz + node.data[CUBA.VELOCITY][2]
                check_ifrcx = check_ifrcx + node.data[CUBA.FORCE][0]
                check_ifrcy = check_ifrcy + node.data[CUBA.FORCE][1]
                check_ifrcz = check_ifrcz + node.data[CUBA.FORCE][2]

        self.assertEqual(check_sum, check_iden)
        self.assertEqual(0.0, check_ivelx)
        self.assertEqual(0.0, check_ively)
        self.assertEqual(check_sum, check_ivelz)
        self.assertEqual(0.0, check_ifrcx)
        self.assertEqual(check_sum, check_ifrcy)
        self.assertEqual(0.0, check_ifrcz)

        # Iterated nodes belong to the solid phase (channel walls)
        # and should not have data related to the given CUBA keywords
        for node in proxy.iter_nodes(np.ndindex(1, self.ny, self.nz)):
            data = node.data
            for key in (CUBA.DENSITY, CUBA.VELOCITY, CUBA.FORCE):
                self.assertNotIn(key, data)

        np_test.assert_array_equal((1*self.h, 2*self.h, 3*self.h),
                                   proxy.get_coordinate((1, 2, 3)))

if __name__ == '__main__':
    unittest.main()
