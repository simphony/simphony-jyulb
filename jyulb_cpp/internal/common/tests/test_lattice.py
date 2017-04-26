"""
    Testing module for jyu-lb proxy-lattice data class.
"""
import unittest
import numpy as np
import numpy.testing as np_test

from simphony.core.cuba import CUBA
from simphony.cuds.lattice import LatticeNode
from simphony.cuds.abstractlattice import ABCLattice
from simphony.testing.abc_check_lattice import ABCCheckLattice
from simphony.core.data_container import DataContainer
from jyulb.internal.common.proxy_lattice import ProxyLattice
from jyulb.internal.common.domain import PyLattice, PyGeometry
from jyulb.internal.common.domain import PyIsothermalData


class TestProxyLattice(ABCCheckLattice, unittest.TestCase):

#    def setUp(self):
#        self.hl_data = DataContainer()
#        self.hl_data[CUBA.STATUS] = 0

    def container_factory(self, name, type_, base_vect, size, origin):
        size_n = len(size)
        orig_n = len(origin)
#        np_size = np.array(size, dtype=np.uint32)
#        np_orig = np.array(origin, dtype=np.float64)
        np_size = np.ones(3, dtype=np.uint32)
        np_orig = np.ones(3, dtype=np.float64)

        np_size[0:len(size)] = size[:]
        np_orig[0:len(origin)] = origin[:]
        
#        if size_n > 0:
#            np_size[0] = size[0]
#            if size_n > 1:
#                np_size[1] = size[1]
#                if size_n > 2:
#                    np_size[2] = size[2]
        
#        if orig_n > 0:
#            np_orig[0] = origin[0]
#            if orig_n > 1:
#                np_orig[1] = origin[1]
#                if orig_n > 2:
#                    np_orig[2] = origin[2]

        self.lat = PyLattice(np_size, np_orig)
        self.geom = PyGeometry(self.lat)

        for index in np.ndindex(tuple(np_size)):
            if index[0] > 0 and index[0] < np_size[0]-1:
                ijk = np.array(index, dtype=np.uint32)
                self.geom.set_material_ijk(ijk, ProxyLattice.FLUID_ENUM)

        ivel = np.array((0.0, 0.0, 0.0), dtype=np.float64)
        ifrc = np.array((0.0, 0.0, 0.0), dtype=np.float64)
        self.fdata = PyIsothermalData(self.lat, 0.0, ivel, ifrc)
        
        return ProxyLattice(name, type_, base_vect, self.geom, self.fdata)

    def supported_cuba(self):
        return set((CUBA.MATERIAL_ID,CUBA.DENSITY,CUBA.VELOCITY,CUBA.FORCE))

    def test_construct_lattice(self):
        """Construction of a lattice."""
        proxy = self.container
        
#        proxy = self.container_factory('test_lat', 'Cubic',
#                                         (0.1, 0.1, 0.1), (10, 20, 30),
#                                         (0, 0, 0))
#        proxy.data = self.hl_data

        self.assertIsInstance(proxy, ABCLattice, "Error: not a ABCLattice!")

        self.assertEqual(proxy._geom, self.geom)
        self.assertEqual(proxy._fdata, self.fdata)
#        self.assertEqual(proxy.data, self.hl_data)

        try:
            proxy.get_node((-1, 0, 0))
        except IndexError:
            pass
        else:
            raise AssertionError('Negative indices incorrectly accepted.')

        test_node = LatticeNode((0, -1, 0))
        try:
            proxy.update_node(test_node)
        except IndexError:
            pass
        else:
            raise AssertionError('Negative indices incorrectly accepted.')

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

        for node in proxy.iter_nodes(np.ndindex(1, self.ny, self.nz)):
            try:
                node.data[CUBA.DENSITY]
                node.data[CUBA.VELOCITY]
                node.data[CUBA.FORCE]
            except:
                pass
            else:
                raise AssertionError('Field data stored for solid nodes.')

        np_test.assert_array_equal((1*self.h, 2*self.h, 3*self.h),
                                   proxy.get_coordinate((1, 2, 3)))


if __name__ == '__main__':
    unittest.main()
