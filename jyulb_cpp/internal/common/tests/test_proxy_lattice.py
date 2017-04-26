"""Testing module for jyu-lb proxy-lattice data class."""
import unittest
import numpy as np

from simphony.core.cuba import CUBA
from simphony.cuds.lattice import LatticeNode
from simphony.core.data_container import DataContainer
from jyulb.internal.common.proxy_lattice import ProxyLattice
from jyulb.internal.common.domain import PyLattice, PyGeometry
from jyulb.internal.common.domain import PyIsothermalData
from simphony.testing.abc_check_lattice import (
    CheckLatticeContainer, CheckLatticeNodeCoordinates)
from jyulb.testing.jyulb_check_proxy_lattice import (
    ProxyLatticeNodeOperations)


def _create_zeroed_lattice(name, primitive_cell, size, origin):
    """ Returns a lattice where the node-data only contains null values

    """
    lat = PyLattice(np.array(size, dtype=np.uint32),
                    np.array(origin, dtype=np.float64))
    geom = PyGeometry(lat)

    for i in xrange(size[0]):
        for j in xrange(size[1]):
            for k in xrange(size[2]):
                geom.set_material_ijk(
                    np.array((i, j, k), dtype=np.uint32),
                    ProxyLattice.FLUID_ENUM)

    ivel = np.array((0.0, 0.0, 0.0), dtype=np.float64)
    ifrc = np.array((0.0, 0.0, 0.0), dtype=np.float64)
    fdata = PyIsothermalData(lat, 0.0, ivel, ifrc)

    return ProxyLattice(name, primitive_cell, geom, fdata)


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
        supported CUBA keys. CUBA.MATERIAL_ID is set to correspond to
        fluid lattice node value.

        """
        return DataContainer(MATERIAL_ID=ProxyLattice.FLUID_ENUM,
                             DENSITY=0.0, FORCE=[0, 0, 0],
                             VELOCITY=[0, 0, 0])

    def test_update_nodes(self):
        """ This test is overridden because ProxyLattice only supports
        updating lattice nodes with CUBA.MATERIAL_ID defined as
        ProxyLattice.FluidEnum
        """
        container = self.container

        indices = ((2, 3, 4), (1, 2, 3))
        nodes = [container.get_node(index) for index in indices]
        for node in nodes:
            node.data = self._create_data_with_zero_values()
        container.update_nodes(nodes)

        for n in xrange(len(indices)):
            index = indices[n]
            new_node = container.get_node(index)
            self.assertEqual(new_node, nodes[n])
            # Check that `new_node` is not the same instance as `node`
            self.assertIsNot(new_node, nodes[n])

    def test_update_nodes_with_extra_keywords(self):
        """ This test is overridden because ProxyLattice only supports
        updating lattice nodes with CUBA.MATERIAL_ID defined as
        ProxyLattice.FluidEnum
        """
        container = self.container

        indices = ((2, 3, 4), (1, 2, 3))
        nodes = [container.get_node(index) for index in indices]
        # Update with DataContainer with non supported keys.
        for node in nodes:
            node.data = self._create_data_with_zero_values()
            node.data = DataContainer(node.data, RADIUS=1)
        container.update_nodes(nodes)

        for n in xrange(len(indices)):
            index = indices[n]
            new_node = container.get_node(index)
            # We expect only the supported CUBA to be stored.
            expected = LatticeNode(
                index=nodes[n].index,
                data=self._create_data_with_zero_values())
            self.assertEqual(new_node, expected)
            # Check that `new_node` is not the same instance as `node`
            self.assertIsNot(new_node, nodes[n])


class TestProxyLatticeNodeCoordinates(CheckLatticeNodeCoordinates,
                                      unittest.TestCase):

    def container_factory(self, name, primitive_cell, size, origin):
        return _create_zeroed_lattice(name, primitive_cell, size, origin)

    def supported_cuba(self):
        return [CUBA.MATERIAL_ID, CUBA.DENSITY, CUBA.VELOCITY, CUBA.FORCE]


if __name__ == '__main__':
    unittest.main()
