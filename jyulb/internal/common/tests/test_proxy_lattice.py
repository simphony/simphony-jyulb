"""Testing module for jyu-lb proxy-lattice data class."""
import unittest
import numpy as np

from simphony.core.cuba import CUBA
from simphony.core.keywords import KEYWORDS
from simphony.cuds.lattice import LatticeNode
from simphony.cuds.abc_lattice import ABCLattice
from simphony.core.data_container import DataContainer
from jyulb.internal.common.proxy_lattice import ProxyLattice
from jyulb.internal.common.domain import PyLattice, PyGeometry
from jyulb.internal.common.domain import PyIsothermalData
from simphony.testing.utils import (create_data_container)
from simphony.testing.abc_check_lattice import (
    CheckLatticeContainer, CheckLatticeNodeOperations,
    CheckLatticeNodeCoordinates)


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

def _create_data_with_zero_values(restricted_cuba):
    """ Return a DataContainer containing ZERO values

    Parameters
    ----------
    restricted_cuba : CUBA
        The cuba keys which should appear

    """
    data = {cuba: _zero_value(cuba) for cuba in restricted_cuba}
    data[CUBA.MATERIAL_ID] = ProxyLattice.FLUID_ENUM
    return DataContainer(data)


def _zero_value(cuba):
    keyword = KEYWORDS[CUBA(cuba).name]
    if np.issubdtype(keyword.dtype, str):
        return keyword.name
    else:
        shape = keyword.shape
        if shape == [1]:
            if np.issubdtype(keyword.dtype, 'float'):
                return float(0)
            if np.issubdtype(keyword.dtype, 'int'):
                return int(0)
        else:
            if np.issubdtype(keyword.dtype, 'float'):
                return np.zeros(shape=shape, dtype=np.float64)
            if np.issubdtype(keyword.dtype, 'int'):
                return np.zeros(shape=shape, dtype=np.int32)


class TestProxyLatticeContainer(CheckLatticeContainer, unittest.TestCase):

    """Test case for ProxyLattice class."""
    def container_factory(self, name, primitive_cell, size, origin):
        return _create_zeroed_lattice(name, primitive_cell, size, origin)

    def supported_cuba(self):
        return set(CUBA)


class TestProxyLatticeNodeOperations(CheckLatticeNodeOperations,
                                     unittest.TestCase):

    def container_factory(self, name, primitive_cell, size, origin):
        return _create_zeroed_lattice(name, primitive_cell, size, origin)

    def supported_cuba(self):
        return [CUBA.MATERIAL_ID, CUBA.DENSITY, CUBA.VELOCITY, CUBA.FORCE]


    def test_get_node(self):
        # TODO https://github.com/simphony/simphony-common/issues/216
        # For some tests in CheckLatticeNodeOperations it is expected
        # that the `container_factory` has produced a lattice with nodes
        # that are empty.  And that then these nodes can be updated on
        # the lattice with the 'supported_cuba'. JYULatticeProxy requires
        # that the supported cuba appear in the `external_node_data` parameter
        # (what CUBA appear on each node.data cannot be changed).  We
        # therefore override some tests to ensure that we are testing
        # the getters and iters using the node-data values provided in
        # the `container_factory`.
        container = self.container

        index = 2, 3, 4
        node = container.get_node(index)

        expected = LatticeNode(index, data=_create_data_with_zero_values(
            self.supported_cuba()))
        self.assertEqual(node, expected)

        # check that mutating the node does not change internal info
        node.data = create_data_container()
        self.assertNotEqual(container.get_node(index), node)

    def test_iter_nodes(self):
        # TODO Overriding this method (see test_get_node above).
        # https://github.com/simphony/simphony-common/issues/216
        container = self.container

        # number of nodes
        number_of_nodes = sum(1 for node in container.iter_nodes())
        self.assertEqual(number_of_nodes, np.prod(self.size))

        # data
        for node in container.iter_nodes():
            self.assertEqual(node.data, _create_data_with_zero_values(
                self.supported_cuba()))

        # indexes
        x, y, z = np.meshgrid(
            range(self.size[0]), range(self.size[1]), range(self.size[2]))
        expected = set(zip(x.flat, y.flat, z.flat))
        indexes = {node.index for node in container.iter_nodes()}
        self.assertEqual(indexes, expected)

    def test_iter_nodes_subset(self):
        # TODO Overriding this method (see test_get_node above).
        # https://github.com/simphony/simphony-common/issues/216
        container = self.container

        x, y, z = np.meshgrid(
            range(2, self.size[0]),
            range(self.size[1]-4),
            range(3, self.size[2], 2))
        expected = set(zip(x.flat, y.flat, z.flat))

        # data
        for node in container.iter_nodes(expected):
            self.assertEqual(node.data, _create_data_with_zero_values(
                self.supported_cuba()))

        # indexes
        indexes = {node.index for node in container.iter_nodes(expected)}
        self.assertEqual(indexes, expected)


    def test_update_nodes(self):
        """ ProxyLattice only supports updating lattice nodes with
        CUBA.MATERIAL_ID defined as ProxyLattice.FluidEnum
        """
        container = self.container

        indices = ((2, 3, 4), (1, 2, 3))
        nodes = [container.get_node(index) for index in indices]
        for node in nodes:
            node.data = _create_data_with_zero_values(self.supported_cuba())
        container.update_nodes(nodes)

        for n in xrange(len(indices)):
            index = indices[n]
            new_node = container.get_node(index)
            self.assertEqual(new_node, nodes[n])
            # Check that `new_node` is not the same instance as `node`
            self.assertIsNot(new_node, nodes[n])

    def test_update_nodes_with_extra_keywords(self):
        """ ProxyLattice only supports updating lattice nodes with
        CUBA.MATERIAL_ID defined as ProxyLattice.FluidEnum
        """
        container = self.container

        indices = ((2, 3, 4), (1, 2, 3))
        nodes = [container.get_node(index) for index in indices]
        # Update with full DataContainer.
        for node in nodes:
            node.data = _create_data_with_zero_values(set(CUBA))
        container.update_nodes(nodes)

        for n in xrange(len(indices)):
            index = indices[n]
            new_node = container.get_node(index)
            # We expect only the supported CUBA to be stored.
            expected = LatticeNode(
                index=nodes[n].index,
                data=_create_data_with_zero_values(self.supported_cuba()))
            self.assertEqual(new_node, expected)
            # Check that `new_node` is not the same instance as `node`
            self.assertIsNot(new_node, nodes[n])



if __name__ == '__main__':
    unittest.main()
