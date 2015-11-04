"""Testing module for jyu-lb proxy-lattice data class."""
import unittest
import numpy as np

from simphony.core.cuba import CUBA
from simphony.core.keywords import KEYWORDS
from simphony.core.data_container import DataContainer
from simphony.cuds.lattice import LatticeNode
from simphony.testing.utils import (create_data_container)
from jyulb.fileio.common.jyu_lattice_proxy import JYULatticeProxy
from simphony.testing.abc_check_lattice import (
    CheckLatticeContainer, CheckLatticeNodeOperations,
    CheckLatticeNodeCoordinates)


def _create_zeroed_lattice(name, primitive_cell, size, origin):
    """ Returns a lattice where the node-data only contains null values

    """
    ext_ndata = DataContainer()
    ext_ndata[CUBA.MATERIAL_ID] = np.zeros(size[::-1], dtype=np.uint8)
    ext_ndata[CUBA.DENSITY] = np.zeros(size[::-1], dtype=np.float64)
    ext_ndata[CUBA.VELOCITY] = np.zeros(size[::-1] + (3,), dtype=np.float64)
    ext_ndata[CUBA.FORCE] = np.zeros(size[::-1] + (3,), dtype=np.float64)
    return JYULatticeProxy(name, primitive_cell, size, origin, ext_ndata)


def _create_data_with_zero_values(restricted_cuba):
    """ Return a DataContainer containing ZERO values

    Parameters
    ----------
    restricted_cuba : CUBA
        The cuba keys which should appear

    """
    data = {cuba: _zero_value(cuba) for cuba in restricted_cuba}
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


class TestJYULatticeProxyContainer(CheckLatticeContainer, unittest.TestCase):

    """Test case for JYULatticeProxy class."""
    def container_factory(self, name, primitive_cell, size, origin):
        return _create_zeroed_lattice(name, primitive_cell, size, origin)

    def supported_cuba(self):
        return set(CUBA)


class TestJYULatticeProxyNodeOperations(CheckLatticeNodeOperations,
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


class TestJYULatticeProxyNodeCoordinates(CheckLatticeNodeCoordinates,
                                         unittest.TestCase):

    def container_factory(self, name, primitive_cell, size, origin):
        return _create_zeroed_lattice(name, primitive_cell, size, origin)

    def supported_cuba(self):
        return [CUBA.MATERIAL_ID, CUBA.DENSITY, CUBA.VELOCITY, CUBA.FORCE]

if __name__ == '__main__':
    unittest.main()
