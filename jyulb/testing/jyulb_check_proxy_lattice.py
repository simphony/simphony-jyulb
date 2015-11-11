import numpy as np
from simphony.cuds.lattice import LatticeNode
from simphony.testing.utils import (create_data_container)
from simphony.testing.abc_check_lattice import CheckLatticeNodeOperations


class ProxyLatticeNodeOperations(CheckLatticeNodeOperations):
    """
    TODO https://github.com/simphony/simphony-common/issues/216

    For some tests in CheckLatticeNodeOperations it is expected
    that the `container_factory` has produced a lattice with nodes
    that are empty.  And that then these nodes can be updated on
    the lattice with the 'supported_cuba'. ProxyLattice and
    JYULatticeProxy require that the supported cuba appear in their
    internal variables and always return LatticeNode object with these
    CUBA keys defined. What CUBA appear on each node.data cannot be
    changed. We therefore override some tests to ensure that we are testing
    the getters and iters using the node-data values provided in
    the `container_factory`.

    """
    def test_get_node(self):
        # TODO Overriding this method (see above).
        # https://github.com/simphony/simphony-common/issues/216
        container = self.container

        index = 2, 3, 4
        node = container.get_node(index)

        expected = LatticeNode(index,
                               data=self._create_data_with_zero_values())

        self.assertEqual(node, expected)

        # check that mutating the node does not change internal info
        node.data = create_data_container()
        self.assertNotEqual(container.get_node(index), node)

    def test_iter_nodes(self):
        # TODO Overriding this method (see above).
        # https://github.com/simphony/simphony-common/issues/216
        container = self.container

        # number of nodes
        number_of_nodes = sum(1 for node in container.iter_nodes())
        self.assertEqual(number_of_nodes, np.prod(self.size))

        # data
        for node in container.iter_nodes():
            self.assertEqual(node.data,
                             self._create_data_with_zero_values())

        # indexes
        x, y, z = np.meshgrid(
            range(self.size[0]), range(self.size[1]), range(self.size[2]))
        expected = set(zip(x.flat, y.flat, z.flat))
        indexes = {node.index for node in container.iter_nodes()}
        self.assertEqual(indexes, expected)

    def test_iter_nodes_subset(self):
        # TODO Overriding this method (see above).
        # https://github.com/simphony/simphony-common/issues/216
        container = self.container

        x, y, z = np.meshgrid(
            range(2, self.size[0]),
            range(self.size[1]-4),
            range(3, self.size[2], 2))
        expected = set(zip(x.flat, y.flat, z.flat))

        # data
        for node in container.iter_nodes(expected):
            self.assertEqual(node.data,
                             self._create_data_with_zero_values())

        # indexes
        indexes = {node.index for node in container.iter_nodes(expected)}
        self.assertEqual(indexes, expected)
