from simphony.cuds.abstractlattice import ABCLattice
from simphony.cuds.lattice import LatticeNode
import numpy as np


class JYULatticeProxy(ABCLattice):
    """
    A Bravais lattice;
    stores references to data containers (node related data).

    Parameters
    ----------
    name : str
    type : str
        Bravais lattice type (should agree with the base_vect below).
    base_vect : D x float
        defines a Bravais lattice (an alternative for primitive vectors).
    size : tuple of D x size
        number of lattice nodes (in the direction of each axis).
    origin : D x float

    """
    def __init__(self, name, type, base_vect, size, origin, data):
        self.name = name
        self._type = type
        self._base_vect = np.array(base_vect, dtype=np.float64)
        self._size = tuple(size)
        self._origin = np.array(origin, dtype=np.float64)
        self._data = data

    @property
    def type(self):
        return self._type

    @property
    def base_vect(self):
        return self._base_vect

    @property
    def size(self):
        return self._size

    @property
    def origin(self):
        return self._origin

    def get_node(self, index):
        """Get a copy of the node corresponding to the given index.

        Parameters
        ----------
        index: tuple of D x int (node index coordinate)

        Returns
        -------
        A reference to a LatticeNode object

        """
        ti = tuple(index)
        rev_ti = tuple(ti[::-1])

        node = LatticeNode(ti)
        for key in self._data:
            node.data[key] = self._data[key][rev_ti]

        return node

    def update_node(self, lat_node):
        """Update the corresponding lattice node (data copied).

        Parameters
        ----------
        lat_node: reference to a LatticeNode object
            data copied from the given node

        """
        ind = lat_node.index
        rev_ind = tuple(ind[::-1])

        for key in self._data:
            if key in lat_node.data:
                self._data[key][rev_ind] = lat_node.data[key]

    def iter_nodes(self, indices=None):
        """Get an iterator over the LatticeNodes described by the indices.

        Parameters
        ----------
        indices: iterable set of D x int, optional
            node index coordinates

        Returns
        -------
        A generator for LatticeNode objects

        """
        if indices is None:
            for index in np.ndindex(self._size):
                yield self.get_node(index)
        else:
            for index in indices:
                yield self.get_node(index)

    def get_coordinate(self, index):
        """Get coordinate of the given index coordinate.

        Parameters
        ----------
        index : D x int (node index coordinate)

        Returns
        -------
        D x float

        """
        return self.origin + self.base_vect*np.array(index)
