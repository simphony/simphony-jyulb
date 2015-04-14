from simphony.core.data_container import DataContainer
from simphony.cuds.abstractlattice import ABCLattice
from simphony.cuds.lattice import LatticeNode
import numpy as np


class JYULatticeProxy(ABCLattice):

    """
    A proxy lattice for accessing state data of a JYU-LB modeling engine.

    Updates and queries of node data are relayed to external data storages.
    Acknowledges only those CUBA keywords which are prescribed at the
    initialization.

    Enumeration of MATERIAL ID values for SOLID and FLUID lattice nodes.

    Attributes
    ----------
    name : str
    type : str
        Bravais lattice type (should agree with the base_vect below).
    base_vect : D x float
        defines a Bravais lattice (an alternative for primitive vectors).
    size : tuple of D x size
        number of lattice nodes (in the direction of each axis).
    origin : D x float
    data : DataContainer
        high level CUBA data assigned to lattice
    external_node_data : dictionary
        references (value) to external data storages (multidimensional
        arrays) for each prescribed CUBA keyword (key)
    """

    # Enumeration of material IDs
    SOLID_ENUM = 0
    FLUID_ENUM = 255

    def __init__(self, name, type, base_vect, size, origin, ext_node_data):
        self.name = name
        self._type = type
        self._base_vect = np.array(base_vect, dtype=np.float64)
        self._size = tuple(size)
        self._origin = np.array(origin, dtype=np.float64)
        self._data = DataContainer()
        self._external_node_data = ext_node_data

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

    @property
    def data(self):
        return DataContainer(self._data)

    @data.setter
    def data(self, value):
        self._data = DataContainer(value)

    def get_node(self, index):
        """Get a copy of the node corresponding to the given index.

        Parameters
        ----------
        index : tuple of D x int (node index coordinate)

        Returns
        -------
        A reference to a LatticeNode object

        Raises
        ------
        IndexError
           if the given index includes negative components.
        """
        # JYU-LB modeling engines assume a specific memory order;
        # the indexes must be reversed in a data access
        ti = index
        rev_ti = ti[::-1]

        if any(value < 0 for value in rev_ti):
            raise IndexError('invalid index: {}'.format(rev_ti))

        node = LatticeNode(ti)
        for key in self._external_node_data:
            node.data[key] = self._external_node_data[key][rev_ti]

        return node

    def update_node(self, lat_node):
        """Update the corresponding lattice node (data copied).

        Parameters
        ----------
        lat_node : reference to a LatticeNode object
            data copied from the given node

        Raises
        ------
        IndexError
           if the index of the given node includes negative components.
        """
        # JYU-LB modeling engines assume a specific memory order;
        # the indexes must be reversed in a data access
        ind = lat_node.index
        rev_ind = ind[::-1]

        if any(value < 0 for value in rev_ind):
            raise IndexError('invalid index: {}'.format(rev_ind))

        for key in self._external_node_data:
            if key in lat_node.data:
                self._external_node_data[key][rev_ind] = lat_node.data[key]

    def iter_nodes(self, indices=None):
        """Get an iterator over the LatticeNodes described by the indices.

        Parameters
        ----------
        indices : iterable set of D x int, optional
            node index coordinates

        Yields
        -------
        LatticeNode
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

        Raises
        ------
        NotImplementedError
           if the lattice type is 'Hexagonal'.
        """
        if self._type == 'Hexagonal':
            raise NotImplementedError("""Get_coordinate for
                Hexagonal system not implemented!""")

        return self.origin + self.base_vect*np.array(index)
