import numpy as np
from simphony.cuds.abc_lattice import ABCLattice
from simphony.core.cuds_item import CUDSItem
from simphony.core.data_container import DataContainer
from simphony.cuds.primitive_cell import PrimitiveCell
from simphony.cuds.lattice import LatticeNode



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
        name of lattice
    primitive_cell : PrimitiveCell
        primitive cell specifying the 3D Bravais lattice
    size : int[3]
        lattice dimensions
    origin : float[3]
        lattice origin
    data : DataContainer
        high level CUBA data assigned to lattice
    external_node_data : dictionary
        references (value) to external data storages (multidimensional
        arrays) for each prescribed CUBA keyword (key)
    """

    # Enumeration of material IDs
    SOLID_ENUM = 0
    FLUID_ENUM = 255

    def __init__(self, name, primitive_cell, size, origin, ext_node_data):
        self.name = name
        self._primitive_cell = primitive_cell
        self._size = size[0], size[1], size[2]
        self._origin = np.array((origin[0], origin[1], origin[2]),
                                dtype=np.float)
        self._data = DataContainer()
        self._items_count = {
            CUDSItem.NODE: lambda: self._size
        }
        self._external_node_data = ext_node_data

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
        # the indices must be reversed in a data access
        rev_ind = tuple(index[::-1])

        if any(value < 0 for value in rev_ind):
            raise IndexError('invalid index: {}'.format(rev_ind))

        node = LatticeNode(tuple(index))
        for key in self._external_node_data:
            node.data[key] = self._external_node_data[key][rev_ind]
        return node

    def update_nodes(self, nodes):
        """Update the corresponding lattice nodes (data copied).

        Parameters
        ----------
        nodes : iterable of LatticeNode objects
            reference to LatticeNode objects from where the data is copied
            to the JYULatticeProxy

        Raises
        ------
        IndexError
           if the index of the given node includes negative components.
        """
        for node in nodes:
            # JYU-LB modeling engines assume a specific memory order;
            # the indexes must be reversed in data access
            ind = node.index
            rev_ind = ind[::-1]
            if any(value < 0 for value in rev_ind):
                raise IndexError('invalid index: {}'.format(rev_ind))

            for key in self._external_node_data:
                if key in node.data:
                    self._external_node_data[key][rev_ind] = node.data[key]

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

    def count_of(self, item_type):
        """ Return the count of item_type in the container.

        Parameters
        ----------
        item_type : CUDSItem
            The CUDSItem enum of the type of the items to return
            the count of.

        Returns
        -------
        count : int
            The number of items of item_type in the container.

        Raises
        ------
        ValueError :
            If the type of the item is not supported in the current
            container.

        """
        try:
            return np.prod(self._items_count[item_type]())
        except KeyError:
            error_str = "Trying to obtain count a of non-supported item: {}"
            raise ValueError(error_str.format(item_type))
