import numpy as np
from simphony.core.cuba import CUBA
from simphony.core.data_container import DataContainer
from simphony.cuds.abstractlattice import ABCLattice
from simphony.cuds.lattice import LatticeNode
from jyulb.internal.common import domain


class ProxyLattice(ABCLattice):

    """
    A proxy lattice for isothermal state data of a JYU-LB modelling engine.

    Updates and queries of node data are relayed to external data storages.
    Only CUBA keywords DENSITY, VELOCITY, and FORCE are acknowledged.

    Enumeration of MATERIAL ID values for SOLID and FLUID lattice nodes.

    Attributes
    ----------
    name : str
    primitive_cell : PrimitiveCell
        primitive cell speficyinf the 3D Bravais lattice
    size : int[3]
        number of lattice nodes (in the direction of primitive vector).
    origin : float[3]
        lattice origin
    data : DataContainer
        high level CUBA data assigned to lattice.
    geom : PyGeom
        lattice structure info and MATERIAL_ID data for lattice nodes.
    fdata : PyIsothermalData
        isothermal field data for fluid lattice nodes.
    """

    # Enumeration of material IDs
    SOLID_ENUM = domain.SOLID_ENUM
    FLUID_ENUM = domain.FLUID_ENUM

    def __init__(self, name, primitive_cell, geom, fdata):
        self.name = name
        self._primitive_cell = primitive_cell
        self._data = DataContainer()
        self._geom = geom
        self._fdata = fdata

        size = np.zeros(3, dtype=np.uint32)
        self._geom.get_lattice().get_size(size)
        self._size = tuple(size)

        origin = np.zeros(3, dtype=np.float64)
        self._geom.get_lattice().get_origin(origin)
        self._origin = tuple(origin)

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
        if any(value < 0 for value in index):
            raise IndexError('invalid index: {}'.format(index))

        node = LatticeNode(index)

        ijk = np.array(index, dtype=np.uint32)
        node.data[CUBA.MATERIAL_ID] = self._geom.get_material_ijk(ijk)

        if node.data[CUBA.MATERIAL_ID] == ProxyLattice.FLUID_ENUM:
            vel = np.zeros(3, dtype=np.float64)
            frc = np.zeros(3, dtype=np.float64)
            self._fdata.get_vel_ijk(ijk, vel)
            self._fdata.get_frc_ijk(ijk, frc)

            node.data[CUBA.DENSITY] = self._fdata.get_den_ijk(ijk)
            node.data[CUBA.VELOCITY] = vel
            node.data[CUBA.FORCE] = frc

        return node

    def update_nodes(self, nodes):
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
        for node in nodes:
            if any(value < 0 for value in node.index):
                raise IndexError('invalid index: {}'.format(node.index))

            ijk = np.array(node.index, dtype=np.uint32)

            if self._geom.get_material_ijk(ijk) == ProxyLattice.FLUID_ENUM:
                if CUBA.DENSITY in node.data:
                    self._fdata.set_den_ijk(ijk, node.data[CUBA.DENSITY])
                if CUBA.VELOCITY in node.data:
                    vel = np.array(node.data[CUBA.VELOCITY], dtype=np.float64)
                    self._fdata.set_vel_ijk(ijk, vel)
                if CUBA.FORCE in node.data:
                    frc = np.array(node.data[CUBA.FORCE], dtype=np.float64)
                    self._fdata.set_frc_ijk(ijk, frc)

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
