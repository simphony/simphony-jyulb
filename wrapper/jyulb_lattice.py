"""A proxy lattice for accessing simulation data of a JYU-LB engine."""
import uuid
import numpy as np

from simphony.core.data_container import DataContainer
from simphony.core.cuba import CUBA

from simphony.cuds.abc_lattice import ABCLattice
from simphony.cuds.lattice_items import LatticeNode
from simphony.cuds.primitive_cell import PrimitiveCell

from jyulb.defs import X, Y, Z, SOLID


class ProxyLattice(ABCLattice):

    """Access to the simulation data of a JYU-LB engine.

    Lattice node data
    =================
    Provided by get
    ---------------
      MATERIAL, PRESSURE, DENSITY, VELOCITY, ACCELERATION

    Acknowledged by update
    ----------------------
      PRESSURE (or alternatively DENSITY),
      VELOCITY (or alternatively MOMENTUM),
      ACCELERATION (or alternatively FORCE)
    """

    cuba_key = CUBA.LATTICE

    def __init__(self, src_lat, jyulb_wrapper):
        """Constructor."""
        # Copy data from the source lattice
        self.name = src_lat.name

        pc = src_lat._primitive_cell
        self._primitive_cell = PrimitiveCell(pc.p1, pc.p2, pc.p3,
                                             pc.bravais_lattice)

        self._size = np.array(src_lat.size, dtype=np.int32)
        self._origin = np.array(src_lat.origin, dtype=np.float)

        self._data = DataContainer(src_lat.data)

        self._items_count = {
            CUBA.NODE: lambda: self._size
        }
        self._uid = uuid.uuid4()

        # References to the simulation data
        self._phase = jyulb_wrapper._phase
        self._den = jyulb_wrapper._den
        self._vel = jyulb_wrapper._vel
        self._acc = jyulb_wrapper._acc

        # Parameters for unit transformations
        # (SI <-> reduced unit system)
        self._ref_den = jyulb_wrapper._ref_den
        self._inv_ref_den = 1.0/self._ref_den

        self._dr = jyulb_wrapper._dr
        self._inv_dr = 1.0/self._dr

        self._dt = jyulb_wrapper._dt
        self._inv_dt = 1.0/self._dt

        self._cr = self._dr/self._dt
        self._inv_cr = 1.0/self._cr

        # References to the fluid and solid materials
        self._fluid_mat = jyulb_wrapper._fluid_mat
        self._solid_mat = jyulb_wrapper._solid_mat

        # Reference to a JYU-LB engine
        self._solver = jyulb_wrapper._solver

    @property
    def uid(self):
        """Universal identifier."""
        return self._uid

    def count_of(self, item_type):
        """ Return the count of item_type in the container.

        Parameters
        ----------
        item_type : CUBA
            The CUBA enum of the type of the items to return
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

    @property
    def size(self):
        """Access (get)."""
        return self._size

    @property
    def origin(self):
        """Access (get)."""
        return self._origin

    @property
    def data(self):
        """Access (get)."""
        return self._data

    @data.setter
    def data(self, value):
        """Access (set)."""
        self._data = DataContainer(value)

    def _get_node(self, index):
        """Get a copy of the node corresponding to the given index.

        Parameters
        ----------
        index : int[3]
            node index coordinate

        Returns
        -------
        A reference to a LatticeNode object

        """
        ijk = tuple(index)
        ijkx = ijk + (X,)
        ijky = ijk + (Y,)
        ijkz = ijk + (Z,)

        if ijk[0] < 0 or ijk[1] < 0 or ijk[2] < 0:
            raise IndexError('invalid index: {}'.format(ijk))

        node = LatticeNode(index)

        if self._phase[ijk] == SOLID:
            node.data[CUBA.MATERIAL] = self._solid_mat
        else:
            node.data[CUBA.MATERIAL] = self._fluid_mat

        # Transformation from the reduced unit system to SI-units
        den = self._den[ijk]                  # reduced units
        p = self._solver.eos_den2p(den)       # reduced units
        p *= self._ref_den*self._cr*self._cr  # SI
        node.data[CUBA.PRESSURE] = p
        node.data[CUBA.DENSITY] = self._ref_den*den

        vx = self._vel[ijkx]  # reduced units
        vy = self._vel[ijky]  # reduced units
        vz = self._vel[ijkz]  # reduced units
        vx *= self._cr  # SI
        vy *= self._cr  # SI
        vz *= self._cr  # SI
        node.data[CUBA.VELOCITY] = (vx, vy, vz)

        ax = self._acc[ijkx]  # reduced units
        ay = self._acc[ijky]  # reduced units
        az = self._acc[ijkz]  # reduced units
        ax *= self._cr*self._inv_dt  # SI
        ay *= self._cr*self._inv_dt  # SI
        az *= self._cr*self._inv_dt  # SI
        node.data[CUBA.ACCELERATION] = (ax, ay, az)

        return node

    def _update_nodes(self, nodes):
        """Update the corresponding lattice nodes (data copied).

        Parameters
        ----------
        nodes : iterable of LatticeNode objects
            reference to LatticeNode objects from where the data is copied
            to the Lattice
        """
        for node in nodes:
            ijk = tuple(node.index)
            ijkx = ijk + (X,)
            ijky = ijk + (Y,)
            ijkz = ijk + (Z,)

            if any(value < 0 for value in ijk):
                raise IndexError('invalid index: {}'.format(ijk))

            # Transformation from SI-units to the reduced unit system
            if CUBA.PRESSURE in node.data:
                p = node.data[CUBA.PRESSURE]  # in SI-units
                p *= self._inv_ref_den*self._inv_cr*self._inv_cr
                self._den[ijk] = self._solver.eos_p2den(p)
            elif CUBA.DENSITY in node.data:
                self._den[ijk] = self._inv_ref_den*node.data[CUBA.DENSITY]

            if CUBA.VELOCITY in node.data:
                self._vel[ijkx] = self._inv_cr*node.data[CUBA.VELOCITY][X]
                self._vel[ijky] = self._inv_cr*node.data[CUBA.VELOCITY][Y]
                self._vel[ijkz] = self._inv_cr*node.data[CUBA.VELOCITY][Z]
            elif CUBA.MOMENTUM in node.data:
                # interpreted as momentum density
                den = self._ref_den*self._den[ijk]        # in SI-units
                inv_den = 1.0/den                         # in SI-units
                vx = inv_den*node.data[CUBA.MOMENTUM][X]  # in SI-units
                vy = inv_den*node.data[CUBA.MOMENTUM][Y]  # in SI-units
                vz = inv_den*node.data[CUBA.MOMENTUM][Z]  # in SI-units
                self._vel[ijkx] = self._inv_cr*vx
                self._vel[ijky] = self._inv_cr*vy
                self._vel[ijkz] = self._inv_cr*vz

            if CUBA.ACCELERATION in node.data:
                ax = node.data[CUBA.ACCELERATION][X]  # in SI-units
                ay = node.data[CUBA.ACCELERATION][Y]  # in SI-units
                az = node.data[CUBA.ACCELERATION][Z]  # in SI-units
                self._acc[ijkx] = self._inv_cr*self._dt*ax
                self._acc[ijky] = self._inv_cr*self._dt*ay
                self._acc[ijkz] = self._inv_cr*self._dt*az
            elif CUBA.FORCE in node.data:
                # interpreted as body force (i.e. force density)
                den = self._ref_den*self._den[ijk]     # in SI-units
                inv_den = 1.0/den                      # in SI-units
                ax = inv_den*node.data[CUBA.FORCE][X]  # in SI-units
                ay = inv_den*node.data[CUBA.FORCE][Y]  # in SI-units
                az = inv_den*node.data[CUBA.FORCE][Z]  # in SI-units
                self._acc[ijkx] = self._inv_cr*self._dt*ax
                self._acc[ijky] = self._inv_cr*self._dt*ay
                self._acc[ijkz] = self._inv_cr*self._dt*az

    def _iter_nodes(self, indices=None):
        """Get an iterator over the LatticeNodes described by the indices.

        Parameters
        ----------
        indices : iterable set of int[3], optional
            When indices (i.e. node index coordinates) are provided, then
            nodes are returned in the same order of the provided indices.
            If indices is None, there is no restriction on the order of the
            returned nodes.

        Returns
        -------
        A generator for LatticeNode objects

        """
        sx = self.size[X]
        sy = self.size[Y]
        sz = self.size[Z]

        if indices is None:
            for index in np.ndindex((sx, sy, sz)):
                yield self.get(index)
        else:
            for index in indices:
                yield self.get(index)
