from jyulb.internal.isothermal import solver
from jyulb.internal.common.proxy_lattice import ProxyLattice
from jyulb.internal.common.domain import PyLattice, PyGeometry
from simphony.cuds.abc_lattice import ABCLattice
from simphony.cuds.primitive_cell import BravaisLattice
from simphony.cuds.abc_modeling_engine import ABCModelingEngine
from jyulb.cuba_extension import CUBAExtension
from simphony.core.cuba import CUBA
import numpy as np


class JYULBEngine(ABCModelingEngine):

    """Internal wrapper for the JYU-LB Isothermal 3D flow modeling engine.

    Only a single lattice can be added for configuring the simulation;
    the wrapper does not accept meshes or particle containers.

    The following CUBA keywords are acknowledged in lattice node data:
    MATERIAL_ID, DENSITY, VELOCITY, and FORCE.

    Some values for the configuration parameters are enumerated:
    STOKES_FLOW/LAMINAR_FLOW/TURBULENT_FLOW for FLOW_TYPE,
    and BGK/TRT/MRT/REG for COLLISION_OPERATOR.

    Attributes
    ----------
    BC : dict
        container of attributes related to the boundary conditions.
    CM : dict
        container of attributes related to the computational method:
        collision operator, number of time steps, and time step.
    SP : dict
        container of attributes related to the system parameters/conditions:
        kinematic viscosity, reference density, gravity, flow type, and
        enforcement of external forcing.
    """

    # Enumeration of some values for the configuration parameters
    STOKES_FLOW_ENUM = solver.STOKES_FLOW_ENUM
    LAMINAR_FLOW_ENUM = solver.LAMINAR_FLOW_ENUM
    TURBULENT_FLOW_ENUM = solver.TURBULENT_FLOW_ENUM

    BGK_ENUM = solver.BGK_ENUM
    TRT_ENUM = solver.TRT_ENUM
    MRT_ENUM = solver.MRT_ENUM
    REG_ENUM = solver.REG_ENUM

    def __init__(self):
        """Initialize and set default parameters for CM, SP, and BC."""
        # Definition of CM, SP, and BC data components
        self._proxy_lattice = None
        self._solver = None
        self._pygeom = None
        self._prms = solver.PyFlowParams()
        self._is_fdata_initialized = False

        self.CM = {}
        self.SP = {}
        self.BC = {}

        # Default Computational Method data
        self.CM[CUBAExtension.COLLISION_OPERATOR] = JYULBEngine.TRT_ENUM
        self.CM[CUBA.NUMBER_OF_TIME_STEPS] = 10000
        self.CM[CUBA.TIME_STEP] = 1.0
        self.CM[CUBA.NAME] = 'default'

        # Default System Parameters data
        self.SP[CUBA.KINEMATIC_VISCOSITY] = 0.1
        self.SP[CUBAExtension.REFERENCE_DENSITY] = 1.0
        self.SP[CUBAExtension.GRAVITY] = (0.0, 0.0, 0.0)
        self.SP[CUBAExtension.FLOW_TYPE] = JYULBEngine.STOKES_FLOW_ENUM

        # Default Boundary Condition data
        self.BC[CUBA.VELOCITY] = {'open': 'periodic',
                                  'wall': 'noSlip'}

        self.BC[CUBA.DENSITY] = {'open': 'periodic',
                                 'wall': 'noFlux'}

    def run(self):
        """Run the modeling engine using the configured settings.

        Raises
        ------
        RuntimeError
           if a lattice has not been added.
        """
        if self._proxy_lattice is None:
            message = 'A lattice is not added before run in JYULBEngine'
            raise RuntimeError(message)

        if not self._is_fdata_initialized:
            self._solver.init_field_data()
            self._is_fdata_initialized = True

        self._solver.evolve(self.CM[CUBA.NUMBER_OF_TIME_STEPS])

    def add_dataset(self, container):
        """Add a CUDS Lattice container

        Parameters
        ----------
        container : {ABCLattice}
            The CUDS Lattice container to add to the engine.

        Raises
        ------
        TypeError:
            If the container type is not supported by the engine.
        ValueError:
            If there is already a dataset with the given name.
        ValueError:
            If the lattice type of the container is not cubic.

        """
        if not isinstance(container, ABCLattice):
            message = 'Only lattice containers are supported in JYULBEngine'
            raise TypeError(message)
        if bool(self._proxy_lattice):
            message = 'A lattice container already exists in JYULBEngine'
            raise ValueError(message)
        lat_type = container.primitive_cell.bravais_lattice
        if lat_type is not BravaisLattice.CUBIC:
            message = 'Lattice type is not cubic'
            raise ValueError(message)

        pylat = PyLattice(np.array(container.size, dtype=np.uint32),
                          np.array(container.origin, dtype=np.float64))
        self._pygeom = PyGeometry(pylat)

        for node in container.iter_nodes():
            ijk = np.array(node.index, dtype=np.uint32)
            self._pygeom.set_material_ijk(ijk, node.data[CUBA.MATERIAL_ID])

        self._prms.time_step = self.CM[CUBA.TIME_STEP]

        rden = self.SP[CUBAExtension.REFERENCE_DENSITY]
        self._prms.reference_density = rden

        self._prms.kinematic_viscosity = self.SP[CUBA.KINEMATIC_VISCOSITY]

        grav = np.array(self.SP[CUBAExtension.GRAVITY], dtype=np.float64)
        self._prms.gravity = grav

        self._prms.flow_type = self.SP[CUBAExtension.FLOW_TYPE]

        coll_op = self.CM[CUBAExtension.COLLISION_OPERATOR]
        self._prms.collision_operator = coll_op

        self._prms.external_forcing = True

        self._solver = solver.PySolver(self._pygeom, self._prms)
        pyfdata = self._solver.get_field_data()

        name = container.name
        pc = container.primitive_cell
        self._proxy_lattice = ProxyLattice(name, pc, self._pygeom, pyfdata)

        self._proxy_lattice.update_nodes(container.iter_nodes())

        self._proxy_lattice.data = container.data

    def remove_dataset(self, name):
        """Delete a lattice.

        Parameters
        ----------
        name : str
            name of the lattice to be deleted.

        Raises
        ------
        ValueError
            if name is not equal to the ProxyLattice name
            if no lattices are added in JYULBEngine

        """
        if self._proxy_lattice is not None:
            if self._proxy_lattice.name is not name:
                message = 'Container does not exist in JYULBEngine'
                raise ValueError(message)
            else:
                self._proxy_lattice = None
                self._is_fdata_initialized = False
                self._solver = None
        else:
            message = 'No lattices added in JYULBEngine'
            raise ValueError(message)

    def get_dataset(self, name):
        """ Get a lattice.

        The returned lattice can be used to query and update the state of the
        lattice inside the modeling engine.

        Parameters
        ----------
        name : str
            name of the lattice to be retrieved.

        Returns
        -------
        ABCLattice

        Raises
        ------
        ValueError
            if any one of the names is not equal to the ProxyLattice name
            if no lattices are added in JYULBEngine

        """
        if self._proxy_lattice is not None:
            if self._proxy_lattice.name is not name:
                message = 'Container does not exists in JYULBEngine'
                raise ValueError(message)
            else:
                return self._proxy_lattice
        else:
            message = 'No lattices added in JYULBEngine'
            raise ValueError(message)

    def get_dataset_names(self):
        """ Returns the names of the all the datasets in the engine workspace.

        """
        if self._proxy_lattice is not None:
            return [self._proxy_lattice.name]
        else:
            return []

    def iter_datasets(self, names=None):
        """Iterate over a subset or all of the lattices.

        Parameters
        ----------
        names : sequence of str, optional
            names of specific lattices to be iterated over. If names is not
            given, then all lattices will be iterated over.

        Returns
        -------
        A generator of ABCLattice objects

        Raises
        ------
        ValueError
            if any one of the names is not equal to the ProxyLattice name
            if no lattices are added in JYULBEngine
        """
        if self._proxy_lattice is not None:
            if names is None:
                yield self._proxy_lattice
            else:
                for name in names:
                    if self._proxy_lattice.name is not name:
                        message = 'State data does not contain requested item'
                        raise ValueError(message)
                    yield self._proxy_lattice
        else:
            message = 'No lattices added in JYULBEngine'
            raise ValueError(message)
