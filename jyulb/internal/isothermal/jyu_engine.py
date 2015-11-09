from jyulb.internal.isothermal import solver
from jyulb.internal.common.proxy_lattice import ProxyLattice
from jyulb.internal.common.domain import PyLattice, PyGeometry
from simphony.cuds.abc_lattice import ABCLattice
from simphony.cuds.primitive_cell import BravaisLattice
from simphony.cuds.abc_modeling_engine import ABCModelingEngine
from jyulb.cuba_extension import CUBAExtension
from simphony.core.cuba import CUBA
import numpy as np


class JYUEngine(ABCModelingEngine):

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
    SD : dict
        container of state data
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
        """Initialize and set default parameters for CM, BC, SP, and SD."""
        # Definition of CM, SP, BC, and SD data components
        self._lattice_proxy = None
        self._solver = None
        self._pygeom = None
        self._prms = solver.PyFlowParams()
        self._is_fdata_initialized = False

        self.CM = {}
        self.SP = {}
        self._SD = {}
        self.BC = {}

        # Default Computational Method data
        self.CM[CUBAExtension.COLLISION_OPERATOR] = JYUEngine.TRT_ENUM
        self.CM[CUBA.NUMBER_OF_TIME_STEPS] = 10000
        self.CM[CUBA.TIME_STEP] = 1.0
        self.CM[CUBA.NAME] = 'default'

        # Default System Parameters data
        self.SP[CUBA.KINEMATIC_VISCOSITY] = 0.1
        self.SP[CUBAExtension.REFERENCE_DENSITY] = 1.0
        self.SP[CUBAExtension.GRAVITY] = (0.0, 0.0, 0.0)
        self.SP[CUBAExtension.FLOW_TYPE] = JYUEngine.STOKES_FLOW_ENUM

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
        if self._lattice_proxy is None:
            message = 'A lattice is not added before run in JYUEngine'
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
        self._update_dataset_names()
        if not isinstance(container, ABCLattice):
            message = 'Only lattice containers are supported in JYUEngine'
            raise TypeError(message)
        if bool(self._SD):
            message = 'A lattice container already exists in JYUEngine'
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
        self._lattice_proxy = ProxyLattice(name, pc, self._pygeom, pyfdata)

        self._lattice_proxy.update_nodes(container.iter_nodes())

        self._lattice_proxy.data = container.data

        self._SD[name] = self._lattice_proxy

    def remove_dataset(self, name):
        """ Remove a dataset from the internal

        Parameters
        ----------
        name: str
            name of CUDS container to be deleted

        Raises
        ------
        ValueError:
            If there is no dataset with the given name

        """
        self._update_dataset_names()
        if name not in self._SD:
            message = 'Container does not exists in JYUEngine'
            raise ValueError(message)
        else:
            del self._SD[name]
            self._lattice_proxy = None
            self._is_fdata_initialized = False
            self._solver = None

    def get_dataset(self, name):
        """ Get the dataset

        Parameters
        ----------
        name: str
            name of CUDS container to be retrieved.

        Returns
        -------
        container :
            A proxy of the dataset named ``name`` that is stored
            internally in the Engine.

        Raises
        ------
        ValueError:
            If there is no dataset with the given name

        """
        self._update_dataset_names()
        if name not in self._SD:
            message = 'Container does not exists in JYUEngine'
            raise ValueError(message)
        return self._lattice_proxy

    def get_dataset_names(self):
        """ Returns the names of the all the datasets in the engine workspace.

        """
        self._update_dataset_names()
        return self._SD.keys()

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
            if any one of the names is not in SD
        """
        self._update_dataset_names()
        if names is None:
            for name in self._SD.keys():
                yield self._SD[name]
        else:
            for name in names:
                if name not in self._SD.keys():
                    message = 'State data does not contain requested item'
                    raise ValueError(message)
                yield self._SD[name]

    def _update_dataset_names(self):
        """ Go through dataset names and update them to self._SD dictionary

        """
        for name in self._SD.keys():
            container = self._SD[name]
            if container.name is not name:
                self._SD[container.name] = self._SD.pop(name)
