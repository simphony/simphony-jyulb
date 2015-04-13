import solver
from jyulb.internal.common.proxy_lattice import ProxyLattice
from jyulb.internal.common.domain import PyLattice, PyGeometry
from simphony.cuds.abc_modeling_engine import ABCModelingEngine
from simphony.core.data_container import DataContainer
from jyulb.cuba_extension import CUBAExtension
from simphony.core.cuba import CUBA
import numpy as np


class JYUEngine(ABCModelingEngine):

    """File-IO wrapper for the JYU-LB Isothermal 3D flow modeling engine.

    Only a single lattice can be added for configuring the simulation;
    the wrapper does not accept meshes or particle containers.

    The following CUBA keywords are acknowledged in lattice node data:
    MATERIAL_ID, DENSITY, VELOCITY, and FORCE.

    Some values for the configuration parameters are enumerated:
    STOKES_FLOW/LAMINAR_FLOW/TURBULENT_FLOW for FLOW_TYPE,
    and BGK/TRT/MRT/REG for COLLISION_OPERATOR.

    Attributes
    ----------
    BC : DataContainer
        container of attributes related to the boundary conditions.
    CM : DataContainer
        container of attributes related to the computational method:
        collision operator, number of time steps, and time step.
    SP : DataContainer
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
        """Initialize and set default parameters for CM, BC, SP, and SD."""
        # Definition of CM, SP, BC, and SD data components
        self._lat = None
        self._solver = None
        self._prms = solver.PyFlowParams()
        self._is_fdata_initialized = False

        self.CM = {}
        self.SP = {}
        self.BC = DataContainer()

        # Default Computational Method data
        self.CM[CUBAExtension.COLLISION_OPERATOR] = JYUEngine.TRT_ENUM
        self.CM[CUBA.NUMBER_OF_TIME_STEPS] = 10000
        self.CM[CUBA.TIME_STEP] = 1.0

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
        if self._lat is None:
            message = 'A lattice is not added before run in JYUEngine'
            raise RuntimeError(message)

        if not self._is_fdata_initialized:
            self._solver.init_field_data()
            self._is_fdata_initialized = True

        self._solver.evolve(self.CM[CUBA.NUMBER_OF_TIME_STEPS])

    def add_lattice(self, lattice):
        """Add a lattice to the modeling engine.

        The following CUBA keywords are acknowledged in node data:
        MATERIAL_ID, DENSITY, VELOCITY, and FORCE.

        Parameters
        ----------
        lattice : ABCLattice
            lattice to be added.

        Returns
        -------
        proxy : ABCLattice
            A lattice to be used to update/query the internal representation
            stored inside the modeling-engine. See get_lattice for more
            information.

        Raises
        ------
        RuntimeError
           if a second lattice is added.
        """
        if self._lat is not None:
            message = 'Not possible to add a second lattice in JYUEngine'
            raise RuntimeError(message)

        pylat = PyLattice(np.array(lattice.size, dtype=np.uint32),
                          np.array(lattice.origin, dtype=np.float64))
        pygeom = PyGeometry(pylat)

        for node in lattice.iter_nodes():
            ijk = np.array(node.index, dtype=np.uint32)
            pygeom.set_material_ijk(ijk, node.data[CUBA.MATERIAL_ID])

        self._prms.lattice_spacing = lattice.base_vect[0]
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

        self._solver = solver.PySolver(pygeom, self._prms)
        pyfdata = self._solver.get_field_data()

        self._lat = ProxyLattice(lattice.name, lattice.type,
                                 lattice.base_vect, pygeom, pyfdata)

        for node in lattice.iter_nodes():
            self._lat.update_node(node)

        return self._lat

    def add_mesh(self, mesh):
        """Add a mesh to the modeling engine.

        Parameters
        ----------
        mesh : ABCMesh
            mesh to be added.

        Returns
        -------
        proxy : ABCMesh
            A proxy mesh to be used to update/query the internal representation
            stored inside the modeling-engine. See get_mesh for more
            information.

        Raises
        ------
        NotImplementedError
           always.
        """
        message = 'JYUEngine does not handle meshes'
        raise NotImplementedError(message)

    def add_particles(self, particles):
        """Add a particle container to the modeling engine.

        Parameters
        ----------
        particles : ABCParticles
            particle container to be added.

        Returns
        -------
        ABCParticles
            A particle container to be used to update/query the internal
            representation stored inside the modeling-engine. See
            get_particles for more information.

        Raises
        ------
        NotImplementedError
           always.
        """
        message = 'JYUEngine does not handle particle containers'
        raise NotImplementedError(message)

    def delete_lattice(self, name):
        """Delete a lattice.

        Parameters
        ----------
        name : str
            name of the lattice to be deleted.
        """
        self._is_fdata_initialized = False
        self._solver = None
        self._lat = None

    def delete_mesh(self, name):
        """Delete a mesh.

        Parameters
        ----------
        name : str
            name of the mesh to be deleted.

        Raises
        ------
        NotImplementedError
           always.
        """
        message = 'JYUEngine does not handle meshes'
        raise NotImplementedError(message)

    def delete_particles(self, name):
        """Delete a particle container.

        Parameters
        ----------
        name : str
            name of the particle container to be deleted.

        Raises
        ------
        NotImplementedError
           always.
        """
        message = 'JYUEngine does not handle particle containers'
        raise NotImplementedError(message)

    def get_lattice(self, name):
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
        """
        return self._lat

    def get_mesh(self, name):
        """Get a mesh.

        The returned mesh can be used to query and update the state of the
        mesh inside the modeling engine.

        Parameters
        ----------
        name : str
            name of the mesh to be retrieved.

        Returns
        -------
        ABCMesh

        Raises
        ------
        NotImplementedError
           always.
        """
        message = 'JYUEngine does not handle meshes'
        raise NotImplementedError(message)

    def get_particles(self, name):
        """Get a particle container.

        The returned particle container can be used to query and update the
        state of the particle container inside the modeling engine.

        Parameters
        ----------
        name : str
            name of the particle container to be retrieved.

        Returns
        -------
        ABCParticles

        Raises
        ------
        NotImplementedError
           always.
        """
        message = 'JYUEngine does not handle particle containers'
        raise NotImplementedError(message)

    def iter_lattices(self, names=None):
        """Iterate over a subset or all of the lattices.

        Parameters
        ----------
        names : sequence of str, optional
            names of specific lattices to be iterated over. If names is not
            given, then all lattices will be iterated over.

        Yields
        -------
        ABCLattice
        """
        yield self._lat

    def iter_meshes(self, names=None):
        """Iterate over a subset or all of the meshes.

        Parameters
        ----------
        names : sequence of str, optional
            names of specific meshes to be iterated over. If names is not
            given, then all meshes will be iterated over.

        Yields
        -------
        ABCMesh

        Raises
        ------
        NotImplementedError
           always.
        """
        message = 'JYUEngine does not handle meshes'
        raise NotImplementedError(message)

    def iter_particles(self, names=None):
        """Iterate over a subset or all of the particle containers.

        Parameters
        ----------
        names : sequence of str, optional
            names of specific particle containers to be iterated over.
            If names is not given, then all particle containers will
            be iterated over.

        Yields
        -------
        ABCParticles

        Raises
        ------
        NotImplementedError
           always.
        """
        message = 'JYUEngine does not handle particle containers'
        raise NotImplementedError(message)
