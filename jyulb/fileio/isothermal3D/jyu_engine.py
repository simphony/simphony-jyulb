from jyulb.fileio.common.jyu_lattice_proxy import JYULatticeProxy
from simphony.cuds.abc_modeling_engine import ABCModelingEngine
from simphony.core.data_container import DataContainer
from simphony.core.cuba import CUBA
import numpy as np
import subprocess
import os


class JYUEngine(ABCModelingEngine):

    """File-IO wrapper for the JYU-LB Isothermal 3D flow modeling engine.

    Only a single lattice can be added for configuring the simulation;
    the wrapper does not accept meshes or particle containers.

    The following CUBA keywords are acknowledged in lattice node data:
    MATERIAL_ID, DENSITY, VELOCITY, and FORCE.

    Some values for the configuration parameters are enumerated: SOLID/FLUID
    for MATERIAL_ID, STOKES_FLOW/LAMINAR_FLOW/TURBULENT_FLOW for FLOW_TYPE,
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
    SOLID_ENUM = 0
    FLUID_ENUM = 255

    STOKES_FLOW_ENUM = 0
    LAMINAR_FLOW_ENUM = 1
    TURBULENT_FLOW_ENUM = 2

    BGK_ENUM = 0
    TRT_ENUM = 1
    MRT_ENUM = 2
    REG_ENUM = 3

    def __init__(self):
        """Initialize and set default parameters for CM, BC, SP, and SD."""
        # Definition of CM, SP, BC, and SD data components
        self._data = {}
        self._lattice = None
        self.CM = DataContainer()
        self.SP = DataContainer()
        self.BC = DataContainer()

        # Default Computational Method data
        self.CM_CUBA_COLLISION_OPERATOR = JYUEngine.TRT_ENUM
        self.CM[CUBA.NUMBER_OF_TIME_STEPS] = 10000
        self.CM[CUBA.TIME_STEP] = 1.0
        self.base_fname = 'jyu_engine'

        # Default System Parameters data
        self.SP[CUBA.KINEMATIC_VISCOSITY] = 0.1
        self.SP_CUBA_REFERENCE_DENSITY = 1.0
        self.SP_CUBA_GRAVITY = (0.0, 0.0, 0.0)
        self.SP_CUBA_FLOW_TYPE = JYUEngine.STOKES_FLOW_ENUM
        self.SP_CUBA_EXTERNAL_FORCING = False

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
           if a lattice has not been added or
           if execution of the modeling engine fails.
        """
        if self._lattice is None:
            message = 'A lattice is not added before run in JYUEngine'
            raise RuntimeError(message)

        # Define input/output file names
        input_script_fname = self.base_fname + '.in'
        evolut_info_fname = self.base_fname + '.evol'
        geom_write_fname = self.base_fname + '.geom.in.raw'
        den_write_fname = self.base_fname + '.den.in.raw'
        vel_write_fname = self.base_fname + '.vel.in.raw'
        frc_write_fname = self.base_fname + '.frc.in.raw'
        den_read_fname = self.base_fname + '.den.out.raw'
        vel_read_fname = self.base_fname + '.vel.out.raw'

        # Write configuration data to files
        self._data[CUBA.MATERIAL_ID].tofile(geom_write_fname)
        self._data[CUBA.DENSITY].tofile(den_write_fname)
        self._data[CUBA.VELOCITY].tofile(vel_write_fname)
        if self.SP_CUBA_EXTERNAL_FORCING:
            self._data[CUBA.FORCE].tofile(frc_write_fname)
        self._write_input_script(input_script_fname)

        # Run the modeling engine
        run_command = './jyu_lb_isothermal3D.exe ' + input_script_fname

        p = subprocess.Popen(run_command, shell=True,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT)

        output = p.communicate()
        print output[0]

        if p.returncode < 0:
            message = 'Execution of the JYU-LB simulator failed'
            raise RuntimeError(message)

        # Read the simulated flow field
        nx = self._lattice.size[0]
        ny = self._lattice.size[1]
        nz = self._lattice.size[2]

        den_data = self._data[CUBA.DENSITY]
        vel_data = self._data[CUBA.VELOCITY]

        # Copy data from the files into the given arrays
        den_data[:] = np.fromfile(den_read_fname,
                                  dtype=np.float64).reshape((nz, ny, nx))
        vel_data[:] = np.fromfile(vel_read_fname,
                                  dtype=np.float64).reshape((nz, ny, nx, 3))

        # Clean up
        os.remove(input_script_fname)
        os.remove(evolut_info_fname)
        os.remove(geom_write_fname)
        os.remove(den_write_fname)
        os.remove(vel_write_fname)
        os.remove(den_read_fname)
        os.remove(vel_read_fname)
        if self.SP_CUBA_EXTERNAL_FORCING:
            os.remove(frc_write_fname)

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
        if self._lattice is not None:
            message = 'Not possible to add a second lattice in JYUEngine'
            raise RuntimeError(message)

        # Copy lattice attributes
        name = lattice.name
        type = lattice.type
        nx = lattice.size[0]
        ny = lattice.size[1]
        nz = lattice.size[2]
        org = lattice.origin
        bvect = lattice.base_vect

        # Allocate arrays for lattice data
        geom = np.zeros((nz, ny, nx), dtype=np.uint8)
        den = np.zeros((nz, ny, nx), dtype=np.float64)
        vel = np.zeros((nz, ny, nx, 3), dtype=np.float64)
        frc = np.zeros((nz, ny, nx, 3), dtype=np.float64)

        self._data[CUBA.MATERIAL_ID] = geom
        self._data[CUBA.DENSITY] = den
        self._data[CUBA.VELOCITY] = vel
        self._data[CUBA.FORCE] = frc

        # Create a proxy lattice
        self._lattice = JYULatticeProxy(name, type, bvect, (nx, ny, nz),
                                        org, self._data)

        for node in lattice.iter_nodes():
            self._lattice.update_node(node)

        return self._lattice

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

    def add_particle_container(self, particle_container):
        """Add a particle container to the modeling engine.

        Parameters
        ----------
        particle_container : ABCParticleContainer
            particle container to be added.

        Returns
        -------
        ABCParticleContainer
            A particle container to be used to update/query the internal
            representation stored inside the modeling-engine. See
            get_particle_container for more information.

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
        self._data = {}
        self._lattice = None

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

    def delete_particle_container(self, name):
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
        return self._lattice

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

    def get_particle_container(self, name):
        """Get a particle container.

        The returned particle container can be used to query and update the
        state of the particle container inside the modeling engine.

        Parameters
        ----------
        name : str
            name of the particle container to be retrieved.

        Returns
        -------
        ABCParticleContainer

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
        yield self._lattice

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

    def iter_particle_containers(self, names=None):
        """Iterate over a subset or all of the particle containers.

        Parameters
        ----------
        names : sequence of str, optional
            names of specific particle containers to be iterated over.
            If names is not given, then all particle containers will
            be iterated over.

        Yields
        -------
        ABCParticleContainer

        Raises
        ------
        NotImplementedError
           always.
        """
        message = 'JYUEngine does not handle particle containers'
        raise NotImplementedError(message)

    def _write_input_script(self, fname):
        """Write an input script file for the modeling engine.

        Parameters
        ----------
        fname : str
            name of the script file.
        """
        f = open(fname, 'w')
        f.write('# Base name of the I/O data files (raw-format files)\n')
        f.write(self.base_fname + '\n')
        f.write('# Lattice size (x y z)\n')
        f.write('%d %d %d\n' % self._lattice.size)
        f.write('# Lattice spacing\n')
        f.write('%e\n' % self._lattice.base_vect[0])
        f.write('# Discrete time step\n')
        f.write('%e\n' % self.CM[CUBA.TIME_STEP])
        f.write('# Reference density\n')
        f.write('%e\n' % self.SP_CUBA_REFERENCE_DENSITY)
        f.write('# Kinematic viscosity\n')
        f.write('%e\n' % self.SP[CUBA.KINEMATIC_VISCOSITY])
        f.write('# Gravity (gx gy gz)\n')
        f.write('%e %e %e\n' % self.SP_CUBA_GRAVITY)

        f.write('# Additional external forcing (No=0, Yes=1)\n')
        if self.SP_CUBA_EXTERNAL_FORCING:
            f.write('1\n')
        else:
            f.write('0\n')

        f.write('# Flow type (Stokes=0, Laminar=1, Turbulent=2)\n')
        f.write('%d\n' % self.SP_CUBA_FLOW_TYPE)
        f.write('# Collision operator (BGK=0, TRT=1, MRT=2, Regul.=3)\n')
        f.write('%d\n' % self.CM_CUBA_COLLISION_OPERATOR)
        f.write('# Number of discrete time steps\n')
        f.write('%d\n' % self.CM[CUBA.NUMBER_OF_TIME_STEPS])

        evol_info_interval = int(self.CM[CUBA.NUMBER_OF_TIME_STEPS]/20)
        if evol_info_interval < 1:
            evol_info_interval = 1

        f.write('# Evolution info output interval (0=No output)\n')
        f.write('%d\n' % evol_info_interval)
        f.write('# Execution info\n')
        f.write('Running JYU-LB software via SimPhoNy file-I/O wrapper\n')

        f.close()
