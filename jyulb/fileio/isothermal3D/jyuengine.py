from jyulb.fileio.common.jyu_lattice_proxy import JYULatticeProxy
from simphony.cuds.abc_modeling_engine import ABCModelingEngine
from simphony.core.data_container import DataContainer
from simphony.core.cuba import CUBA
import numpy as np
import subprocess
import math


class JYUEngine(ABCModelingEngine):
    """Abstract base class for modeling engines in SimPhoNy.

    Through this interface, the user controls and interacts with the
    simulation/calculation (which is being performed by the modeling
    engine).

    Attributes
    ----------
    BC : DataContainer
        container of attributes related to the boundary conditions
    CM : DataContainer
        container of attributes related to the computational method
    SP : DataContainer
        container of attributes related to the system parameters/conditions

    """
    # Enumeration of CM and SP data values
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
        # Definition of CM, SP, BC, and SD data components
        self.CM = DataContainer()
        self.SP = DataContainer()
        self.BC = DataContainer()
        self._lattice = None
        self._data = {}

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
        
    def run(self):
        """Run the modeling engine

        Run the modeling engine using the configured settings (e.g. CM, BC,
        and SP) and the configured state data (e.g. particle, mesh and
        lattice data).

        """
        if self._lattice is None:
            message = 'A lattice is not added before run in JYUEngine'
            raise RuntimeError(message)
            
        input_script_fname = self.base_fname + '.in'
        evolut_info_fname = self.base_fname + '.evol'
        geom_write_fname = self.base_fname + '.geom.in.raw'
        den_write_fname = self.base_fname + '.den.in.raw'
        vel_write_fname = self.base_fname + '.vel.in.raw'
        frc_write_fname = self.base_fname + '.frc.in.raw'
        den_read_fname = self.base_fname + '.den.out.raw'
        vel_read_fname = self.base_fname + '.vel.out.raw'
                
#        self._data[CUBA.MATERIAL_ID].tofile(geom_write_fname)
#        self._data[CUBA.DENSITY].tofile(den_write_fname)
#        self._data[CUBA.VELOCITY].tofile(vel_write_fname)

#        if self.SP_CUBA_EXTERNAL_FORCING:
#            self._data[CUBA.FORCE].tofile(frc_write_fname)

        self._write_input_script(input_script_fname)
        run_command = './jyu_lb_isothermal3D.exe ' + input_script_fname
        
#        p = subprocess.Popen(run_command, shell=True,
#                             stdout=subprocess.PIPE,
#                             stderr=subprocess.STDOUT)
#        retval = p.wait()

        nx = self._lattice.size[0]        
        ny = self._lattice.size[1]        
        nz = self._lattice.size[2]

        den_data = self._data[CUBA.DENSITY]
        vel_data = self._data[CUBA.VELOCITY]
        
#        den_data[:] = np.fromfile(den_read_fname,
#                                  dtype=np.float64).reshape((nz, ny, nx))
#        vel_data[:] = np.fromfile(vel_read_fname,
#                                  dtype=np.float64).reshape((nz, ny, nx, 3))
                                  
#        os.remove(input_script_fname)                           
#        os.remove(evolut_info_fname)                           
#        os.remove(geom_write_fname)                           
#        os.remove(den_write_fname)                           
#        os.remove(vel_write_fname)                           
#        os.remove(den_read_fname)                           
#        os.remove(vel_read_fname)                           

#        if self.SP_CUBA_EXTERNAL_FORCING:
#            os.remove(frc_write_fname)                           
        
    def add_lattice(self, lattice):
        """Add lattice to the modeling engine

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

        """
        if self._lattice is not None:
            message = 'Not possible to add a second lattice in JYUEngine'
            raise RuntimeError(message)
        
        name = lattice.name
        type = lattice.type
        nx = lattice.size[0]        
        ny = lattice.size[1]        
        nz = lattice.size[2]
        org = lattice.origin
        bvect = lattice.base_vect

        geom = np.zeros((nz, ny, nx), dtype=np.uint8)
        den = np.zeros((nz, ny, nx), dtype=np.float64)
        vel = np.zeros((nz, ny, nx, 3), dtype=np.float64)
        frc = np.zeros((nz, ny, nx, 3), dtype=np.float64)

        self._data[CUBA.MATERIAL_ID] = geom
        self._data[CUBA.DENSITY] = den
        self._data[CUBA.VELOCITY] = vel
        self._data[CUBA.FORCE] = frc
        
        self._lattice = JYULatticeProxy(name, type, bvect, (nx, ny, nz),
                                        org, self._data)
        
        for node in lattice.iter_nodes():
            self._lattice.update_node(node)
            
        return self._lattice

    def add_mesh(self, mesh):
        """Add mesh to the modeling engine

        Parameters
        ----------
        mesh: ABCMesh
            mesh to be added.

        Returns
        -------
        proxy : ABCMesh
            A proxy mesh to be used to update/query the internal representation
            stored inside the modeling-engine. See get_mesh for more
            information.
        """
        message = 'JYUEngine does not handle meshes'
        raise NotImplementedError(message)
        
    def add_particle_container(self, particle_container):
        """Add particle container to the modeling engine

        Parameters
        ----------
        particle_container: ABCParticleContainer
            particle container to be added.

        Returns
        -------
        ABCParticleContainer
            A particle container to be used to update/query the internal
            representation stored inside the modeling-engine. See
            get_particle_container for more information.

        """
        message = 'JYUEngine does not handle particle containers'
        raise NotImplementedError(message)
        
    def delete_lattice(self, name):
        """Delete a lattice

        Parameters
        ----------
        name: str
            name of lattice to be deleted

        """
        self._data = {}
        self._lattice = None

    def delete_mesh(self, name):
        """Delete a mesh

        Parameters
        ----------
        name: str
            name of mesh to be deleted

        """
        message = 'JYUEngine does not handle meshes'
        raise NotImplementedError(message)

    def delete_particle_container(self, name):
        """Delete a particle container

        Parameters
        ----------
        name: str
            name of particle container to be deleted

        """
        message = 'JYUEngine does not handle particle containers'
        raise NotImplementedError(message)

    def get_lattice(self, name):
        """ Get lattice

        The returned lattice can be used to query and update the state of the
        lattice inside the modeling engine.

        Returns
        -------
        ABCLattice

        """
        return self._lattice

    def get_mesh(self, name):
        """ Get mesh

        The returned mesh can be used to query and update the state of the
        mesh inside the modeling engine.

        Returns
        -------
        ABCMesh

        """
        message = 'JYUEngine does not handle meshes'
        raise NotImplementedError(message)
        
    def get_particle_container(self, name):
        """ Get particle container

        The returned particle container can be used to query and update the
        state of the particle container inside the modeling engine.

        Returns
        -------
        ABCParticleContainer

        """
        message = 'JYUEngine does not handle particle containers'
        raise NotImplementedError(message)
        
    def iter_lattices(self, names=None):
        """ Returns an iterator over a subset or all of the lattices.

        Parameters
        ----------
        names : sequence of str, optional
            names of specific lattices to be iterated over. If names is not
            given, then all lattices will be iterated over.

        """
        yield self._lattice

    def iter_meshes(self, names=None):
        """ Returns an iterator over a subset or all of the meshes.

        Parameters
        ----------
        names : sequence of str, optional
            names of specific meshes to be iterated over. If names is not
            given, then all meshes will be iterated over.

        """
        message = 'JYUEngine does not handle meshes'
        raise NotImplementedError(message)
        
    def iter_particle_containers(self, names=None):
        """ Returns an iterator over a subset or all of the particle containers.

        Parameters
        ----------
        names : sequence of str, optional
            names of specific particle containers to be iterated over.
            If names is not given, then all particle containers will
            be iterated over.

        """
        message = 'JYUEngine does not handle particle containers'
        raise NotImplementedError(message)

    def _write_input_script(self, fname):
        """ Returns an iterator over a subset or all of the particle containers.

        Parameters
        ----------
        names : sequence of str, optional
            names of specific particle containers to be iterated over.
            If names is not given, then all particle containers will
            be iterated over.

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

#        evol_info_interval = math.floor(self.CM[CUBA.NUMBER_OF_TIME_STEPS]/20)
        evol_info_interval = int(self.CM[CUBA.NUMBER_OF_TIME_STEPS]/20)
        if evol_info_interval < 1:
            evol_info_interval = 1

        f.write('# Evolution info output interval (0=No output)\n')
        f.write('%d\n' % evol_info_interval)
        f.write('# Execution info\n')
        f.write('Running JYU-LB software via SimPhoNy file-I/O wrapper\n')

        f.close()