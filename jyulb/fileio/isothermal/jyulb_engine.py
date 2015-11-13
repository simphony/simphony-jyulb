from jyulb.fileio.common.proxy_lattice import ProxyLattice
from simphony.cuds.abc_lattice import ABCLattice
from simphony.cuds.primitive_cell import BravaisLattice
from simphony.cuds.abc_modeling_engine import ABCModelingEngine
from jyulb.cuba_extension import CUBAExtension
from simphony.core.cuba import CUBA
import numpy as np
import subprocess
import os


class JYULBEngine(ABCModelingEngine):

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
    STOKES_FLOW_ENUM = 0
    LAMINAR_FLOW_ENUM = 1
    TURBULENT_FLOW_ENUM = 2

    BGK_ENUM = 0
    TRT_ENUM = 1
    MRT_ENUM = 2
    REG_ENUM = 3

    def __init__(self):
        """Initialize and set default parameters for CM, SP, and BC."""
        # Definition of CM, SP, and BC data components
        self._data = {}
        self._proxy_lattice = None
        self.CM = {}
        self.SP = {}
        self.BC = {}

        # Default Computational Method data
        self.CM[CUBAExtension.COLLISION_OPERATOR] = JYULBEngine.TRT_ENUM
        self.CM[CUBA.NUMBER_OF_TIME_STEPS] = 10000
        self.CM[CUBA.TIME_STEP] = 1.0
        self.base_fname = 'jyu_engine'

        # Default System Parameters data
        self.SP[CUBA.KINEMATIC_VISCOSITY] = 0.1
        self.SP[CUBAExtension.REFERENCE_DENSITY] = 1.0
        self.SP[CUBAExtension.GRAVITY] = (0.0, 0.0, 0.0)
        self.SP[CUBAExtension.FLOW_TYPE] = JYULBEngine.STOKES_FLOW_ENUM
        self.SP[CUBAExtension.EXTERNAL_FORCING] = False

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
        if self._proxy_lattice is None:
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
        if self.SP[CUBAExtension.EXTERNAL_FORCING]:
            self._data[CUBA.FORCE].tofile(frc_write_fname)
        self._write_input_script(input_script_fname)

        # Run the modeling engine
        run_command = 'jyu_lb_isothermal.exe ' + input_script_fname

        p = subprocess.Popen(run_command, shell=True,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT)

        output = p.communicate()
        print output[0]

        if p.returncode < 0:
            message = 'Execution of the JYU-LB simulator failed'
            raise RuntimeError(message)

        # Read the simulated flow field
        nx = self._proxy_lattice.size[0]
        ny = self._proxy_lattice.size[1]
        nz = self._proxy_lattice.size[2]

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
        if self.SP[CUBAExtension.EXTERNAL_FORCING]:
            os.remove(frc_write_fname)

    def add_dataset(self, container):
        """Add a CUDS Lattice container

        The following CUBA keywords are acknowledged in node data:
        MATERIAL_ID, DENSITY, VELOCITY, and FORCE.

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

        # Copy lattice attributes
        name = container.name
        pc = container.primitive_cell
        nx = container.size[0]
        ny = container.size[1]
        nz = container.size[2]
        org = container.origin

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
        self._proxy_lattice = ProxyLattice(name, pc, (nx, ny, nz),
                                           org, self._data)

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
                self._data = {}
                self._proxy_lattice._data = None
                self._proxy_lattice = None
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
        f.write('%d %d %d\n' % self._proxy_lattice.size)
        f.write('# Lattice spacing\n')
        p1 = self._proxy_lattice.primitive_cell.p1
        f.write('%e\n' % np.sqrt(np.dot(p1, p1)))
        f.write('# Discrete time step\n')
        f.write('%e\n' % self.CM[CUBA.TIME_STEP])
        f.write('# Reference density\n')
        f.write('%e\n' % self.SP[CUBAExtension.REFERENCE_DENSITY])
        f.write('# Kinematic viscosity\n')
        f.write('%e\n' % self.SP[CUBA.KINEMATIC_VISCOSITY])
        f.write('# Gravity (gx gy gz)\n')
        f.write('%e %e %e\n' % self.SP[CUBAExtension.GRAVITY])

        f.write('# Additional external forcing (No=0, Yes=1)\n')
        if self.SP[CUBAExtension.EXTERNAL_FORCING]:
            f.write('1\n')
        else:
            f.write('0\n')

        f.write('# Flow type (Stokes=0, Laminar=1, Turbulent=2)\n')
        f.write('%d\n' % self.SP[CUBAExtension.FLOW_TYPE])
        f.write('# Collision operator (BGK=0, TRT=1, MRT=2, Regul.=3)\n')
        f.write('%d\n' % self.CM[CUBAExtension.COLLISION_OPERATOR])
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
