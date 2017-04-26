"""JYU-LB wrapper for the SimPhoNy framework.

The wrapper relies on the jyulb-extension module
(i.e. an internal interface wrapper).
"""
import numpy as np

from simphony.cuds.meta import api
from simphony.core.cuba import CUBA
from simphony.cuds.abc_lattice import ABCLattice
from simphony.cuds.primitive_cell import BravaisLattice
from simphony.cuds.abc_modeling_engine import ABCModelingEngine

from jyulb.isothermal import Solver
from jyulb.defs import FACE_1X, FACE_X1, FACE_1Y, FACE_Y1, FACE_1Z, FACE_Z1
from jyulb.defs import PERIODIC, FIXED_DEN, FIXED_VEL
from jyulb.defs import BGK, TRT, REGULARIZATION
from jyulb.defs import FLUID, SOLID
from jyulb.defs import X, Y, Z

from .jyulb_lattice import ProxyLattice


class Wrapper(ABCModelingEngine):

    """JYU-LB wrapper (internal interface) for the SimPhoNy framework.

    CUDS configuration
    ==================
    Physics model
    -------------
      - Expecting CFD (exactly one) with submodels
        incompressible,
        single phase,
        isothermal,
        newtonian,
        laminar
      - Acceleration due to gravity is read from the Gravity submodel of CFD
        (i.e. Gravity models not associated with the CFD are ignored);
        acceleration due to gravity is updated before each run
      - The electrostatic submodel of CFD is ignored

    Materials
    ---------
      - Assume 1 (fluid) or 2 (fluid and solid/wall) materials
      - If 2 materials, one material must be named 'solid' or 'wall'
      - The non-solid/non-wall material is interpreted as fluid with data
        DENSITY or REFERENCE_DENSITY and
        KINEMATIC_VISCOSITY or DYNAMIC_VISCOSITY or
        VISCOSITY (interpreted as DYNAMIC_VISCOSITY)

    Dataset
    -------
      - Cubic Bravais lattice (exactly one),
        length of a primitive vector defines the lattice spacing
      - Lattice node data acknowledged includes
        MATERIAL (the default material is fluid if not specified),
        PRESSURE (or alternatively DENSITY),
        VELOCITY (or alternatively MOMENTUM),
        ACCELERATION (or alternatively FORCE);
        see also the documentation of ProxyLattice

    Solver parameters
    -----------------
      - TimeIntegration (exactly one),
        the attribute size defines the discrete time step;
        attributes final and current are updated by the wrapper
        after each run
      - The collision operator can be specified by providing
        a SolverParameter with data COLLISION_OPERATOR
        (value 'bgk' or 'trt' or 'regularization');
        e.g. the TimeIntegration can be used for this purpose

    Boundaries
    ----------
      - In order to configure inlet/outlet/etc. boundary conditions,
        exterior/outer boundary faces of the simulation domain
        must be specified
      - Exterior boundary faces of the simulation domain are specified
        by providing Boundaries with particular unit outward normals,
        i.e. expecting Boundaries with data DIRECTION = (-1,0,0) or
        (1,0,0) or (0,-1,0) or (0,1,0) or (0,0,-1) or (0,0,1)
      - Two or more Boundaries with equal DIRECTION value (see above)
        is not allowed

    Conditions
    ----------
      - Each Boundary (see above) can be associated with conditions
        (i.e. conditions not associated with Boundaries are ignored)
      - The default boundary condition is periodic
      - Inlet boundary condition (non-periodic) can be assigned for
        exactly one exterior boundary face, and the opposite boundary face
        is then expected to be the outlet (with non-periodic condition)
      - Allowed inlet/outlet boundary conditions are
        Dirichlet with data VARIABLE = PRESSURE
        (where also the in-plane velocity components are fixed) or
        Dirichlet with data VARIABLE = VELOCITY
      - Note that
        Dirichlet conditions must be assigned for the fluid material,
        only one Dirichlet condition per boundary is acknowledged, and
        the boundary condition data is read directly from the flow field
        (i.e. from the data of associated lattice nodes)
    """

    def __init__(self, **kwargs):
        """Constructor."""
        # Solver and simulation data
        self._solver = None
        self._proxy_lat = None
        self._is_init_finalized = False

        self._phase = None
        self._den = None
        self._vel = None
        self._acc = None

        # CFD model
        self._is_incompressible = True
        self._is_singlephase = True
        self._is_isothermal = True
        self._is_newtonian = True
        self._is_laminar = True
        self._cfd = None

        # Gravity
        self._g = np.zeros((3), dtype=np.float64)

        # Materials
        self._fluid_mat = None
        self._solid_mat = None
        self._material_cnt = 0
        self._kvisc = 1.0/6.0
        self._ref_den = 1.0
        self._inv_ref_den = 1.0

        # Lattice
        self._lat = None
        self._lat_size = np.zeros((3), dtype=np.int32)
        self._lat_fnode_cnt = 0  # number of fluid lattice sites
        self._dr = 0             # lattice spacing
        self._inv_dr = 0
        self._V = 0              # volume of the lattice cell
        self._inv_V = 0

        # Solver parameters
        self._time_integr = None
        self._coll_oper = BGK
        self._dt = 0              # discrete time step
        self._inv_dt = 0
        self._cr = 0              # lattice reference velocity cr = dr/dt
        self._inv_cr = 0.0

        # Boundary conditions
        self._bc = np.full((6), PERIODIC, dtype=np.int32)

        # Call the base class in order to load CUDS
        super(Wrapper, self).__init__(**kwargs)

    def _update_gravity(self):
        """Read acceleration due to gravity from the CFD model."""
        self._g[0] = self._cfd.gravity_model.acceleration[0]
        self._g[1] = self._cfd.gravity_model.acceleration[1]
        self._g[2] = self._cfd.gravity_model.acceleration[2]

    def _check_cuds_model(self, cuds):
        """Check physics model in the given cuds."""
        print 'Physics model'
        print '-------------'

        # General CFD model
        cfd_cnt = cuds.count_of(CUBA.CFD)
        self._cfd = next(cuds.iter(item_type=CUBA.CFD), None)

        if cfd_cnt != 1:
            err_str = 'Expecting exactly one CFD physics model '
            err_str += '(not {0:d})!'.format(cfd_cnt)
            raise ValueError(err_str)

        # Electrostatic submodel (ignored)
        print '* Ignoring CFD electrostatic model'

        # Gravity submodel
        self._update_gravity()

        # Compressibility submodel
        self._is_incompressible = True
        if not isinstance(self._cfd.compressibility_model,
                          api.IncompressibleFluidModel):
            err_str = 'Compressibility models are not yet supported!'
            raise NotImplementedError(err_str)

        # Multiphase submodel
        self._is_singlephase = True
        if not isinstance(self._cfd.multiphase_model,
                          api.SinglePhaseModel):
            err_str = 'Multiphase models are not yet supported!'
            raise NotImplementedError(err_str)

        # Thermal submodel
        self._is_isothermal = True
        if not isinstance(self._cfd.thermal_model,
                          api.IsothermalModel):
            err_str = 'Thermal models are not yet supported!'
            raise NotImplementedError(err_str)

        # Rheology submodel
        self._is_newtonian = True
        if not isinstance(self._cfd.rheology_model,
                          api.NewtonianFluidModel):
            err_str = 'Rheology models are not yet supported!'
            raise NotImplementedError(err_str)

        # Turbulence submodel
        self._is_laminar = True
        if not isinstance(self._cfd.turbulence_model,
                          api.LaminarFlowModel):
            err_str = 'Turbulence models are not yet supported!'
            raise NotImplementedError(err_str)

        print '* Single phase, incompressible, isothermal, '\
              'laminar, newtonian'
        print '* Gravity: ({0:.4e}, {1:.4e}, {2:.4e})'.format(self._g[X],
                                                              self._g[Y],
                                                              self._g[Z])

    def _check_cuds_material(self, cuds):
        """Check materials in the given cuds."""
        print '---------'
        print 'Materials'
        print '---------'

        # Fluid and solid (or wall) material
        self._fluid_mat = None
        self._solid_mat = None
        self._material_cnt = cuds.count_of(CUBA.MATERIAL)

        if self._material_cnt < 1 or self._material_cnt > 2:
            err_str = 'Expecting exactly 1 or 2 materials '
            err_str += '(not {0:d})!'.format(self._material_cnt)
            raise ValueError(err_str)

        mat_ind = 1
        for imat in cuds.iter(item_type=CUBA.MATERIAL):
            mat_name = imat.name.lower()
            if mat_name == 'solid' or mat_name == 'wall':
                self._solid_mat = imat
            else:
                self._fluid_mat = imat

            print '* Material {0:d}: {1:s}'.format(mat_ind, mat_name)
            mat_ind += 1

        if self._fluid_mat is None:
            err_str = 'Expecting one material which is '
            err_str += 'not named as solid or wall!'
            raise ValueError(err_str)

        if self._material_cnt == 2 and self._solid_mat is None:
            err_str = 'With 2 materials the other must be '
            err_str += 'named as solid or wall!'
            raise ValueError(err_str)

        # Density and kinematic viscosity for the fluid material
        if CUBA.DENSITY in self._fluid_mat.data:
            self._ref_den = self._fluid_mat.data[CUBA.DENSITY]
        else:  # try
            self._ref_den = self._fluid_mat.data[CUBA.REFERENCE_DENSITY]

        self._inv_ref_den = 1.0/self._ref_den

        if CUBA.KINEMATIC_VISCOSITY in self._fluid_mat.data:
            self._kvisc = self._fluid_mat.data[CUBA.KINEMATIC_VISCOSITY]
        elif CUBA.DYNAMIC_VISCOSITY in self._fluid_mat.data:
            visc = self._fluid_mat.data[CUBA.DYNAMIC_VISCOSITY]
            self._kvisc = visc/self._ref_den
        else:  # try
            print '* Trying to interpret VISCOSITY as dynamic viscosity'

            visc = self._fluid_mat.data[CUBA.VISCOSITY]
            self._kvisc = visc/self._ref_den

        print '* Fluid density: {0:.4e}'.format(self._ref_den)
        print '* Fluid kinematic viscosity: {0:.4e}'.format(self._kvisc)

    def _check_cuds_dataset(self, cuds):
        """Check dataset in the given cuds."""
        print '-------'
        print 'Dataset'
        print '-------'

        # Cubic Bravais lattice
        cubic_lat_cnt = 0
        self._lat = None

        for ic in cuds.iter():
            if isinstance(ic, ABCLattice):
                lat_type = ic.primitive_cell.bravais_lattice
                if lat_type is BravaisLattice.CUBIC:
                    cubic_lat_cnt += 1
                    self._lat = ic

        if cubic_lat_cnt != 1:
            err_str = 'Expecting exactly one cubic Bravais lattice dataset '
            err_str += '(not {0:d})!'.format(cubic_lat_cnt)
            raise ValueError(err_str)

        self._lat_size[0] = self._lat.size[0]
        self._lat_size[1] = self._lat.size[1]
        self._lat_size[2] = self._lat.size[2]

        p1x = self._lat.primitive_cell.p1[X]
        p1y = self._lat.primitive_cell.p1[Y]
        p1z = self._lat.primitive_cell.p1[Z]

        dr = np.sqrt(p1x*p1x + p1y*p1y + p1z*p1z)
        self._dr = dr
        self._inv_dr = 1.0/self._dr

        self._V = self._lat.primitive_cell.volume
        self._inv_V = 1.0/self._V

        print '* Cubic lattice: {0:s}'.format(self._lat.name)
        print '* Lattice size: {0:}'.format(self._lat.size)
        print '* Lattice spacing: {0:.4e}'.format(self._dr)

    def _check_cuds_solver_parameter(self, cuds):
        """Check solver parameters in the given cuds."""
        print '-----------------'
        print 'Solver parameters'
        print '-----------------'

        # Time integration
        tint_cnt = cuds.count_of(CUBA.INTEGRATION_TIME)
        self._time_integr = next(cuds.iter(item_type=CUBA.INTEGRATION_TIME),
                                 None)

        if tint_cnt != 1:
            err_str = 'JYU-LB wrapper (load cuds): must configure exactly '
            err_str += 'one time integration (not {0:d})!'.format(tint_cnt)
            raise ValueError(err_str)

        self._dt = self._time_integr.size
        self._inv_dt = 1.0/self._dt
        self._cr = self._dr/self._dt
        self._inv_cr = 1.0/self._cr

        print '* Discrete time step: {0:.4e}'.format(self._dt)
        print '* Lattice reference velocity: {0:.4e}'.format(self._cr)

        # Check collision operator;
        # related to time integration from the CFD point of view,
        # i.e. here NOT part of the (kinetic) model equation specification
        coll_oper_cnt = 0
        coll_oper_str = 'bgk'
        self._coll_oper = BGK  # default

        for sp in cuds.iter(item_type=CUBA.SOLVER_PARAMETER):
            if CUBA.COLLISION_OPERATOR in sp.data:
                aux_str = sp.data[CUBA.COLLISION_OPERATOR].lower()
                coll_oper_cnt += 1

                if aux_str == 'bgk':
                    self._coll_oper = BGK
                    coll_oper_str = aux_str
                elif aux_str == 'trt':
                    self._coll_oper = TRT
                    coll_oper_str = aux_str
                elif aux_str == 'regularization':
                    self._coll_oper = REGULARIZATION
                    coll_oper_str = aux_str
                else:
                    wrn_str = 'Ignoring the unknown collision operator '
                    wrn_str += '({0:s})'.format(aux_str)
                    print wrn_str

        if coll_oper_cnt > 1:
            wrn_str = '  Several ({0:d}) '.format(coll_oper_cnt)
            wrn_str += 'collision operators specified'
            print wrn_str

        print '* Using collision operator: {0:s}'.format(coll_oper_str)

    def _check_cuds_boundary_condition(self, cuds):
        """Check boundary conditions in the given cuds."""
        print '-------------------'
        print 'Boundary conditions'
        print '-------------------'

        # Unit outward normals for
        # the exterior boundary faces of the simulation domain
        bface_uon_str = np.empty((6), dtype="S10")
        bface_uon_str[FACE_1X] = '(-1,0,0)'
        bface_uon_str[FACE_X1] = '(1,0,0)'
        bface_uon_str[FACE_1Y] = '(0,-1,0)'
        bface_uon_str[FACE_Y1] = '(0,1,0)'
        bface_uon_str[FACE_1Z] = '(0,0,-1)'
        bface_uon_str[FACE_Z1] = '(0,0,1)'

        bface_cnt = np.zeros((6), dtype=np.int32)
        bface = np.full((6), None, dtype=object)

        # Check the boundaries
        for b in cuds.iter(item_type=CUBA.BOUNDARY):
            # Check the unit outward normal (aligned with cart.coord.axes)
            if CUBA.DIRECTION in b.data:
                nx = b.data[CUBA.DIRECTION][0]
                ny = b.data[CUBA.DIRECTION][1]
                nz = b.data[CUBA.DIRECTION][2]
                nx2 = nx*nx
                ny2 = ny*ny
                nz2 = nz*nz
                n2 = nx2 + ny2 + nz2

                # Not a unit vector
                if np.fabs(n2-1.0) > 1e-06:
                    continue

                # Not aligned with cart.coord.axes
                if (np.fabs(nx2-1.0) > 1e-06 and np.fabs(ny2-1.0) > 1e-06 and
                   np.fabs(nz2-1.0) > 1e-06):
                    continue

                if nx < 0.0:
                    bface_cnt[FACE_1X] += 1
                    bface[FACE_1X] = b
                elif nx > 0.0:
                    bface_cnt[FACE_X1] += 1
                    bface[FACE_X1] = b
                elif ny < 0.0:
                    bface_cnt[FACE_1Y] += 1
                    bface[FACE_1Y] = b
                elif ny > 0.0:
                    bface_cnt[FACE_Y1] += 1
                    bface[FACE_Y1] = b
                elif nz < 0.0:
                    bface_cnt[FACE_1Z] += 1
                    bface[FACE_1Z] = b
                else:  # nz > 0.0
                    bface_cnt[FACE_Z1] += 1
                    bface[FACE_Z1] = b

        for f in range(6):
            msg_str = '* Exterior boundary face with unit outward normal '
            msg_str += '{0:s}'.format(bface_uon_str[f])
            print msg_str

            bc_str = 'Periodic'
            bc = PERIODIC

            if bface_cnt[f] == 0:
                print '  Not configured'
            elif bface_cnt[f] > 1:
                err_str = 'Found more than one ({0:d})!'.format(bface_cnt[f])
                raise ValueError(err_str)
            else:
                b = bface[f]
                bc_cnt = len(b.condition)

                if bc_cnt == 0:
                    print '  No conditions configured'
                elif bc_cnt > 1:
                    wrn_str = '  Expected 0 or 1 conditions '
                    wrn_str += '(not {0:d})'.format(bc_cnt)
                    print wrn_str

                for c in b.condition:
                    if isinstance(c, api.Periodic):
                        bc_str = 'Periodic'
                        bc = PERIODIC
                    elif isinstance(c, api.Dirichlet):
                        c_mat = c.data[CUBA.MATERIAL]
                        if c_mat != self._fluid_mat:
                            err_str = 'Dirichlet condition must be specified'
                            err_str += ' for the fluid material (not for '
                            err_str += '{0:s})!'.format(c_mat.name)
                            raise ValueError(err_str)

                        if c.data[CUBA.VARIABLE] == CUBA.PRESSURE:
                            bc = FIXED_DEN
                            bc_str = 'Dirichlet (pressure '
                            bc_str += '+ in-plane velocity components)'
                        elif c.data[CUBA.VARIABLE] == CUBA.VELOCITY:
                            bc = FIXED_VEL
                            bc_str = 'Dirichlet (velocity)'
                        else:
                            wrn_str = '  Dirichlet bc supported only for '
                            wrn_str += 'pressure or velocity (not for '
                            wrn_str += '{0:})'.format(c.data[CUBA.VARIABLE])
                            print wrn_str
                            wrn_str = '  Ignoring Dirichlet bc ('
                            wrn_str += 'name: {0:s})'.format(c.name)
                            print wrn_str
                    else:
                        wrn_str = '  Found unsupported bc ('
                        wrn_str += 'name: {0:s}, '.format(c.name)
                        wrn_str += 'cuba key: {0:}), '.format(c.cuba_key)
                        wrn_str += 'ignored'
                        print wrn_str

            self._bc[f] = bc
            print '  Using {0:s} bc'.format(bc_str)

    def _load_cuds(self):
        """Load CUDS data for the configuration of the JYU-LB engine."""
        print 'JYU-LB wrapper, load cuds'
        print '-------------------------'

        cuds = self.get_cuds()
        if cuds is None:
            raise ValueError('Unable to load cuds!')

        # Check the cuds for consistency and extract simulation parameters
        self._check_cuds_model(cuds)
        self._check_cuds_material(cuds)
        self._check_cuds_dataset(cuds)
        self._check_cuds_solver_parameter(cuds)
        self._check_cuds_boundary_condition(cuds)

        # Set up the solver
        print '-------------'
        print 'Solver set-up'
        print '-------------'
        kvisc_redu = self._inv_dr*self._inv_cr*self._kvisc
        g_redu = np.zeros((3), dtype=np.float64)
        g_redu[X] = self._dt*self._inv_cr*self._g[X]
        g_redu[Y] = self._dt*self._inv_cr*self._g[Y]
        g_redu[Z] = self._dt*self._inv_cr*self._g[Z]

        msg_str = '* Kinematic viscosity (reduced units): '
        msg_str += '{0:.4e}'.format(kvisc_redu)
        print msg_str
        msg_str = '* Gravity (reduced units): '
        msg_str += '({0:.4e}, {1:.4e}, {2:.4e})'.format(g_redu[X],
                                                        g_redu[Y],
                                                        g_redu[Z])
        print msg_str

        self._solver = Solver(self._lat_size, self._coll_oper,
                              kvisc_redu, g_redu, self._bc)

        if self._solver.error_occured():
            raise RuntimeError(self._solver.error_msg)

        # Simulation data
        self._den = self._solver.get_den()  # density field updated by solver
        self._vel = self._solver.get_vel()  # vel.field updated by solver
        self._acc = self._solver.get_acc()  # accel.field utilizied by solver
        self._phase = self._solver.get_phase()  # phase field util. by solver

        # Initialize simulation data for the lattice nodes
        # First initialize the phase information (geometry)
        fluid_node_cnt = 0
        for node in self._lat.iter():
            ijk = node.index

            if CUBA.MATERIAL in node.data:
                mat = node.data[CUBA.MATERIAL]
                if mat.uid == self._fluid_mat.uid:
                    self._phase[ijk] = FLUID
                    fluid_node_cnt += 1
                elif mat.uid == self._solid_mat.uid:
                    self._phase[ijk] = SOLID
                else:  # unidentified material
                    err_str = 'JYU-LB wrapper (load cuds): unidentified '
                    err_str += 'material (name {0:s}) '.format(mat.name)
                    err_str += 'associated with the lattice node '
                    err_str += '{0}!'.format(ijk)
                    raise ValueError(err_str)
            else:  # default material is fluid
                self._phase[ijk] = FLUID
                fluid_node_cnt += 1

        print '* Fluid lattice nodes: {0:d}'.format(fluid_node_cnt)

        self._proxy_lat = ProxyLattice(self._lat, self)
        self._proxy_lat.update(self._lat.iter())  # Then init. the flow field

        # Replace the original input lattice
        # with the newly created proxy lattice
        cuds.remove([self._lat.uid])
        cuds.add([self._proxy_lat])

        self._is_init_finalized = False
        self._lat = None  # the input dataset/lattice not needed anymore

    def run(self):
        """Run modeling engine using the configured settings."""
        # Solver initialization finalized before the first run
        if not self._is_init_finalized:
            self._lat_fnode_cnt = self._solver.finalize_init()
            self._is_init_finalized = True

        # Update particuler simulation parameters before each run
        # (bc data updated internally by the solver)
        self._update_gravity()
        gx_redu = self._dt*self._inv_cr*self._g[X]
        gy_redu = self._dt*self._inv_cr*self._g[Y]
        gz_redu = self._dt*self._inv_cr*self._g[Z]
        self._solver.gravity = (gx_redu, gy_redu, gz_redu)

        tperiod = self._time_integr.final - self._time_integr.current
        tsteps = int(tperiod/self._dt)
        tsimulated = tsteps*self._dt

        self._solver.evolve(tsteps)

        # Advance the integration times after each run
        self._time_integr.current += tsimulated
        self._time_integr.final += tsimulated

    def add_dataset(self, container):
        """Add a CUDS Lattice container.

        Parameters
        ----------
        container : {ABCLattice}
            The CUDS Lattice container to add to the engine.

        Raises
        ------
        NotImplementedError:
            Allways (no support for this feature).
        """
        err_str = 'JYU-LB wrapper (add dataset): adding datasets '
        err_str += 'directly to the wrapper not supported!'
        raise NotImplementedError(err_str)

    def remove_dataset(self, name):
        """Delete a lattice.

        Parameters
        ----------
        name : str
            name of the lattice to be deleted.

        Raises
        ------
        NotImplementedError:
            Allways (no support for this feature).
        """
        err_str = 'JYU-LB wrapper (add dataset): removing datasets '
        err_str += 'from the wrapper not supported!'
        raise NotImplementedError(err_str)

    def get_dataset(self, name):
        """ Get a lattice.

        The returned lattice can be used to query and update the state of the
        lattice inside the modeling engine.

        Parameters
        ----------
        name : str
            name of the lattice to be retrieved
            (not used, the proxy lattice is always returned).

        Returns
        -------
        ABCLattice
        """
        return self._proxy_lattice

    def get_dataset_names(self):
        """ Return the names of all datasets in the engine workspace."""
        return [self._proxy_lattice.name]

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
        NotImplementedError:
            Allways (no support for this feature; forces to use get_dataset).
        """
        err_str = 'JYU-LB wrapper (iter.dataset): iterating over datasets is'
        err_str += ' not supported (use get_dataset to access the single '
        err_str += 'dataset)!'
        raise NotImplementedError(err_str)
