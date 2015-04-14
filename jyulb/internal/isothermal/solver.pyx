import cython
import numpy as np
cimport numpy as cnp
from jyulb.internal.common.domain cimport PyAbstractIsothermalData


STOKES_FLOW_ENUM = STOKES_FLOW
LAMINAR_FLOW_ENUM = LAMINAR_FLOW
TURBULENT_FLOW_ENUM = TURBULENT_FLOW

BGK_ENUM = BGK
TRT_ENUM = TRT
MRT_ENUM = MRT
REG_ENUM = REG


cdef class PyFlowParams:

    def __cinit__(self):
        self.vals.dr = 1.0
        self.vals.dt = 1.0
        self.vals.ref_den = 1.0
        self.vals.kvisc = 1.0/6.0
        self.vals.gx = 0.0
        self.vals.gy = 0.0
        self.vals.gz = 0.0
        self.vals.flow_type = STOKES_FLOW_ENUM
        self.vals.collision_operator = TRT_ENUM
        self.vals.external_forcing = False

    property lattice_spacing:
        def __get__(self):
            return self.vals.dr
        def __set__(self, double dr):
            self.vals.dr = dr

    property time_step:
        def __get__(self):
            return self.vals.dt
        def __set__(self, double dt):
            self.vals.dt = dt

    property reference_density:
        def __get__(self):
            return self.vals.ref_den
        def __set__(self, double ref_den):
            self.vals.ref_den = ref_den
        
    property kinematic_viscosity:
        def __get__(self):
            return self.vals.kvisc
        def __set__(self, double kvisc):
            self.vals.kvisc = kvisc

    property gravity:
        def __get__(self):
            cdef cnp.ndarray[double, ndim=1, mode="c"] grav = np.empty(3, dtype=np.float64)
            grav[0] = self.vals.gx
            grav[1] = self.vals.gy
            grav[2] = self.vals.gz
            return grav
        def __set__(self, cnp.ndarray[double, ndim=1, mode="c"] grav):
            self.vals.gx = grav[0] 
            self.vals.gy = grav[1]
            self.vals.gz = grav[2]

    property flow_type:
        def __get__(self):
            return self.vals.flow_type
        def __set__(self, unsigned char flow_type):
            self.vals.flow_type = flow_type

    property collision_operator:
        def __get__(self):
            return self.vals.collision_operator
        def __set__(self, unsigned char coll_oper):
            self.vals.collision_operator = coll_oper

    property external_forcing:
        def __get__(self):
            return self.vals.external_forcing
        def __set__(self, bool ext_frc):
            self.vals.external_forcing = ext_frc


cdef class PyIsothermalSolverData(PyAbstractIsothermalData):

    def __cinit__(self, PySolver solver):
        self.thisptr = solver.thisptr.get_fluid_node_data()

        
cdef class PySolver:

    def __cinit__(self, PyGeometry geom, PyFlowParams params):
        self.thisptr = new IsothermalSolver(geom.thisptr, &params.vals)

    def __dealloc__(self):
        del self.thisptr

    def init_field_data(self):
        self.thisptr.init_field_data()
        
    def evolve(self, unsigned int time_steps):
        self.thisptr.evolve(time_steps)

    def get_field_data(self):
        return PyIsothermalSolverData(self)
        