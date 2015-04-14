import cython
import numpy as np
cimport numpy as cnp
from libcpp cimport bool
from jyulb.internal.common.domain cimport Geometry, IsothermalNodeData, PyGeometry


cdef extern from "const.h":
    cdef unsigned char STOKES_FLOW
    cdef unsigned char LAMINAR_FLOW
    cdef unsigned char TURBULENT_FLOW

    cdef unsigned char BGK
    cdef unsigned char TRT
    cdef unsigned char MRT
    cdef unsigned char REG


cdef extern from "kernel.h":
    cdef cppclass IsothermalFlowParams:
        double dr, dt, ref_den, kvisc, gx, gy, gz
        unsigned char flow_type, collision_operator
        bool external_forcing


cdef extern from "solver.h":
    cdef cppclass IsothermalSolver:
        IsothermalSolver(Geometry *geom, IsothermalFlowParams *params)
      
        set_kvisc(double kvisc)
        set_gravity(double gx, double gy, double gz)
        void init_field_data() except +
        void evolve(unsigned int tsteps) except +

        IsothermalNodeData *get_fluid_node_data() const
        
        
cdef class PyFlowParams:
    cdef IsothermalFlowParams vals

    
cdef class PySolver:
    cdef IsothermalSolver *thisptr
