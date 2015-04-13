import cython
import numpy as np
cimport numpy as np
from domain cimport PyAbstractIsothermalData


SOLID_ENUM = SOLID_NODE # material id for solid nodes
FLUID_ENUM = FLUID_NODE # material id for fluid nodes


cdef class PyNodeSet:

    def get_n(self, np.ndarray[np.uint32_t, ndim=1, mode="c"] ijk):
        cdef int ijk_size
        ijk_size = ijk.shape[0]

        if ijk_size < 3:
            raise RuntimeError('Invalid array size: {}!'.format(ijk_size))

        return self.thisptr.get_n(&(ijk[0]))

    def get_ijk(self, unsigned int n, np.ndarray[np.uint32_t, ndim=1, mode="c"] ijk):
        cdef int ijk_size
        ijk_size = ijk.shape[0]

        if ijk_size < 3:
            raise RuntimeError('Invalid array size: {}!'.format(ijk_size))

        self.thisptr.get_ijk(n, &(ijk[0]))

    def get_node_count(self):
        return self.thisptr.get_node_count()

        
cdef class PyLattice(PyNodeSet):

    def __cinit__(self, np.ndarray[np.uint32_t, ndim=1, mode="c"] size,
                  np.ndarray[double, ndim=1, mode="c"] origin):
        self.origin[0] = origin[0]
        self.origin[1] = origin[1]
        self.origin[2] = origin[2]

        self.size[0] = size[0]
        self.size[1] = size[1]
        self.size[2] = size[2]

        self._lat_ptr = new Lattice(&(size[0]), &(origin[0]))
        self.thisptr = self._lat_ptr
        
    def __dealloc__(self):
        del self.thisptr

    def get_size(self, np.ndarray[np.uint32_t, ndim=1, mode="c"] size):
        size[0] = self.size[0]
        size[1] = self.size[1]
        size[2] = self.size[2]
        
    def get_origin(self, np.ndarray[double, ndim=1, mode="c"] origin):
        origin[0] = self.origin[0]
        origin[1] = self.origin[1]
        origin[2] = self.origin[2]

    @classmethod
    def fromlattice(cls, PyLattice lat):
        lat_size = np.zeros(3, dtype=np.uint32)
        lat_origin = np.zeros(3, dtype=np.float64)

        lat.get_size(lat_size)
        lat.get_origin(lat_origin)

        return cls(lat_size, lat_origin)

        
cdef class PyGeometry:

    def __cinit__(self, PyLattice lattice):
        self.thisptr = new Geometry(lattice._lat_ptr)
        self._lattice = lattice

    def __dealloc__(self):
        del self.thisptr

    def set_material_n(self, unsigned int n, unsigned char value):
        self.thisptr.get_phase().set_val_n(n, value)
        
    def set_material_ijk(self, np.ndarray[np.uint32_t, ndim=1, mode="c"] ijk, unsigned char value):
        cdef unsigned int n
        n = self._lattice.get_n(ijk)
        self.thisptr.get_phase().set_val_n(n, value)

    def get_material_n(self, unsigned int n):
        return self.thisptr.get_phase().get_val_n(n)

    def get_material_ijk(self, np.ndarray[np.uint32_t, ndim=1, mode="c"] ijk):
        cdef unsigned int n
        n = self._lattice.get_n(ijk)
        return self.thisptr.get_phase().get_val_n(n)

    def get_lattice(self):
        return self._lattice


cdef class PyAbstractIsothermalData:

    def set_den_n(self, unsigned int n, double den):
        self.thisptr.den().set_val_n(n, den)

    def set_vel_n(self, unsigned int n, np.ndarray[double, ndim=1, mode="c"] vel):
        self.thisptr.velx().set_val_n(n, vel[0])
        self.thisptr.vely().set_val_n(n, vel[1])
        self.thisptr.velz().set_val_n(n, vel[2])

    def set_frc_n(self, unsigned int n, np.ndarray[double, ndim=1, mode="c"] frc):
        self.thisptr.frcx().set_val_n(n, frc[0])
        self.thisptr.frcy().set_val_n(n, frc[1])
        self.thisptr.frcz().set_val_n(n, frc[2])

    def set_den_ijk(self, np.ndarray[np.uint32_t, ndim=1, mode="c"] ijk, double den):
        cdef unsigned int n
        n = self.thisptr.get_nodeset().get_n(&(ijk[0]))
#        n = self._nodes.get_n(ijk)
        self.set_den_n(n, den)

    def set_vel_ijk(self, np.ndarray[np.uint32_t, ndim=1, mode="c"] ijk,
                    np.ndarray[double, ndim=1, mode="c"] vel):
        cdef unsigned int n
        n = self.thisptr.get_nodeset().get_n(&(ijk[0]))
#        n = self._nodes.get_n(ijk)
        self.set_vel_n(n, vel)

    def set_frc_ijk(self, np.ndarray[np.uint32_t, ndim=1, mode="c"] ijk,
                    np.ndarray[double, ndim=1, mode="c"] frc):
        cdef unsigned int n
        n = self.thisptr.get_nodeset().get_n(&(ijk[0]))
#        n = self._nodes.get_n(ijk)
        self.set_frc_n(n, frc)

    def get_den_n(self, unsigned int n):
        return self.thisptr.den().get_val_n(n)

    def get_vel_n(self, unsigned int n, np.ndarray[double, ndim=1, mode="c"] vel):
        vel[0] = self.thisptr.velx().get_val_n(n)
        vel[1] = self.thisptr.vely().get_val_n(n)
        vel[2] = self.thisptr.velz().get_val_n(n)

    def get_frc_n(self, unsigned int n, np.ndarray[double, ndim=1, mode="c"] frc):
        frc[0] = self.thisptr.frcx().get_val_n(n)
        frc[1] = self.thisptr.frcy().get_val_n(n)
        frc[2] = self.thisptr.frcz().get_val_n(n)

    def get_den_ijk(self, np.ndarray[np.uint32_t, ndim=1, mode="c"] ijk):
        cdef unsigned int n
        n = self.thisptr.get_nodeset().get_n(&(ijk[0]))
#        n = self._nodes.get_n(ijk)
        return self.get_den_n(n)

    def get_vel_ijk(self, np.ndarray[np.uint32_t, ndim=1, mode="c"] ijk, np.ndarray[double, ndim=1, mode="c"] vel):
        cdef unsigned int n
        n = self.thisptr.get_nodeset().get_n(&(ijk[0]))
#        n = self._nodes.get_n(ijk)
        self.get_vel_n(n, vel)
        
    def get_frc_ijk(self, np.ndarray[np.uint32_t, ndim=1, mode="c"] ijk, np.ndarray[double, ndim=1, mode="c"] frc):
        cdef unsigned int n
        n = self.thisptr.get_nodeset().get_n(&(ijk[0]))
#        n = self._nodes.get_n(ijk)
        self.get_frc_n(n, frc)

#    def get_node_set(self):
#        return self._nodes

    def get_n(self, np.ndarray[np.uint32_t, ndim=1, mode="c"] ijk):
        return self.thisptr.get_nodeset().get_n(&(ijk[0]))

    def get_ijk(self, unsigned int n, np.ndarray[np.uint32_t, ndim=1, mode="c"] ijk):
        self.thisptr.get_nodeset().get_ijk(n, &(ijk[0]))

    def get_data_count(self):
        return self.thisptr.get_nodeset().get_node_count()
        
cdef class PyIsothermalData(PyAbstractIsothermalData):

    def __cinit__(self, PyNodeSet nodes, double iden, np.ndarray[double, ndim=1, mode="c"] ivel,
                  np.ndarray[double, ndim=1, mode="c"] ifrc):
        cdef int ivel_n, ifrc_n
        ivel_n, ifrc_n = ivel.shape[0], ifrc.shape[0]

        if ivel_n < 3 or ifrc_n < 3:
            raise RuntimeError('Invalid array size: {} {}!'.format(ivel_n, ifrc_n))

        self.thisptr = new IsothermalNodeData(nodes.thisptr, iden, &(ivel[0]), &(ifrc[0]))
#        self._nodes = nodes
        
    def __dealloc__(self):
        del self.thisptr
