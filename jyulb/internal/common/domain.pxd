import cython
import numpy as np
cimport numpy as np


cdef extern from "const.h":
    cdef unsigned char SOLID_NODE # constant
    cdef unsigned char FLUID_NODE # constant

    
cdef extern from "node.h":
    cdef cppclass NodeSet:
        unsigned int get_n(const unsigned int ijk[3]) except +
        void get_ijk(unsigned int n, unsigned int ijk[3]) except +

        unsigned int get_node_count() const

    cdef cppclass Lattice(NodeSet):
        Lattice(const unsigned int size[3], const double origin[3]) except +

        void get_size(unsigned int size[3]) const
        void get_origin(double origin[3]) const


cdef extern from "data.h":
    cdef cppclass GeomData:
        GeomData(NodeSet *nodes, unsigned char ivalue) except +

        void set_val_n(unsigned int n, unsigned char value) except +
        void set_val_ijk(const unsigned int ijk[3], unsigned char value) except +

        unsigned char get_val_n(unsigned int n) except +
        unsigned char get_val_ijk(const unsigned int ijk[3]) except +
        
    cdef cppclass FieldData:
        FieldData(NodeSet *nodes, double ivalue) except +

        void set_val_n(unsigned int n, double value) except +
        void set_val_ijk(const unsigned int ijk[3], double value) except +

        double get_val_n(unsigned int n) except +
        double get_val_ijk(const unsigned int ijk[3]) except +
        
    cdef cppclass Geometry:
        Geometry(Lattice *lat)

        Lattice *get_lattice() const 
        GeomData *get_phase() const

    cdef cppclass IsothermalNodeData:
        IsothermalNodeData(NodeSet *nodes, double iden, double ivel[3], double ifrc[3])

        FieldData *den() const
        FieldData *velx() const
        FieldData *vely() const
        FieldData *velz() const
        FieldData *frcx() const
        FieldData *frcy() const
        FieldData *frcz() const

        NodeSet *get_nodeset() const
        
        
cdef class PyNodeSet:
    cdef NodeSet *thisptr
    
    
cdef class PyLattice(PyNodeSet):
    cdef unsigned int size[3]
    cdef double origin[3]
    cdef Lattice *_lat_ptr


cdef class PyGeometry:
    cdef Geometry *thisptr
    cdef PyLattice _lattice

    
cdef class PyAbstractIsothermalData:
    cdef IsothermalNodeData *thisptr
#    cdef PyNodeSet _nodes

#cdef class PyIsothermalData(PyAbstractIsothermalData):
      
