import domain
import time
import numpy as np

print 'Solid material id = {}'.format(domain.SOLID_ENUM)
print 'Fluid material id = {}'.format(domain.FLUID_ENUM)

lat1 = domain.PyLattice(np.array((4, 2, 3), dtype=np.uint32),
                        np.array((0.0, 1.0, 2.0), dtype=np.float64))

size1 = np.zeros(3, dtype=np.uint32)
origin1 = np.zeros(3, dtype=np.float64)

lat1.get_size(size1)
lat1.get_origin(origin1)
ncount1 = lat1.get_node_count()

print type(size1), size1, ncount1
print type(origin1), origin1 
                        
lat2 = domain.PyLattice.fromlattice(lat1)

size2 = np.zeros(3, dtype=np.uint32)
origin2 = np.zeros(3, dtype=np.float64)

lat2.get_size(size2)
lat2.get_origin(origin2)
ncount2 = lat2.get_node_count()

print type(size2), size2, ncount2
print type(origin2), origin2
                        
for index in np.ndindex(tuple(size1)):
    print index, lat1.get_n(np.array(index, dtype=np.uint32))
                        
for n in np.arange(ncount1):
    ijk = np.zeros(3, dtype=np.uint32)
    lat1.get_ijk(n, ijk)
    print n, ijk
            
geom1 = domain.PyGeometry(lat1)
            
for index in np.ndindex(tuple(size1)):
    if index[0] > 0 and index[0] < size1[0]-1:
        geom1.set_material_ijk(np.array(index, dtype=np.uint32), domain.FLUID_ENUM)

fluid_ncount = 0
solid_ncount = 0        
for index in np.ndindex(tuple(size1)):
    material = geom1.get_material_ijk(np.array(index, dtype=np.uint32))
    print 'Node {} material id = {}'.format(index, material)

    if material == domain.SOLID_ENUM:
        solid_ncount = solid_ncount + 1
    else:
        fluid_ncount = fluid_ncount + 1

print 'Solid node count = {}'.format(solid_ncount)
print 'Fluid node count = {}'.format(fluid_ncount)

ivel = np.array((0.0, 1.0, 2.0), dtype=np.float64)
ifrc = np.array((0.0, 0.0, 0.0), dtype=np.float64)
field_data = domain.PyIsothermalData(lat1, 1.0, ivel, ifrc)

#data_count = field_data.get_node_set().get_node_count()
data_count = field_data.get_data_count()
print 'Isothermal data count = {}'.format(data_count)

vel = np.zeros(3, dtype=np.float64)
frc = np.ones(3, dtype=np.float64)
ijk = np.zeros(3, dtype=np.uint32)

for n in np.arange(data_count):
#    field_data.get_node_set().get_ijk(n, ijk)
    field_data.get_ijk(n, ijk)

    vel[0] = n*1.0
    vel[1] = 0.0
    vel[2] = 0.0
    
    frc[0] = 0.0
    frc[1] = 0.0
    frc[2] = n*1.0

    field_data.set_vel_ijk(ijk, vel)
    field_data.set_frc_ijk(ijk, frc)
    
    den = field_data.get_den_n(n)
    field_data.get_vel_n(n, vel)
    field_data.get_frc_n(n, frc)
    
    print n, ijk, den, vel, frc

