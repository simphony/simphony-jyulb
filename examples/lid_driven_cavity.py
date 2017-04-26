"""Lid-driven cavity flow (2D setup).

Details
-------
- 3D simulation (2D setup, i.e. periodic in one direction)
- flow driven by a moving wall
"""
import time
import numpy as np
from simphony.cuds.meta import api
from simphony.core.cuba import CUBA
from simphony.api import CUDS, Simulation
from simphony.engine import EngineInterface
from simphony.cuds.lattice import make_cubic_lattice

# CUDS
cuds = CUDS(name='lid-driven cavity')

# Physics model
cfd = api.Cfd(name='fluid flow')

# Submodels (in fact, these are the default submodels in CFD)
cfd.thermal_model = api.IsothermalModel(name='isothermal')
cfd.turbulence_model = api.LaminarFlowModel(name='laminar')
cfd.rheology_model = api.NewtonianFluidModel(name='newtonian')
cfd.multiphase_model = api.SinglePhaseModel(name='singlephase')
cfd.compressibility_model = api.IncompressibleFluidModel(name='incompressible')
cuds.add([cfd])

# Materials
fluid = api.Material(name='fluid')
solid = api.Material(name='wall')
fluid.data[CUBA.DENSITY] = 1.0
fluid.data[CUBA.KINEMATIC_VISCOSITY] = 1.0/120.0
cuds.add([fluid, solid])

# Dataset (lattice)
H = 100  # Distance between the lower (not moving) and upper (moving) wall
W = 100  # Distance between the side walls (not moving)
lat_size = (W+2, 1, H+2)  # number of lattice nodes
lat_h = 1.0               # lattice spacing
lat = make_cubic_lattice('cavity', lat_h, lat_size)

# Create a cavity geometry, set the initial flow field
# The moving wall is located at the upper exterior boundary face (z-dir.)
# No walls in the "dummy" direction (y-dir.)
mwall_vel = 5.0e-2
for node in lat.iter():
    ijk = node.index

    if ijk[2] == (lat_size[2]-1):
        node.data[CUBA.MATERIAL] = solid
        node.data[CUBA.VELOCITY] = (mwall_vel, 0, 0)
    elif ijk[0] == 0 or ijk[0] == (lat_size[0]-1) or ijk[2] == 0:
        node.data[CUBA.MATERIAL] = solid
        node.data[CUBA.VELOCITY] = (0, 0, 0)
    else:
        node.data[CUBA.MATERIAL] = fluid
        node.data[CUBA.VELOCITY] = (0, 0, 0)

    lat.update([node])

cuds.add([lat])

# Reynolds number
Re = H*mwall_vel/fluid.data[CUBA.KINEMATIC_VISCOSITY]

# Boundary conditions
# for the exterior boundary faces of the simulation domain
# (the directions refer to unit outward normals)
face1y = api.Boundary(name='face1y', condition=[api.Periodic()])
face1y.data[CUBA.DIRECTION] = (0, -1, 0)

facey1 = api.Boundary(name='facey1', condition=[api.Periodic()])
facey1.data[CUBA.DIRECTION] = (0, 1, 0)

face1z = api.Boundary(name='face1z', condition=[api.Periodic()])
face1z.data[CUBA.DIRECTION] = (0, 0, -1)

facez1 = api.Boundary(name='facez1', condition=[api.Periodic()])
facez1.data[CUBA.DIRECTION] = (0, 0, 1)

cuds.add([face1y, facey1, face1z, facez1])

# Solver parameters (time integration, collision operator)
ti = api.IntegrationTime(name='simulation time',
                         current=0.0, final=5000.0, size=1.0)

ti.data[CUBA.COLLISION_OPERATOR] = 'Regularization'
cuds.add([ti])

# Simulate
print '-'*77
sim = Simulation(cuds, 'JYU-LB', engine_interface=EngineInterface.Internal)
lat = cuds.get_by_name(lat.name)

# Total fluid mass
tot_fmass = 0.0
for node in lat.iter():
    if node.data[CUBA.MATERIAL] == fluid:
        tot_fmass += node.data[CUBA.DENSITY]

print '-'*77
print 'Lid-Driven Cavity flow simulation (Re = {0:f})'.format(Re)
print '-'*77
print "Total fluid mass before simulation: {0:.4e}".format(tot_fmass)
print '-'*77

start = time.time()

for e in range(25):
    sim.run()

    # Total velocity
    tot_vx = 0.0
    tot_vy = 0.0
    tot_vz = 0.0
    for node in lat.iter():
        if node.data[CUBA.MATERIAL] == fluid:
            tot_vx += node.data[CUBA.VELOCITY][0]
            tot_vy += node.data[CUBA.VELOCITY][1]
            tot_vz += node.data[CUBA.VELOCITY][2]

    sf = 'Time = {0:f}, tot.vel = ({1:11.4e}, {2:11.4e}, {3:11.4e})'
    print sf.format(ti.current, tot_vx, tot_vy, tot_vz)

end = time.time()

# Total fluid mass
tot_fmass = 0.0
for node in lat.iter():
    if node.data[CUBA.MATERIAL] == fluid:
        tot_fmass += node.data[CUBA.DENSITY]

print '-'*77
print "Total fluid mass after simulation: {0:.4e}".format(tot_fmass)
print '-'*77
print "Time spend in run (s): ", end-start
print '-'*77

# ---------------------------------------------------------------------------
# Post-processing: output
# ---------------------------------------------------------------------------
# Matlab: 2D flow field
# ---------------------------------------------------------------------------
coord = np.zeros(3, dtype=np.int32)
coord[0] = int(lat_size[0]/2)
coord[1] = int(lat_size[1]/2)
coord[2] = int(lat_size[2]/2)

with open('flow_2D.field', 'w') as f:
    for x1 in range(lat_size[0]):
        for x2 in range(lat_size[2]):
            coord[0] = x1
            coord[2] = x2

            node = lat.get(coord)
            if node.data[CUBA.MATERIAL] == solid:
                continue

            p = node.data[CUBA.PRESSURE]
            dn = node.data[CUBA.DENSITY]
            v1 = node.data[CUBA.VELOCITY][0]
            v2 = node.data[CUBA.VELOCITY][2]

            sf = '{0:d} {1:d} {2:e} {3:e} {4:e}\n'
            f.write(sf.format(x1, x2, dn, v1, v2))

# ---------------------------------------------------------------------------
