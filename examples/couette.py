"""Coeutte fluid flow simulation (flow between two moving, parallel plates).

Details
-------
- 3D simulation (with two "dummy" direction)
- flow driven by moving walls
- periodic bc in the "dummy" directions
"""
import time
import numpy as np
from simphony.cuds.meta import api
from simphony.core.cuba import CUBA
from simphony.api import CUDS, Simulation
from simphony.engine import EngineInterface
from simphony.cuds.lattice import make_cubic_lattice
from jyulb.flow_field import couette_plate_steady_state

# CUDS
cuds = CUDS(name='couette')

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
fluid.data[CUBA.KINEMATIC_VISCOSITY] = 1.0/6.0
cuds.add([fluid, solid])

# Dataset (lattice)
lat_size = (10, 1, 10)  # number of lattice nodes
lat_h = 1.0             # lattice spacing
lat = make_cubic_lattice('channel', lat_h, lat_size)

# Create a channel geometry, set the initial flow field
lower_wall_vel1 = -1.0e-3
lower_wall_vel2 = 0.0e-3
upper_wall_vel1 = 1.0e-3
upper_wall_vel2 = -2.0e-3

for node in lat.iter():
    ijk = node.index
    # The two plates are separated in the z-direction
    if ijk[2] == 0:
        node.data[CUBA.MATERIAL] = solid
        node.data[CUBA.VELOCITY] = (lower_wall_vel1, lower_wall_vel2, 0)
    elif ijk[2] == (lat_size[2]-1):
        node.data[CUBA.MATERIAL] = solid
        node.data[CUBA.VELOCITY] = (upper_wall_vel1, upper_wall_vel2, 0)
    else:
        node.data[CUBA.MATERIAL] = fluid
        node.data[CUBA.VELOCITY] = (0, 0, 0)

    lat.update([node])

cuds.add([lat])

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
                         current=0.0, final=1000.0, size=1.0)

ti.data[CUBA.COLLISION_OPERATOR] = 'TRT'
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

print '='*77
print 'Couette flow simulation'
print '='*77
print "Total fluid mass before simulation: {0:.4e}".format(tot_fmass)
print '-'*77

start = time.time()

for e in range(10):
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
# Post-processing:
#   Analytical solution for the Couette flow, error, and output
# ---------------------------------------------------------------------------
# GLE: flow speed in y-data format
# ---------------------------------------------------------------------------
H = lat_size[2]-2  # plates separated in the z-direction
centrel_dist = np.zeros((H), dtype=np.float64)
ana_vel = np.zeros((H, 2), dtype=np.float64)

couette_plate_steady_state(H, lower_wall_vel1, lower_wall_vel2,
                           upper_wall_vel1, upper_wall_vel2,
                           centrel_dist, ana_vel)

sim_loc_vel = np.zeros(3, dtype=np.float64)
coord = np.zeros(3, dtype=np.int32)
coord[0] = int(lat_size[0]/2)
coord[1] = int(lat_size[1]/2)
coord[2] = int(lat_size[2]/2)

ana2 = 0.0
sim_ana_diff2 = 0.0
sim_loc_vel = np.zeros((3), dtype=np.float64)

with open('flow_prof.txt', 'w') as f:
    for h in range(H):
        coord[2] = h+1
        node = lat.get(coord)
        dist = centrel_dist[h]

        if node.data[CUBA.MATERIAL] == solid:
            continue

        ana_vel1 = ana_vel[h, 0]
        ana_vel2 = ana_vel[h, 1]

        sim_loc_p = node.data[CUBA.PRESSURE]
        sim_loc_den = node.data[CUBA.DENSITY]
        sim_loc_vel[0] = node.data[CUBA.VELOCITY][0]
        sim_loc_vel[1] = node.data[CUBA.VELOCITY][1]
        sim_loc_vel[2] = node.data[CUBA.VELOCITY][2]

        sim_loc_speed = np.sqrt(sim_loc_vel[0]*sim_loc_vel[0] +
                                sim_loc_vel[1]*sim_loc_vel[1] +
                                sim_loc_vel[2]*sim_loc_vel[2])

        ana_loc_speed = np.sqrt(ana_vel1*ana_vel1 + ana_vel2*ana_vel2)

        sim_ana_diff = sim_loc_speed - ana_loc_speed
        sim_ana_diff2 += sim_ana_diff*sim_ana_diff
        ana2 += ana_loc_speed*ana_loc_speed

        f.write('{0:f} {1:e} {2:e} {3:e} {4:e} {5:e} {6:e}\n'.format(dist,
                sim_loc_den, sim_loc_vel[0], sim_loc_vel[1], sim_loc_vel[2],
                sim_loc_speed, ana_loc_speed))

rel_l2_err_norm_vel = np.sqrt(sim_ana_diff2/ana2)

print 'Relative L2-error norm: {0:.4e}'.format(rel_l2_err_norm_vel)
print '-'*77

# ---------------------------------------------------------------------------
# Matlab: 2D flow field
# ---------------------------------------------------------------------------
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
            v2 = node.data[CUBA.VELOCITY][1]

            sf = '{0:d} {1:d} {2:e} {3:e} {4:e}\n'
            f.write(sf.format(x1, x2, dn, v1, v2))

# ---------------------------------------------------------------------------
