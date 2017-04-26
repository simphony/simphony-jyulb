"""Testing module for a file-io based wrapper for JYU-LB modeling engine."""
import math
import numpy as np
from simphony.core.cuba import CUBA
from simphony.cuds.lattice import make_cubic_lattice
from jyulb.fileio.isothermal3D.jyu_engine import JYUEngine

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import ticker

dr = 1.0
nx = 10
ny = 5
nz = 20

coll_oper = JYUEngine.TRT_ENUM
dt = 1.0
tsteps = 10000

gx = 0.0
gy = 0.0
gz = 1.0e-5
kvisc = 0.1
rden = 1.0
flow_type = JYUEngine.STOKES_FLOW_ENUM
ext_frc = False

channel_h = 0.5*(nx-2.0)
max_vel = 0.5*gz*channel_h*channel_h/kvisc

def calc_poiseuille_vel(index):
    wall_dist = (float(index-1) + 0.5)
    centerl = (wall_dist) - channel_h
    d = (centerl/channel_h)*(centerl/channel_h)
    return max_vel*(1.0 - d)

engine = JYUEngine()

# Computational Method data
engine.CM_CUBA_COLLISION_OPERATOR = coll_oper
engine.CM[CUBA.TIME_STEP] = dt
engine.CM[CUBA.NUMBER_OF_TIME_STEPS] = tsteps

# System Parameters data
engine.SP_CUBA_REFERENCE_DENSITY = rden
engine.SP[CUBA.KINEMATIC_VISCOSITY] = kvisc
engine.SP_CUBA_GRAVITY = (gx, gy, gz)
engine.SP_CUBA_FLOW_TYPE = flow_type
engine.SP_CUBA_EXTERNAL_FORCING = ext_frc

# Boundary Conditions data
engine.BC[CUBA.VELOCITY] = {'open': 'periodic',
                            'wall': 'noSlip'}

engine.BC[CUBA.DENSITY] = {'open': 'periodic',
                           'wall': 'noFlux'}

# Configure a lattice
lat = make_cubic_lattice("lattice1", dr, (nx, ny, nz))

# Set geometry for a Poiseuille channel
for node in lat.iter_nodes():
    if node.index[0] == 0 or node.index[0] == nx-1:
        node.data[CUBA.MATERIAL_ID] = engine.SOLID_ENUM
    else:
        node.data[CUBA.MATERIAL_ID] = engine.FLUID_ENUM
    lat.update_node(node)

# Initialize flow variables at fluid lattice nodes
for node in lat.iter_nodes():
    if node.data[CUBA.MATERIAL_ID] == engine.FLUID_ENUM:
        node.data[CUBA.VELOCITY] = (0, 0, 0)
        node.data[CUBA.DENSITY] = 1.0
    lat.update_node(node)

# Add lattice to the engine
engine.add_lattice(lat)

# Run the case
engine.run()

# Analyse the results
proxy_lat = engine.get_lattice(lat.name)

# Compute the relative L2-error norm
tot_diff2 = 0.0
tot_ana2 = 0.0
tot_ux = 0.0
tot_uy = 0.0

sim_uzs = np.zeros(nx, dtype=np.float64)
ana_uzs = np.zeros(nx, dtype=np.float64)

for node in proxy_lat.iter_nodes():
    if node.data[CUBA.MATERIAL_ID] == engine.FLUID_ENUM:
        sim_ux = node.data[CUBA.VELOCITY][0]
        sim_uy = node.data[CUBA.VELOCITY][1]
        sim_uz = node.data[CUBA.VELOCITY][2]
        ana_uz = calc_poiseuille_vel(node.index[0])
        diff = ana_uz - sim_uz
        tot_diff2 = tot_diff2 + diff*diff
        tot_ana2 = tot_ana2 + ana_uz*ana_uz
        tot_ux = tot_ux + sim_ux
        tot_uy = tot_uy + sim_uy

        if node.index[1] == 1 and node.index[2] == 1:
            sim_uzs[node.index[0]] = sim_uz
            ana_uzs[node.index[0]] = ana_uz
        
rel_l2_error = math.sqrt(tot_diff2/tot_ana2)
print ('Relative L2-error norm = %e\n' % (rel_l2_error))

xvals = np.arange(1, nx-1, dtype=np.int)

#mpl.rcParams.update({'font.size': 12, 'text.usetex': True})
fig, ax = plt.subplots()

sim_line, = plt.plot(xvals, sim_uzs[1:-1], color='blue', markeredgewidth=1, marker='o', markersize=15, fillstyle=u'none', label='Simulated')
ana_line, = plt.plot(xvals, ana_uzs[1:-1], color='red', marker='*', markersize=12, fillstyle=u'full', label='Analytical')

formatter = ticker.ScalarFormatter(useMathText=True)
formatter.set_scientific(True) 
formatter.set_powerlimits((-1,1)) 

ax.yaxis.set_major_formatter(formatter) 
plt.axis([0.85, 8.15, 1.5e-4, 8.1e-4])

plt.legend(handles=[sim_line, ana_line])
plt.legend(loc=10)

plt.title('Poiseuille channel flow')
plt.ylabel('Velocity in the z-direction')
plt.xlabel('x-coordinate')
plt.show()
plt.savefig('comp_vel_prof.png')