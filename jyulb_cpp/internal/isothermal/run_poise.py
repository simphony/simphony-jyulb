"""A script for simulating simple Poiseuille flow
with the isothermal JYU-LB modeling engine."""
import time
import math
import numpy as np
from simphony.core.cuba import CUBA
from jyulb.cuba_extension import CUBAExtension
from simphony.cuds.lattice import make_cubic_lattice
from simphony.engine import jyulb_internal_isothermal as lb
from jyulb.internal.common.proxy_lattice import ProxyLattice

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import ticker

# Simulation parameters
dr = 1.0
nx = 20
ny = 5
nz = 10

coll_oper = lb.JYULBEngine.TRT_ENUM
dt = 1.0
run_tsteps = 500
number_of_runs = 20

gx = 0.0
gy = 0.0
gz = 1.0e-5
kvisc = 0.1
rden = 1.0
flow_type = lb.JYULBEngine.STOKES_FLOW_ENUM

# Auxiliarry variables for comparing the simulated flow profile
# against the analytical expression
channel_h = 0.5*(nx-2.0)
max_vel = 0.5*gz*channel_h*channel_h/kvisc

sim_uzs = np.zeros(nx, dtype=np.float64)
ana_uzs = np.zeros(nx, dtype=np.float64)

# Analytical expression for the velocity in a Poiseuille flow
def calc_poiseuille_vel(index):
    wall_dist = (float(index-1) + 0.5)
    centerl = (wall_dist) - channel_h
    d = (centerl/channel_h)*(centerl/channel_h)
    return max_vel*(1.0 - d)

# Calculate statistics from an evolving flow field
def calc_stats(tstep, lat):
    tot_diff2 = 0.0
    tot_ana2 = 0.0
    tot_den = 0.0
    tot_ux = 0.0
    tot_uy = 0.0
    tot_uz = 0.0

    for node in proxy_lat.iter_nodes():
        if node.data[CUBA.MATERIAL_ID] == ProxyLattice.FLUID_ENUM:
            sim_den = node.data[CUBA.DENSITY]
            sim_ux = node.data[CUBA.VELOCITY][0]
            sim_uy = node.data[CUBA.VELOCITY][1]
            sim_uz = node.data[CUBA.VELOCITY][2]
            ana_uz = calc_poiseuille_vel(node.index[0])

            diff = ana_uz - sim_uz
            tot_diff2 = tot_diff2 + diff*diff
            tot_ana2 = tot_ana2 + ana_uz*ana_uz

            tot_den = tot_den + sim_den
            tot_ux = tot_ux + sim_ux
            tot_uy = tot_uy + sim_uy
            tot_uz = tot_uz + sim_uz

            if node.index[1] == 1 and node.index[2] == 1:
                sim_uzs[node.index[0]] = sim_uz
                ana_uzs[node.index[0]] = ana_uz

    rel_l2_error = math.sqrt(tot_diff2/tot_ana2)
    print 'Tstep = {}, tot.mass = {}, tot.vel. = ({}, {}, {}), rel.L2-error = {}'.format(tstep, tot_den, tot_ux, tot_uy, tot_uz, rel_l2_error)        

# Create a wrapper for the JYU-LB solver
engine = lb.JYULBEngine()

# Computational Method data
engine.CM[CUBAExtension.COLLISION_OPERATOR] = coll_oper
engine.CM[CUBA.TIME_STEP] = dt
engine.CM[CUBA.NUMBER_OF_TIME_STEPS] = run_tsteps

# System Parameters data
engine.SP[CUBAExtension.REFERENCE_DENSITY] = rden
engine.SP[CUBA.KINEMATIC_VISCOSITY] = kvisc
engine.SP[CUBAExtension.GRAVITY] = (gx, gy, gz)
engine.SP[CUBAExtension.FLOW_TYPE] = flow_type

# Boundary Conditions data
engine.BC[CUBA.VELOCITY] = {'open': 'periodic',
                            'wall': 'noSlip'}

engine.BC[CUBA.DENSITY] = {'open': 'periodic',
                           'wall': 'noFlux'}

# Configure a lattice
lat = make_cubic_lattice("lattice1", dr, (nx, ny, nz))

# Set geometry (flow between two plates)
# and initialize the flow field
for node in lat.iter_nodes():
    if node.index[0] > 0 and node.index[0] < nx-1:
        node.data[CUBA.MATERIAL_ID] = ProxyLattice.FLUID_ENUM
        node.data[CUBA.VELOCITY] = (0, 0, 0)
        node.data[CUBA.DENSITY] = rden
    else:
        node.data[CUBA.MATERIAL_ID] = ProxyLattice.SOLID_ENUM
    lat.update_nodes([node])

# Add lattice to the engine
engine.add_dataset(lat)
proxy_lat = engine.get_dataset(lat.name)

print '====================================================================='
print 'Poiseuille channel flow'
print '====================================================================='
calc_stats(0, proxy_lat)

# Run the case (and measure the execution time)
start_time = time.time()

for r in np.arange(1,number_of_runs+1):
    engine.run()
    calc_stats(r*run_tsteps, proxy_lat)

end_time = time.time()
comp_time = end_time - start_time
MFLUP = (nx-2)*ny*nz*run_tsteps*number_of_runs/1e6

print '====================================================================='
print 'Comp.time (s) = {}, MFLUPS = {}'.format(comp_time, MFLUP/comp_time)
print '====================================================================='

# Plot the flow profile
xvals = np.arange(1, nx-1, dtype=np.int)
fig, ax = plt.subplots()

sim_line, = plt.plot(xvals, sim_uzs[1:-1], color='blue', markeredgewidth=1, marker='o', markersize=15, fillstyle=u'none', label='Simulated')
ana_line, = plt.plot(xvals, ana_uzs[1:-1], color='red', marker='*', markersize=12, fillstyle=u'full', label='Analytical')

formatter = ticker.ScalarFormatter(useMathText=True)
formatter.set_scientific(True) 
formatter.set_powerlimits((-1,1)) 

ax.yaxis.set_major_formatter(formatter) 
#plt.axis([0.85, 8.15, 1.5e-4, 8.1e-4])

plt.legend(handles=[sim_line, ana_line])
plt.legend(loc=10)

plt.title('Poiseuille channel flow')
plt.ylabel('Velocity in the z-direction')
plt.xlabel('x-coordinate')
plt.show()
plt.savefig('comp_vel_prof.png')