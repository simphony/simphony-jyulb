"""Testing module for the JYU-LB wrapper."""
import os
import shutil
import tempfile
import unittest
import numpy as np

from simphony.cuds.meta import api
from simphony.core.cuba import CUBA
from simphony.api import CUDS, Simulation
from simphony.engine import EngineInterface
from simphony.cuds.lattice import make_cubic_lattice
from jyulb.flow_field import poise_plate_steady_state_pgrad


class JYULBEngineTestCaseSimulate(unittest.TestCase):

    """Simulation test for checking that installation is ok."""

    def setUp(self):
        """Initialize the test case."""
        self.temp_dir = tempfile.mkdtemp()
        self.saved_path = os.getcwd()

        os.chdir(self.temp_dir)
        self.addCleanup(self.cleanup)

    def cleanup(self):
        """Clean files after the test case."""
        os.chdir(self.saved_path)
        shutil.rmtree(self.temp_dir)

    def cuds_factory(self):
        """Create CUDS for poiseuille fluid flow simulation."""
        # CUDS
        cuds = CUDS(name='poiseuille')

        # Physics model
        cfd = api.Cfd(name='fluid flow')

        # Submodels (in fact, these are the default submodels in CFD)
        cfd.thermal_model = api.IsothermalModel(name='isothermal')
        cfd.turbulence_model = api.LaminarFlowModel(name='laminar')
        cfd.rheology_model = api.NewtonianFluidModel(name='newtonian')
        cfd.multiphase_model = api.SinglePhaseModel(name='singlephase')
        compr_model = api.IncompressibleFluidModel(name='incompressible')
        cfd.compressibility_model = compr_model

        # Gravity in the z-direction
        cfd.gravity_model.acceleration = (0, 0, 1e-5)
        cuds.add([cfd])

        # Materials
        fluid = api.Material(name='fluid')
        solid = api.Material(name='wall')
        fluid.data[CUBA.DENSITY] = 1.0
        fluid.data[CUBA.KINEMATIC_VISCOSITY] = 2.0/6.0
        cuds.add([fluid, solid])

        # Dataset (lattice)
        lat_size = (5, 1, 1)  # number of lattice nodes
        lat_h = 1.0             # lattice spacing
        lat = make_cubic_lattice('channel', lat_h, lat_size)

        # Create a channel geometry, set the initial flow field
        for node in lat.iter():
            ijk = node.index
            # The two plates are separated in the x-direction
            if ijk[0] == 0 or ijk[0] == (lat_size[0]-1):
                node.data[CUBA.MATERIAL] = solid
                node.data[CUBA.VELOCITY] = (0, 0, 0)
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

        self._fluid = fluid
        self._solid = solid
        self._cuds = cuds
        self._cfd = cfd
        self._lat = lat
        self._ti = ti

        return cuds

    def analyse_results(self):
        """Compute the relative L2-error norm."""
        H = self._lat.size[0]-2  # plates separated in the x-direction

        den = self._fluid.data[CUBA.DENSITY]
        dvisc = den*self._fluid.data[CUBA.KINEMATIC_VISCOSITY]
        eff_pgrad = -den*self._cfd.gravity_model.acceleration[2]

        centrel_dist = np.zeros((H), dtype=np.float64)
        ana_vel = np.zeros((H), dtype=np.float64)

        poise_plate_steady_state_pgrad(H, dvisc, eff_pgrad,
                                       centrel_dist, ana_vel)

        sim_loc_vel = np.zeros(3, dtype=np.float64)
        coord = np.zeros(3, dtype=np.int32)
        coord[0] = int(self._lat.size[0]/2)
        coord[1] = int(self._lat.size[1]/2)
        coord[2] = int(self._lat.size[2]/2)

        ana2 = 0.0
        sim_ana_diff2 = 0.0
        sim_loc_vel = np.zeros((3), dtype=np.float64)

        for h in range(H):
            coord[0] = h+1
            node = self._lat.get(coord)

            if node.data[CUBA.MATERIAL] == self._solid:
                continue

            sim_loc_vel[0] = node.data[CUBA.VELOCITY][0]
            sim_loc_vel[1] = node.data[CUBA.VELOCITY][1]
            sim_loc_vel[2] = node.data[CUBA.VELOCITY][2]

            sim_loc_speed = np.sqrt(sim_loc_vel[0]*sim_loc_vel[0] +
                                    sim_loc_vel[1]*sim_loc_vel[1] +
                                    sim_loc_vel[2]*sim_loc_vel[2])

            ana_loc_vel = ana_vel[h]

            sim_ana_diff = sim_loc_speed - ana_loc_vel
            sim_ana_diff2 += sim_ana_diff*sim_ana_diff
            ana2 += ana_loc_vel*ana_loc_vel

        rel_l2_err_norm_vel = np.sqrt(sim_ana_diff2/ana2)

        print 'Relative L2-error norm: {0:.4e}'.format(rel_l2_err_norm_vel)
        print '-'*77

        return rel_l2_err_norm_vel

    def test_run_engine(self):
        """Running the JYU-LB engine."""
        cuds = self.cuds_factory()

        sim = Simulation(cuds, 'JYU-LB',
                         engine_interface=EngineInterface.Internal)

        lat = cuds.get_by_name(self._lat.name)
        self._lat = lat

        # Total fluid mass before simulation
        tot_fmass = 0.0
        for node in lat.iter():
            if node.data[CUBA.MATERIAL] == self._fluid:
                tot_fmass += node.data[CUBA.DENSITY]

        print '='*77
        print 'Poiseuille flow simulation'
        print '='*77
        print "Total fluid mass before simulation: {0:.4e}".format(tot_fmass)
        print '-'*77

        # Simulate
        sim.run()

        # Total velocity after simulation
        tot_vx = 0.0
        tot_vy = 0.0
        tot_vz = 0.0
        for node in lat.iter():
            if node.data[CUBA.MATERIAL] == self._fluid:
                tot_vx += node.data[CUBA.VELOCITY][0]
                tot_vy += node.data[CUBA.VELOCITY][1]
                tot_vz += node.data[CUBA.VELOCITY][2]

        sf = 'Time = {0:f}, tot.vel = ({1:11.4e}, {2:11.4e}, {3:11.4e})'
        print sf.format(self._ti.current, tot_vx, tot_vy, tot_vz)

        # Total fluid mass after simulation
        tot_fmass = 0.0
        for node in lat.iter():
            if node.data[CUBA.MATERIAL] == self._fluid:
                tot_fmass += node.data[CUBA.DENSITY]

        print '-'*77
        print "Total fluid mass after simulation: {0:.4e}".format(tot_fmass)
        print '-'*77

        rel_l2_error = self.analyse_results()

        self.assertTrue(rel_l2_error < 1.0e-10)


if __name__ == '__main__':
    unittest.main()
