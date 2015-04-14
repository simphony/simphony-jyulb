"""Testing module for an internal wrapper for JYU-LB modeling engine."""
import time
import math
import os
import tempfile
import shutil
import unittest
from simphony.core.cuba import CUBA
from jyulb.cuba_extension import CUBAExtension
from simphony.cuds.lattice import make_cubic_lattice
from simphony.engine import jyulb_internal_isothermal as lb
from jyulb.internal.common.proxy_lattice import ProxyLattice
from simphony.cuds.abc_modeling_engine import ABCModelingEngine


class JYUEngineTestCase(unittest.TestCase):

    """Test case for JYUEngine class."""

    def setUp(self):
        self.dr = 1.0
        self.nx = 5
        self.ny = 3
        self.nz = 4

        self.coll_oper = lb.JYUEngine.TRT_ENUM
        self.dt = 1.0
        self.tsteps = 1000

        self.gx = 0.0
        self.gy = 0.0
        self.gz = 1.0e-5
        self.kvisc = 0.1
        self.rden = 1.0
        self.flow_type = lb.JYUEngine.STOKES_FLOW_ENUM
        self.ext_frc = False

        self.channel_h = 0.5*(self.nx-2.0)
        self.max_vel = 0.5*self.gz*self.channel_h*self.channel_h/self.kvisc

        self.temp_dir = tempfile.mkdtemp()
        self.saved_path = os.getcwd()
        os.chdir(self.temp_dir)
        self.addCleanup(self.cleanup)

    def cleanup(self):
        os.chdir(self.saved_path)
        shutil.rmtree(self.temp_dir)

    def test_run_engine(self):
        """Running the jyu-lb modeling engine."""
        engine = lb.JYUEngine()

        self.assertIsInstance(engine, ABCModelingEngine,
                              "Error: not a ABCModelingEngine!")

        # Computational Method data
        engine.CM[CUBAExtension.COLLISION_OPERATOR] = self.coll_oper
        engine.CM[CUBA.TIME_STEP] = self.dt
        engine.CM[CUBA.NUMBER_OF_TIME_STEPS] = self.tsteps

        # System Parameters data
        engine.SP[CUBAExtension.REFERENCE_DENSITY] = self.rden
        engine.SP[CUBA.KINEMATIC_VISCOSITY] = self.kvisc
        engine.SP[CUBAExtension.GRAVITY] = (self.gx, self.gy, self.gz)
        engine.SP[CUBAExtension.FLOW_TYPE] = self.flow_type
        engine.SP[CUBAExtension.EXTERNAL_FORCING] = self.ext_frc

        # Boundary Conditions data
        engine.BC[CUBA.VELOCITY] = {'open': 'periodic',
                                    'wall': 'noSlip'}

        engine.BC[CUBA.DENSITY] = {'open': 'periodic',
                                   'wall': 'noFlux'}

        # Configure a lattice
        lat = make_cubic_lattice("lattice1", self.dr,
                                 (self.nx, self.ny, self.nz))

        # Set geometry for a Poiseuille channel
        for node in lat.iter_nodes():
            if node.index[0] == 0 or node.index[0] == self.nx-1:
                node.data[CUBA.MATERIAL_ID] = ProxyLattice.SOLID_ENUM
            else:
                node.data[CUBA.MATERIAL_ID] = ProxyLattice.FLUID_ENUM
            lat.update_node(node)

        # Initialize flow variables at fluid lattice nodes
        for node in lat.iter_nodes():
            if node.data[CUBA.MATERIAL_ID] == ProxyLattice.FLUID_ENUM:
                node.data[CUBA.VELOCITY] = (0, 0, 0)
                node.data[CUBA.DENSITY] = 1.0
            lat.update_node(node)

        # Add lattice to the engine
        engine.add_lattice(lat)

        # Run the case
        start_time = time.time()

        engine.run()

        end_time = time.time()
        comp_time = end_time - start_time
        MFLUP = (self.nx-2)*self.ny*self.nz*self.tsteps/1e6

        # Analyse the results
        proxy_lat = engine.get_lattice(lat.name)

        # Compute the relative L2-error norm
        tot_diff2 = 0.0
        tot_ana2 = 0.0
        tot_ux = 0.0
        tot_uy = 0.0
        for node in proxy_lat.iter_nodes():
            if node.data[CUBA.MATERIAL_ID] == ProxyLattice.FLUID_ENUM:
                sim_ux = node.data[CUBA.VELOCITY][0]
                sim_uy = node.data[CUBA.VELOCITY][1]
                sim_uz = node.data[CUBA.VELOCITY][2]
                ana_uz = self._calc_poiseuille_vel(node.index[0])
                diff = ana_uz - sim_uz
                tot_diff2 = tot_diff2 + diff*diff
                tot_ana2 = tot_ana2 + ana_uz*ana_uz
                tot_ux = tot_ux + sim_ux
                tot_uy = tot_uy + sim_uy

        rel_l2_error = math.sqrt(tot_diff2/tot_ana2)
        print ('\nRelative L2-error norm = %e' % (rel_l2_error))
        print 'Comp.time (s) = {}, MFLUPS = {}'.format(comp_time,
                                                       MFLUP/comp_time)

        self.assertTrue(rel_l2_error < 1.0e-10)
        self.assertTrue(math.fabs(tot_ux) < 1.0e-10)
        self.assertTrue(math.fabs(tot_uy) < 1.0e-10)

        # Test iteration and removal of lattices
        for lat in engine.iter_lattices():
            self.assertEqual(lat, proxy_lat)

        engine.delete_lattice(proxy_lat.name)
        none_lat = engine.get_lattice(proxy_lat.name)

        self.assertEqual(none_lat, None)

    def _calc_poiseuille_vel(self, index):
        wall_dist = (float(index-1) + 0.5)
        centerl = (wall_dist) - self.channel_h
        d = (centerl/self.channel_h)*(centerl/self.channel_h)
        return self.max_vel*(1.0 - d)

if __name__ == '__main__':
    unittest.main()
