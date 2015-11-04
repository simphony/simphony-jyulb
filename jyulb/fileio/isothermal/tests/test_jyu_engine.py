"""Testing module for a file-io based wrapper for JYU-LB modeling engine."""
import math
import os
import tempfile
import shutil
import unittest
from simphony.core.cuba import CUBA
from jyulb.cuba_extension import CUBAExtension
from simphony.cuds.lattice import make_cubic_lattice
from simphony.engine import jyulb_fileio_isothermal as lb
from simphony.cuds.abc_lattice import ABCLattice
from simphony.cuds.abc_modeling_engine import ABCModelingEngine
from jyulb.fileio.common.jyu_lattice_proxy import JYULatticeProxy
from simphony.testing.abc_check_engine import LatticeEngineCheck


class JYUEngineTestCase(LatticeEngineCheck, unittest.TestCase):

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

    def engine_factory(self):
        return lb.JYUEngine()

    def create_dataset(self, name):
        """ This method is overriden, because JyuLB requires that certain
        CUBA keys are always defined on Proxy Lattice objects.

        """
        lat = make_cubic_lattice(name, 1.0, (2, 3, 4))
        for node in lat.iter_nodes():
            node.data = {CUBA.MATERIAL_ID: 0, CUBA.DENSITY: 0,
                         CUBA.VELOCITY: [0, 0, 0], CUBA.FORCE: [0, 0, 0]}
            lat.update_nodes([node])
        return lat

    def create_dataset_items(self):
        """ Not applicable to JYU-LB
        """
        pass

    def check_instance_of_dataset(self, ds):
        self.assertIsInstance(ds, ABCLattice,
                              "Error: Dataset must be ABCLattice!")

    def test_delete_dataset(self):
        """ JYU-LB does not support adding multiple datasets, therefore
        this test is overriden
        """
        engine = self.engine_factory()
        engine.add_dataset(self.create_dataset("test"))
        for ds in engine.iter_datasets():
            engine.remove_dataset(ds.name)
            with self.assertRaises(ValueError):
                engine.get_dataset("test")

    def test_dataset_rename(self):
        """ JYU-LB does not support adding multiple datasets, therefore
        this test is overriden
        """
        engine = self.engine_factory()
        engine.add_dataset(self.create_dataset(name='foo'))
        ds = engine.get_dataset("foo")
        ds.name = "bar"
        self.assertEqual(ds.name, "bar")

        # we should not be able to use the old name "foo"
        with self.assertRaises(ValueError):
            engine.get_dataset("foo")
        with self.assertRaises(ValueError):
            engine.remove_dataset("foo")
        with self.assertRaises(ValueError):
            [_ for _ in engine.iter_datasets(names=["foo"])]

        # we should be able to access using the new "bar" name
        ds_bar = engine.get_dataset("bar")
        self.assertEqual(ds_bar.name, "bar")

        # and we should be able to use the no-longer used
        # "foo" name when adding another dataset
        # remove the other dataset first
        engine.remove_dataset("bar")
        ds = engine.add_dataset(self.create_dataset(name='foo'))

    def test_add_dataset_data_copy(self):
        """ JYU-LB does not support multiple datasets, therefore
        this test is overridden
        """
        pass

    def test_get_dataset_names(self):
        """ JYU-LB does not support multiple datasets, therefore
        this test is overridden
        """

        engine = self.engine_factory()
        # add a few empty datasets
        ds_names = ["test"]

        engine.add_dataset(self.create_dataset("test"))

        # test that we are getting all the names
        names = [
            n for n in engine.get_dataset_names()]
        self.assertEqual(names, ds_names)

    def test_iter_dataset(self):
        """ JYU-LB does not support multiple datasets, therefore
        this test is overridden
        """
        engine = self.engine_factory()

        ds_names = ["test"]
        engine.add_dataset(self.create_dataset("test"))

        # test iterating over all
        names = [ds.name for ds in engine.iter_datasets()]
        self.assertEqual(names, ds_names)

        for ds in engine.iter_datasets(ds_names):
            self.check_instance_of_dataset(ds)


    def test_run_engine(self):
        """Running the jyu-lb modeling engine."""
        engine = lb.JYUEngine()

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
                node.data[CUBA.MATERIAL_ID] = JYULatticeProxy.SOLID_ENUM
            else:
                node.data[CUBA.MATERIAL_ID] = JYULatticeProxy.FLUID_ENUM
            lat.update_nodes([node])

        # Initialize flow variables at fluid lattice nodes
        for node in lat.iter_nodes():
            if node.data[CUBA.MATERIAL_ID] == JYULatticeProxy.FLUID_ENUM:
                node.data[CUBA.VELOCITY] = (0, 0, 0)
                node.data[CUBA.DENSITY] = 1.0
            lat.update_nodes([node])

        # Add lattice to the engine
        engine.add_dataset(lat)

        # Run the case
        engine.run()

        # Analyse the results
        proxy_lat = engine.get_dataset(lat.name)

        # Compute the relative L2-error norm
        tot_diff2 = 0.0
        tot_ana2 = 0.0
        tot_ux = 0.0
        tot_uy = 0.0
        for node in proxy_lat.iter_nodes():
            if node.data[CUBA.MATERIAL_ID] == JYULatticeProxy.FLUID_ENUM:
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
        print ('Relative L2-error norm = %e\n' % (rel_l2_error))

        self.assertTrue(rel_l2_error < 1.0e-10)
        self.assertTrue(math.fabs(tot_ux) < 1.0e-10)
        self.assertTrue(math.fabs(tot_uy) < 1.0e-10)

    def _calc_poiseuille_vel(self, index):
        wall_dist = (float(index-1) + 0.5)
        centerl = (wall_dist) - self.channel_h
        d = (centerl/self.channel_h)*(centerl/self.channel_h)
        return self.max_vel*(1.0 - d)

if __name__ == '__main__':
    unittest.main()
