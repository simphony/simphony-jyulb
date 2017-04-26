"""Testing module for a file-io based wrapper for JYU-LB modeling engine."""
import time
import os
import tempfile
import shutil
import unittest
from simphony.engine import jyulb_fileio_isothermal as lb
from jyulb.fileio.common.proxy_lattice import ProxyLattice
from jyulb.testing.jyulb_check_engine import JYULBEngineCheck
from simphony.testing.abc_check_engine import LatticeEngineCheck


class JYULBEngineTestCase(JYULBEngineCheck, unittest.TestCase):
    """Test case for JYULBEngine class."""

    def setUp(self):
        LatticeEngineCheck.setUp(self, number_datasets_used_in_testing=1)

        self.fluid_enum = ProxyLattice.FLUID_ENUM
        self.solid_enum = ProxyLattice.SOLID_ENUM
        self.temp_dir = tempfile.mkdtemp()
        self.saved_path = os.getcwd()
        os.chdir(self.temp_dir)
        self.addCleanup(self.cleanup)

    def cleanup(self):
        os.chdir(self.saved_path)
        shutil.rmtree(self.temp_dir)

    def engine_factory(self):
        return lb.JYULBEngine()

    def test_run_engine(self):
        """Running the JYU-LB modeling engine."""
        engine = lb.JYULBEngine()

        # Set engine dependent parameters
        self.coll_oper = lb.JYULBEngine.TRT_ENUM
        self.flow_type = lb.JYULBEngine.STOKES_FLOW_ENUM

        # Set other problem parameters
        self._setup_test_problem(engine)

        # Run the case and measure performance
        start_time = time.time()
        engine.run()
        end_time = time.time()
        comp_time = end_time - start_time
        MFLUP = (self.nx-2)*self.ny*self.nz*self.tsteps/1e6
        print 'Comp.time (s) = {}, MFLUPS = {}'.format(comp_time,
                                                       MFLUP/comp_time)

        # Analyse the results
        self._analyse_test_problem_results(engine)


if __name__ == '__main__':
    unittest.main()
