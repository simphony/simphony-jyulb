import unittest


class TestPluginIntegration(unittest.TestCase):

    """Test case for JYUEngine class."""

    def test_plugin_integration(self):
        """Test to run JYU-LB modeling engine."""
        # Assert that we can import the jyulb plugin
        from simphony.engine import jyulb_internal_isothermal as lb

        # Check that the expected top level objects are available
        self.assertTrue(hasattr(lb, 'JYUEngine'))


if __name__ == '__main__':
    unittest.main()
