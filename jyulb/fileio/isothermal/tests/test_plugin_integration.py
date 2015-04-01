import unittest


class TestPluginIntegration(unittest.TestCase):

    """Test case for JYUEngine class."""

    def test_plugin_integration(self):
        """Test to run JYU-LB modeling engine."""
        # Assert that we can import the jyulb plugin
        from simphony.engine import jyulb

        # Check that the expected top level objects are available
        self.assertTrue(hasattr(jyulb, 'jyu_engine'))


if __name__ == '__main__':
    unittest.main()
