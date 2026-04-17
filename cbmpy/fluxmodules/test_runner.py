"""Test module to verify CBMPy can be imported and basic functionality works."""
import unittest


class TestCBMPyConfig(unittest.TestCase):
    """Test CBMPy configuration and imports."""

    def test_import_cbmpy(self):
        """Test that cbmpy can be imported."""
        import cbmpy
        self.assertTrue(hasattr(cbmpy, 'CBModel'))

    def test_version(self):
        """Test that version is available."""
        from cbmpy.CBConfig import current_version
        version = current_version()
        self.assertTrue(len(version) > 0)

    def test_config(self):
        """Test that configuration works."""
        from cbmpy.CBConfig import __CBCONFIG__
        self.assertIn('VERSION', __CBCONFIG__)
        self.assertIn('SOLVER_PREF', __CBCONFIG__)


if __name__ == '__main__':
    unittest.main()
