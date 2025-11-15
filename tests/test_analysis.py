"""
This module contains unit tests for the analysis functions.
"""
import unittest
import os
import sys

# Add the src directory to the Python path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../src')))

from analysis import perform_gamma_analysis
from file_handlers import DicomFileHandler, MCCFileHandler

class TestGammaAnalysis(unittest.TestCase):
    """Test suite for the gamma analysis functions."""

    def setUp(self):
        """Set up the test case by loading the data files."""
        self.dicom_path = 'data/case2/1G240_2cm.dcm'
        self.mcc_path = 'data/case2/1G240_2cm.mcc'

        self.dicom_handler = DicomFileHandler()
        self.mcc_handler = MCCFileHandler()

        dicom_success, dicom_error = self.dicom_handler.open_file(self.dicom_path)
        self.assertTrue(dicom_success, f"Failed to load DICOM file: {dicom_error}")

        mcc_success, mcc_error = self.mcc_handler.open_file(self.mcc_path)
        self.assertTrue(mcc_success, f"Failed to load MCC file: {mcc_error}")

    def test_perform_gamma_analysis(self):
        """
        Test the perform_gamma_analysis function with example data.
        """
        result = perform_gamma_analysis(
            reference_handler=self.mcc_handler,
            evaluation_handler=self.dicom_handler,
            dose_percent_threshold=3,
            distance_mm_threshold=3
        )
        # Unpack all 11 values (8 original + 3 interpolated maps)
        (gamma_map, gamma_stats, phys_extent, _, dd_map, dta_map, dd_stats, dta_stats,
         gamma_map_interp, dd_map_interp, dta_map_interp) = result

        # Check that the gamma analysis produced results
        self.assertIsNotNone(gamma_map)
        self.assertIsNotNone(gamma_stats)
        self.assertIsNotNone(phys_extent)
        self.assertIsNotNone(dd_map)
        self.assertIsNotNone(dta_map)
        self.assertIsNotNone(dd_stats)
        self.assertIsNotNone(dta_stats)

        # Check the statistics
        self.assertIn('pass_rate', gamma_stats)
        self.assertIn('mean', gamma_stats)
        self.assertIn('max', gamma_stats)
        self.assertIn('min', gamma_stats)
        self.assertIn('total_points', gamma_stats)

        # Check that the pass rate is a float between 0 and 100
        self.assertIsInstance(gamma_stats['pass_rate'], float)
        self.assertGreaterEqual(gamma_stats['pass_rate'], 0)
        self.assertLessEqual(gamma_stats['pass_rate'], 100)

        # Check that the number of points is greater than 0
        self.assertGreater(gamma_stats['total_points'], 0)

        # Check the DD statistics
        self.assertIn('mean', dd_stats)
        self.assertIn('max', dd_stats)
        self.assertIn('min', dd_stats)
        self.assertIn('std', dd_stats)
        self.assertIn('total_points', dd_stats)

        # Check the DTA statistics
        self.assertIn('mean', dta_stats)
        self.assertIn('max', dta_stats)
        self.assertIn('min', dta_stats)
        self.assertIn('std', dta_stats)
        self.assertIn('total_points', dta_stats)

    def test_interpolated_analysis_maps(self):
        """
        Test that gamma, DD, and DTA maps are interpolated to DICOM grid.
        This ensures the maps have no gaps when visualized.
        """
        result = perform_gamma_analysis(
            reference_handler=self.mcc_handler,
            evaluation_handler=self.dicom_handler,
            dose_percent_threshold=3,
            distance_mm_threshold=3
        )

        # Unpack all return values including interpolated maps
        (gamma_map, gamma_stats, phys_extent, mcc_interp_data,
         dd_map, dta_map, dd_stats, dta_stats,
         gamma_map_interp, dd_map_interp, dta_map_interp) = result

        # Get DICOM grid shape
        dicom_data = self.dicom_handler.get_pixel_data()
        expected_shape = dicom_data.shape

        # Check that interpolated maps are returned
        self.assertIsNotNone(gamma_map_interp, "Gamma map interpolated should be returned")
        self.assertIsNotNone(dd_map_interp, "DD map interpolated should be returned")
        self.assertIsNotNone(dta_map_interp, "DTA map interpolated should be returned")

        # Check that interpolated maps have the same shape as DICOM grid
        self.assertEqual(gamma_map_interp.shape, expected_shape,
                        f"Gamma map interpolated should have shape {expected_shape}")
        self.assertEqual(dd_map_interp.shape, expected_shape,
                        f"DD map interpolated should have shape {expected_shape}")
        self.assertEqual(dta_map_interp.shape, expected_shape,
                        f"DTA map interpolated should have shape {expected_shape}")

        # Check that interpolated maps contain valid data (not all NaN or zeros)
        import numpy as np
        self.assertTrue(np.any(~np.isnan(gamma_map_interp)),
                       "Gamma map interpolated should contain non-NaN values")
        self.assertTrue(np.any(~np.isnan(dd_map_interp)),
                       "DD map interpolated should contain non-NaN values")
        self.assertTrue(np.any(~np.isnan(dta_map_interp)),
                       "DTA map interpolated should contain non-NaN values")

        # Check that interpolated maps have more non-NaN values than sparse maps
        # (since they fill in the gaps)
        gamma_sparse_valid = np.sum(~np.isnan(gamma_map))
        gamma_interp_valid = np.sum(~np.isnan(gamma_map_interp))
        self.assertGreater(gamma_interp_valid, gamma_sparse_valid,
                          "Interpolated gamma map should have more valid values than sparse map")


if __name__ == '__main__':
    unittest.main()
