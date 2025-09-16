
"""
This module provides functions for analyzing and comparing DICOM and MCC data,
including profile extraction and gamma analysis.
"""
import numpy as np
from scipy.interpolate import interpn
from utils import logger, find_nearest_index

def extract_profile_data(direction, fixed_position, dicom_handler, mcc_handler=None):
    """
    Extracts profile data from DICOM and MCC datasets.

    Args:
        direction (str): "vertical" or "horizontal".
        fixed_position (float): The fixed position in mm.
        dicom_handler: Handler object containing DICOM image and coordinate information.
        mcc_handler: Handler object containing MCC image and coordinate information (optional).

    Returns:
        dict: A dictionary containing the profile data.
    """
    profile_data = {'type': direction, 'fixed_pos': fixed_position}
    
    dicom_image = dicom_handler.get_pixel_data()
    if dicom_image is None:
        logger.warning("DICOM pixel data not available for profile extraction.")
        return profile_data
    
    try:
        # Define axis-dependent variables based on direction
        if direction == "vertical":
            # Vertical profile: x is fixed, y is the profile axis
            fixed_axis_coords_dicom = dicom_handler.phys_x_mesh[0, :]
            profile_axis_mesh_dicom = dicom_handler.phys_y_mesh
            slicer_dicom = lambda idx: (slice(None), idx)
            if mcc_handler:
                mcc_fixed_axis_coords = mcc_handler.phys_x_mesh[0, :]
                mcc_profile_axis_mesh = mcc_handler.phys_y_mesh
                slicer_mcc = lambda idx: (slice(None), idx)
                sort_required = True  # MCC y-axis is inverted and needs sorting for interpolation
        else:  # "horizontal"
            # Horizontal profile: y is fixed, x is the profile axis
            fixed_axis_coords_dicom = dicom_handler.phys_y_mesh[:, 0]
            profile_axis_mesh_dicom = dicom_handler.phys_x_mesh
            slicer_dicom = lambda idx: (idx, slice(None))
            if mcc_handler:
                mcc_fixed_axis_coords = mcc_handler.phys_y_mesh[:, 0]
                mcc_profile_axis_mesh = mcc_handler.phys_x_mesh
                slicer_mcc = lambda idx: (idx, slice(None))
                sort_required = False # MCC x-axis is already sorted

        # --- Common Logic for both directions ---

        # 1. Extract DICOM profile
        closest_idx_dicom = find_nearest_index(fixed_axis_coords_dicom, fixed_position)
        profile_coords_dicom = profile_axis_mesh_dicom[slicer_dicom(closest_idx_dicom)]
        dicom_values = dicom_image[slicer_dicom(closest_idx_dicom)]
        
        profile_data['phys_coords'] = profile_coords_dicom
        profile_data['dicom_values'] = dicom_values

        # 2. Process MCC data if available
        if mcc_handler and mcc_handler.get_matrix_data() is not None:
            mcc_image = mcc_handler.get_matrix_data()
            
            closest_idx_mcc = find_nearest_index(mcc_fixed_axis_coords, fixed_position)
            mcc_line_values = mcc_image[slicer_mcc(closest_idx_mcc)]
            
            valid_indices = np.where(mcc_line_values >= 0)[0]
            
            if len(valid_indices) > 0:
                # Get physical coordinates and dose values for valid MCC points
                if direction == "vertical":
                    mcc_phys_coords = mcc_profile_axis_mesh[valid_indices, closest_idx_mcc]
                else: # horizontal
                    mcc_phys_coords = mcc_profile_axis_mesh[closest_idx_mcc, valid_indices]
                
                mcc_values = mcc_line_values[valid_indices]
                
                # Find corresponding DICOM values at MCC measurement points
                dicom_at_mcc_positions = np.array([
                    dicom_values[find_nearest_index(profile_coords_dicom, pos)] for pos in mcc_phys_coords
                ])
                
                # Interpolate MCC data for smooth plotting
                if len(mcc_values) > 1:
                    interp_coords = mcc_phys_coords
                    interp_values = mcc_values
                    if sort_required:
                        sort_indices = np.argsort(interp_coords)
                        interp_coords = interp_coords[sort_indices]
                        interp_values = interp_values[sort_indices]
                    
                    mcc_interp = np.interp(
                        profile_coords_dicom,
                        interp_coords,
                        interp_values,
                        left=np.nan, right=np.nan
                    )
                    profile_data['mcc_interp'] = mcc_interp

                profile_data['mcc_phys_coords'] = mcc_phys_coords
                profile_data['mcc_values'] = mcc_values
                profile_data['dicom_at_mcc'] = dicom_at_mcc_positions

        return profile_data
        
    except Exception as e:
        logger.error(f"프로파일 데이터 추출 오류: {str(e)}", exc_info=True)
        return profile_data

import pymedphys
from scipy.interpolate import griddata

def perform_gamma_analysis(reference_handler, evaluation_handler,
                           dose_percent_threshold, distance_mm_threshold,
                           global_normalisation=True, threshold=10, max_gamma=None):
    """
    Performs gamma analysis using pymedphys.gamma.

    The reference data (MCC) is sparse, and the evaluation data (DICOM) is a grid.
    This function uses the sparse MCC measurement points, filtered by a dose threshold,
    directly as the reference for the gamma calculation against the DICOM dose grid.

    Args:
        reference_handler: Handler for the reference data (MCC).
        evaluation_handler: Handler for the evaluation data (DICOM).
        dose_percent_threshold (float): Dose difference threshold in percent.
        distance_mm_threshold (float): Distance to agreement threshold in mm.
        global_normalisation (bool, optional): Whether to use global normalization. Defaults to True.
        threshold (int, optional): Lower dose cutoff threshold as a percentage of the normalization dose. Defaults to 10.
        max_gamma (float, optional): Maximum gamma value to consider. Defaults to None.

    Returns:
        tuple: A tuple containing the gamma map for display, gamma statistics, physical extent,
               and the interpolated MCC dose grid for visualization.
    """
    try:
        # --- Step 1: Extract and filter reference data (MCC) ---
        mcc_dose_data = reference_handler.get_matrix_data()
        if mcc_dose_data is None:
            raise ValueError("MCC data not found in reference_handler.")

        # Get all valid (>=0) measurement points and their coordinates
        all_valid_indices = np.where(mcc_dose_data >= 0)
        if all_valid_indices[0].size == 0:
            raise ValueError("No valid measurement data (>= 0) in MCC file.")

        # Vectorized calculation of physical coordinates from pixel coordinates for clarity and robustness
        y_pix, x_pix = all_valid_indices
        handler = reference_handler
        full_grid_px = x_pix + handler.crop_pixel_offset[0]
        full_grid_py = y_pix + handler.crop_pixel_offset[1]
        phys_x_all = (full_grid_px - handler.mcc_origin_x) * handler.mcc_spacing_x
        phys_y_all = -(full_grid_py - handler.mcc_origin_y) * handler.mcc_spacing_y
        # Reconstruct the (N, 2) array for griddata, which expects points as pairs
        all_mcc_coords_phys = np.vstack((phys_x_all, phys_y_all)).T

        all_mcc_dose_values = mcc_dose_data[all_valid_indices]

        # Determine normalization dose for thresholding
        norm_dose = np.max(all_mcc_dose_values)
        if norm_dose == 0:
            raise ValueError("Cannot determine normalization dose (max reference dose is zero).")

        # Filter points based on the lower dose cutoff threshold
        threshold_dose = (threshold / 100.0) * norm_dose
        analysis_mask = all_mcc_dose_values >= threshold_dose

        # --- Step 2: Extract evaluation data (DICOM) ---
        dicom_dose_grid = evaluation_handler.get_pixel_data()
        if dicom_dose_grid is None:
            raise ValueError("DICOM data not found in evaluation_handler.")

        dicom_phys_x_mesh, dicom_phys_y_mesh = evaluation_handler.phys_x_mesh, evaluation_handler.phys_y_mesh
        axes_eval_grid = (dicom_phys_y_mesh[:, 0], dicom_phys_x_mesh[0, :])
        phys_extent = reference_handler.get_physical_extent()

        # If no points are left after filtering, return early.
        if not np.any(analysis_mask):
            logger.warning(f"No MCC data points above the {threshold}% dose threshold ({threshold_dose:.2f} Gy). Gamma analysis will be skipped.")
            gamma_stats = {'pass_rate': 100, 'mean': 0, 'max': 0, 'min': 0, 'total_points': 0}
            gamma_map_for_display = np.full_like(mcc_dose_data, np.nan)
            # Still create interpolated data for the report
            mcc_interp_data = griddata(
                all_mcc_coords_phys, all_mcc_dose_values,
                (dicom_phys_x_mesh, dicom_phys_y_mesh),
                method='linear', fill_value=0
            )
            return gamma_map_for_display, gamma_stats, phys_extent, mcc_interp_data

        # These are the points that will be used in the gamma calculation
        mcc_dose_for_gamma = all_mcc_dose_values[analysis_mask]
        y_coords_for_gamma = phys_y_all[analysis_mask]
        x_coords_for_gamma = phys_x_all[analysis_mask]

        # Get the original pixel indices of the analyzed points to place results back into the map
        original_indices_for_gamma = (all_valid_indices[0][analysis_mask], all_valid_indices[1][analysis_mask])

        # --- Step 3: Perform gamma analysis using sparse reference points ---
        # The points have already been filtered, so set lower_percent_dose_cutoff to 0.
        gamma_values = pymedphys.gamma(
            axes_reference=(y_coords_for_gamma, x_coords_for_gamma),
            dose_reference=mcc_dose_for_gamma,
            axes_evaluation=axes_eval_grid,
            dose_evaluation=dicom_dose_grid,
            dose_percent_threshold=dose_percent_threshold,
            distance_mm_threshold=distance_mm_threshold,
            lower_percent_dose_cutoff=0,  # Points are already filtered
            global_normalisation=norm_dose if global_normalisation else None,
            interp_algo='scipy',
        )

        # --- Step 4: Create gamma map and calculate statistics ---
        gamma_map_for_display = np.full_like(mcc_dose_data, np.nan)
        gamma_map_for_display[original_indices_for_gamma] = gamma_values

        gamma_stats = {}
        valid_gamma = gamma_values[~np.isnan(gamma_values)]

        if len(valid_gamma) > 0:
            passed = valid_gamma <= 1
            gamma_stats['pass_rate'] = 100 * np.sum(passed) / len(valid_gamma)
            gamma_stats['mean'] = np.mean(valid_gamma)
            gamma_stats['max'] = np.max(valid_gamma)
            gamma_stats['min'] = np.min(valid_gamma)
            gamma_stats['total_points'] = len(valid_gamma)
        else:
            gamma_stats.update({'pass_rate': 0, 'mean': 0, 'max': 0, 'min': 0, 'total_points': 0})

        logger.info(f"Gamma analysis complete: {gamma_stats.get('total_points', 0)} points analyzed, pass rate {gamma_stats.get('pass_rate', 0):.1f}%")

        # --- Step 5: Create interpolated MCC data for visualization purposes only ---
        mcc_interp_data = griddata(
            all_mcc_coords_phys,
            all_mcc_dose_values,
            (dicom_phys_x_mesh, dicom_phys_y_mesh),
            method='linear',
            fill_value=0
        )

        return gamma_map_for_display, gamma_stats, phys_extent, mcc_interp_data

    except Exception as e:
        logger.error(f"Error during gamma analysis: {str(e)}")
        raise
