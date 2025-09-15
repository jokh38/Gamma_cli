import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle

def generate_report(
    output_path,
    dicom_handler,
    mcc_handler,
    gamma_map,
    gamma_stats,
    dta,
    dd,
    suppression_level,
    ver_profile_data,
    hor_profile_data,
    bounds=None,
    mcc_interp_data=None
):
    """Generates a report with all the analysis information."""
    fig = plt.figure(figsize=(12, 18))
    gs = fig.add_gridspec(3, 2, height_ratios=[1, 1, 1])

    # Patient Info Header
    institution, patient_id, patient_name = dicom_handler.get_patient_info()
    fig.suptitle(f'Institution: {institution} | Patient: {patient_name} ({patient_id})', fontsize=16)

    # 1. 2D Dose Plots
    # DICOM Dose
    ax_dicom = fig.add_subplot(gs[0, 0])
    dicom_crop_data, dicom_crop_extent = dicom_handler.get_cropped_data(bounds)
    if dicom_crop_data is not None:
        im_dicom = ax_dicom.imshow(dicom_crop_data, cmap='jet', extent=dicom_crop_extent, aspect='equal')
        fig.colorbar(im_dicom, ax=ax_dicom, label='Dose (Gy)')
    ax_dicom.set_title('DICOM RT Dose')
    ax_dicom.set_xlabel('Position (mm)')
    ax_dicom.set_ylabel('Position (mm)')

    # MCC Dose
    ax_mcc = fig.add_subplot(gs[0, 1])
    if mcc_interp_data is not None and bounds is not None:
        # Since mcc_interp_data is on the DICOM grid, we can use dicom_handler's methods
        min_px, min_py = dicom_handler.physical_to_pixel_coord(bounds['min_x'], bounds['min_y'])
        max_px, max_py = dicom_handler.physical_to_pixel_coord(bounds['max_x'], bounds['max_y'])

        # Swap if coordinates are descending
        if min_py > max_py: min_py, max_py = max_py, min_py
        if min_px > max_px: min_px, max_px = max_px, min_px

        # Ensure indices are within array bounds
        min_py = max(0, min_py)
        min_px = max(0, min_px)
        max_py = min(mcc_interp_data.shape[0], max_py)
        max_px = min(mcc_interp_data.shape[1], max_px)

        mcc_crop_data = mcc_interp_data[min_py:max_py, min_px:max_px]

        im_mcc = ax_mcc.imshow(mcc_crop_data, cmap='jet', extent=dicom_crop_extent, aspect='equal')
        fig.colorbar(im_mcc, ax=ax_mcc, label='Dose')

    ax_mcc.set_title('MCC Dose (Interpolated)')
    ax_mcc.set_xlabel('Position (mm)')
    ax_mcc.set_ylabel('Position (mm)')


    # 2. Profile Plots
    # Horizontal Profile
    ax_hor_profile = fig.add_subplot(gs[1, 0])
    if hor_profile_data:
        ax_hor_profile.plot(hor_profile_data['phys_coords'], hor_profile_data['dicom_values'], 'b-', label='RT dose')
        if 'mcc_interp' in hor_profile_data:
            ax_hor_profile.plot(hor_profile_data['phys_coords'], hor_profile_data['mcc_interp'], 'r-', label='mcc dose')
        if 'mcc_values' in hor_profile_data and 'mcc_phys_coords' in hor_profile_data:
            ax_hor_profile.plot(hor_profile_data['mcc_phys_coords'], hor_profile_data['mcc_values'], 'r.', markersize=3)
        ax_hor_profile.set_xlabel('Position (mm)')
        ax_hor_profile.set_ylabel('Dose (Gy)')
        ax_hor_profile.set_title('Left-Right Profile (Horizontal)')
        ax_hor_profile.legend()
        ax_hor_profile.grid(True)
        if bounds:
            ax_hor_profile.set_xlim(bounds['min_x'], bounds['max_x'])

    # Vertical Profile
    ax_ver_profile = fig.add_subplot(gs[1, 1])
    if ver_profile_data:
        ax_ver_profile.plot(ver_profile_data['phys_coords'], ver_profile_data['dicom_values'], 'b-', label='RT dose')
        if 'mcc_interp' in ver_profile_data:
            ax_ver_profile.plot(ver_profile_data['phys_coords'], ver_profile_data['mcc_interp'], 'r-', label='mcc dose')
        if 'mcc_values' in ver_profile_data and 'mcc_phys_coords' in ver_profile_data:
            ax_ver_profile.plot(ver_profile_data['mcc_phys_coords'], ver_profile_data['mcc_values'], 'r.', markersize=3)
        ax_ver_profile.set_xlabel('Position (mm)')
        ax_ver_profile.set_ylabel('Dose (Gy)')
        ax_ver_profile.set_title('In-Out Profile (Vertical)')
        ax_ver_profile.legend()
        ax_ver_profile.grid(True)
        if bounds:
            ax_ver_profile.set_xlim(bounds['min_y'], bounds['max_y'])

    # 3. Gamma Analysis
    # Gamma Map
    ax_gamma = fig.add_subplot(gs[2, 0])
    gamma_crop_data, gamma_crop_extent = mcc_handler.crop_array_by_bounds(gamma_map, bounds)
    if gamma_crop_data is not None:
        im_gamma = ax_gamma.imshow(gamma_crop_data, cmap='coolwarm', extent=gamma_crop_extent, vmin=0, vmax=2, aspect='equal')
        fig.colorbar(im_gamma, ax=ax_gamma, label='Gamma Index')
    ax_gamma.set_title(f'Gamma Analysis (Pass: {gamma_stats.get("pass_rate", 0):.1f}%)')
    ax_gamma.set_xlabel('Position (mm)')
    ax_gamma.set_ylabel('Position (mm)')


    # Gamma Stats Text
    ax_gamma_text = fig.add_subplot(gs[2, 1])
    ax_gamma_text.axis('off')

    pass_rate = gamma_stats.get('pass_rate', 0)
    total_points = gamma_stats.get('total_points', 0)
    passed_points = int(total_points * pass_rate / 100)
    failed_points = total_points - passed_points

    stats_text = (
        f"Gamma Analysis Results\n\n"
        f"Acceptance Criteria:\n"
        f"  DTA: {dta} mm\n"
        f"  Dose Difference: {dd} %\n"
        f"  Threshold: {suppression_level} %\n\n"
        f"Results:\n"
        f"  Passing Rate: {pass_rate:.2f} %\n"
        f"  Analyzed Pixels: {total_points}\n"
        f"  Passed Pixels: {passed_points}\n"
        f"  Failed Pixels: {failed_points}\n"
    )
    ax_gamma_text.text(0.05, 0.95, stats_text, transform=ax_gamma_text.transAxes, fontsize=12,
                     verticalalignment='top', bbox=dict(boxstyle='round,pad=0.5', fc='aliceblue', alpha=0.5))

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig(output_path, format='jpeg')
    plt.close(fig)
