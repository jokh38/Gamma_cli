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
    hor_profile_data
):
    """Generates a report with all the analysis information."""
    fig = plt.figure(figsize=(16, 12))
    gs = fig.add_gridspec(4, 2)

    # Header
    ax_header = fig.add_subplot(gs[0, 0])
    ax_header.axis('off')
    institution, patient_id, patient_name = dicom_handler.get_patient_info()
    header_text = (
        f"Institution: {institution}\n"
        f"Patient ID: {patient_id}\n"
        f"Patient Name: {patient_name}"
    )
    ax_header.text(0.01, 0.5, header_text, transform=ax_header.transAxes, fontsize=12, verticalalignment='center')

    # Parameters
    ax_params = fig.add_subplot(gs[0, 1])
    ax_params.axis('off')
    params_text = (
        f"DTA: {dta} mm\n"
        f"DD: {dd} %\n"
        f"Suppression Level: {suppression_level}%"
    )
    ax_params.text(0.99, 0.5, params_text, transform=ax_params.transAxes, fontsize=12, verticalalignment='center', horizontalalignment='right')

    # DICOM Dose
    ax_dicom = fig.add_subplot(gs[1, 0])
    im_dicom = ax_dicom.imshow(dicom_handler.get_pixel_data(), cmap='jet', extent=dicom_handler.get_physical_extent())
    fig.colorbar(im_dicom, ax=ax_dicom, label='Dose (Gy)')
    ax_dicom.set_title('DICOM RT Dose')
    bounds = dicom_handler.dose_bounds
    if bounds:
        width = bounds['max_x'] - bounds['min_x']
        height = bounds['max_y'] - bounds['min_y']
        rect = Rectangle((bounds['min_x'], bounds['min_y']), width, height, linewidth=1, edgecolor='r', facecolor='none', alpha=0.7)
        ax_dicom.add_patch(rect)

    # MCC Dose
    ax_mcc = fig.add_subplot(gs[1, 1])
    im_mcc = ax_mcc.imshow(mcc_handler.get_interpolated_matrix_data(), cmap='jet', extent=mcc_handler.get_physical_extent())
    fig.colorbar(im_mcc, ax=ax_mcc, label='Dose')
    ax_mcc.set_title('MCC Dose')
    bounds = mcc_handler.dose_bounds
    if bounds:
        width = bounds['max_x'] - bounds['min_x']
        height = bounds['min_y'] - bounds['max_y']
        rect = Rectangle((bounds['min_x'], bounds['max_y']), width, height, linewidth=1, edgecolor='r', facecolor='none', alpha=0.7)
        ax_mcc.add_patch(rect)

    # Gamma Map
    ax_gamma = fig.add_subplot(gs[2, 0])
    im_gamma = ax_gamma.imshow(gamma_map, cmap='coolwarm', extent=mcc_handler.get_physical_extent(), vmin=0, vmax=2)
    fig.colorbar(im_gamma, ax=ax_gamma, label='Gamma Index')
    ax_gamma.set_title(f'Gamma Analysis (Pass: {gamma_stats.get("pass_rate", 0):.1f}%)')

    # Vertical Profile
    ax_ver_profile = fig.add_subplot(gs[3, 0])
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

        bounds = dicom_handler.dose_bounds
        if bounds:
            ax_ver_profile.set_xlim(bounds['min_y'], bounds['max_y'])


    # Horizontal Profile
    ax_hor_profile = fig.add_subplot(gs[3, 1])
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

        bounds = dicom_handler.dose_bounds
        if bounds:
            ax_hor_profile.set_xlim(bounds['min_x'], bounds['max_x'])

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(output_path, format='jpeg')
    plt.close(fig)
