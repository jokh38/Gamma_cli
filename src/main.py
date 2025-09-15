import argparse
import yaml
import sys
import os

from file_handlers import DicomFileHandler, MCCFileHandler
from analysis import perform_gamma_analysis, extract_profile_data
from reporting import generate_report
from utils import logger

def main():
    """Main function to run the gamma analysis from the command line."""
    parser = argparse.ArgumentParser(description="Run gamma analysis on DICOM RT and MCC files.")
    parser.add_argument("--rtplan", type=str, required=True, help="Path to the DICOM RT plan file.")
    parser.add_argument("--measure", type=str, required=True, help="Path to the measurement data file (MCC).")
    args = parser.parse_args()

    logger.info(f"Starting gamma analysis with the following arguments:")
    logger.info(f"  RT Plan: {args.rtplan}")
    logger.info(f"  Measurement: {args.measure}")

    # Load configuration
    try:
        with open("config.yaml", "r") as f:
            config = yaml.safe_load(f)
        dta = config.get("dta", 3)
        dd = config.get("dd", 3)
        suppression_level = config.get("suppression_level", 10)
        logger.info(f"Loaded configuration: dta={dta}, dd={dd}, suppression_level={suppression_level}")
    except FileNotFoundError:
        logger.error("Error: config.yaml not found. Please create it.")
        sys.exit(1)
    except yaml.YAMLError as e:
        logger.error(f"Error: Could not parse config.yaml. Please check its format: {e}")
        sys.exit(1)

    # Initialize file handlers
    dicom_handler = DicomFileHandler()
    mcc_handler = MCCFileHandler()

    # Load files
    if not dicom_handler.open_file(args.rtplan):
        logger.error(f"Failed to load DICOM file: {args.rtplan}")
        sys.exit(1)

    if not mcc_handler.open_file(args.measure):
        logger.error(f"Failed to load MCC file: {args.measure}")
        sys.exit(1)

    # Crop MCC data to match the DICOM's calculated dose bounds
    if dicom_handler.dose_bounds:
        mcc_handler.crop_to_bounds(dicom_handler.dose_bounds)
        logger.info("Applied DICOM bounds to MCC data.")

    # Perform gamma analysis
    try:
        gamma_map, gamma_stats, phys_extent, mcc_interp_data = perform_gamma_analysis(
            mcc_handler, dicom_handler,
            dd, dta,
            global_normalisation=True  # Assuming global gamma, can be made configurable
        )

        if 'pass_rate' in gamma_stats:
            print("Gamma Analysis Results:")
            print(f"  Pass Rate: {gamma_stats['pass_rate']:.2f}%")
            print(f"  Mean Gamma: {gamma_stats['mean']:.3f}")
            print(f"  Max Gamma: {gamma_stats['max']:.3f}")
            print(f"  Min Gamma: {gamma_stats['min']:.3f}")
            print(f"  Total Points: {gamma_stats['total_points']}")

            # Generate profile data for the report
            ver_profile_data = extract_profile_data(
                direction="vertical",
                fixed_position=0,
                dicom_handler=dicom_handler,
                mcc_handler=mcc_handler
            )
            hor_profile_data = extract_profile_data(
                direction="horizontal",
                fixed_position=0,
                dicom_handler=dicom_handler,
                mcc_handler=mcc_handler
            )

            # Generate report
            reports_dir = "reports"
            os.makedirs(reports_dir, exist_ok=True)
            _, patient_id, _ = dicom_handler.get_patient_info()
            output_filename = os.path.join(reports_dir, f"report_{patient_id}.jpg")
            generate_report(
                output_filename,
                dicom_handler,
                mcc_handler,
                gamma_map,
                gamma_stats,
                dta,
                dd,
                suppression_level,
                ver_profile_data,
                hor_profile_data,
                mcc_interp_data=mcc_interp_data
            )
            logger.info(f"Report saved to {output_filename}")

        else:
            print("Gamma analysis did not produce valid results.")

    except Exception as e:
        logger.error(f"An error occurred during gamma analysis: {e}", exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
