import os
import psutil
import logging
import argparse
from pathlib import Path
from lib.safe.safe_s1 import S1SAFEFile
from lib.safe.safe_s2 import S2SAFEFile
from lib.netcdf.netcdf_s1 import S1NetCDFFile
from lib.netcdf.netcdf_s2 import S2NetCDFFile
from lib.utils import load_variable_attributes, load_global_attributes, get_key
#import config.constants as cst

# Configure logging with timestamps
logging.basicConfig(
    level=logging.DEBUG,
    format="%(asctime)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S"
)
logger = logging.getLogger(__name__)

def log_memory_usage(stage):
    process = psutil.Process(os.getpid())
    memory_info = process.memory_info()
    logger.info(f"[{stage}] Memory usage: {memory_info.rss / 1024**2:.2f} MB (RSS), "
                f"{memory_info.vms / 1024**2:.2f} MB (VMS)")

def transform(
        product_name,
        safedir=None,
        netcdfdir=None,
        geotiffdir=None,
        tmpdir = '/home/lukem/Documents/MET/Projects/ESA_NBS/Git_repos/safe_to_netcdf/tmp',
        global_attributes_config = 'config/global_attributes.yaml',
        variable_attributes_config = 'config/variable_attributes.yaml',
        product_id = None,
        target = 'netCDF', # netCDF or geoTIFF
    ):
    log_memory_usage("Before transformation")

    # Validation logic
    if target == 'netCDF':
        if not safedir or not netcdfdir:
            raise ValueError("For target 'netCDF', both 'safedir' and 'netcdfdir' must be provided.")

    if target == 'geoTIFF':
        if not geotiffdir:
            raise ValueError("For target 'geoTIFF', 'geotiffdir' must be provided.")
        if not (netcdfdir or safedir):
            raise ValueError("For target 'geoTIFF', either 'netcdfdir' or 'safedir' must be provided.")

    #TODO: Do different things depending on target and whether netCDF exists

    logger.info("Starting the transformation workflow.")

    platform = product_name[0:3]
    mission = product_name[0:2]

    if safedir:
        safedir = Path(safedir)
    if netcdfdir:
        netcdfdir = Path(netcdfdir)
    if geotiffdir:
        geotiffdir = Path(geotiffdir)
    if tmpdir:
        tmpdir = Path(tmpdir)

    if mission == 'S1':
        safe_file = S1SAFEFile(product=product_name, zipdir=safedir, tmpdir=tmpdir)
        netcdf_file = S1NetCDFFile(product=product_name, directory=netcdfdir)
    elif mission == 'S2':
        safe_file = S2SAFEFile(product=product_name, zipdir=safedir, tmpdir=tmpdir)
        netcdf_file = S2NetCDFFile(product=product_name, directory=netcdfdir)

    # Prepare the SAFEFile instance for use
    logger.info("Preparing SAFEFile instance for use.")
    safe_file.prepare_for_use()

    logger.debug("Loaded variable attributes from config file.")
    variable_attributes = load_variable_attributes(
        variable_attributes_config,
        mission,
    )

    logger.info("Initializing NetCDF file.")
    netcdf_file.initialise()

    logger.info("Creating NetCDF dimensions.")
    netcdf_file.create_dimensions(safe_file)

    logger.info("Creating time variable in NetCDF.")
    if mission == 'S1':
        netcdf_file.create_time(safe_file.globalAttribs["ACQUISITION_START_TIME"])

        logger.info("Initializing latitude and longitude variables.")
        netcdf_file.init_lat_lon()
        logger.info("Generating latitude and longitude grids.")
        lat, lon = safe_file.genLatLon_regGrid()

        logger.info("Writing latitude variable to NetCDF.")
        netcdf_file.write_variable_with_preprocessing(
            'lat', lat,
            process_chunk=None,
            workers=8,
            sync_every=None,
        )
        logger.info("Writing longitude variable to NetCDF.")
        netcdf_file.write_variable_with_preprocessing(
            'lon', lon,
            process_chunk=None,
            workers=8,
            sync_every=None,
        )
        del lat, lon

        logger.info("Adding raw measurement layers to NetCDF.")
        netcdf_file.add_raw_measurement_layers(safe_file)

        logger.info("Writing grid mapping variable to NetCDF.")
        netcdf_file.write_grid_mapping_variable()

        logger.info("Adding calibration layers to NetCDF.")
        netcdf_file.add_calibration_layers(safe_file)

        logger.info("Adding noise layers to NetCDF.")
        netcdf_file.add_noise_layers(safe_file)

        logger.info("Adding subswath layers to NetCDF.")
        netcdf_file.add_subswath_layers(safe_file)

        logger.info("Adding GCP (Ground Control Points) information to NetCDF.")
        netcdf_file.add_gcp_information(safe_file)

        logger.info("Adding product annotation metadata to NetCDF.")
        netcdf_file.add_product_annotation_metadata(safe_file)

    elif mission == 'S2':
        netcdf_file.create_time(safe_file.globalAttribs["PRODUCT_START_TIME"])

        logger.info('Adding projection coordinates')
        netcdf_file.add_projection_coordinates(safe_file)

        logger.info('Adding band_id variable')
        netcdf_file.add_raster_band_id_variable()

        logger.info('Adding frequency bands layers')
        netcdf_file.add_frequency_band_layers(safe_file)

        #! Temporarily updating chunk size to keep preexisting consistency in chunk sizes between variables in previously processed S2 files
        netcdf_file.chunk_size = (1, 91, 99)

        logger.info('Adding to grid mapping variable')
        netcdf_file.write_grid_mapping_variable(safe_file)

        logger.info('Adding specific L1C or L2A layers')
        netcdf_file.add_specific_L1C_or_L2A_layers(safe_file)

        logger.info('Adding sun and view angles in native resolution')
        netcdf_file.add_sun_and_view_angles(safe_file)

        logger.info('Adding satellite orbit specific data')
        netcdf_file.add_orbit_specific_data(safe_file)

        logger.info('Adding integrated information as variables')
        netcdf_file.add_integrated_information_as_variables(safe_file)

    logger.info("Writing variable attributes to NetCDF.")
    netcdf_file.write_variable_attributes(variable_attributes)

    logger.info("Retrieving global attributes from SAFEFile.")
    safe_file.get_global_attributes()

    logger.info("Writing global attributes to NetCDF.")
    global_attributes_from_config = load_global_attributes(global_attributes_config, platform)
    global_attributes = global_attributes_from_config | safe_file.globalAttribs
    if product_id:
        global_attributes['id'] = product_id

    netcdf_file.write_global_attributes(global_attributes)

    logger.info("Closing NetCDF file and deleting unzipped file.")
    netcdf_file.close()
    safe_file.finalize_usage()

    log_memory_usage("After transformation")
    logger.info("End of job")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Transform product data to netCDF or GeoTIFF format.")

    # Required positional argument
    parser.add_argument("product_name", type=str, help="Name of the product.")

    # Optional directory arguments
    parser.add_argument("--safedir", type=str, default=None,
                        help="Directory containing SAFE file for the product (default: %(default)s).")
    parser.add_argument("--netcdfdir", type=str, default=None,
                        help="Directory containing netCDF file for the product (default: %(default)s).")
    parser.add_argument("--geotiffdir", type=str, default=None,
                        help="Directory containing GeoTIFF file for the product (default: %(default)s).")

    # Other optional arguments with default values
    parser.add_argument("--tmpdir", type=str, default='/home/lukem/Documents/MET/Projects/ESA_NBS/Git_repos/safe_to_netcdf/tmp',
                        help="Path to the temporary directory (default: %(default)s).")
    parser.add_argument("--global_attributes_config", type=str, default='config/global_attributes.yaml',
                        help="Path to the global attributes configuration file (default: %(default)s).")
    parser.add_argument("--variable_attributes_config", type=str, default='config/variable_attributes.yaml',
                        help="Path to the variable attributes configuration file (default: %(default)s).")
    parser.add_argument("--product_id", type=str, default=None,
                        help="Optional product ID (default: %(default)s).")
    parser.add_argument("--target", type=str, choices=["netCDF", "geoTIFF"], default="netCDF",
                        help="Target format, either 'netCDF' or 'geoTIFF' (default: %(default)s).")

    # Parse arguments
    args = parser.parse_args()

    # Pass arguments to the transform function
    transform(
        product_name=args.product_name,
        safedir=args.safedir,
        netcdfdir=args.netcdfdir,
        geotiffdir=args.geotiffdir,
        tmpdir=args.tmpdir,
        global_attributes_config=args.global_attributes_config,
        variable_attributes_config=args.variable_attributes_config,
        product_id=args.product_id,
        target=args.target,
    )

# TODO: Redo for S1
# TODO: Develop for S2 and S3 OLCI
# TODO: Check memory usage to see how many can run in parallel
# TODO: Test run and write times with different compressions and chunk sizes
# TODO: Write to geoTIFF
# TODO: Update routine to retrieve global attributes from MMD where possible