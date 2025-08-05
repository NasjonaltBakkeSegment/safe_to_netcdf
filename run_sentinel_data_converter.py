# import os
# import re
# import sys
# import math
# import logging
# import pathlib
# from pathlib import Path
# from collections import defaultdict
# import datetime as dt
# import xml.etree.ElementTree as ET
# import numpy as np
# import netCDF4
# import xarray as xr
# import rasterio
# import scipy.ndimage
# import geopandas as geopd
# import pyproj
# from pyproj import CRS
# import utils as utils
# import constants as cst
# import argparse
# import io
# from io import BytesIO
# from netCDF4 import Dataset
# from rasterio.io import MemoryFile
# import tempfile
# import rioxarray as rio
# from s1_reader_and_NetCDF_converter import Sentinel1_reader_and_NetCDF_converter
# from s2_reader_and_NetCDF_converter import Sentinel2_reader_and_NetCDF_converter
# from s3_olci_l1_reader_and_CF_converter import Sentinel3_olci_reader_and_CF_converter


import os
import sys
import logging
from pathlib import Path
import utils as utils
import argparse
from s1_reader_and_NetCDF_converter import Sentinel1_reader_and_NetCDF_converter
from s2_reader_and_NetCDF_converter import Sentinel2_reader_and_NetCDF_converter
from s3_olci_l1_reader_and_CF_converter import Sentinel3_olci_reader_and_CF_converter

logger = logging.getLogger(__name__)


def parse_args():
    parser = argparse.ArgumentParser(
    formatter_class=argparse.RawDescriptionHelpFormatter,
            description='Script to read and convert Sentinel 2 data to NetCDF or GeoTIFF. Inputs are: \n'+
            'input: either a list with paths to the sentinel data or a path to a directory containing sentinel data\n'+
            'output (not mandatory): either the path only of where to write the parent file or the path + parent filename or the parent filename only\n'+
            'format: either NetCDF or GeoTIFF (if GeoTIFF, output is a directory containing files for each individual raster band)\n'+
            'data_type: either Sentinel1 (S1), Sentinel 2 (S2), or Sentinel 3 (S3)\n'+
            'e.g. S2_reader_and_converter_NetCDF_GeoTIFF -i /path/to/sentinel/data -f NetCDF -o /path/to/output.nc - dt S2\n'+
            'e.g. S2_reader_and_converter_NetCDF_GeoTIFF -i /path/to/sentinel/data -f GeoTIFF -o /path/to/outputfiles -dt S2\n'+
            '...'
            )
    
    parser.add_argument(
        "--input", '-i',
        required=True,
        type=parse_input,
        help="Required: either a .txt file with one /path/to/SAFE.zip per line, or a single valid /path/to/SAFE.zip"
    )

    parser.add_argument(
        '--format', '-f',
        choices=['netcdf', 'geotiff'],
        required=True,
        help="Output format: 'netcdf' for NetCDF file or 'geotiff' for GeoTIFF file."
    )

    parser.add_argument(
        '--output', '-o',
        type=str,
        default=os.getcwd(),
        help="Path to the output directory. If not provided; current directory"
    )

    parser.add_argument(
        '--data_type', '-dt',
        choices=['S1', 'S2', 'S3'],
        help='Type of sentinel data to read and convert.'
    )

    return parser.parse_args()




def parse_input(path):
    if path.endswith(".txt") and path.is_file():
        # Read file and return list of paths
        with path.open() as f:
            return [line.strip() for line in f if line.strip()]
        print('Path is list!')
    elif Path(path).exists():
        # Return single path in a list
        return [str(path)]
    else:
        raise argparse.ArgumentTypeError(f"Invalid path or file: {path}")




def main():

    # Log to console
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    log_info = logging.StreamHandler(sys.stdout)
    log_info.setFormatter(logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s'))
    logger.addHandler(log_info)



    args = parse_args()

    for path in args.input:

        indir = Path(path).parent
        product = str(os.path.splitext(os.path.basename(path))[0])
        outdir = Path(args.output)
        outdir.parent.mkdir(parents=True, exist_ok=True)

        
        if product.startswith("S1") or args.data_type == 'S1':
            print(f'Product: {product}')
            print(f'Indir: {indir}')
            print(f'Outdir: {outdir}')
            conversion_object = Sentinel1_reader_and_NetCDF_converter(
                product=product,
                indir=indir,
                outdir=outdir)

        if product.startswith("S2") or args.data_type == 'S2':
            conversion_object = Sentinel2_reader_and_NetCDF_converter(
            product=product,
            indir=indir,
            outdir=outdir
            )

        if product.startswith("S3") or args.data_type == 'S3':

            conversion_object = Sentinel3_olci_reader_and_CF_converter(
                product=product,
                indir=indir,
                outdir=outdir)
            
        conversion_object.run()

        if args.format == 'geotiff':
        
            utils.write_to_geotiff(conversion_object, indir, product, outdir)                    

        if args.format == 'netcdf':

            if conversion_object.read_ok:
                conversion_object.write_to_NetCDF(outdir, outdir, product, 7)




if __name__ == '__main__':

    main()

