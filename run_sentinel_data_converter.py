import os
import sys
import logging
from pathlib import Path
import utils as utils
import argparse
from s1_reader_and_NetCDF_converter import Sentinel1_reader_and_NetCDF_converter
from s2_reader_and_NetCDF_converter import Sentinel2_reader_and_NetCDF_converter

logger = logging.getLogger(__name__)


def parse_args():
    parser = argparse.ArgumentParser(
    formatter_class=argparse.RawDescriptionHelpFormatter,
            description='Script to read and convert Sentinel data to NetCDF or GeoTIFF. Inputs are: \n'+
            'input: either a list with paths to the sentinel data or a path to a directory containing sentinel data\n'+
            'output (not mandatory): either the path only of where to write the parent file or the path + parent filename or the parent filename only\n'
            )

    parser.add_argument(
        "--input_filepath", '-i',
        required=True,
        help="Required: a single valid /path/to/SAFE.zip"
    )

    parser.add_argument(
        '--output', '-o',
        type=str,
        default=os.getcwd(),
        help="Path to the output directory. If not provided; current directory"
    )

    return parser.parse_args()

def main():
    # Log to console
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    log_info = logging.StreamHandler(sys.stdout)
    log_info.setFormatter(logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s'))
    logger.addHandler(log_info)

    args = parse_args()

    # Get input directory and product name
    indir = os.path.dirname(args.input_filepath)
    product = os.path.splitext(os.path.basename(args.input_filepath))[0]

    outdir = args.output

    print(product)
    print(outdir)
    print(indir)

    conversion_object = Sentinel2_reader_and_NetCDF_converter(
        product=product,
        indir=Path(indir),
        outdir=Path(outdir)
    )

    if conversion_object.read_ok:
        conversion_object.write_to_NetCDF(outdir, 7)


if __name__ == '__main__':
    main()
