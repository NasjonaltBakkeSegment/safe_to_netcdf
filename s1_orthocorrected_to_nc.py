#!/usr/bin/python3

# Name:          XX.py
# Purpose:       TBD
# Author(s):     Trygve Halsne,
# Created:
# Modifications:
# Copyright:     (c) Norwegian Meteorological Institute, 2022
#

import sys
from collections import defaultdict
import netCDF4
import numpy as np
import pathlib
import datetime as dt
#import safe_to_netcdf.utils as utils
#import utils
import logging
import xarray as xa
from pyproj import CRS
import os

#logger = logging.getLogger(__name__)


class S1_orthocorrected_to_nc:
    """ Should do:
            1. read Products (including all bands)
            2. create a netCDF4 file

        Keyword arguments:
    """

    def __init__(self, product, srcdir, outdir):
        self.product = product
        self.srcdir = srcdir
        self.t0 = dt.datetime.now(dt.timezone.utc)
        self.ncout = None  # NetCDF output file

        self.input_files = self.main()
        print(self.input_files)

    def main(self):
        """test
        """
        # 1. read all products
        input_files = []
        for file in os.listdir(self.srcdir):
            if file.endswith(".tiff"):
                input_files.append(self.srcdir +'/'+ file)
        return input_files

    def write_to_NetCDF(self, nc_outpath, compression_level, chunk_size=(1, 32, 32)):

        """ Method writing each input_file to NetCDF product.

        Keyword arguments:
        nc_outpath -- output path where NetCDF file should be stored
        compression_level -- compression level on output NetCDF file (1-9)
        chunk_size -- chunk_size
        """

        logger.info("------------START CONVERSION FROM SAFE TO NETCDF-------------")

        # Status
        logger.info('Creating NetCDF file')
        with (netCDF4.Dataset(out_netcdf, 'w', format='NETCDF4')) as ncout:
            for i, current_file in enumerate(self.input_files):
                da = xa.open_rasterio(current_file)
                varName = current_file.split('/')[-1].split('.')[0]
                if i == 0:
                    ncout.createDimension('time', 0)
                    ncout.createDimension('x', nx)
                    ncout.createDimension('y', ny)

                    utils.create_time(ncout, self.globalAttribs["PRODUCT_START_TIME"])

                    # Add projection coordinates
                    ##########################################################
                    # Status
                    logger.info('Adding projection coordinates')
                    utils.memory_use(self.t0)


                    ncx = ncout.createVariable('x', 'i4', 'x', zlib=True, complevel=compression_level)
                    ncx.units = 'm'
                    ncx.standard_name = 'projection_x_coordinate'
                    ncx[:] = da.x.astype(np.int32)

                    ncy = ncout.createVariable('y', 'i4', 'y', zlib=True, complevel=compression_level)
                    ncy.units = 'm'
                    ncy.standard_name = 'projection_y_coordinate'
                    ncy[:] = da.y.astype(np.int32)

                    # set grid mapping
                    ##########################################################
                    crs_proj4 = da.crs
                    source_crs = CRS.from_proj4(crs_proj4)
                    crs_cf = source_crs.to_cf()

                    nc_crs = ncout.createVariable('crs', 'i1')
                    nc_crs.setncatts(crs_cf)
                    ncout.sync()

                # Iterate over each band
                assert da.band.size == 1, "Tiffs should be single band files, only"

                print('NOTE: MAKE DATA TYPE CONFIG SPECIFIC')
                dataType = np.float32
                varout = ncout.createVariable(varName, DataType, ('time', 'y', 'x'), fill_value=np.nan,
                                              zlib=True, complevel=compression_level, chunksizes=chunk_size)

                print('NOTE: read varible attributes from config file')
                longname='test'
                varout.grid_mapping = "crs"
                varout.long_name = longName
                break
                #varout[0, :] = da.data




#if __name__ == '__main__':
# Log to console
#logger = logging.getLogger()
#logger.setLevel(logging.DEBUG)
#log_info = logging.StreamHandler(sys.stdout)
#log_info.setFormatter(logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s'))
#logger.addHandler(log_info)
