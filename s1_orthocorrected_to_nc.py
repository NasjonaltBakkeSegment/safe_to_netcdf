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
import utils as utils
import constants as cst
import logging
import xarray as xa
from pyproj import CRS
import os
import xarray as xa

logger = logging.getLogger(__name__)

print('TODO:Add following metadata info(?): Add generation time, Add interpolation method?, ADD DEM?,If thermal denoising? pixel spacing')
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

    def main(self):
        """test
        """
        # 1. read all products
        input_files = []
        for file in os.listdir(self.srcdir):
            if file.endswith(".tiff"):
                input_files.append(self.srcdir +'/'+ file)
        return input_files

    def write_to_NetCDF(self, fname, compression_level, chunk_size=(1, 32, 32)):

        """ Method writing each input_file to NetCDF product.

        Keyword arguments:
        nc_outpath -- output path where NetCDF file should be stored
        compression_level -- compression level on output NetCDF file (1-9)
        chunk_size -- chunk_size
        """

        logger.info("------------START CONVERSION FROM SAFE TO NETCDF-------------")

        # Status
        logger.info('Creating NetCDF file')
        with (netCDF4.Dataset(fname, 'w', format='NETCDF4')) as ncout:
            for i, current_file in enumerate(self.input_files):
                da = xa.open_rasterio(current_file)
                varName = current_file.split('/')[-1].split('.')[0]
                if i == 0:
                    ncout.createDimension('time', 0)
                    ncout.createDimension('x', da.x.size)
                    ncout.createDimension('y', da.y.size)


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

                if varName=='mask':
                    dataType = 'u1' #Byte unsigned integer 8
                else:
                    dataType = np.float32
                varout = ncout.createVariable(varName, dataType, ('time', 'y', 'x'), fill_value=0,
                                              zlib=True, complevel=compression_level, chunksizes=chunk_size)

                if varName == 'mask':
                    varout.flag_values = np.arange(4, dtype=np.int8)+1
                    varout.flag_meanings = ' '.join(['Missing','Shadow','Folding','Valid'])
                    longName='Processing mask'
                elif varName=='incidence':
                    varout.standard_name = "angle_of_incidence"
                    varout.units="rad"
                    longName='local incidence angle'
                elif varName=='height':
                    varout.standard_name = "height"
                    varout.units="m"
                    varout.positive='up'
                    longName='Surface height from DEM'
                else:
                    varout.standard_name = "surface_backwards_scattering_coefficient_of_radar_wave"
                    varout.units="1"
                    varout.polarisation = varName.split('_')[-1]
                    longName=varName

                varout.grid_mapping = "crs"
                varout.long_name = longName
                varout[0, :] = da.data

            #Add global attributes
            globalAttribs={'rectified_time':da.TIFFTAG_DATETIME,
                            'software':da.TIFFTAG_SOFTWARE,
                            'title':'Rectified Sentinel-1 GRD data',
                            'file_creation_date': '{}'.format(dt.datetime.now()),
                            'Conventions': 'CF-1.8',
                            'summary': 'Rectified Sentinel-1 C-band SAR GRD product',
                            'keywords': '[Earth Science, Spectral/Engineering, RADAR, RADAR backscatter], [Earth Science, Spectral/Engineering, RADAR, RADAR imagery], [Earth Science, Spectral/Engineering, Microwave, Microwave Imagery]',
                            'keywords_vocabulary': 'GCMD Science Keywords',
                            'institution': 'Norwegian Meteorological Institute',
                            'history': '{} rectified from original product by the NBS team'.format(da.TIFFTAG_DATETIME),
                            'source': 'surface observation'
            }
            ncout.setncatts(globalAttribs)
            ncout.sync()

        # Status
        logger.info('Finished.')
        #utils.memory_use(self.t0)
                #if i ==3:
                #    break




#if __name__ == '__main__':
# Log to console
#logger = logging.getLogger()
#logger.setLevel(logging.DEBUG)
#log_info = logging.StreamHandler(sys.stdout)
#log_info.setFormatter(logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s'))
#logger.addHandler(log_info)
