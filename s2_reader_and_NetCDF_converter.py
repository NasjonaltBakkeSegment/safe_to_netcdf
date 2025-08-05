#!/usr/bin/python3

# Name:          s2_reader_and_NetCDF_converter.py
# Purpose:       Read Sentinel-2 MSI L1C/L2A data from ESA SAFE into a single
#                object with methods for extracting variables, tabels, gml
#                files, etc. as raster layers as well as a method converting
#                the product to netCDF.
#
#                Note: the routine also works for S2 MSI L1C products produced
#                by ESA with Norwegian DEM (*DTERRENGDATA* products).
#
# Author(s):     Trygve Halsne, Elodie Fernandez
# Created:
# Modifications:
# Copyright:     (c) Norwegian Meteorological Institute, 2018
#
# Need to use gdal 2.1.1-> to have support of the SAFE reader

# # Standard library imports
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
# import argparse
# import io
# from io import BytesIO
# import tempfile

# # Third-party library imports
# import numpy as np
# import netCDF4
# from netCDF4 import Dataset
# import xarray as xr
# import rasterio
# from rasterio.io import MemoryFile
# import rioxarray as rio
# import imageio.v3 as iio
# import scipy.ndimage
# import geopandas as geopd
# import pyproj
# from pyproj import CRS
# from pyproj import CRS as PyCRS
# from affine import Affine

# # Custom library imports
# import utils as utils
# import constants as cst


# Standard library imports
import math
import logging
from pathlib import Path
from collections import defaultdict
import datetime as dt

# Third-party library imports
import numpy as np
import netCDF4
import rasterio
import scipy.ndimage
import geopandas as geopd
import pyproj
from pyproj import CRS
from pyproj import CRS as PyCRS


# Custom library imports
import utils as utils
import constants as cst

logger = logging.getLogger(__name__)



class Sentinel2_reader_and_NetCDF_converter:
    ''' Class for reading Sentinel-2 MSI L1C/L2A products from SAFE with methods for
        reading auxilary information as e.g. clouds, solar and view angles.
        In addition, it is possible to convert product into NetCDF4/CF (1.6).

        The implemented methods uses standard python libraries as
        gdal(v. > 2.1.1), numpy, lxml etc.

        Keyword arguments:
        SAFE_file -- absolute path to zipped file
        SAFE_outdir -- output storage location for unzipped SAFE product
        '''

    def __init__(self, product, indir, outdir, colhub_uuid=None):
        

        self.uuid = colhub_uuid
        self.product_id = product
        file_path = indir / (product + '.zip')
        print(f'Indir = {indir}')
        print(f'Product = {product}')
        print(f'Path = {file_path}')
        if file_path.exists():
            self.input_zip = file_path
        else:
            file_path_safe = indir / (product + '.SAFE.zip')
            if file_path_safe.exists():
                self.input_zip = file_path_safe
        self.baseline = self.product_id.split('_')[3]
        self.SAFE_dir = (outdir / self.product_id).with_suffix('.SAFE')
        self.processing_level = 'Level-' + self.product_id.split('_')[1][4:6]
        self.xmlFiles = defaultdict(list)
        self.imageFiles = defaultdict(list)
        self.globalAttribs = {}
        self.src = None
        self.t0 = dt.datetime.now(dt.timezone.utc)
        self.ncout = None  # NetCDF output file
        self.reference_band = None
        self.dterrengdata = False  # variable saying if products is Norwegian DEM L1C
        self.sunAndViewAngles = defaultdict(list)
        self.vectorInformation = defaultdict(list)
        self.image_list_dterreng = []
        self.read_ok = True

        #self.run()

    def run(self):
        """ Main method for traversing and reading key values from SAFE
            directory.
        """

        # 1) unzip SAFE archive
        utils.uncompress(self)

        # 2) Set some of the global __init__ variables
        utils.initializer(self)

        # 3) Read sun and view angles
        logger.info('Read view and sun angles')

        if self.dterrengdata:
            currXml = self.xmlFiles.get('MTD_TL', None)
        elif self.baseline == "N0207":
            currXml = self.xmlFiles['S2_{}_Tile1_Data'.format(self.processing_level)]
        else:
            currXml = self.xmlFiles['S2_{}_Tile1_Metadata'.format(self.processing_level)]
        # Check for both None and empty list
        if currXml is None or not currXml:
            logger.error("xml file not found in SAFE directory. Hence exiting")
            self.read_ok = False
            return False
        self.readSunAndViewAngles(currXml)


    def write_to_NetCDF(self, nc_out, outdir, product, compression_level, chunk_size=(1, 91, 99)):
        logger.info("------------START CONVERSION FROM SAFE TO NETCDF-------------")
        utils.memory_use(self.t0)

        SAFE_root = Path(outdir) / f"{product}.SAFE"

        #print(f'SAFE_root: ------------------- {SAFE_root} --------------------------------------------')

        img_data_dirs = list(SAFE_root.rglob("IMG_DATA"))

        ten_meter_bands = {"B02", "B03", "B04", "B08"}

        # Get reference band (10m resolution)
        reference_band_set = False
        for img_data_dir in img_data_dirs:
            if img_data_dir.is_dir():
                for subdir in sorted(img_data_dir.iterdir()):
                    if subdir.is_file():
                        band_code = subdir.stem[-3:]
                        if band_code in ten_meter_bands:
                            with rasterio.open(subdir) as ref:
                                self.reference_band = ref
                                ny, nx = ref.height, ref.width
                                reference_band_set = True
                            break
                    elif subdir.is_dir() and '10m' in subdir.name:
                        for file_path in sorted(subdir.iterdir()):
                            if file_path.is_file() and not any(x in str(file_path) for x in ['TCI', 'True color image']):
                                with rasterio.open(file_path) as ref:
                                    self.reference_band = ref
                                    ny, nx = ref.height, ref.width
                                    reference_band_set = True
                                break
                if reference_band_set:
                    break

        if not reference_band_set:
            raise RuntimeError("No 10m reference band found")

        nxa, nya = self.sunAndViewAngles[list(self.sunAndViewAngles)[0]].shape
        out_netcdf = (nc_out / self.product_id).with_suffix('.nc')

        with netCDF4.Dataset(out_netcdf, 'w', format='NETCDF4') as ncout:
            # Dimensions
            ncout.createDimension('time', 0)
            ncout.createDimension('x', nx)
            ncout.createDimension('y', ny)
            ncout.createDimension('raster_band_id', len(cst.s2_bands_order.keys()))
            ncout.createDimension('xa', nxa)
            ncout.createDimension('ya', nya)

            utils.create_time(ncout, self.globalAttribs["PRODUCT_START_TIME"])

            logger.info('Adding projection coordinates')
            utils.memory_use(self.t0)

            xnp, ynp = self.genLatLon(nx, ny, latlon=False)

            x_var = ncout.createVariable('x', 'i4', 'x', zlib=True, complevel=compression_level)
            x_var.units = 'm'
            x_var.standard_name = 'projection_x_coordinate'
            x_var[:] = xnp

            y_var = ncout.createVariable('y', 'i4', 'y', zlib=True, complevel=compression_level)
            y_var.units = 'm'
            y_var.standard_name = 'projection_y_coordinate'
            y_var[:] = ynp

            # Raster band ID
            nc_rasterband_id = ncout.createVariable('band_id', 'i4', 'raster_band_id', zlib=True, complevel=compression_level)
            band_keys = list(cst.s2_bands_order.keys())
            nc_rasterband_id[:] = np.array(band_keys)
            nc_rasterband_id.long_name = 'raster band id'
            nc_rasterband_id.flag_values = np.array(band_keys, dtype=np.int8)
            nc_rasterband_id.flag_meanings = ' '.join(cst.s2_bands_order.values())

            # Set CRS from first valid file
            first_image_path = None
            for img_data_dir in img_data_dirs:
                for item in img_data_dir.rglob('*'):
                    if item.is_file() and not any(x in str(item) for x in ['TCI', 'True color image']):
                        first_image_path = item
                        break
                if first_image_path:
                    break

            if first_image_path is None:
                raise RuntimeError("No valid image files found.")
            
            with rasterio.open(first_image_path) as src:
                try:
                    crs = src.crs if src.crs and src.crs.is_valid else CRS.from_epsg(32631)  # UTM zone 31N
                    nc_crs = ncout.createVariable('UTM_projection', np.int32)
                    nc_crs.crs_wkt = crs.to_wkt()
                    nc_crs.grid_mapping_name = crs.to_dict().get('proj', 'unknown')
                    epsg = crs.to_epsg()
                    if epsg:
                        nc_crs.epsg_code = epsg
                except Exception as e:
                    raise RuntimeError(f"Failed to extract CRS: {e}")

            # Process and write all bands (flat or subdirectory layout)
            logger.info('Adding frequency bands layers')
            utils.memory_use(self.t0)

            if self.dterrengdata:
                # For DTERR data, gdal fails to properly do the src.GetSubDatasets()
                # so manually read the list of images created beforehand
                image_list = [[str(i), i.stem] for i in self.image_list_dterreng]
            else:
                #Checking whether image files are structured by resolution in subdirectories 
                image_list = []
                for img_data_dir in img_data_dirs:
                    for subitem in sorted(img_data_dir.iterdir()): # sorting to assert highest resolution is added to nc-file
                        if subitem.is_file():
                            image_list.append(subitem)
                        elif subitem.is_dir():
                            for file in subitem.iterdir():
                                if file.is_file():
                                    image_list.append(file)
                        else:
                            continue

            # Defining the valid band names
            valid_band_names = set(cst.s2_bands_aliases.keys())

            for image_path in image_list:
                if not image_path.is_file() or ("True color image" in str(image_path)) or ('TCI' in str(image_path)):
                #print(image_path)
                    continue

                with rasterio.open(image_path) as src:
                    #print(f'FOUND IMAGE {src} ----------------------------------------------------------------------------------')
                    is_multiband = src.count > 1
                    band_tags = src.tags()
                    if self.dterrengdata:
                        band_name = cst.s2_bands_aliases[image_path[-3::]]
                    for i in range(1, src.count + 1):
                        if self.dterrengdata:
                            band_name = cst.s2_bands_aliases[image_path[-3::]]
                        if subitem.is_file():
                            band_name = image_path.stem[-3:]
                            if band_name in valid_band_names:
                                if band_name.startswith('B0'):
                                    band_name = cst.s2_bands_aliases[f'{band_name}']
                            else:
                                continue
                        if subitem.is_dir():
                            band_name = image_path.stem[-7:-4]
                            if band_name in valid_band_names:
                                #print(band_name)
                                if band_name.startswith('B'):
                                    band_name = cst.s2_bands_aliases[f'{band_name}']
                            else:
                                continue

                        if not band_name:
                            continue

                        # === Write raster band ===
                        data = src.read(i)
                        if data.shape != (ny, nx):
                            scale_y = ny / data.shape[0]
                            scale_x = nx / data.shape[1]
                            data = scipy.ndimage.zoom(data, (scale_y, scale_x), order=0)
                            logger.debug(f"Resampled {band_name} from {data.shape} to {(ny, nx)}")
                                                    # Check final shape before writing
                            assert data.shape == (ny, nx), f"{band_name}: expected {(ny, nx)}, got {src.shape}"
                        
                        if band_name in ncout.variables:
                            logger.warning(f"Higher resolution of '{band_name}' already exists in NetCDF. Skipping.")
                            continue
                        varout = ncout.createVariable(
                            band_name, np.uint16, ('time', 'y', 'x'),
                            zlib=True, complevel=compression_level, fill_value=0
                        )
                        varout.units = "1"
                        varout.grid_mapping = "UTM_projection"
                        varout.long_name = f"Reflectance in band {band_name}"
                        varout._Unsigned = "true"
                        varout.standard_name = (
                            'surface_bidirectional_reflectance'
                            if self.processing_level == 'Level-2A'
                            else 'toa_bidirectional_reflectance'
                        )

                        # Adding if not present
                        for key in ['BANDWIDTH', 'BANDWIDTH_UNIT', 'WAVELENGTH', 'WAVELENGTH_UNIT',
                                    'SOLAR_IRRADIANCE', 'SOLAR_IRRADIANCE_UNIT']:
                            val = band_tags.get(key)
                            if val and val.lower() != 'none':
                                setattr(varout, key.lower(), val)
                        
                        varout[0, :, :] = data

            # Set grid mapping
            ##########################################################
            cf_proj_map = {
                "utm": "transverse_mercator",
                "longlat": "latitude_longitude",
                # Add more as needed
            }

            # # Dummy dimension to allow the variable
            ncout.createDimension("crs", 1)
            #nc_crs = ncout.createVariable('UTM_projection', np.int32)
            rio_crs = src.crs  
            pyproj_crs = PyCRS.from_wkt(rio_crs.to_wkt())  # convert to pyproj CRS
            # Projection strings
            proj_name = pyproj_crs.to_dict().get("proj", "") # hardcoding name of grid_projection to comply with CF
            nc_crs.grid_mapping_name = cf_proj_map.get(proj_name, proj_name)
            nc_crs.crs_wkt = pyproj_crs.to_wkt()
            nc_crs.proj4_string = pyproj_crs.to_proj4()

            # Semi-axes 
            ellipsoid = pyproj_crs.ellipsoid
            nc_crs.semi_major_axis = ellipsoid.semi_major_metre
            nc_crs.semi_minor_axis = ellipsoid.semi_minor_metre

            # Projection parameters (with safe fallback)
            params = pyproj_crs.to_dict()
            nc_crs.false_easting = float(params.get("x_0", 0.0))
            nc_crs.false_northing = float(params.get("y_0", 0.0))
            nc_crs.latitude_of_projection_origin = float(params.get("lat_0", 0.0))
            nc_crs.longitude_of_central_meridian = float(params.get("lon_0", 0.0))
            nc_crs.scale_factor_at_central_meridian = float(params.get("k", 1.0))
            
            # Add sun and view angles
            ##########################################################

            # Status
            logger.info('Adding sun and view angles in native resolution')
            utils.memory_use(self.t0)

            varout_view_azimuth = ncout.createVariable('view_azimuth', np.float32, ('time','raster_band_id', 'ya', 'xa'), fill_value=netCDF4.default_fillvals['f4'],
                                                  zlib=True, complevel=compression_level)
            varout_view_azimuth.units = 'degree'
            varout_view_azimuth.long_name = 'Viewing incidence azimuth angle'
            varout_view_azimuth.comment = 'Original 22x22 pixel resolution'

            varout_view_zenith = ncout.createVariable('view_zenith', np.float32, ('time','raster_band_id', 'ya', 'xa'), fill_value=netCDF4.default_fillvals['f4'],
                                                  zlib=True, complevel=compression_level)
            varout_view_zenith.units = 'degree'
            varout_view_zenith.long_name = 'Viewing incidence zenith angle'
            varout_view_zenith.comment = 'Original 22x22 pixel resolution'

            counter = 1
            for k, v in list(self.sunAndViewAngles.items()):
                logger.debug(("Handeling %i of %i" % (counter, len(self.sunAndViewAngles))))

                if 'sun' in k:
                    varout = ncout.createVariable(k, np.float32, ('time', 'ya', 'xa'), fill_value=netCDF4.default_fillvals['f4'],
                                                  zlib=True, complevel=compression_level)
                    varout.units = 'degree'
                    varout.long_name = 'Solar %s angle' % k.split('_')[-1]
                    varout.comment = 'Original 22x22 pixel resolution'
                    varout[0, :, :] = v
                elif 'zenith' in k :
                    band_id = k.split('_')[-1]
                    varout_view_zenith[0,utils.get_key(cst.s2_bands_order,band_id), :, :] = v
                elif 'azimuth' in k :
                    band_id = k.split('_')[-1]
                    varout_view_azimuth[0,utils.get_key(cst.s2_bands_order,band_id), :, :] = v

                counter += 1


            # Add specific Level-1C or Level-2A layers
            ##########################################################

            # Mapping GDAL types to NumPy types (rasterio returns np arrays, so you only need dtype)
            gdal_nc_data_types = {'Byte': 'u1', 'UInt16': 'u2'}

            # Choose layers to include
            specific_layers = {}
            if self.processing_level == 'Level-1C':
                logger.info('Adding Level-1C specific layers')
                for layer_key in cst.s2_l1c_layers:
                    for k, v in self.imageFiles.items():
                        if layer_key in k or layer_key in str(v):
                            specific_layers[k] = cst.s2_l1c_layers[layer_key]
            elif self.processing_level == 'Level-2A':
                logger.info('Adding Level-2A specific layers')
                for layer_key in cst.s2_l2a_layers:
                    for k, v in self.imageFiles.items():
                        if layer_key in k or layer_key in str(v):
                            specific_layers[k] = cst.s2_l2a_layers[layer_key]

            # Process each auxiliary layer
            for k, v in specific_layers.items():
                image_path = Path(self.imageFiles[k])
                var_name, long_name = v.split(',')

                with rasterio.open(image_path) as src:
                    nb_bands = src.count
                    data_type_name = src.dtypes[0]  # e.g. 'uint16', 'uint8'
                    data_type = data_type_name if data_type_name in ['u1', 'u2'] else 'u2'  # fallback

                    if nb_bands > 1:
                        logger.info(f"Raster data contains more than one layer")
                    
                    for i in range(1, nb_bands + 1):
                        if nb_bands > 1:
                            var_name = v.split(',')[0].split()[i - 1]
                            long_name = v.split(',')[1].split('-')[i - 1]

                        data = src.read(i)
                        NDV = src.nodata or 0

                        # Resample if not 10m resolution
                        if data.shape != (ny, nx):
                            scale_y = ny / data.shape[0]
                            scale_x = nx / data.shape[1]
                            data = scipy.ndimage.zoom(data, (scale_y, scale_x), order=0)
                            logger.debug(f"Resampled {var_name} to shape {(ny, nx)}")
                        
                        if var_name in ncout.variables:
                            logger.warning(f"Skipping duplicate variable '{var_name}'")
                            continue

                        varout = ncout.createVariable(
                            var_name, data.dtype.str[1:], ('time', 'y', 'x'),
                            zlib=True, complevel=compression_level, fill_value=NDV,
                            chunksizes=chunk_size
                        )
                        varout.grid_mapping = "UTM_projection"
                        varout.long_name = long_name

                        # Add classification flags (e.g. for SCL)
                        if var_name == "SCL":
                            varout.flag_values = np.array(
                                list(cst.s2_scene_classification_flags.values()), dtype=np.int8)
                            varout.flag_meanings = ' '.join(cst.s2_scene_classification_flags.keys())

                        varout[0, :, :] = data


            
            # Add orbit specific data
            ##########################################################
            # Status
            logger.info('Adding satellite orbit specific data')
            utils.memory_use(self.t0)

            root = utils.xml_read(self.mainXML)
            if not self.dterrengdata:
                self.globalAttribs['orbitNumber'] = root.find('.//safe:orbitNumber',
                                                              namespaces=root.nsmap).text
            else:
                self.globalAttribs['orbitNumber'] = root.find('.//SENSING_ORBIT_NUMBER').text


            ncout.createDimension('orbit_dim', 3)
            nc_orb = ncout.createVariable('orbit_data', np.int32, ('time', 'orbit_dim'))
            nc_orb[0, :] = [
                int(self.globalAttribs['DATATAKE_1_SENSING_ORBIT_NUMBER']),
                int(self.globalAttribs['orbitNumber']),
                cst.platform_id[self.globalAttribs['DATATAKE_1_SPACECRAFT_NAME']]
            ]

            # Add global attributes
            logger.info('Adding global attributes')
            utils.memory_use(self.t0)
            utils.get_global_attributes(self)
            ncout.setncatts(self.globalAttribs)
            ncout.sync()

        logger.info("Finished conversion.")

        return out_netcdf.is_file()
    

    def readSunAndViewAngles(self, xmlfile):
        """ Method for reading sun and view angles from Sentinel-2
            annotation files.
        """

        root = utils.xml_read(xmlfile)

        angles_view_list = root.findall('.//Tile_Angles')[0]
        angle_step = float(root.findall('.//COL_STEP')[0].text)  # m
        col_step = float(root.findall('.//ROW_STEP')[0].text)  # m
        nx = int(root.xpath(str(
            '//n1:{}_Tile_ID/n1:Geometric_Info/Tile_Geocoding/Size[@resolution=10]/NROWS'.format(
                self.processing_level)), namespaces=root.nsmap)[0].text)  # nb of rows
        ny = int(root.xpath(str(
            '//n1:{}_Tile_ID/n1:Geometric_Info/Tile_Geocoding/Size[@resolution=10]/NCOLS'.format(
                self.processing_level)), namespaces=root.nsmap)[0].text)  # nb of columns
        spatial_resolution = 10

        angle_len = int(math.ceil(nx * spatial_resolution / angle_step))
        sun_zenith = np.zeros((angle_len, angle_len), dtype=np.float32)
        sun_azimuth = np.zeros((angle_len, angle_len), dtype=np.float32)
        angle_step = int(math.ceil(nx / float(angle_len)))
        incidence_angles_list = angles_view_list.findall('Viewing_Incidence_Angles_Grids')

        # Sun angles
        for angle in angles_view_list.find('Sun_Angles_Grid'):
            # print('\t',angle.tag)
            counter_entry = 0
            values_list = angle.find('Values_List')
            for value_entry in values_list[0:-1]:
                if angle.tag == 'Zenith':
                    tmp_sun = np.array([float(i) for i in value_entry.text.split()])[0:-1]
                    sun_zenith[counter_entry, :] = tmp_sun
                    counter_entry += 1
                if angle.tag == 'Azimuth':
                    tmp_sun = np.array([float(i) for i in value_entry.text.split()])[0:-1]
                    sun_azimuth[counter_entry, :] = tmp_sun
                    counter_entry += 1

        self.sunAndViewAngles['sun_zenith'] = sun_zenith
        self.sunAndViewAngles['sun_azimuth'] = sun_azimuth

        # View angles
        counter_angle = 0
        for BANDID in np.array(list(cst.s2_bands_order.keys())):
            tmp_view_zenith = np.zeros((angle_len, angle_len), dtype=np.float32)
            tmp_view_azimuth = np.zeros((angle_len, angle_len), dtype=np.float32)
            tmp_view_zenith[:] = np.nan
            tmp_view_azimuth[:] = np.nan
            for incidence_angles in incidence_angles_list:
                if int(incidence_angles.attrib['bandId']) == BANDID:
                    # print('\t',incidence_angles.attrib)
                    for angle in incidence_angles:
                        values_list = angle.find('Values_List')
                        counter_entry = 0
                        for value_entry in values_list[0:-1]:
                            if angle.tag == 'Zenith':
                                tmp_angle = np.array([float(i) for i in value_entry.text.split()])[
                                            0:-1]
                                tmp_view_zenith[counter_entry, np.isnan(tmp_angle) == False] = \
                                tmp_angle[np.isnan(tmp_angle) == False]
                                counter_entry += 1
                            if angle.tag == 'Azimuth':
                                tmp_angle = np.array([float(i) for i in value_entry.text.split()])[
                                            0:-1]
                                tmp_view_azimuth[counter_entry, np.isnan(tmp_angle) == False] = \
                                tmp_angle[np.isnan(tmp_angle) == False]
                                counter_entry += 1
                    counter_angle += 1
                self.sunAndViewAngles[
                    str('view_zenith_' + cst.s2_bands_order[BANDID])] = tmp_view_zenith
                self.sunAndViewAngles[
                    str('view_azimuth_' + cst.s2_bands_order[BANDID])] = tmp_view_azimuth

    def resample_angles(self, angles, new_dim, angles_length, angles_height, step, type=np.float32):
        ''' Resample angles to get 1-1 with original output.
            angles: numpy array
            new_dim: new dimension (one number, assumes quadratic)
            angles_length: nb. columns in angles array
            angles_height: nb. rows in angles array
            step: stepsize for new dimension
            type: numpy dtype. float32 default
            '''
        angles_resampled = np.zeros((new_dim, new_dim), dtype=type)
        for i in range(angles_length):
            for j in range(angles_height):
                if not i == angles_length - 1 and not j == angles_height - 1:
                    angles_resampled[i * step:i * step + step, j * step:j * step + step] = angles[
                        i, j]
                else:
                    angles_resampled[i * step:new_dim, j * step:new_dim] = angles[i, j]
        return angles_resampled

    def genLatLon(self, nx, ny, latlon=True):
        """Generate latitude/longitude or projection coordinates using rasterio + pyproj."""
        
        # Get affine transform from rasterio
        transform = self.reference_band.transform  # Affine transform
        ulx, xres, xskew, uly, yskew, yres = (
            transform.c,
            transform.a,
            transform.b,
            transform.f,
            transform.d,
            transform.e
        )

        # UTM coordinates (projected)
        xnp = np.arange(nx) * xres + ulx
        ynp = np.arange(ny) * yres + uly

        if not latlon:
            return xnp, ynp

        # Create coordinate grid in UTM
        xx, yy = np.meshgrid(xnp + xres / 2, ynp + yres / 2)  # center of pixels

        # Reproject UTM to WGS84 lat/lon
        transformer = pyproj.Transformer.from_crs(
            self.reference_band.crs,
            "EPSG:4326",  # WGS84
            always_xy=True
        )
        lon, lat = transformer.transform(xx, yy)

        return lat, lon


    def write_vector(self, vectorfile, ncfile):
        """

        Write content of a shapefile in netcdf - following CF 1.8 convention
        Input file can be: any vector-based spatial data (including shape, gml, geojson, ...).
        Input data can be: polygon(s) with no holes.

        Args:
            vectorfile [pathlib.Path]: input vector file
            ncfile     [pathlib.Path]: output nc file, already existing and open for writing
        Returns: N/A

        """

        # -------------------------------
        #    Read input information
        # -------------------------------

        # Read data from vector file. Exits if no data found.
        try:
            src = geopd.read_file(vectorfile)
        except ValueError:
            print(f'No data in {vectorfile}')
            return

        # - nb of shapes within file
        nGeometries = src.shape[0]

        # - list of nodes for each shape
        nodes = src['geometry'].apply(lambda x: x.exterior.coords)

        # - nb of nodes for each shape
        nNodesPerGeom = [len(x) for x in nodes]

        # - x, y, z coordinates for each node of each shape
        flat_nodes = [y for x in nodes for y in x]
        try:
            x, y = zip(*flat_nodes)
        except ValueError:
            x, y, z = zip(*flat_nodes)

        # Variable name in output nc file
        name = vectorfile.stem
        if name == "MSK_CLOUDS_B00":
            name = 'Clouds'

        # Read actual vector information
        flags = src.gml_id
        if name == "Clouds":
            geom_values = flags.apply(lambda x: int(x.split('.')[-1]) + 1)
        else:
            geom_values = flags.apply(lambda x: int(x[-1]) + 1)

        # -------------------------------
        #    Write to nc file
        # -------------------------------

        # Write shape data in specific group
        #shapes = ncfile.createGroup('masks')

        # Add new dimensions
        ncfile.createDimension(f'instance_{name}', nGeometries)
        ncfile.createDimension(f'node_{name}', sum(nNodesPerGeom))

        # Add new variables:

        # - geometry container
        geom = ncfile.createVariable(f'geometry_container_{name}', 'i4')
        geom.geometry_type = "polygon"
        geom.node_count = f'node_count_{name}'
        geom.node_coordinates = f'x_node_{name} y_node_{name}'    # variables containing spatial
        geom.grid_mapping = "UTM_projection"

        # - node count
        nodecount = ncfile.createVariable(f'node_count_{name}', 'i4', f'instance_{name}')
        nodecount.long_name = "count of coordinates in each instance geometry"
        nodecount[:] = nNodesPerGeom

        # - y coordinates
        ncynode = ncfile.createVariable(f'y_node_{name}', 'i4', f'node_{name}', zlib=True)
        ncynode.units = 'm'
        ncynode.standard_name = 'projection_y_coordinate'
        ncynode.axis = 'Y'
        ncynode[:] = y

        # - x coordinates
        #todo: check what standard_name should those have?
        # is it ok to have several x y coordinates?
        ncxnode = ncfile.createVariable(f'x_node_{name}', 'i4', f'node_{name}', zlib=True)
        ncxnode.units = 'm'
        ncxnode.standard_name = 'projection_x_coordinate'
        ncxnode.axis = 'X'
        ncxnode[:] = x

        # - vector information
        varout = ncfile.createVariable(name, 'byte', ('time', f'instance_{name}'))
        # todo: add a comment? before it was "Rasterized vector information."
        varout.long_name = f"{name} mask 10m resolution"
        varout.grid_mapping = "UTM_projection"
        varout.geometry = f'geometry_container_{name}'
        varout.flag_values = geom_values.values.astype('b')
        varout.flag_meanings = ' '.join(flags)
        varout[0, :] = geom_values

        return
    

