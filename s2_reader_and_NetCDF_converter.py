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

import pathlib
import sys
import math
from collections import defaultdict
import datetime as dt
import netCDF4
import numpy as np
import osgeo.osr as osr
import pyproj
import scipy.ndimage
from osgeo import gdal
import geopandas as geopd
import safe_to_netcdf.utils as utils
import safe_to_netcdf.constants as cst
import logging
gdal.UseExceptions()


logger = logging.getLogger(__name__)


class Sentinel2_reader_and_NetCDF_converter:
    ''' Class for reading Sentinel-2 MSI L1C/L2A products from SAFE with methods for
        reading auxilary information as e.g. clouds, solar and view angles.
        In addition, it is possible to convert product into NetCDF4/CF (1.6).

        The implemented methods uses standard python libraries as
        gdal(v. > 2.1.1), numpy, lxml etc.

        Keyword arguments:
        SAFE_file -- absolute path to zipped file
        SAFE_outpath -- output storage location for unzipped SAFE product
        '''

    def __init__(self, product, indir, outdir, colhub_uuid=None):
        self.uuid = colhub_uuid
        self.product_id = product
        self.input_zip = (indir / product).with_suffix('.zip')
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

        self.main()

    def main(self):
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

    def write_to_NetCDF(self, nc_outpath, compression_level, chunk_size=(1, 32, 32)):
        """ Method writing output NetCDF product.

        Keyword arguments:
        nc_outpath -- output path where NetCDF file should be stored
        compression_level -- compression level on output NetCDF file (1-9)
        chunk_size -- chunk_size
        """

        logger.info("------------START CONVERSION FROM SAFE TO NETCDF-------------")

        # Status
        logger.info('Creating NetCDF file')
        utils.memory_use(self.t0)

        # Deciding a reference band
        #todo dterreng warning coming from here?
        # yes -> self.src.GetSubDatasets() ok but the gdal.Open does not work
        # add break? how to remove warning from dterr?
        for k, v in self.src.GetSubDatasets():
            if v.find('10m') > 0:
                self.reference_band = gdal.Open(k)

        # frequency bands
        nx = self.reference_band.RasterXSize  # number of pixels for 10m spatial resolution
        ny = self.reference_band.RasterYSize  # number of pixels for 10m spatial resolution
        
        # sun and view angles raster resolution
        nxa, nya = self.sunAndViewAngles[list(self.sunAndViewAngles)[0]].shape  

        # output filename
        out_netcdf = (nc_outpath / self.product_id).with_suffix('.nc')

        with (netCDF4.Dataset(out_netcdf, 'w', format='NETCDF4')) as ncout:
            ncout.createDimension('time', 0)
            ncout.createDimension('x', nx)
            ncout.createDimension('y', ny)
            ncout.createDimension('raster_band_id', len(cst.s2_bands_order.keys()))
            ncout.createDimension('xa', nxa)
            ncout.createDimension('ya', nya)

            utils.create_time(ncout, self.globalAttribs["PRODUCT_START_TIME"])

            # Add projection coordinates
            ##########################################################
            # Status
            logger.info('Adding projection coordinates')
            utils.memory_use(self.t0)

            xnp, ynp = self.genLatLon(nx, ny, latlon=False)  # Assume gcps are on a regular grid

            ncx = ncout.createVariable('x', 'i4', 'x', zlib=True, complevel=compression_level)
            ncx.units = 'm'
            ncx.standard_name = 'projection_x_coordinate'
            ncx[:] = xnp

            ncy = ncout.createVariable('y', 'i4', 'y', zlib=True, complevel=compression_level)
            ncy.units = 'm'
            ncy.standard_name = 'projection_y_coordinate'
            ncy[:] = ynp
            
            # Add projection raster band id variable
            ##########################################################
            nc_rasterband_id = ncout.createVariable('band_id', 'i4', 'raster_band_id', zlib=True, complevel=compression_level)
            nc_rasterband_id.long_name = 'raster band id'
            nc_rasterband_id[:] = np.array(list(cst.s2_bands_order.keys()))
            nc_rasterband_id.flag_values = np.array(list(
                            cst.s2_bands_order.keys()),
                                                      dtype=np.int8)
            nc_rasterband_id.flag_meanings = ' '.join(
                            [value for value in list(cst.s2_bands_order.values())])
            

            # Add raw measurement layers
            # Currently adding TCI
            # NODATA = 0 (ie. fillvalue) from
            # https://sentinel.esa.int/documents/247904/685211/Sentinel-2-Products-Specification
            # -Document
            ##########################################################
            # Status
            logger.info('Adding frequency bands layers')
            utils.memory_use(self.t0)

            if self.dterrengdata:
                # For DTERR data, gdal fails to properly do the src.GetSubDatasets()
                # so manually read the list of images created beforehand
                images = [[str(i), i.stem] for i in self.image_list_dterreng]
            else:
                images = self.src.GetSubDatasets()
            for k, v in images:
                subdataset = gdal.Open(k)
                subdataset_geotransform = subdataset.GetGeoTransform()
                # True color image (8 bit true color image)
                if ("True color image" in v) or ('TCI' in v):
                    continue
                # Reflectance data for each band
                else:
                    for i in range(1, subdataset.RasterCount + 1):
                        current_band = subdataset.GetRasterBand(i)
                        if self.dterrengdata:
                            band_metadata = None
                            varName = cst.s2_bands_aliases[v[-3::]]
                        else:
                            band_metadata = current_band.GetMetadata()
                            varName = band_metadata['BANDNAME']
                        if varName.startswith('B'):
                            varout = ncout.createVariable(varName, np.uint16, ('time', 'y', 'x'), fill_value=0,
                                                          zlib=True, complevel=compression_level)
                            varout.units = "1"
                            varout.grid_mapping = "UTM_projection"
                            if self.processing_level == 'Level-2A':
                                varout.standard_name = 'surface_bidirectional_reflectance'
                            else:
                                varout.standard_name = 'toa_bidirectional_reflectance'
                            varout.long_name = 'Reflectance in band %s' % varName
                            if band_metadata:
                                varout.bandwidth = band_metadata['BANDWIDTH']
                                varout.bandwidth_unit = band_metadata['BANDWIDTH_UNIT']
                                varout.wavelength = band_metadata['WAVELENGTH']
                                varout.wavelength_unit = band_metadata['WAVELENGTH_UNIT']
                                varout.solar_irradiance = band_metadata['SOLAR_IRRADIANCE']
                                varout.solar_irradiance_unit = band_metadata['SOLAR_IRRADIANCE_UNIT']
                            varout._Unsigned = "true"
                            # from DN to reflectance
                            logger.debug((varName, subdataset_geotransform))
                            if subdataset_geotransform[1] != 10:
                                current_size = current_band.XSize
                                band_measurement = scipy.ndimage.zoom(
                                    input=current_band.GetVirtualMemArray(), zoom=nx / current_size,
                                    order=0)
                            else:
                                band_measurement = current_band.GetVirtualMemArray()
                            varout[0, :, :] = band_measurement

            # set grid mapping
            ##########################################################
            source_crs = osr.SpatialReference()
            source_crs.ImportFromWkt(self.reference_band.GetProjection())
            nc_crs = ncout.createVariable('UTM_projection', np.int32)
            nc_crs.latitude_of_projection_origin = source_crs.GetProjParm('latitude_of_origin')
            nc_crs.proj4_string = source_crs.ExportToProj4()
            nc_crs.crs_wkt = source_crs.ExportToWkt()
            nc_crs.semi_major_axis = source_crs.GetSemiMajor()
            nc_crs.scale_factor_at_central_meridian = source_crs.GetProjParm('scale_factor')
            nc_crs.longitude_of_central_meridian = source_crs.GetProjParm('central_meridian')
            nc_crs.grid_mapping_name = source_crs.GetAttrValue('PROJECTION').lower()
            nc_crs.semi_minor_axis = source_crs.GetSemiMinor()
            nc_crs.false_easting = source_crs.GetProjParm('false_easting')
            nc_crs.false_northing = source_crs.GetProjParm('false_northing')
            nc_crs.epsg_code = source_crs.GetAttrValue('AUTHORITY', 1)
            nc_crs.crs_wkt = self.reference_band.GetProjection()

            # Add vector layers
            ##########################################################
            # Status
            logger.info('Adding masks layers')
            utils.memory_use(self.t0)
            gdal_nc_data_types = {'Byte': 'u1', 'UInt16': 'u2'}

            # Add specific Level-1C or Level-2A layers
            ##########################################################
            # Status
            specific_layers_kv = {}
            utils.memory_use(self.t0)
            gdal_nc_data_types = {'Byte': 'u1', 'UInt16': 'u2'}

            if self.processing_level == 'Level-1C':
                logger.info('Adding Level-1C specific layers')
                for layer in cst.s2_l1c_layers:
                    for k, v in self.imageFiles.items():
                        if layer in k:
                            specific_layers_kv[k] = cst.s2_l1c_layers[k]
                        elif layer in str(v):
                            logger.debug((layer, str(v), k))
                            specific_layers_kv[k] = cst.s2_l1c_layers[layer]
                

            elif self.processing_level == 'Level-2A':
                logger.info('Adding Level-2A specific layers')
                for layer in cst.s2_l2a_layers:
                    for k, v in self.imageFiles.items():
                        if layer in k:
                            specific_layers_kv[k] = cst.s2_l2a_layers[k]
                        elif layer in str(v):
                            logger.debug((layer, str(v), k))
                            specific_layers_kv[k] = cst.s2_l2a_layers[layer]

            for k, v in specific_layers_kv.items():
                logger.debug((k, v))
                varName, longName = v.split(',')
                SourceDS = gdal.Open(str(self.imageFiles[k]), gdal.GA_ReadOnly)
                nb_rasterBands =  SourceDS.RasterCount 
                    
                if SourceDS.RasterCount > 1:
                    logger.info("Raster data contains more than one layer")
                         
                for i in range(1,nb_rasterBands+1):
                    if nb_rasterBands>1:
                        varName =  v.split(',')[0].split()[i-1]
                        longName =  v.split(',')[1].split('-')[i-1]
                    NDV = SourceDS.GetRasterBand(i).GetNoDataValue()
                    xsize = SourceDS.RasterXSize
                    ysize = SourceDS.RasterYSize
                    GeoT = SourceDS.GetGeoTransform()
                    logger.info("{}".format(GeoT))
                    DataType = gdal_nc_data_types[
                        gdal.GetDataTypeName(SourceDS.GetRasterBand(i).DataType)]
                    varout = ncout.createVariable(varName, DataType, ('time', 'y', 'x'), fill_value=0,
                                                  zlib=True, complevel=compression_level, chunksizes=chunk_size)
                    varout.grid_mapping = "UTM_projection"
                    varout.long_name = longName
                    if varName == "SCL":
                        varout.flag_values = np.array(list(
                            cst.s2_scene_classification_flags.values()),
                                                      dtype=np.int8)
                        varout.flag_meanings = ' '.join(
                            [key for key in list(cst.s2_scene_classification_flags.keys())])
    
                    if GeoT[1] != 10:
                        raster_data = scipy.ndimage.zoom(input=SourceDS.GetRasterBand(i).GetVirtualMemArray(),
                                                         zoom=nx / xsize, order=0)
                    else:
                        raster_data = SourceDS.GetRasterBand(i).GetVirtualMemArray()
                    varout[0, :] = raster_data

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
            rel_orb_nb = self.globalAttribs['DATATAKE_1_SENSING_ORBIT_NUMBER']
            orb_nb = self.globalAttribs['orbitNumber']
            orb_dir = self.globalAttribs['DATATAKE_1_SENSING_ORBIT_DIRECTION']
            platform = self.globalAttribs['DATATAKE_1_SPACECRAFT_NAME']

            nc_orb.relativeOrbitNumber = rel_orb_nb
            nc_orb.orbitNumber = orb_nb
            nc_orb.orbitDirection = orb_dir
            nc_orb.platform = platform
            nc_orb.description = "Values structured as [relative orbit number, orbit number, " \
                                 "platform]. platform corresponds to 0:Sentinel-2A, 1:Sentinel-2B.."
            nc_orb[0, :] = [int(rel_orb_nb), int(orb_nb), cst.platform_id[platform]]

            # Add global attributes

            logger.info('Adding global attributes')
            utils.memory_use(self.t0)

            utils.get_global_attributes(self)
            ncout.setncatts(self.globalAttribs)
            ncout.sync()

            ### Status
            logger.info('Finished.')
            utils.memory_use(self.t0)

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
        """ Method providing latitude and longitude arrays or projection
            coordinates depending on latlon argument."""

        ulx, xres, xskew, uly, yskew, yres = self.reference_band.GetGeoTransform()  # ulx - upper
        # left x, uly - upper left y

        # x and y in UTM coordinates
        xnp = np.arange(nx) * xres + ulx
        ynp = np.arange(ny) * yres + uly

        if not latlon:
            return xnp, ynp

        # Generate coordinate mesh (UTM) Correct lat lon for center pixel
        indices = np.indices((nx, ny), dtype=np.int32)
        xp = np.int32(ulx + (xres * 0.5)) + indices[1] * np.int32(xres) + indices[1] * np.int32(
            xskew)
        yp = np.int32(uly - (yres * 0.5)) + indices[0] * np.int32(yres) + indices[0] * np.int32(
            yskew)

        current_projection = pyproj.CRS.from_string(self.reference_band.GetProjection())
        target_projection = pyproj.CRS.from_proj4('+proj=longlat +ellps=WGS84')
        longitude, latitude = pyproj.Transformer.from_crs(current_projection, target_projection).transform(xp, yp)

        return latitude, longitude

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


if __name__ == '__main__':

    # Log to console
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    log_info = logging.StreamHandler(sys.stdout)
    log_info.setFormatter(logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s'))
    logger.addHandler(log_info)

    workdir = pathlib.Path('/lustre/storeB/project/NBS2/sentinel/production/NorwAREA/netCDFNBS_work/test_environment/test_s2_N0400_updated')
    workdir = pathlib.Path('/lustre/storeA/users/elodief/NBS_test_data/new_s2_02')
    workdir = pathlib.Path('/home/elodief/Data/NBS/NBS_test_data/merge')

    products = [
            'S2A_MSIL1C_20220321T140851_N0400_R053_T31XFJ_20220321T162746',
            'S2A_MSIL2A_20220321T104731_N0400_R051_T32VNN_20220321T145616'

    ]
    products = ['S2A_MSIL1C_20201022T100051_N0202_R122_T35WPU_20201026T035024_DTERRENGDATA']

    for product in products:

        outdir = workdir / product
        outdir.parent.mkdir(parents=False, exist_ok=True)
        conversion_object = Sentinel2_reader_and_NetCDF_converter(
            product=product,
            indir=outdir,
            outdir=outdir)
        if conversion_object.read_ok:
            conversion_object.write_to_NetCDF(outdir, 7)



