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
# Author(s):     Trygve Halsne
# Created:
# Modifications:
# Copyright:     (c) Norwegian Meteorological Institute, 2018
#
# Need to use gdal 2.1.1-> to have support of the SAFE reader

import pathlib
import math
from collections import defaultdict
from datetime import datetime
import lxml.etree as ET
import netCDF4
import numpy as np
import osgeo.ogr as ogr
import osgeo.osr as osr
import pyproj
import scipy.ndimage
from osgeo import gdal
import geopandas as geopd
import safe_to_netcdf.utils as utils
import safe_to_netcdf.constants as cst
import os
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

    def __init__(self, product, indir, outdir):
        self.product_id = product
        self.input_zip = (indir / product).with_suffix('.zip')
        self.SAFE_dir = (outdir / self.product_id).with_suffix('.SAFE')
        self.processing_level = 'Level-' + self.product_id.split('_')[1][4:6]
        self.xmlFiles = defaultdict(list)
        self.imageFiles = defaultdict(list)
        self.globalAttribs = {}
        self.src = None
        self.t0 = datetime.now()
        self.ncout = None  # NetCDF output file
        self.reference_band = None
        self.dterrengdata = False  # variable saying if products is Norwegian DEM L1C
        self.sunAndViewAngles = defaultdict(list)
        self.vectorInformation = defaultdict(list)
        self.SAFE_structure = None
        self.image_list_dterreng = []

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
        if not self.dterrengdata:
            currXml = self.xmlFiles['S2_{}_Tile1_Metadata'.format(self.processing_level)]
        else:
            currXml = self.xmlFiles['MTD_TL']
        self.readSunAndViewAngles(currXml)

        # 5) Retrieve SAFE product structure
        # much difficulty afterwards to be able to save this to netCDF
        ##self.SAFE_structure = zipfile.ZipFile(self.input_zip).namelist()
        self.SAFE_structure = self.list_product_structure()

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

        nx = self.reference_band.RasterXSize  # number of pixels for 10m spatial resolution
        # frequency bands
        ny = self.reference_band.RasterYSize  # number of pixels for 10m spatial resolution
        # frequency bands

        # output filename
        out_netcdf = (nc_outpath / self.product_id).with_suffix('.nc')

        with (netCDF4.Dataset(out_netcdf, 'w', format='NETCDF4')) as ncout:
            ncout.createDimension('time', 0)
            ncout.createDimension('x', nx)
            ncout.createDimension('y', ny)

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
                    ncout.createDimension('dimension_rgb', subdataset.RasterCount)
                    varout = ncout.createVariable('TCI', 'u1', ('time', 'dimension_rgb', 'y', 'x'),
                                                  fill_value=0, zlib=True, complevel=compression_level)
                    varout.units = "1"
                    varout.grid_mapping = "UTM_projection"
                    varout.long_name = 'TCI RGB from B4, B3 and B2'
                    varout._Unsigned = "true"
                    for i in range(1, subdataset.RasterCount + 1):
                        current_band = subdataset.GetRasterBand(i)
                        band_measurement = current_band.GetVirtualMemArray()
                        varout[0, i - 1, :, :] = band_measurement
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
            nc_crs = ncout.createVariable('UTM_projection', np.int32, ('time'))
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
            logger.info('Adding vector layers')
            utils.memory_use(self.t0)

            for gmlfile in self.xmlFiles.values():
                if gmlfile and gmlfile.suffix == '.gml':
                    self.write_vector(gmlfile, ncout)

            # Add Level-2A layers
            ##########################################################
            # Status
            if self.processing_level == 'Level-2A':
                logger.info('Adding Level-2A specific layers')
                utils.memory_use(self.t0)
                gdal_nc_data_types = {'Byte': 'u1', 'UInt16': 'u2'}
                l2a_kv = {}
                for layer in cst.s2_l2a_layers:
                    for k, v in self.imageFiles.items():
                        if layer in k:
                            l2a_kv[k] = cst.s2_l2a_layers[k]
                        elif layer in str(v):
                            logger.debug((layer, str(v), k))
                            l2a_kv[k] = cst.s2_l2a_layers[layer]

                for k, v in l2a_kv.items():
                    logger.debug((k, v))
                    varName, longName = v.split(',')
                    SourceDS = gdal.Open(str(self.imageFiles[k]), gdal.GA_ReadOnly)
                    if SourceDS.RasterCount > 1:
                        logger.info("Raster data contains more than one layer")
                    NDV = SourceDS.GetRasterBand(1).GetNoDataValue()
                    xsize = SourceDS.RasterXSize
                    ysize = SourceDS.RasterYSize
                    GeoT = SourceDS.GetGeoTransform()
                    DataType = gdal_nc_data_types[
                        gdal.GetDataTypeName(SourceDS.GetRasterBand(1).DataType)]
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
                        raster_data = scipy.ndimage.zoom(input=SourceDS.GetVirtualMemArray(),
                                                         zoom=nx / xsize, order=0)
                    else:
                        raster_data = SourceDS.GetVirtualMemArray()
                    varout[0, :] = raster_data

            # Add sun and view angles
            ##########################################################
            # Status
            logger.info('Adding sun and view angles')
            utils.memory_use(self.t0)

            counter = 1
            for k, v in list(self.sunAndViewAngles.items()):
                logger.debug(("Handeling %i of %i" % (counter, len(self.sunAndViewAngles))))
                angle_step = int(math.ceil(nx / float(v.shape[0])))

                resampled_angles = self.resample_angles(v, nx, v.shape[0], v.shape[1], angle_step,
                                                        type=np.float32)

                varout = ncout.createVariable(k, np.float32, ('time', 'y', 'x'), fill_value=netCDF4.default_fillvals['f4'],
                                                  zlib=True, complevel=compression_level)
                varout.units = 'degree'
                if 'sun' in k:
                    varout.long_name = 'Solar %s angle' % k.split('_')[-1]
                else:
                    varout.long_name = 'Viewing incidence %s angle' % k.split('_')[1]
                varout.grid_mapping = "UTM_projection"
                varout.comment = '1 to 1 with original 22x22 resolution'
                varout[0, :, :] = resampled_angles
                counter += 1

            # Add xml files as character values see:
            # https://stackoverflow.com/questions/37079883/string-handling-in-python-netcdf4
            ##########################################################
            # Status
            logger.info('Adding XML files as character variables')
            utils.memory_use(self.t0)

            for k, xmlfile in self.xmlFiles.items():
                if xmlfile and xmlfile.suffix == '.xml':
                        xmlString = self.xmlToString(xmlfile)

                        if xmlString:
                            dim_name = str('dimension_' + k.replace('-', '_'))
                            ncout.createDimension(dim_name, len(xmlString))
                            msg_var = ncout.createVariable(k.replace('-', '_'), 'S1', dim_name)
                            msg_var.long_name = str("SAFE xml file: " + k)
                            msg_var.comment = "Original SAFE xml file added as character values."
                            # todo DeprecationWarning: tostring() is deprecated. Use tobytes()
                            # instead.
                            msg_var[:] = netCDF4.stringtochar(np.array([xmlString], 'S'))

            # Add SAFE product structure as character values
            ##########################################################
            # Status
            logger.info('Adding SAFE product structure as character variable')
            if self.SAFE_structure:
                dim_name = str('dimension_SAFE_structure')
                ncout.createDimension(dim_name, len(self.SAFE_structure))
                msg_var = ncout.createVariable("SAFE_structure", 'S1', dim_name)
                msg_var.comment = "Original SAFE product structure xml file as character values."
                msg_var.long_name = "Original SAFE product structure."
                msg_var[:] = netCDF4.stringtochar(np.array([self.SAFE_structure], 'S'))

            # Add orbit specific data
            ##########################################################
            # Status
            print('\nAdding satellite orbit specific data')
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
            ##########################################################
            # Status
            logger.info('Adding global attributes')
            utils.memory_use(self.t0)

            nowstr = self.t0.strftime("%Y-%m-%dT%H:%M:%SZ")
            ncout.title = 'Sentinel-2 {} data'.format(self.processing_level)
            ncout.netcdf4_version_id = netCDF4.__netcdf4libversion__
            ncout.file_creation_date = nowstr

            self.globalAttribs[
                'summary'] = 'Sentinel-2 Multi-Spectral Instrument {} product.'.format(
                self.processing_level)
            self.globalAttribs[
                'keywords'] = '[Earth Science, Atmosphere, Atmospheric radiation, Reflectance]'
            self.globalAttribs['keywords_vocabulary'] = "GCMD Science Keywords"
            self.globalAttribs['institution'] = "Norwegian Meteorological Institute"
            self.globalAttribs['history'] = nowstr + ". Converted from SAFE to NetCDF by NBS team."
            self.globalAttribs['source'] = "surface observation"
            self.globalAttribs['relativeOrbitNumber'] = self.globalAttribs.pop(
                'DATATAKE_1_SENSING_ORBIT_NUMBER')
            self.globalAttribs['Conventions'] = "CF-1.8"
            ncout.setncatts(self.globalAttribs)
            ncout.sync()

            # Status
            logger.info('Finished.')
            utils.memory_use(self.t0)

        return out_netcdf.is_file()

    def xmlToString(self, xmlfile):
        """ Method for reading XML files returning the entire file as single
            string.
        """
        if not xmlfile.is_file():
            logger.error(('Error: Can\'t find xmlfile %s' % (xmlfile)))
            return False
        try:
            parser = ET.XMLParser(recover=True)
            tree = ET.parse(str(xmlfile), parser)
            return ET.tostring(tree)
        except:
            logger.info(("Could not parse %s as xmlFile. Try to open regularly." % xmlfile))
            with open(xmlfile, 'r') as infile:
                text = infile.read()
            if text:
                return text
            else:
                logger.error(("Could not parse %s. Something wrong with file." % xmlfile))
                return False

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

        source = osr.SpatialReference()
        source.ImportFromWkt(self.reference_band.GetProjection())
        target = osr.SpatialReference()
        # target.ImportFromEPSG(4326)
        target.ImportFromProj4('+proj=longlat +ellps=WGS84')

        current_projection = pyproj.Proj(source.ExportToProj4())
        target_projection = pyproj.Proj(target.ExportToProj4())
        target2 = current_projection.to_latlong()

        longitude, latitude = pyproj.transform(current_projection, target_projection, xp, yp)

        return latitude, longitude

    def list_product_structure(self):
        """ Traverse SAFE file structure (or any file structure) and
            creates a xml file containing the file structure.

            Returns a string representation of the xml file."""

        startpath = str(self.SAFE_dir)

        root_path = startpath.split('/')[-1]
        ET_root = ET.Element(root_path)
        # Create dictionary that contains all elements
        elements = {root_path: ET_root}

        # Iterate throgh the file structure
        for root, dirs, files in os.walk(startpath):
            level = root.replace(startpath, '').count(os.sep)
            current_xpath = str('/' + root_path + root.replace(startpath, ''))
            current_path = root.replace(startpath, '')

            # Create all folder elements
            if dirs:
                for dir in dirs:
                    element = ET.SubElement(elements[current_xpath.strip('/')], dir)
                    elements[str(current_xpath.strip('/') + '/' + dir)] = element
            # Add all files in the correct folder
            for f in files:
                sub_element = ET.SubElement(elements[current_xpath.strip('/')], 'file')
                sub_element.text = f

        return ET.tostring(ET_root)

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
            values = flags.apply(lambda x: int(x.split('.')[-1]) + 1)
        else:
            values = flags.apply(lambda x: int(x[-1]) + 1)

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
        varout = ncfile.createVariable(name, 'i1', ('time', f'instance_{name}'))
        # todo: update long_name and comment
        #varout.long_name = f"{name} mask 10m resolution"
        varout.grid_mapping = "UTM_projection"
        varout.geometry = f'geometry_container_{name}'
        varout.flag_values = values
        varout.flag_meanings = ' '.join(flags)
        varout[0, :] = values

        return


if __name__ == '__main__':

    workdir = pathlib.Path('/home/elodief/Data/NBS')


    products = ['S2A_MSIL1C_20201028T102141_N0209_R065_T34WDA_20201028T104239',
                'S2A_MSIL1C_20201022T100051_N0202_R122_T35WPU_20201026T035024_DTERRENGDATA',
                'S2A_MSIL2A_20210714T105031_N0301_R051_T32VMK_20210714T135226']

    products = ['S2A_MSIL2A_20210714T105031_N0301_R051_T32VMK_20210714T135226']

    for product in products:

        outdir = workdir / 'NBS_test_data' / 'cf18_04' / product
        outdir.parent.mkdir(parents=False, exist_ok=True)
        conversion_object = Sentinel2_reader_and_NetCDF_converter(
            product=product,
            indir=workdir / 'NBS_reference_data' / 'reference_datain_local',
            #indir=workdir / 'NBS_test_data' / 'cf18_01',
            outdir=outdir)
        conversion_object.write_to_NetCDF(outdir, 7)
