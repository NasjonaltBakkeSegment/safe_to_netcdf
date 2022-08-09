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
import sys
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
import safe_to_netcdf.utils as utils
import safe_to_netcdf.constants as cst
import os
import logging
import rasterio
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
        elif "_N0207_" in self.product_id:
            currXml = self.xmlFiles['S2_{}_Tile1_Data'.format(self.processing_level)]
        else:
            currXml = self.xmlFiles['S2_{}_Tile1_Metadata'.format(self.processing_level)]
        # Check for both None and empty list
        if currXml is None or not currXml:
            logger.error("xml file not found in SAFE directory. Hence exiting")
            self.read_ok = False
            return False
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
            if self.processing_level == 'Level-2A':
                ncout.createDimension('time', 0)
            else:
                ncout.createDimension('time', 1)
            ncout.createDimension('x', nx)
            ncout.createDimension('y', ny)

            utils.create_time(ncout, self.globalAttribs["PRODUCT_START_TIME"])

            if not self.processing_level == 'Level-2A':
                nclat = ncout.createVariable('lat', 'f4', ('y', 'x',), zlib=True, complevel=compression_level)
                nclon = ncout.createVariable('lon', 'f4', ('y', 'x',), zlib=True, complevel=compression_level)
                lat,lon = self.genLatLon(nx, ny) #Assume gcps are on a regular grid
                nclat.long_name = 'latitude'
                nclat.units = 'degrees_north'
                nclat.standard_name = 'latitude'
                nclat[:,:]=lat

                nclon.long_name = 'longitude'
                nclon.units = 'degrees_east'
                nclon.standard_name = 'longitude'
                nclon[:,:]=lon

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
                    varout = ncout.createVariable('TCI', 'u1', ('dimension_rgb', 'y', 'x'),
                                                  fill_value=0, zlib=True, complevel=compression_level)
                    varout.units = "1"
                    if not self.processing_level == 'Level-2A':
                        varout.coordinates = "lat lon"
                    varout.grid_mapping = "UTM_projection"
                    varout.long_name = 'TCI RGB from B4, B3 and B2'
                    varout._Unsigned = "true"
                    for i in range(1, subdataset.RasterCount + 1):
                        current_band = subdataset.GetRasterBand(i)
                        band_measurement = current_band.GetVirtualMemArray()
                        varout[i - 1, :, :] = band_measurement
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
                        logger.debug(varName)
                        if varName.startswith('B'):
                            if self.processing_level == 'Level-2A':
                                varout = ncout.createVariable(varName, np.uint16, ('time', 'y', 'x'), fill_value=0,
                                                              zlib=True, complevel=compression_level, chunksizes=chunk_size)
                            else:
                                varout = ncout.createVariable(varName, np.uint16, ('time', 'y', 'x'), fill_value=0,
                                                          zlib=True, complevel=compression_level)
                            varout.units = "1"
                            varout.grid_mapping = "UTM_projection"
                            if self.processing_level == 'Level-2A':
                                varout.standard_name = 'surface_bidirectional_reflectance'
                            else:
                                varout.coordinates = "lat lon"
                                varout.standard_name = 'toa_bidirectional_reflectance'
                            varout.long_name = 'Reflectance in band %s' % varName
                            if band_metadata:
                                #try:
                                varout.bandwidth = band_metadata['BANDWIDTH']
                                varout.bandwidth_unit = band_metadata['BANDWIDTH_UNIT']
                                varout.wavelength = band_metadata['WAVELENGTH']
                                varout.wavelength_unit = band_metadata['WAVELENGTH_UNIT']
                                varout.solar_irradiance = band_metadata['SOLAR_IRRADIANCE']
                                varout.solar_irradiance_unit = band_metadata['SOLAR_IRRADIANCE_UNIT']
                                #except KeyError as e:
                                #    print(f'Metadata not found: {e}')
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
                            #print(band_measurement.shape)
                            varout[0, :, :] = band_measurement

            # set grid mapping
            ##########################################################
            source_crs = osr.SpatialReference()
            source_crs.ImportFromWkt(self.reference_band.GetProjection())
            nc_crs = ncout.createVariable('UTM_projection', np.int32)
            nc_crs.latitude_of_projection_origin = source_crs.GetProjParm('latitude_of_origin')
            nc_crs.proj4_string = source_crs.ExportToProj4()
            nc_crs.semi_major_axis = source_crs.GetSemiMajor()
            nc_crs.scale_factor_at_central_meridian = source_crs.GetProjParm('scale_factor')
            nc_crs.longitude_of_central_meridian = source_crs.GetProjParm('central_meridian')
            nc_crs.grid_mapping_name = source_crs.GetAttrValue('PROJECTION').lower()
            nc_crs.semi_minor_axis = source_crs.GetSemiMinor()
            nc_crs.false_easting = source_crs.GetProjParm('false_easting')
            nc_crs.false_northing = source_crs.GetProjParm('false_northing')
            nc_crs.epsg_code = source_crs.GetAttrValue('AUTHORITY', 1)

            # Add vector layers
            ##########################################################
            # Status
            logger.info('Adding masks layers')
            utils.memory_use(self.t0)
            gdal_nc_data_types = {'Byte': 'u1', 'UInt16': 'u2'}

            # Old masks format
            for gmlfile in self.xmlFiles.values():
                if gmlfile and gmlfile.suffix == '.gml':
                    layer = gmlfile.stem
                    rasterized_ok, layer_mask, mask = self.rasterizeVectorLayers(nx, ny, gmlfile)
                    # build transformer, assuming matching coordinate systems.
                    if rasterized_ok:
                        if layer == "MSK_CLOUDS_B00":
                            layer_name = 'Clouds'
                            comment_name = 'cloud'
                        else:
                            layer_name = layer
                            comment_name = 'vector'
                        varout = ncout.createVariable(layer_name, 'i1', ('time', 'y', 'x'), fill_value=-1,
                                                      zlib=True, complevel=compression_level)
                        varout.long_name = f"{layer_name} mask 10m resolution"
                        varout.comment = f"Rasterized {comment_name} information."
                        if not self.processing_level == 'Level-2A':
                             varout.coordinates = "lat lon"
                        varout.grid_mapping = "UTM_projection"
                        varout.flag_values = np.array(list(layer_mask.values()), dtype=np.int8)
                        varout.flag_meanings = ' '.join(
                            [key.replace('-', '_') for key in list(layer_mask.keys())])
                        varout[0, :] = mask

            # Add mask information, new jp2 format
            for k, v in self.imageFiles.items():
                if v.suffix == '.jp2' and (v.stem.startswith('MSK_DETFOO') or v.stem.startswith('MSK_CLASSI')):
                    utils.memory_use(self.t0)
                    logger.debug((k, v))
                    # Read and resample image to 10m resolution
                    src = rasterio.open(v)
                    upscale_factor = nx / src.shape[1]
                    img = src.read(out_shape=(src.count,
                        int(src.height * upscale_factor), int(src.width * upscale_factor)),
                        resampling=rasterio.enums.Resampling.bilinear)
                    # Several masks per file
                    for i in src.indexes:
                        dataType = src.dtypes[i-1]
                        if v.stem.startswith('MSK_CLASSI'):
                            masks = ['OPAQUE', 'CIRRUS', 'SNOW']
                            comment = ['Bit is 1 when pixel is OPAQUE', 'Bit is 1 when pixel is CIRRUS', 'Bit is 1 when pixel is SNOW/ICE']
                        elif v.stem.startswith('MSK_DETFOO'):
                            masks = ['DETECTOR FOOTPRINT']
                            comment = ['4 bits to encode the 12 detectors']
                        elif v.stem.startswith('MSK_QUALIT'):
                            masks = ['ANC_LOST', 'ANC_DEG', 'MSI_LOST', 'MSI_DEG', 'QT_DEFECTIVE_PIXELS', 'QT_NODATA_PIXELS', 'QT_PARTIALLY_CORRECTED_PIXELS', 'QT_SATURATED_PIXELS']
                            comment = ['Bit is 1 when pixel is ANC data lost', 'Bit is 1 when pixel is ANC data degraded',
                                       'Bit is 1 when pixel is MSI data lost', 'Bit is 1 when pixel is MSI data degraded',
                                       'Bit is 1 when pixel is defective', 'Bit is 1 when pixel is NO_DATA',
                                       'Bit is 1 when pixel is PARTIALLY_CORRECTED', 'Bit is 1 when pixel is saturated at L1A or L1B']
                        else:
                           continue
                        varout = ncout.createVariable('_'.join([v.stem, masks[i-1]]), dataType, ('time', 'y', 'x'), fill_value=0,
                                                  zlib=True, complevel=compression_level, chunksizes=chunk_size)
                        # varout.coordinates = "lat lon" ;
                        varout.grid_mapping = "UTM_projection"
                        varout.long_name = f'{masks[i-1]} from {k}'
                        varout.comment = comment[i-1]
                        varout[0, :, :] = img[i-1, :, :]
                        if not self.processing_level == 'Level-2A':
                            varout.coordinates = "lat lon"

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
                    # print(NDV, xsize, ysize, GeoT, DataType)

                    varout = ncout.createVariable(varName, DataType, ('time', 'y', 'x'), fill_value=0,
                                                  zlib=True, complevel=compression_level, chunksizes=chunk_size)
                    # varout.coordinates = "lat lon" ;
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

                if self.processing_level == 'Level-2A':
                    varout = ncout.createVariable(k, np.float32, ('time', 'y', 'x'), fill_value=netCDF4.default_fillvals['f4'],
                                              zlib=True, complevel=compression_level, chunksizes=chunk_size)
                else:
                    varout = ncout.createVariable(k, np.float32, ('time', 'y', 'x'), fill_value=netCDF4.default_fillvals['f4'],
                                                  zlib=True, complevel=compression_level)
                varout.units = 'degree'
                if 'sun' in k:
                    varout.long_name = 'Solar %s angle' % k.split('_')[-1]
                else:
                    varout.long_name = 'Viewing incidence %s angle' % k.split('_')[1]
                varout.coordinates = 'lat lon'
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

            if self.processing_level == 'Level-2A':
                # Add orbit specific data
                ##########################################################
                # Status
                logger.info('Adding satellite orbit specific data')

                platform_id = {"Sentinel-2A":0, "Sentinel-2B":1,
                               "Sentinel-2C":2,"Sentinel-2D":3,}
                #orb_dir_id = {"DESCENDING":0, "":1,
                #print(self.xmlFiles)

                #if self.xmlFiles['manifest']:
                #tree = ET.parse(self.SAFE_dir+'/manifest.safe')
                tree = ET.parse(str(self.SAFE_dir/'manifest.safe'))
                root = tree.getroot()
                self.globalAttribs['orbitNumber'] = root.find('.//safe:orbitNumber',namespaces=root.nsmap).text

                dim_orb = ncout.createDimension('orbit_dim',3)
                nc_orb = ncout.createVariable('orbit_data',np.int32,('time','orbit_dim'))
                rel_orb_nb = self.globalAttribs['DATATAKE_1_SENSING_ORBIT_NUMBER']
                orb_nb = root.find('.//safe:orbitNumber',namespaces=root.nsmap).text
                orb_dir = self.globalAttribs['DATATAKE_1_SENSING_ORBIT_DIRECTION']
                platform = self.globalAttribs['DATATAKE_1_SPACECRAFT_NAME']

                nc_orb.relativeOrbitNumber = rel_orb_nb
                nc_orb.orbitNumber = orb_nb
                nc_orb.orbitDirection = orb_dir
                nc_orb.platform = platform
                nc_orb.description = "Values structured as [relative orbit number, orbit number, platform]. platform corresponds to 0:Sentinel-2A, 1:Sentinel-2B.."

                nc_orb[0,:] = [int(rel_orb_nb),int(orb_nb),platform_id[platform]]


            # Add global attributes
            ##########################################################
            # Status
            logger.info('Adding global attributes')
            utils.memory_use(self.t0)

            nowstr = self.t0.strftime("%Y-%m-%dT%H:%M:%SZ")
            ncout.title = 'Sentinel-2 {} data'.format(self.processing_level)
            ncout.netcdf4_version_id = netCDF4.__netcdf4libversion__
            ncout.file_creation_date = nowstr

            self.globalAttribs['Conventions'] = "CF-1.6"
            self.globalAttribs[
                'summary'] = 'Sentinel-2 Multi-Spectral Instrument {} product.'.format(
                self.processing_level)
            self.globalAttribs[
                'keywords'] = '[Earth Science, Atmosphere, Atmospheric radiation, Reflectance]'
            self.globalAttribs['keywords_vocabulary'] = "GCMD Science Keywords"
            self.globalAttribs['institution'] = "Norwegian Meteorological Institute"
            self.globalAttribs['history'] = nowstr + ". Converted from SAFE to NetCDF by NBS team."
            self.globalAttribs['source'] = "surface observation"
            root = utils.xml_read(self.mainXML)
            if not self.dterrengdata:
                self.globalAttribs['orbitNumber'] = root.find('.//safe:orbitNumber',
                                                              namespaces=root.nsmap).text
            # Commented out to be stricly identical to older SAFE2NC version in production
            #else:
            #    self.globalAttribs['orbitNumber'] = root.find('.//SENSING_ORBIT_NUMBER').text

            self.globalAttribs['relativeOrbitNumber'] = self.globalAttribs.pop(
                'DATATAKE_1_SENSING_ORBIT_NUMBER')
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

    def rasterizeVectorLayers(self, nx, ny, gmlfile):

        # Open the data source and read in the extent
        NoData_value = 0

        source_ds = ogr.Open(str(gmlfile))
        source_layer = source_ds.GetLayer()
        # Check that gml file contains features
        if source_layer is None:
            return False, None, None

        geotransform = self.reference_band.GetGeoTransform()
        dst_ds = gdal.GetDriverByName('MEM').Create('', nx, ny, 1, gdal.GDT_Byte)
        dst_rb = dst_ds.GetRasterBand(1)
        dst_rb.SetNoDataValue(NoData_value)
        dst_ds.SetGeoTransform(geotransform)

        flag_values = []
        flag_meanings = []
        for i in range(0, source_layer.GetFeatureCount()):
            feat = source_layer.GetFeature(i)
            # Some SAFE files have their features starting at index 1 instead of 0
            # (it's the the case for the QT_PARTIALLY_CORRECTED_PIXELS for S2 L2A products)
            # So try and get the next feature if one is None, if still None, exit
            if feat is None:
                feat = source_layer.GetFeature(i+1)
                if feat is None:
                    return False, None, None
            name = feat.GetField('gml_id')
            source_layer.SetAttributeFilter("gml_id = \'%s\'" % name)
            if 'OPAQUE' in name:
                value = int(name.split('.')[-1]) + 1
            else:
                value = int(name[-1]) + 1
            gdal.RasterizeLayer(dst_ds, [1], source_layer, burn_values=[value])
            flag_values.append(value)
            flag_meanings.append(name)

        dst_ds.FlushCache()
        mask_arr = dst_ds.GetRasterBand(1).ReadAsArray()
        dst_ds = None

        # order flags for easier comparison
        sorted_pairs = sorted(zip(flag_meanings, flag_values))
        layer_mask = dict(sorted_pairs)

        return True, layer_mask, mask_arr

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


if __name__ == '__main__':

    # Log to console
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    log_info = logging.StreamHandler(sys.stdout)
    log_info.setFormatter(logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s'))
    logger.addHandler(log_info)

    workdir = pathlib.Path('/lustre/storeB/project/NBS2/sentinel/production/NorwAREA/netCDFNBS_work/test_environment/test_s2_N0400_updated_prod')
    #workdir = pathlib.Path('/lustre/storeA/users/elodief/NBS_test_data/fix_s2_11')
    workdir = pathlib.Path('/home/elodief/Data/NBS/NBS_test_data/processing_errors')

    products = ['S2A_MSIL2A_20220321T104731_N0400_R051_T32VNN_20220321T145616']

    for product in products:

        outdir = workdir / product
        outdir.parent.mkdir(parents=False, exist_ok=True)
        conversion_object = Sentinel2_reader_and_NetCDF_converter(
            product=product,
            indir=workdir / product,
            outdir=outdir)
        if conversion_object.read_ok:
            conversion_object.write_to_NetCDF(outdir, 7)



