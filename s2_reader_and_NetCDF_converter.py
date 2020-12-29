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
import subprocess
import time
from collections import defaultdict
from datetime import datetime
import sys
import gdal
import lxml.etree as ET
import netCDF4
import numpy as np
import osgeo.ogr as ogr
import osgeo.osr as osr
import pyproj
import scipy.ndimage
import safe_to_netcdf.utils as utils


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
        self.input_zip = indir / product.with_suffix('.zip')
        self.SAFE_dir = outdir / self.product_id.with_suffix('.SAFE')
        self.processing_level = None
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
        self.bands_alias_bandID = {0: 'B1', 1: 'B2', 2: 'B3', 3: 'B4', 4: 'B5',
                                   5: 'B6', 6: 'B7', 7: 'B8', 8: 'B8A', 9: 'B9', 10: 'B10',
                                   11: 'B11',
                                   12: 'B12'}
        self.SAFE_structure = None
        self.image_list_dterreng = []

        self.main()

    def main(self):
        """ Main method for traversing and reading key values from SAFE
            directory.
        """

        # 1) Fetch main file
        self.xmlFiles['manifest'] = self.uncompress('manifest.safe')
        print((self.xmlFiles['manifest']))
        if not self.xmlFiles['manifest']:
            print("\nNo manifest.safe file. Most likely S2 L1C Norwegian DEM product")
            self.dterrengdata = True
            self.xmlFiles['mainXML'] = self.uncompress('MTD*.xml')

        # 2) Set some of the global __init__ variables
        if not self.dterrengdata:
            initializer_ok = self.initializer(self.xmlFiles['manifest'])
        else:
            initializer_ok = self.initializer_dterr(self.xmlFiles['mainXML'])

        # 3) Read sun and view angles
        print('\nRead view and sun angles')
        if not self.dterrengdata:
            self.readSunAndViewAngles(self.xmlFiles['S2_{}_Tile1_Metadata'.format(
                self.processing_level)])

            # self.readSunAndViewAngles(self.SAFE_path + self.xmlFiles[
            # 'S2_Level-2A_Tile1_Metadata'])
        else:
            self.readSunAndViewAngles(self.xmlFiles['MTD_TL'])

        # 4) Read vector information
        print('\nRead vector information')

        for gmlfile in self.xmlFiles.values():
            if gmlfile and gmlfile.suffix == '.gml':
                    vectorID, vectorPath = self.readVectorInformation(gmlfile)
                    self.vectorInformation[vectorID] = vectorPath

        # 5) Retrieve SAFE product structure
        ##self.SAFE_structure = self.list_product_structure()

    def uncompress(self, xmlReturn):
        """ Uncompress SAFE zip file. Return xmlReturn file

            Keyword arguments:
            xmlReturn -- string specifying name of xmlFile to return
        """
        """ Uncompress SAFE zip file. Return manifest file """

        # If zip not extracted yet
        if not self.SAFE_dir.is_dir():
            cmd = f'/usr/bin/unzip {self.input_zip} -d {self.SAFE_dir.parent}'
            subprocess.call(cmd, shell=True)
            if not self.SAFE_dir.is_dir():
                print(f'Error unzipping file {self.input_zip}')
                sys.exit()

        # todo: will not work for dterreng data as it's a regular expression
        xmlFile = self.SAFE_dir / xmlReturn
        if not xmlFile.is_file():
            print(f'Main file not available {xmlFile}')
            return False

        return xmlFile

    def initializer(self, mainXML):
        """ Traverse manifest file for setting additional variables
            in __init__ """

        root = utils.xml_read(mainXML)

        # Set xml-files
        dataObjectSection = root.find('./dataObjectSection')
        for dataObject in dataObjectSection.findall('./'):
            repID = dataObject.attrib['ID']
            ftype = None
            href = None
            for element in dataObject.iter():
                attrib = element.attrib
                if 'mimeType' in attrib:
                    ftype = attrib['mimeType']
                if 'href' in attrib:
                    href = attrib['href'][1:]
            if (ftype == 'text/xml' or ftype == 'application/xml') and href:
                self.xmlFiles[repID] = self.SAFE_dir / href[1:]
            elif ftype == 'application/octet-stream':
                self.imageFiles[repID] = self.SAFE_dir / href[1:]

        # Set processing level
        for k in list(self.xmlFiles.keys()):
            if "Level-1C" in k:
                self.processing_level = "Level-1C"
                break
            elif "Level-2A" in k:
                self.processing_level = "Level-2A"
                break

        if not self.processing_level:
            print("Unknown processing level. Hence terminating.")
            sys.exit(1)

        # Set gdal object
        # str(self.xmlFiles['manifest'])
        self.src = gdal.Open(
            str(self.xmlFiles['S2_{}_Product_Metadata'.format(self.processing_level)]))
        print((self.src))

        # Set global metadata attributes from gdal
        self.globalAttribs = self.src.GetMetadata()

        return True

##    def initializer_dterr(self, mainXML):
##        """Set additional variables in __init__ for dterrengdata products"""
##        SAFE_outpath = self.SAFE_outpath
##        SAFE_file = self.SAFE_file
##        SAFE_id = self.SAFE_id
##        fdirName = '%s%s.SAFE' % (SAFE_outpath, SAFE_id)
##
##        # Adding xml files
##        for dirName, subdirList, fileList in os.walk(fdirName):
##            for fname in fileList:
##                if fname.endswith('.xml') or fname.endswith('.gml'):
##                    fID = fname.split('.')[0]
##                    self.xmlFiles[fID] = '/'.join((dirName, fname))
##
##        # Set gdal object
##        self.src = gdal.Open(self.xmlFiles['MTD_MSIL1C'])
##
##        # Read relative image path (since gdal can't open all these products..)
##        tree = ET.parse(self.xmlFiles['MTD_MSIL1C'])
##        root = tree.getroot()
##
##        for element in root.findall('.//IMAGE_FILE'):
##            img = self.SAFE_path + '/' + element.text
##            self.image_list_dterreng.append(img)
##
##        # Set global metadata attributes from gdal
##        self.globalAttribs = self.src.GetMetadata()
##
    def write_to_NetCDF(self, nc_outpath, compression_level, chunk_size=(1, 32, 32)):
        """ Method writing output NetCDF product.

        Keyword arguments:
        nc_outpath -- output path where NetCDF file should be stored
        compression_level -- compression level on output NetCDF file (1-9)
        chunk_size -- chunk_size
        """

        print("------------START CONVERSION FROM SAFE TO NETCDF-------------")
        print("------------DEBUG-------------")

        # Status
        print('\nCreating NetCDF file')
        utils.memory_use(self.t0)

        # Deciding a reference band
        for k, v in self.src.GetSubDatasets():
            if v.find('10m') > 0:
                self.reference_band = gdal.Open(k)

        nx = self.reference_band.RasterXSize  # number of pixels for 10m spatial resolution
        # frequency bands
        ny = self.reference_band.RasterYSize  # number of pixels for 10m spatial resolution
        # frequency bands

        # output filename
        out_netcdf = nc_outpath / self.product_id.with_suffix('.nc')

        with (netCDF4.Dataset(out_netcdf, 'w', format='NETCDF4')) as ncout:
            ncout.createDimension('time', 0)
            ncout.createDimension('x', nx)
            ncout.createDimension('y', ny)

            nctime = ncout.createVariable('time', 'i4', ('time',))
            # nclat = ncout.createVariable('lat','f4',('y','x',),zlib=True,
            # complevel=compression_level, chunksizes=chunk_size[1:])
            # nclon = ncout.createVariable('lon','f4',('y','x',),zlib=True,
            # complevel=compression_level, chunksizes=chunk_size[1:])

            # Set time value
            ##########################################################
            nctime.long_name = 'reference time of satellite image'
            nctime.units = 'seconds since 1981-01-01 00:00:00'
            nctime.calendar = 'gregorian'
            nctime[:] = self.getReferenceTime()

            # Add latitude and longitude layers
            ##########################################################
            # Status
            """
            print('\nCreating latitude longitude')
            self.memoryUsage()
            print(datetime.now() - self.t0)

            lat,lon = self.genLatLon(nx,ny) #Assume gcps are on a regular grid
            nclat.long_name = 'latitude'
            nclat.units = 'degrees_north'
            nclat.standard_name = 'latitude'
            nclat[:,:]=lat

            nclon.long_name = 'longitude'
            nclon.units = 'degrees_east'
            nclon.standard_name = 'longitude'
            nclon[:,:]=lon
            """
            # Add projection coordinates
            ##########################################################
            # Status
            print('\nAdding projection coordinates')
            utils.memory_use(self.t0)

            xnp, ynp = self.genLatLon(nx, ny, latlon=False)  # Assume gcps are on a regular grid

            ncx = ncout.createVariable('x', 'i4', 'x', zlib=True)
            ncx.units = 'm'
            ncx.standard_name = 'projection_x_coordinate'
            ncx[:] = xnp

            ncy = ncout.createVariable('y', 'i4', 'y', zlib=True)
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
            print('\nAdding frequency bands layers')
            utils.memory_use(self.t0)

            if not self.dterrengdata:
                for k, v in self.src.GetSubDatasets():
                    subdataset = gdal.Open(k)
                    subdataset_geotransform = subdataset.GetGeoTransform()
                    if not "True color image" in v:
                        for i in range(1, subdataset.RasterCount + 1):
                            current_band = subdataset.GetRasterBand(i)
                            band_metadata = current_band.GetMetadata()

                            varName = band_metadata['BANDNAME']
                            varout = ncout.createVariable(varName, np.uint16,
                                                          ('time', 'y', 'x'), fill_value=0,
                                                          zlib=True, complevel=compression_level,
                                                          chunksizes=chunk_size)
                            varout.units = "1"
                            # varout.coordinates = "lat lon" ;
                            varout.grid_mapping = "UTM_projection"
                            if not self.processing_level == 'Level-2A':
                                varout.standard_name = 'toa_bidirectional_reflectance'
                            else:
                                varout.standard_name = 'surface_bidirectional_reflectance'
                            varout.long_name = 'Reflectance in band %s' % band_metadata['BANDNAME']
                            varout.bandwidth = band_metadata['BANDWIDTH']
                            varout.bandwidth_unit = band_metadata['BANDWIDTH_UNIT']
                            varout.wavelength = band_metadata['WAVELENGTH']
                            varout.wavelength_unit = band_metadata['WAVELENGTH_UNIT']
                            varout.solar_irradiance = band_metadata['SOLAR_IRRADIANCE']
                            varout.solar_irradiance_unit = band_metadata['SOLAR_IRRADIANCE_UNIT']
                            varout._Unsigned = "true"
                            # varout.scale_factor = 0.0001 # 1/(quanitification value) converts
                            # from DN to reflectance
                            print((varName, subdataset_geotransform))
                            if subdataset_geotransform[1] != 10:
                                current_size = current_band.XSize
                                band_measurement = scipy.ndimage.zoom(
                                    input=current_band.GetVirtualMemArray(), zoom=nx / current_size,
                                    order=0)
                            else:
                                band_measurement = current_band.GetVirtualMemArray()
                            print((band_measurement.shape))

                            # varout[0,:,:] = band_measurement

                    else:  # 8 bit true color image ("u1")
                        # create new rgb dimension
                        ncout.createDimension('dimension_rgb', subdataset.RasterCount)
                        varout = ncout.createVariable('TCI', 'u1',
                                                      ('time', 'dimension_rgb', 'y', 'x'),
                                                      fill_value=0, zlib=True,
                                                      complevel=compression_level,
                                                      chunksizes=(1,) + chunk_size)

                        varout.units = "1"
                        # varout.coordinates = "lat lon" ;
                        varout.grid_mapping = "UTM_projection"
                        varout.long_name = 'TCI RGB from B4, B3 and B2'
                        varout._Unsigned = "true"
                        for i in range(1, subdataset.RasterCount + 1):
                            current_band = subdataset.GetRasterBand(i)
                            band_metadata = current_band.GetMetadata()
                            band_measurement = current_band.GetVirtualMemArray()
                            varout[0, i - 1, :, :] = band_measurement
            else:
                bands_alias = {'B01': 'B1', 'B02': 'B2', 'B03': 'B3', 'B04': 'B4',
                               'B05': 'B5', 'B06': 'B6', 'B07': 'B7', 'B08': 'B8',
                               'B8A': 'B8A', 'B09': 'B9', 'B10': 'B10',
                               'B11': 'B11', 'B12': 'B12'}
                for raster_path in self.image_list_dterreng:
                    subdataset = gdal.Open(raster_path)
                    subdataset_geotransform = subdataset.GetGeoTransform()
                    if not "TCI" in raster_path:
                        for i in range(1, subdataset.RasterCount + 1):
                            current_band = subdataset.GetRasterBand(i)
                            band_metadata = current_band.GetMetadata()

                            varName = bands_alias[raster_path.split('.')[-2][-3::]]
                            varout = ncout.createVariable(varName, np.int16,
                                                          ('time', 'y', 'x'), fill_value=0,
                                                          zlib=True, complevel=compression_level,
                                                          chunksizes=chunk_size)
                            varout.units = "1"
                            # varout.coordinates = "lat lon" ;
                            varout.grid_mapping = "UTM_projection"
                            varout.standard_name = 'toa_bidirectional_reflectance'
                            varout.long_name = 'Reflectance in band %s' % bands_alias[
                                raster_path.split('.')[-2][-3::]]
                            # varout.bandwidth = band_metadata['BANDWIDTH']
                            # varout.bandwidth_unit = band_metadata['BANDWIDTH_UNIT']
                            # varout.wavelength = band_metadata['WAVELENGTH']
                            # varout.wavelength_unit = band_metadata['WAVELENGTH_UNIT']
                            # varout.solar_irradiance = band_metadata['SOLAR_IRRADIANCE']
                            # varout.solar_irradiance_unit = band_metadata['SOLAR_IRRADIANCE_UNIT']
                            varout._Unsigned = "true"
                            # varout.scale_factor = 0.0001 # 1/(quanitification value) converts
                            # from DN to reflectance
                            if subdataset_geotransform[1] != 10:
                                current_size = current_band.XSize
                                band_measurement = scipy.ndimage.zoom(
                                    input=current_band.GetVirtualMemArray(), zoom=nx / current_size,
                                    order=0)
                            else:
                                band_measurement = current_band.GetVirtualMemArray()
                            varout[:] = band_measurement

                    else:  # 8 bit true color image ("u1")
                        # create new rgb dimension
                        ncout.createDimension('dimension_rgb', subdataset.RasterCount)
                        varout = ncout.createVariable('TCI', 'u1',
                                                      ('time', 'dimension_rgb', 'y', 'x'),
                                                      fill_value=0, zlib=True,
                                                      complevel=compression_level,
                                                      chunksizes=(1,) + chunk_size)

                        varout.units = "1"
                        # varout.coordinates = "lat lon" ;
                        varout.grid_mapping = "UTM_projection"
                        varout.long_name = 'TCI RGB from B4, B3 and B2'
                        varout._Unsigned = "true"
                        for i in range(1, subdataset.RasterCount + 1):
                            current_band = subdataset.GetRasterBand(i)
                            band_metadata = current_band.GetMetadata()
                            band_measurement = current_band.GetVirtualMemArray()
                            varout[i - 1, :, :] = band_measurement

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
            print('\nAdding vector layers')
            utils.memory_use(self.t0)

            for layer, path in self.vectorInformation.items():
                if path:
                    output_file = (self.SAFE_dir / 'tmp' / layer).with_suffix('.tiff')
                    rasterized_ok, layer_mask = self.rasterizeVectorLayers(nx, ny, path,
                                                                           output_file)
                    # Warning 1: Failed to fetch spatial reference on layer MSK_CLOUDS_B00 to
                    # build transformer, assuming matching coordinate systems.
                    if rasterized_ok:
                        # todo: fusionner les 2 cas
                        if layer == "MSK_CLOUDS_B00":
                            varout = ncout.createVariable('Clouds', 'i1', ('time', 'y', 'x'),
                                                          fill_value=-1, zlib=True,
                                                          chunksizes=chunk_size)
                            varout.long_name = "Cloud mask 10m resolution"
                            varout.comment = "Rasterized cloud information."
                        else:
                            varout = ncout.createVariable(layer, 'i1', ('time', 'y', 'x'),
                                                          fill_value=-1, zlib=True,
                                                          chunksizes=chunk_size)
                            varout.long_name = "%s mask 10m resolution" % layer
                            varout.comment = "Rasterized vector information."

                        # varout.coordinates = "lat lon"
                        varout.grid_mapping = "UTM_projection"
                        varout.flag_values = np.array(list(layer_mask.values()), dtype=np.int8)
                        # varout.flag_values = ', '.join([str(i) for i in layer_mask.values()])
                        # varout.flag_meanings = ' '.join(layer_mask.keys())
                        varout.flag_meanings = ' '.join(
                            [key.replace('-', '_') for key in list(layer_mask.keys())])
                        vector_band = gdal.Open(str(output_file))
                        varout[0, :] = vector_band.GetVirtualMemArray()
                        break

            # Add Level-2A layers
            ##########################################################
            # Status
            if self.processing_level == 'Level-2A':
                print('\nAdding Level-2A specific layers')
                utils.memory_use(self.t0)
                l2a_layers = {"MSK_CLDPRB_20m": "MSK_CLDPRB, Cloud Probabilities",
                              "MSK_SNWPRB_20m": "MSK_SNWPRB, Snow Probabilities",
                              'IMG_DATA_Band_AOT_10m_Tile1_Data': "AOT, Aerosol Optical Thickness",
                              'IMG_DATA_Band_WVP_10m_Tile1_Data': "WVP, Water Vapour",
                              'IMG_DATA_Band_SCL_20m_Tile1_Data': "SCL, Scene Classification"}
                scene_classifcation_flags = {'NODATA': 0, 'SATURATED_DEFECTIVE': 1,
                                             'DARK_FEATURE_SHADOW': 2,
                                             'CLOUD_SHADOW': 3,
                                             'VEGETATION': 4,
                                             'NOT_VEGETATED': 5,
                                             'WATER': 6,
                                             'UNCLASSIFIED': 7,
                                             'CLOUD_MEDIUM_PROBA': 8,
                                             'CLOUD_HIGH_PROBA': 9,
                                             'THIN_CIRRUS': 10,
                                             'SNOW_ICE': 11}
                gdal_nc_data_types = {'Byte': 'u1', 'UInt16': 'u2'}

                l2a_kv = {}
                for layer in list(l2a_layers.keys()):
                    for k, v in list(self.imageFiles.items()):
                        if layer in k:
                            l2a_kv[k] = l2a_layers[k]
                        elif layer in v:
                            print((layer, v, k))
                            l2a_kv[k] = l2a_layers[layer]

                for k, v in list(l2a_kv.items()):
                    print((k, v))
                    varName, longName = v.split(',')
                    SourceDS = gdal.Open(self.SAFE_path + self.imageFiles[k], gdal.GA_ReadOnly)
                    if SourceDS.RasterCount > 1:
                        print("Raster data contains more than one layer")
                    NDV = SourceDS.GetRasterBand(1).GetNoDataValue()
                    xsize = SourceDS.RasterXSize
                    ysize = SourceDS.RasterYSize
                    GeoT = SourceDS.GetGeoTransform()
                    DataType = gdal_nc_data_types[
                        gdal.GetDataTypeName(SourceDS.GetRasterBand(1).DataType)]
                    # print(NDV, xsize, ysize, GeoT, DataType)

                    varout = ncout.createVariable(varName, DataType,
                                                  ('time', 'y', 'x'), fill_value=0, zlib=True,
                                                  complevel=compression_level,
                                                  chunksizes=chunk_size)
                    # varout.coordinates = "lat lon" ;
                    varout.grid_mapping = "UTM_projection"
                    varout.long_name = longName
                    if varName == "SCL":
                        varout.flag_values = np.array(list(scene_classifcation_flags.values()),
                                                      dtype=np.int8)
                        varout.flag_meanings = ' '.join(
                            [key for key in list(scene_classifcation_flags.keys())])

                    if GeoT[1] != 10:
                        raster_data = scipy.ndimage.zoom(input=SourceDS.GetVirtualMemArray(),
                                                         zoom=nx / xsize, order=0)
                    else:
                        raster_data = SourceDS.GetVirtualMemArray()
                    varout[0, :] = raster_data

            # Add sun and view angles
            ##########################################################
            # Status
            print('\nAdding sun and view angles')
            utils.memory_use(self.t0)

            counter = 1
            for k, v in list(self.sunAndViewAngles.items()):
                print(("\tHandeling %i of %i" % (counter, len(self.sunAndViewAngles))))
                angle_step = int(math.ceil(nx / float(v.shape[0])))

                resampled_angles = self.resample_angles(v, nx, v.shape[0], v.shape[1], angle_step,
                                                        type=np.float32)

                varout = ncout.createVariable(k, np.float32, ('time', 'y', 'x'),
                                              fill_value=netCDF4.default_fillvals['f4'], zlib=True,
                                              chunksizes=chunk_size)
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
            print('\nAdding XML files as character variables')
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
            print('\nAdding SAFE product structure as character variable')
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

            platform_id = {"Sentinel-2A": 0, "Sentinel-2B": 1,
                           "Sentinel-2C": 2, "Sentinel-2D": 3, }
            # orb_dir_id = {"DESCENDING":0, "":1,

            if self.xmlFiles['manifest']:
                root = utils.xml_read(self.xmlFiles['manifest'])
                self.globalAttribs['orbitNumber'] = root.find('.//safe:orbitNumber',
                                                              namespaces=root.nsmap).text

            ncout.createDimension('orbit_dim', 3)
            nc_orb = ncout.createVariable('orbit_data', np.int32, ('time', 'orbit_dim'))
            rel_orb_nb = self.globalAttribs['DATATAKE_1_SENSING_ORBIT_NUMBER']
            orb_nb = root.find('.//safe:orbitNumber', namespaces=root.nsmap).text
            orb_dir = self.globalAttribs['DATATAKE_1_SENSING_ORBIT_DIRECTION']
            platform = self.globalAttribs['DATATAKE_1_SPACECRAFT_NAME']

            nc_orb.relativeOrbitNumber = rel_orb_nb
            nc_orb.orbitNumber = orb_nb
            nc_orb.orbitDirection = orb_dir
            nc_orb.platform = platform
            nc_orb.description = "Values structured as [relative orbit number, orbit number, " \
                                 "platform]. platform corresponds to 0:Sentinel-2A, 1:Sentinel-2B.."

            nc_orb[0, :] = [int(rel_orb_nb), int(orb_nb), platform_id[platform]]

            # Add global attributes
            ##########################################################
            # Status
            print('\nAdding global attributes')
            utils.memory_use(self.t0)

            # nowstr = self.t0.strftime("%Y-%m-%dT%H:%M:%SZ")
            # todo: use ref time
            nowstr = datetime.strftime(
                datetime(self.t0.year, self.t0.month, self.t0.day, self.t0.hour, self.t0.minute,
                         self.t0.second), "%Y-%m-%dT%H:%M:%SZ")

            ncout.title = 'Sentinel-2 {} data'.format(self.processing_level)
            ncout.netcdf4_version_id = netCDF4.__netcdf4libversion__
            ncout.file_creation_date = nowstr

            self.globalAttribs['Conventions'] = "CF-1.8"
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
            ncout.setncatts(self.globalAttribs)
            ncout.sync()

            # Status
            print('\nFinished.')
            utils.memory_use(self.t0)

        return out_netcdf.is_file()

    def xmlToString(self, xmlfile):
        """ Method for reading XML files returning the entire file as single
            string.
        """
        if not xmlfile.is_file():
            print(('Error: Can\'t find xmlfile %s' % (xmlfile)))
            return False
        try:
            parser = ET.XMLParser(recover=True)
            tree = ET.parse(str(xmlfile), parser)
            return ET.tostring(tree)
        except:
            print(("Could not parse %s as xmlFile. Try to open regularly." % xmlfile))
            with open(xmlfile, 'r') as infile:
                text = infile.read()
            if text:
                return text
            else:
                print(("Could not parse %s. Something wrong with file." % xmlfile))
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
        for BANDID in np.array(list(self.bands_alias_bandID.keys())):
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
                    str('view_zenith_' + self.bands_alias_bandID[BANDID])] = tmp_view_zenith
                self.sunAndViewAngles[
                    str('view_azimuth_' + self.bands_alias_bandID[BANDID])] = tmp_view_azimuth

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

    def readVectorInformation(self, gmlfile):
        """ Reading vector information from Sentinel-2 .gml files.
        """

        gmlID = gmlfile.stem

        if not gmlfile.is_file():
            print(('Error: Can\'t find gmlfile %s' % (gmlfile)))
            return gmlID, False

        # Create directories for output
        output_dir = self.SAFE_dir / 'tmp'
        # if python < 3.5
        ##if not output_dir.is_dir():
        ##    output_dir.mkdir()
        output_dir.mkdir(exist_ok=True)

        destName = (output_dir / gmlID).with_suffix('.shp')

        vto = gdal.VectorTranslateOptions(format='ESRI Shapefile')
        vt = gdal.VectorTranslate(destNameOrDestDS=str(destName), srcDS=str(gmlfile), options=vto)
        vt.FlushCache()
        vt = None

        if destName.is_file():
            return gmlID, destName
        else:
            return gmlID, False

    def rasterizeVectorLayers(self, nx, ny, shapefile,
                              output_file):  # ,output_dir=str(self.SAFE_path + '/tmp/'),
        # maskType=None):
        """ Supports all vector layers. Rasterized by means of gdal tools."""
        ulx, xres, xskew, uly, yskew, yres = self.reference_band.GetGeoTransform()  # ulx - upper
        # left x, uly - upper left y

        """ For the future....?
        if shapefile:
            if maskType:
                ro =  gdal.RasterizeOptions(xRes=xres, yRes=yres,
                    outputBounds=[ulx, uly+(yres*ny),ulx+(xres*nx), uly], width=nx,
                    height=ny,format="GTiff", where="maskType = 'OPAQUE'", burnValues=1,
                    layers=[layer_name], noData=0)
            else:
                ro =  gdal.RasterizeOptions(xRes=xres, yRes=yres,
                    outputBounds=[ulx, uly+(yres*ny),ulx+(xres*nx), uly], width=nx,
                    height=ny,format="GTiff", burnValues=1,
                    layers=[layer_name], noData=0)
            tmp = gdal.Rasterize(destNameOrDestDS=str(output_dir + layer_name + '.tif'),
            srcDS=shapefile, options=ro)
            tmp.FlushCache()
            tmp=None
        vector_fn = src
        nx = ny = 10980
        # Define pixel_size and NoData value of new raster
        pixel_size = 10
        NoData_value = 255
        """

        # Open the data source and read in the extent
        NoData_value = 0
        source_ds = ogr.Open(str(shapefile))
        source_layer = source_ds.GetLayer()
        source_srs = source_layer.GetSpatialRef()
        x_min, x_max, y_min, y_max = source_layer.GetExtent()
        target_ds = gdal.GetDriverByName('GTiff').Create(str(output_file), nx, ny, 1, gdal.GDT_Byte)
        # target_ds.SetGeoTransform((x_min, xres, 0, y_max, 0, yres))
        target_ds.SetGeoTransform(self.reference_band.GetGeoTransform())
        target_ds.SetProjection(self.reference_band.GetProjection())
        band = target_ds.GetRasterBand(1)
        band.SetNoDataValue(NoData_value)

        # Find features:
        gml_id = []
        maskType = []
        for i in range(0, source_layer.GetFeatureCount()):
            f = source_layer.GetFeature(i)
            gml_id.append(f.GetField('gml_id'))
            maskType.append(f.GetField('maskType'))

        unique_gml_id = list(set(gml_id))
        unique_maskType = list(set(maskType))

        layer_mask = {}
        if ('CIRRUS' or 'OPAQUE') in unique_maskType:
            for i, filt in enumerate(unique_maskType):  # 1 - OPAQUE, 2 - CIRRUS
                source_layer.SetAttributeFilter("maskType = \'%s\'" % filt)

                # Rasterize
                gdal.RasterizeLayer(target_ds, [1], source_layer, burn_values=[i + 1])
                layer_mask[filt] = i + 1
            target_ds.FlushCache()
            target_ds = None
        else:
            for i, filt in enumerate(unique_gml_id):
                source_layer.SetAttributeFilter("gml_id = \'%s\'" % filt)

                # Rasterize
                gdal.RasterizeLayer(target_ds, [1], source_layer, burn_values=[i + 1])
                layer_mask[filt] = i + 1

            target_ds.FlushCache()
            target_ds = None
        if output_file.is_file():
            return True, layer_mask
        else:
            return False, []

    def getReferenceTime(self):
        """ Method for retrieving Geo Location Point parameter from xml file."""
        dt = datetime.strptime(self.globalAttribs["PRODUCT_START_TIME"].split('.')[0],
                               '%Y-%m-%dT%H:%M:%S')
        AQ_START = time.mktime((dt.year, dt.month, dt.day, dt.hour, dt.minute, dt.second, 0, 0, 0))
        dt1981 = time.mktime((1981, 1, 1, 0, 0, 0, 0, 0, 0))
        delta_t = int(AQ_START - dt1981)
        return delta_t

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

# todo
##    def list_product_structure(self, startpath=None):
##        """ Traverse SAFE file structure (or any file structure) and
##            creates a xml file containing the file structure.
##
##            Returns a string representation of the xml file."""
##
##        if not startpath:
##            startpath = self.SAFE_path
##
##        root_path = startpath.split('/')[-1]
##        ET_root = ET.Element(root_path)
##        # Create dictionary that contains all elements
##        elements = {root_path: ET_root}
##
##        # Iterate throgh the file structure
##        for root, dirs, files in os.walk(startpath):
##            level = root.replace(startpath, '').count(os.sep)
##            current_xpath = str('/' + root_path + root.replace(startpath, ''))
##            current_path = root.replace(startpath, '')
##
##            # Create all folder elements
##            if dirs:
##                for dir in dirs:
##                    element = ET.SubElement(elements[current_xpath.strip('/')], dir)
##                    elements[str(current_xpath.strip('/') + '/' + dir)] = element
##            # Add all files in the correct folder
##            for f in files:
##                sub_element = ET.SubElement(elements[current_xpath.strip('/')], 'file')
##                sub_element.text = f
##
##        return ET.tostring(ET_root)
##

if __name__ == '__main__':

    workdir = pathlib.Path('/home/elodief/Data/NBS')

    products = ['S2A_MSIL1C_20201028T102141_N0209_R065_T34WDA_20201028T104239']

    for product in products:

        conversion_object = Sentinel2_reader_and_NetCDF_converter(
            product=pathlib.Path(product),
            indir=workdir / 'NBS_reference_data' / 'reference_datain',
            outdir=workdir / 'NBS_test_data' / 'safe2nc_latest_local_02' / product)

        nc_outpath = workdir / 'NBS_test_data' / 'safe2nc_latest_local_02' / product
        cl = 7
        conversion_object.write_to_NetCDF(nc_outpath, cl)
