import numpy as np
import logging
import osgeo.osr as osr
import scipy.ndimage
import netCDF4
from osgeo import gdal

from lib.utils import get_key
import config.constants as cst
from lib.netcdf.netcdf_base import NetCDFFile

logger = logging.getLogger(__name__)

class S2NetCDFFile(NetCDFFile):
    """
    Subclass for working with Sentinel-2 SAFE files
    """
    def __init__(self, product, directory, compression_level=7, chunk_size=(1, 2745, 2745)):

        super().__init__(product, directory, compression_level, chunk_size)

    def create_dimensions(self, safe_file):
        self.ncout.createDimension('time', 0)
        self.ncout.createDimension('x', safe_file.xSize)
        self.ncout.createDimension('y', safe_file.ySize)
        self.ncout.createDimension('raster_band_id', len(cst.s2_bands_order.keys()))
        self.ncout.createDimension('xa', safe_file.xaSize)
        self.ncout.createDimension('ya', safe_file.yaSize)

    def add_projection_coordinates(self, safe_file):
        xnp, ynp = safe_file.genLatLon(safe_file.xSize, safe_file.ySize, latlon=False)

        for varName, dim, data in [('x', 'x', xnp), ('y', 'y', ynp)]:
            self.variables[varName] = self.ncout.createVariable(
                varName, 'i4', dim,
                zlib=True, complevel=self.compression_level
            )
            self.variables[varName][:] = data

    def add_raster_band_id_variable(self):
        varName = 'band_id'
        self.variables[varName] = self.ncout.createVariable(varName, 'i4', 'raster_band_id', zlib=True, complevel=self.compression_level)
        band_keys = list(cst.s2_bands_order.keys())
        self.variables[varName][:] = np.array(band_keys)
        self.variables[varName].flag_values = np.array(band_keys, dtype=np.int8)
        self.variables[varName].flag_meanings = ' '.join(cst.s2_bands_order.values())

    def add_frequency_band_layers(self, safe_file):

        if safe_file.dterrengdata:
            # For DTERR data, gdal fails to properly do the src.GetSubDatasets()
            # so manually read the list of images created beforehand
            images = [[str(i), i.stem] for i in safe_file.image_list_dterreng]
        else:
            images = safe_file.src.GetSubDatasets()
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
                    if safe_file.dterrengdata:
                        band_metadata = None
                        varName = cst.s2_bands_aliases[v[-3::]]
                    else:
                        band_metadata = current_band.GetMetadata()
                        varName = band_metadata['BANDNAME']
                    if varName.startswith('B'):
                        self.variables[varName] = self.ncout.createVariable(
                            varName, np.uint16, ('time', 'y', 'x'), fill_value=0,
                            zlib=True, complevel=self.compression_level,
                            chunksizes=self.chunk_size
                        )
                        if safe_file.processing_level == 'Level-2A':
                            self.variables[varName].standard_name = 'surface_bidirectional_reflectance'
                        else:
                            self.variables[varName].standard_name = 'toa_bidirectional_reflectance'
                        self.variables[varName].long_name = 'Reflectance in band %s' % varName
                        if band_metadata:
                            self.variables[varName].bandwidth = band_metadata['BANDWIDTH']
                            self.variables[varName].bandwidth_unit = band_metadata['BANDWIDTH_UNIT']
                            self.variables[varName].wavelength = band_metadata['WAVELENGTH']
                            self.variables[varName].wavelength_unit = band_metadata['WAVELENGTH_UNIT']
                            self.variables[varName].solar_irradiance = band_metadata['SOLAR_IRRADIANCE']
                            self.variables[varName].solar_irradiance_unit = band_metadata['SOLAR_IRRADIANCE_UNIT']
                        # from DN to reflectance
                        logger.debug((varName, subdataset_geotransform))
                        if subdataset_geotransform[1] != 10:
                            current_size = current_band.XSize
                            band_measurement = scipy.ndimage.zoom(
                                input=current_band.GetVirtualMemArray(), zoom=safe_file.xSize / current_size,
                                order=0)
                        else:
                            band_measurement = current_band.GetVirtualMemArray()

                        self.write_variable_with_preprocessing(
                            varName, band_measurement,
                            process_chunk=None,
                            workers=8,
                            sync_every=None
                        )

    def write_variable_attributes(self, variable_attributes):

        super().write_variable_attributes(variable_attributes)

        for variable in self.variables:
            if variable.startswith('B'):
                for attribute, value in variable_attributes['B'].items():
                    setattr(self.variables[variable], attribute, value)

    def write_grid_mapping_variable(self, safe_file):

        # set grid mapping
        ##########################################################
        source_crs = osr.SpatialReference()
        source_crs.ImportFromWkt(safe_file.reference_band.GetProjection())
        varName = 'UTM_projection'
        self.variables[varName] = self.ncout.createVariable(varName, np.int32)
        self.variables[varName].latitude_of_projection_origin = source_crs.GetProjParm('latitude_of_origin')
        self.variables[varName].proj4_string = source_crs.ExportToProj4()
        self.variables[varName].crs_wkt = source_crs.ExportToWkt()
        self.variables[varName].semi_major_axis = source_crs.GetSemiMajor()
        self.variables[varName].scale_factor_at_central_meridian = source_crs.GetProjParm('scale_factor')
        self.variables[varName].longitude_of_central_meridian = source_crs.GetProjParm('central_meridian')
        self.variables[varName].grid_mapping_name = source_crs.GetAttrValue('PROJECTION').lower()
        self.variables[varName].semi_minor_axis = source_crs.GetSemiMinor()
        self.variables[varName].false_easting = source_crs.GetProjParm('false_easting')
        self.variables[varName].false_northing = source_crs.GetProjParm('false_northing')
        self.variables[varName].epsg_code = source_crs.GetAttrValue('AUTHORITY', 1)
        self.variables[varName].crs_wkt = safe_file.reference_band.GetProjection()

    def add_sun_and_view_angles(self, safe_file):
        # Add sun and view angles
        ##########################################################
        for varName in ['view_azimuth', 'view_zenith']:
            self.variables[varName] = self.ncout.createVariable(
                varName, np.float32, ('time','raster_band_id', 'ya', 'xa'),
                fill_value=netCDF4.default_fillvals['f4'],
                zlib=True, complevel=self.compression_level
            )

        counter = 1
        for k, v in list(safe_file.sunAndViewAngles.items()):
            logger.debug(("Handling %i of %i" % (counter, len(safe_file.sunAndViewAngles))))

            if 'sun' in k:
                self.variables[k] = self.ncout.createVariable(k, np.float32, ('time', 'ya', 'xa'), fill_value=netCDF4.default_fillvals['f4'],
                                                zlib=True, complevel=self.compression_level)
                self.variables[k].long_name = 'Solar %s angle' % k.split('_')[-1]
                self.write_variable_with_preprocessing(
                    k, v,
                    process_chunk=None,
                    workers=8,
                    sync_every=None
                )

            elif 'zenith' in k :
                band_id = k.split('_')[-1]
                self.variables['view_zenith'][0,get_key(cst.s2_bands_order,band_id), :, :] = v
            elif 'azimuth' in k :
                band_id = k.split('_')[-1]
                self.variables['view_azimuth'][0,get_key(cst.s2_bands_order,band_id), :, :] = v

            counter += 1

    def add_specific_L1C_or_L2A_layers(self, safe_file):
        specific_layers_kv = {}
        gdal_nc_data_types = {'Byte': 'u1', 'UInt16': 'u2'}

        if safe_file.processing_level == 'Level-1C':
            logger.info('Adding Level-1C specific layers')
            for layer in cst.s2_l1c_layers:
                for k, v in safe_file.imageFiles.items():
                    if layer in k:
                        specific_layers_kv[k] = cst.s2_l1c_layers[k]
                    elif layer in str(v):
                        logger.debug((layer, str(v), k))
                        specific_layers_kv[k] = cst.s2_l1c_layers[layer]

        elif safe_file.processing_level == 'Level-2A':
            logger.info('Adding Level-2A specific layers')
            for layer in cst.s2_l2a_layers:
                for k, v in safe_file.imageFiles.items():
                    if layer in k:
                        specific_layers_kv[k] = cst.s2_l2a_layers[k]
                    elif layer in str(v):
                        logger.debug((layer, str(v), k))
                        specific_layers_kv[k] = cst.s2_l2a_layers[layer]

        for k, v in specific_layers_kv.items():
            logger.debug((k, v))
            varName, longName = v.split(',')
            imageFile = safe_file.SAFE_dir / str(safe_file.imageFiles[k]).lstrip("/")
            SourceDS = gdal.Open(imageFile, gdal.GA_ReadOnly)
            nb_rasterBands =  SourceDS.RasterCount

            if SourceDS.RasterCount > 1:
                logger.info("Raster data contains more than one layer")

            for i in range(1,nb_rasterBands+1):
                if nb_rasterBands>1:
                    varName =  v.split(',')[0].split()[i-1]
                    longName =  v.split(',')[1].split('-')[i-1]
                #NDV = SourceDS.GetRasterBand(i).GetNoDataValue()
                xsize = SourceDS.RasterXSize
                GeoT = SourceDS.GetGeoTransform()
                logger.info("{}".format(GeoT))
                DataType = gdal_nc_data_types[
                    gdal.GetDataTypeName(SourceDS.GetRasterBand(i).DataType)]
                self.variables[varName] = self.ncout.createVariable(
                    varName, DataType, ('time', 'y', 'x'), fill_value=0,
                    zlib=True, complevel=self.compression_level, chunksizes=self.chunk_size
                )

                self.variables[varName].grid_mapping = "UTM_projection"
                self.variables[varName].long_name = longName
                if varName == "SCL":
                    self.variables[varName].flag_values = np.array(list(
                        cst.s2_scene_classification_flags.values()),
                                                    dtype=np.int8)
                    self.variables[varName].flag_meanings = ' '.join(
                        [key for key in list(cst.s2_scene_classification_flags.keys())])

                if GeoT[1] != 10:
                    raster_data = scipy.ndimage.zoom(input=SourceDS.GetRasterBand(i).GetVirtualMemArray(),
                                                        zoom=safe_file.xSize / xsize, order=0)
                else:
                    raster_data = SourceDS.GetRasterBand(i).GetVirtualMemArray()
                self.write_variable_with_preprocessing(
                    varName, raster_data,
                    process_chunk=None,
                    workers=8,
                    sync_every=None
                )

    def add_orbit_specific_data(self, safe_file):
        if not safe_file.dterrengdata:
            safe_file.globalAttribs['orbitNumber'] = safe_file.root.find('.//safe:orbitNumber',
                                                            namespaces=safe_file.root.nsmap).text
        else:
            safe_file.globalAttribs['orbitNumber'] = safe_file.root.find('.//SENSING_ORBIT_NUMBER').text

        self.ncout.createDimension('orbit_dim', 3)
        varName = 'orbit_data'
        self.variables[varName] = self.ncout.createVariable(varName, np.int32, ('time', 'orbit_dim'))
        self.variables[varName][0, :] = [
            int(safe_file.globalAttribs['DATATAKE_1_SENSING_ORBIT_NUMBER']),
            int(safe_file.globalAttribs['orbitNumber']),
            cst.platform_id[safe_file.globalAttribs['DATATAKE_1_SPACECRAFT_NAME']]
        ]

        rel_orb_nb = safe_file.globalAttribs['DATATAKE_1_SENSING_ORBIT_NUMBER']
        orb_nb = safe_file.globalAttribs['orbitNumber']
        orb_dir = safe_file.globalAttribs['DATATAKE_1_SENSING_ORBIT_DIRECTION']
        platform = safe_file.globalAttribs['DATATAKE_1_SPACECRAFT_NAME']

        self.variables[varName].relativeOrbitNumber = rel_orb_nb
        self.variables[varName].orbitNumber = orb_nb
        self.variables[varName].orbitDirection = orb_dir
        self.variables[varName].platform = platform
        self.variables[varName][0, :] = [int(rel_orb_nb), int(orb_nb), cst.platform_id[platform]]
    
    def add_integrated_information_as_variables(self, safe_file):
        varName = 'cloud_coverage'
        self.variables[varName] = self.ncout.createVariable(varName, np.float32, ('time'))
        self.variables[varName][0] = safe_file.globalAttribs['CLOUD_COVERAGE_ASSESSMENT']