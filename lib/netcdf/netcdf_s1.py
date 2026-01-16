import numpy as np
import datetime as dt
try:
    from lib.netcdf.netcdf_base import NetCDFFile
except ModuleNotFoundError:
    from netcdf_base import NetCDFFile

class S1NetCDFFile(NetCDFFile):
    """
    Subclass for working with Sentinel-1 SAFE files
    """
    def __init__(self, product, directory, compression_level=7, chunk_size=(1, 91, 99)):

        super().__init__(product, directory, compression_level, chunk_size)

    def create_dimensions(self, xSize, ySize):

        super().create_dimensions()
        self.ncout.createDimension('x', xSize)
        self.ncout.createDimension('y', ySize)

    def init_lat_lon(self):
        self.variables['lat'] = self.ncout.createVariable('lat', 'f4', ('y', 'x',), zlib=True,
        complevel=self.compression_level, chunksizes=self.chunk_size[1:])
        self.variables['lon'] = self.ncout.createVariable('lon', 'f4', ('y', 'x',), zlib=True,
        complevel=self.compression_level, chunksizes=self.chunk_size[1:])

    def add_raw_measurement_layers(self, safe_file):

        for i in range(1, safe_file.src.count + 1):
            band_data = safe_file.src.read(i)
            band_metadata = safe_file.src.tags(i)
            try:
                varName = 'Amplitude_%s' % band_metadata['POLARISATION']
            except:
                varName = 'Amplitude_%s' % band_metadata['POLARIZATION']
            self.variables[varName] = self.ncout.createVariable(
                varName, 'u2', ('time', 'y', 'x',),
                fill_value=0, zlib=True, complevel=self.compression_level,
                chunksizes=self.chunk_size
            )
            try:
                self.variables[varName].long_name = 'Amplitude %s-polarisation' % band_metadata['POLARIZATION']
            except:
                self.variables[varName].long_name = 'Amplitude %s-polarisation' % band_metadata['POLARISATION']
            try:
                self.variables[varName].polarisation = "%s" % band_metadata['POLARIZATION']
            except:
                self.variables[varName].polarisation = "%s" % band_metadata['POLARISATION']

            self.write_variable_with_preprocessing(
                varName, band_data,
                process_chunk=None,
                workers=8,
                sync_every=None,
            )

    def add_calibration_layers(self, safe_file):
        for calibration in safe_file.xmlCalLUTs:
            varName = str(calibration)
            current_polarisation = calibration.split('_')[-1]
            pixels, lines = safe_file.xmlCalPixelLines[current_polarisation]
            calibration_LUT = safe_file.xmlCalLUTs[calibration]
            resampled_calibration = safe_file.getCalLayer(pixels, lines, calibration_LUT)

            self.variables[varName] = self.ncout.createVariable(
                varName, 'f4', ('time', 'y', 'x',),
                zlib=True, complevel=self.compression_level,
                chunksizes=self.chunk_size
            )
            self.variables[varName].long_name = '%s calibration table' % calibration
            self.variables[varName].polarisation = "%s" % current_polarisation

            self.write_variable_with_preprocessing(
                varName, resampled_calibration,
                process_chunk=None,
                workers=8,
                sync_every=None
            )

    def write_grid_mapping_variable(self):
        self.variables['crsWGS84'] = self.ncout.createVariable('crsWGS84', np.int32)

    def add_noise_layers(self, safe_file):

        for polarisation in safe_file.polarisation:
            noiseCorrectionMatrix = safe_file.getNoiseCorrectionMatrix(safe_file.noiseVectors[polarisation],
                                                                  polarisation)
            varName = str('noiseCorrectionMatrix_' + polarisation)
            print(varName)
            self.variables[varName] = self.ncout.createVariable(
                varName, 'f4',
                ('time', 'y', 'x',),
                zlib=True, complevel=self.compression_level,
                chunksizes=self.chunk_size
            )

            self.variables[varName].polarisation = "%s" % polarisation

            self.write_variable_with_preprocessing(
                varName, noiseCorrectionMatrix,
                process_chunk=None,
                workers=8,
                sync_every=None
            )

    def add_subswath_layers(self, safe_file):
        for polarisation in safe_file.polarisation:
            swathLayer, flags = safe_file.getSwathList(polarisation)
            flag_values = np.array(sorted(flags.values()), dtype=np.int8)
            flags_meanings = ""
            for key in sorted(flags.keys()):
                flags_meanings += str(key + ' ')

            varName = 'swathList'

            self.variables[varName] = self.ncout.createVariable(
                varName, 'i1', ('y', 'x',), fill_value=0,
                zlib=True, complevel=self.compression_level, chunksizes=self.chunk_size[1:]
            )

            self.variables[varName].flag_values = flag_values
            self.variables[varName].valid_range = np.array([flag_values.min(), flag_values.max()])
            self.variables[varName].flag_meanings = flags_meanings.strip()

            self.write_variable_with_preprocessing(
                varName, swathLayer,
                process_chunk=None,
                workers=8,
                sync_every=None
            )
            break

    def add_gcp_information(self, safe_file):
        self.ncout.createDimension('gcp_index', len(safe_file.gcps))
        for key, value in safe_file.xmlGCPs.items():
            current_variable = key.split('_')[0]
            varName = str('GCP_%s' % key)
            if current_variable == 'azimuthTime':
                self.variables[varName] = self.ncout.createVariable(
                    varName, 'f4', ('gcp_index'), zlib=True
                )

                dates = np.array([dt.datetime.strptime(t, '%Y-%m-%dT%H:%M:%S.%f') for t in value])
                ref_date = dates.min()
                value = np.array([td.total_seconds() for td in dates - ref_date])
                self.variables[varName].comment = 'Seconds since %s' % ref_date.strftime('%Y-%m-%dT%H:%M:%S.%f')
            else:
                self.variables[varName] = self.ncout.createVariable(
                    varName, value.dtype, ('gcp_index'), zlib=True
                )

            self.variables[varName][:] = value

    def add_product_annotation_metadata(self, safe_file):

        for polarisation in safe_file.productMetadata:
            varName = str('s1Level1ProductSchema_' + polarisation)
            productMetadata = safe_file.productMetadata[polarisation]
            self.variables[varName] = self.ncout.createVariable(varName, 'i1')
            self.variables[varName].setncatts(productMetadata)

        productMetadataListComment = {
            'swathMergeList': 'index:{swath:[firstAzimuthLine, firstRangeSample, lastAzimuthLine, '
                              'lastRangeSample, azimuthTime]}',
            'orbitList': 'time:[frame, position (x,y,z), velocity (x,y,z)]',
            'coordinateConversionList': 'index:[azimuthTime, slantRangeTime, sr0, '
                                        'srgrCoefficients, gr0, grsrCoefficients ]',
            'antennaPatternList': 'index:[swath, azimuthTime, slantRangeTime, elevationAngle, '
                                  'elevationPattern, incidenceAngle, terrainHeight, roll]'
            }

        for polarisation in safe_file.productMetadataList:
            for subkey in safe_file.productMetadataList[polarisation]:
                varName = str(subkey + '_' + polarisation)
                if True:
                    productMetadataList = safe_file.productMetadataList[polarisation][subkey]
                    tmp_dict = {}
                    for k, v in productMetadataList.items():
                        tmp_dict[str(k)] = str(v)
                    self.variables[varName]= self.ncout.createVariable(varName, 'i1')
                    self.variables[varName].comment = productMetadataListComment[subkey]
                    self.variables[varName].setncatts(tmp_dict)