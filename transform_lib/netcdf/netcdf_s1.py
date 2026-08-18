import numpy as np
import datetime as dt
try:
    from transform_lib.netcdf.netcdf_base import NetCDFFile
except ModuleNotFoundError:
    from netcdf_base import NetCDFFile

class S1NetCDFFile(NetCDFFile):
    """
    Subclass for working with Sentinel-1 NetCDF files
    """
    def __init__(self, product, directory, compression_level=7, chunk_size=(1, 91, 99)):

        super().__init__(product, directory, compression_level, chunk_size)

    def create_dimensions(self, safe_file):

        self.ncout.createDimension('time', 1)
        self.ncout.createDimension('x', safe_file.xSize)
        self.ncout.createDimension('y', safe_file.ySize)

    def init_lat_lon(self):
        self.variables['lat'] = self.ncout.createVariable('lat', 'f4', ('y', 'x',), zlib=True,
        complevel=self.compression_level, chunksizes=self.chunk_size[1:])
        self.variables['lon'] = self.ncout.createVariable('lon', 'f4', ('y', 'x',), zlib=True,
        complevel=self.compression_level, chunksizes=self.chunk_size[1:])

    def add_raw_measurement_layers(self, safe_file):
        for i in range(1, safe_file.src.RasterCount + 1):
            band = safe_file.src.GetRasterBand(i)
            band_metadata = band.GetMetadata()
            try:
                baseName = 'Amplitude_%s' % band_metadata['POLARISATION']
            except:
                baseName = 'Amplitude_%s' % band_metadata['POLARIZATION']

            # Read array robustly
            arr = band.ReadAsArray()

            # If complex (SLC), store real and imaginary parts separately
            if np.iscomplexobj(arr):
                real_name = f"{baseName}_real"
                imag_name = f"{baseName}_imag"

                self.variables[real_name] = self.ncout.createVariable(
                    real_name, 'f4', ('time', 'y', 'x',),
                    fill_value=0.0, zlib=True, complevel=self.compression_level,
                    chunksizes=self.chunk_size
                )
                self.variables[imag_name] = self.ncout.createVariable(
                    imag_name, 'f4', ('time', 'y', 'x',),
                    fill_value=0.0, zlib=True, complevel=self.compression_level,
                    chunksizes=self.chunk_size
                )

                try:
                    longpol = band_metadata['POLARIZATION']
                except Exception:
                    longpol = band_metadata.get('POLARISATION', '')

                self.variables[real_name].long_name = f"Real part of complex samples ({longpol})"
                self.variables[imag_name].long_name = f"Imaginary part of complex samples ({longpol})"
                self.variables[real_name].polarisation = self.variables[imag_name].polarisation = longpol

                # Write real and imaginary parts
                self.write_variable_with_preprocessing(real_name, np.real(arr).astype(np.float32),
                                                       process_chunk=None, workers=8, sync_every=None)
                self.write_variable_with_preprocessing(imag_name, np.imag(arr).astype(np.float32),
                                                       process_chunk=None, workers=8, sync_every=None)

            else:
                # Non-complex (e.g., GRD). Preserve reasonable dtype mapping and behavior
                try:
                    pol = band_metadata['POLARIZATION']
                except Exception:
                    pol = band_metadata.get('POLARISATION', '')
                varName = f"{baseName}"

                # Determine netCDF dtype and ensure array is suitable
                if arr.dtype == np.uint16:
                    nc_dtype = 'u2'
                elif arr.dtype == np.int16:
                    nc_dtype = 'i2'
                elif arr.dtype == np.uint8:
                    nc_dtype = 'u1'
                elif arr.dtype == np.int32:
                    nc_dtype = 'i4'
                else:
                    nc_dtype = 'f4'
                    arr = arr.astype(np.float32)

                self.variables[varName] = self.ncout.createVariable(
                    varName, nc_dtype, ('time', 'y', 'x',),
                    fill_value=0, zlib=True, complevel=self.compression_level,
                    chunksizes=self.chunk_size
                )
                try:
                    self.variables[varName].long_name = f"Amplitude {pol}-polarisation"
                except:
                    self.variables[varName].long_name = f"Amplitude {pol}-polarisation"
                self.variables[varName].polarisation = pol

                self.write_variable_with_preprocessing(varName, arr,
                                                       process_chunk=None, workers=8, sync_every=None)

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