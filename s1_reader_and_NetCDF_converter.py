#!/usr/bin/python3

# Name:          Sentinel1_reader_and_NetCDF_converter.py
# Purpose:       Read Sentinel-1 data from ESA SAFE and convert to
#                to netCDF
# Author(s):     Trygve Halsne
# Created:
# Modifications:
# Copyright:     (c) Norwegian Meteorological Institute, 2018
#
# Need to use gdal 2.1.1-> to have support of the SAFE reader

import subprocess
import sys
from collections import defaultdict
from datetime import datetime
from datetime import timedelta
import gdal
import netCDF4
import numpy as np
from scipy import interpolate
import pathlib
import safe_to_netcdf.utils as utils


class Sentinel1_reader_and_NetCDF_converter:
    """
        Class for reading Sentinel-1 products from SAFE with methods for
        creating variables for calibration, noise correction etc. In addition,
        it is possible to convert product into NetCDF4/CF (1.6).

        The implemented methods uses standard python libraries as
        gdal(v. > 2.1.1), numpy, lxml etc.

        Keyword arguments:
        SAFE_file -- absolute path to zipped file
        SAFE_outpath -- output storage location for unzipped SAFE product
    """

    def __init__(self, product, indir, outdir):
        self.product_id = product
        self.input_zip = (indir / product).with_suffix('.zip')
        self.SAFE_dir = (outdir / self.product_id).with_suffix('.SAFE')
        self.gcps = []  # GCPs from gdal used for generation of lat lon
        self.polarisation = []
        self.xSize = None
        self.ySize = None
        self.xmlFiles = defaultdict(list)
        self.globalAttribs = {}
        self.src = None
        self.t0 = datetime.now()
        self.ncout = None  # NetCDF output file
        self.xmlCalPixelLines = defaultdict(list)
        self.xmlCalLUTs = defaultdict(list)
        self.xmlGCPs = defaultdict(list)
        self.imageAnnotation = defaultdict(dict)
        self.noiseVectors = defaultdict(list)
        self.productMetadata = defaultdict(dict)  # list of values from image annotation files
        self.productMetadataList = defaultdict(dict)  # list of lists from image annotation files
        self.main()

    def main(self):
        """ Main method for traversing and reading key parameters from SAFE
            directory.
        """

        # 1) Fetch manifest.xml file
        self.xmlFiles['manifest'] = self.uncompress()

        # 2) Set some of the gloal parameters
        utils.initializer(self, self.xmlFiles['manifest'])

        gcps_ok = self.getGCPs()

        # Calibration tables
        calibrationTables = ['sigmaNought', 'betaNought', 'gamma', 'dn']
        for calXmlFile in self.xmlFiles['s1Level1CalibrationSchema']:
            #calibrationXmlFile = self.SAFE_dir / calXmlFile

            # Retrieve pixels and lines
            polarisation, cal_pixels, cal_lines = self.readPixelsLines(calXmlFile)
            self.xmlCalPixelLines[polarisation] = [np.array(cal_pixels, np.int16),
                                                   np.array(cal_lines, np.int16)]

            # Retrieve Look Up Tables
            for ct in calibrationTables:
                self.xmlCalLUTs[str(ct + '_' + polarisation)] = np.array(self.getCalTable(calXmlFile, ct), np.float32)

        # Retrieve thermal noise vectors
        for nXmlFile in self.xmlFiles['s1Level1NoiseSchema']:
            noiseVector, polarisation = self.readNoiseData(nXmlFile)
            self.noiseVectors[str(polarisation)] = noiseVector

        # Retrieve GCP parameters
        gcp_parameters = ['azimuthTime', 'slantRangeTime', 'line', 'pixel',
                          'latitude', 'longitude', 'height', 'incidenceAngle', 'elevationAngle']
        for xmlFile in self.xmlFiles['s1Level1ProductSchema']:

            for parameter in gcp_parameters:
                polarisation, values = self.getGCPValues(xmlFile, parameter)

                if not parameter == 'azimuthTime':
                    self.xmlGCPs[str(parameter + '_' + polarisation)] = np.array(values, np.float32)
                else:
                    self.xmlGCPs[str(parameter + '_' + polarisation)] = np.array(values, np.str)

        # retrieve product metadata from image annotation files
        productMetadata_parameters = [
            'missionId', 'productType', 'polarisation', 'mode', 'startTime',
            'stopTime', 'absoluteOrbitNumber', 'missionDataTakeId', 'imageNumber',
            'productQualityIndex', 'iInputDataMean', 'qInputDataMean',
            'inputDataMeanOutsideNominalRangeFlag',
            'iInputDataStdDev', 'qInputDataStdDev', 'inputDataStDevOutsideNominalRangeFlag',
            'numDownlinkInputDataGaps', 'downlinkGapsInInputDataSignificantFlag',
            'numDownlinkInputMissingLines', 'downlinkMissingLinesSignificantFlag',
            'numInstrumentInputDataGaps', 'instrumentGapsInInputDataSignificantFlag',
            'numInstrumentInputMissingLines', 'instrumentMissingLinesSignificantFlag',
            'numSsbErrorInputDataGaps', 'ssbErrorGapsInInputDataSignificantFlag',
            'numSsbErrorInputMissingLines', 'ssbErrorMissingLinesSignificantFlag',
            'chirpSourceUsed', 'pgSourceUsed', 'rrfSpectrumUsed', 'replicaReconstructionFailedFlag',
            'meanPgProductAmplitude', 'stdDevPgProductAmplitude', 'meanPgProductPhase',
            'stdDevPgProductPhase', 'pgProductDerivationFailedFlag',
            'invalidDownlinkParamsFlag', 'iBiasSignificanceFlag',
            'qBiasSignificanceFlag', 'iqGainSignificanceFlag',
            'iqQuadratureDepartureSignificanceFlag',
            'platformHeading', 'projection', 'rangeSamplingRate',
            'radarFrequency', 'azimuthSteeringRate', 'rangePixelSpacing',
            'azimuthPixelSpacing', 'azimuthTimeInterval', 'azimuthFrequency', 'numberOfSamples',
            'numberOfLines', 'zeroDopMinusAcqTime', 'incidenceAngleMidSwath',
            'rawDataAnalysisUsed', 'orbitDataFileUsed', 'attitudeDataFileUsed',
            'rxVariationCorrectionApplied', 'antennaElevationPatternApplied',
            'antennaAzimuthPatternApplied', 'antennaAzimuthElementPatternApplied',
            'rangeSpreadingLossCompensationApplied', 'srgrConversionApplied',
            'detectionPerformed', 'thermalNoiseCorrectionPerformed',
            'referenceRange', 'ellipsoidName', 'ellipsoidSemiMajorAxis',
            'ellipsoidSemiMinorAxis', 'bistaticDelayCorrectionApplied',
            'topsFilterConvention']

        productMetadataList_parameters = [
            'orbitList', 'attitudeList', 'noiseList',
            'terrainHeightList', 'azimuthFmRateList', 'sliceList',
            'inputDimensionsList', 'dcEstimateList', 'antennaPatternList',
            'coordinateConversionList', 'swathMergeList']

        for xmlFile in self.xmlFiles['s1Level1ProductSchema']:
            root = utils.xml_read(xmlFile)
            polarisation = root.find('.//polarisation').text

            for pm in productMetadata_parameters:
                variable = root.find(str('.//' + pm))
                self.productMetadata[polarisation][pm] = variable.text

            for pml in productMetadataList_parameters:
                variable = root.find(str('.//' + pml))
                self.extractProductMetadataList(variable, polarisation)

    def extractProductMetadataList(self, mother_element, polarisation):
        """ Write the input mother_element from the product xml annotation file
            to the extractProductMetadataList variable.

            Could/should be extended depending on peoples need.
        """
        listType = mother_element.tag
        if listType == 'orbitList':
            # {time1:[frame, position (x,y,z), velocity (x,y,z)], time2: ...}
            orbitDict = {}
            for orbit in mother_element.getchildren():
                time = orbit.find('.//time').text
                values = []
                position = (orbit.find('.//position/x').text,
                            orbit.find('.//position/y').text,
                            orbit.find('.//position/z').text)
                position = ' '.join(position)
                velocity = (orbit.find('.//velocity/x').text,
                            orbit.find('.//velocity/y').text,
                            orbit.find('.//velocity/z').text)
                velocity = ' '.join(velocity)
                values.append(orbit.find('.//frame').text)
                values.append(position)
                values.append(velocity)
                orbitDict[time] = values

            self.productMetadataList[polarisation][listType] = orbitDict
        elif listType == 'attitudeList':  # List of attitude and angular velocity annotation
            # records.
            pass
        elif listType == 'noiseList':  # List of noise packet records.
            pass
        elif listType == 'terrainHeightList':  # Terrain height list. This element is a list of
            # terrainHeight records that contain the average terrain height at the given zero
            # Doppler azimuth time. The actual terrain heights used by the IPF may represent
            # bilinearly interpolated values from this list. The list contains an entry for each
            # terrain height update made along azimuth.
            pass
        elif listType == 'azimuthFmRateList':  # Azimuth Frequency Modulation (FM) rate list.
            # This element is a list of azimuthFmRate records that contain the parameters needed
            # to calculate the azimuth FM rate. The list contains an entry for each azimuth FM
            # rate update made along azimuth.
            pass
        elif listType == 'sliceList':  # List of annotations for all slices in segment.  The
            # total size of the list represents the number of slices in the segment.  If product
            # composition type is Individual or Assembled, the total size of this list is 0.
            pass
        elif listType == 'inputDimensionsList':  # Input dimensions list. This element contains a
            # list of inputDimensions records which describe the number of input range samples
            # and azimuth lines.
            pass
        elif listType == 'dcEstimateList':  # List of Doppler centroid estimates that have been
            # calculated by the IPF during image processing. The list contains an entry for each
            # Doppler centroid estimate made along azimuth.
            pass
        elif listType == 'antennaPatternList':  # Antenna pattern list. This element is a list of
            # antennaPattern records that describe the antenna elevation pattern as it is updated
            # in azimuth. The list contains an entry for each AEP update made along azimuth.
            # {index1:['swath','azimuthTime', 'slantRangeTime', 'elevationAngle',
            # 'elevationPattern','incidenceAngle','terrainHeight', 'roll' ], index2: ...}
            antennaPatternDict = {}
            variables = ['swath', 'azimuthTime', 'slantRangeTime', 'elevationAngle',
                         'elevationPattern',
                         'incidenceAngle', 'terrainHeight', 'roll']
            index = 0
            for ap in mother_element.iterchildren():
                values = []
                for v in variables:
                    values.append(ap.find(str('.//' + v)).text)
                antennaPatternDict[index] = values
                index += 1

            # self.productMetadataList[polarisation][listType] = antennaPatternDict

        elif listType == 'coordinateConversionList':
            # {index1:[azimuthTime, slantRangeTime, sr0, srgrCoefficients, gr0,grsrCoefficients
            # ], index2: ...}
            coordinateConversionDict = {}
            variables = ['azimuthTime', 'slantRangeTime', 'sr0', 'srgrCoefficients', 'gr0',
                         'grsrCoefficients']
            index = 0
            for cc in mother_element.getchildren():
                values = []
                for v in variables:
                    values.append(cc.find(str('.//' + v)).text)
                coordinateConversionDict[index] = values
                index += 1

            self.productMetadataList[polarisation][listType] = coordinateConversionDict
        elif listType == 'swathMergeList':  # Merge information for IW and EW GRD products. This
            # list contains one record per swath.
            # {index1:{(sub)swath:['firstAzimuthLine','firstRangeSample','lastAzimuthLine',
            # 'lastRangeSample', 'azimuthTime']}, index2:{}}
            swathBoundsList = defaultdict(dict)
            counter = 0
            for swath in mother_element.findall('.//swathMerge'):
                current_swath = swath.find('.//swath').text
                for burst in swath.findall('.//swathBounds'):
                    values = []
                    values.append(burst.find('./firstAzimuthLine').text)
                    values.append(burst.find('./firstRangeSample').text)
                    values.append(burst.find('./lastAzimuthLine').text)
                    values.append(burst.find('./lastRangeSample').text)
                    values.append(burst.find('./azimuthTime').text)

                    swathBoundsList[counter][current_swath] = values
                    swathBoundsList[counter][current_swath] = values
                    counter += 1
            self.productMetadataList[polarisation][listType] = swathBoundsList
        else:
            print("Extraction of %s is not implemented" % listType)

    def uncompress(self):
        """ Uncompress SAFE zip file. Return manifest file """

        # If zip not extracted yet
        if not self.SAFE_dir.is_dir():
            cmd = f'/usr/bin/unzip {self.input_zip} -d {self.SAFE_dir.parent}'
            subprocess.call(cmd, shell=True)
            if not self.SAFE_dir.is_dir():
                print(f'Error unzipping file {self.input_zip}')
                sys.exit()

        xmlFile = self.SAFE_dir / 'manifest.safe'
        if not xmlFile.is_file():
            print(f'Manifest file not available {xmlFile}')
            sys.exit()

        return xmlFile

    def initializer(self, manifest):
        """ Traverse manifest file for setting additional parameters
        in __init__ """
        root = utils.xml_read(manifest)

        # Set xml-files
        dataObjectSection = root.find('./dataObjectSection')
        for dataObject in dataObjectSection.findall('./'):
            repID = dataObject.attrib['repID']
            ftype = None
            href = None
            for element in dataObject.iter():
                attrib = element.attrib
                if 'mimeType' in attrib:
                    ftype = attrib['mimeType']
                if 'href' in attrib:
                    href = attrib['href'][1:]
            if ftype == 'text/xml' and href:
                self.xmlFiles[repID].append(self.SAFE_dir / href[1:])

        # Set gdal object
        self.src = gdal.Open(str(self.xmlFiles['manifest']))

        # Set raster size parameters
        self.xSize = self.src.RasterXSize
        self.ySize = self.src.RasterYSize

        # Set global metadata attributes
        self.globalAttribs = self.src.GetMetadata()

        polarisations = root.findall('.//s1sarl1:transmitterReceiverPolarisation',
                                     namespaces=root.nsmap)
        for polarisation in polarisations:
            self.polarisation.append(polarisation.text)
        self.globalAttribs['polarisation'] = self.polarisation

        self.globalAttribs['ProductTimelinessCategory'] = root.find(
            './/s1sarl1:productTimelinessCategory', namespaces=root.nsmap).text

        return True

    def getGCPs(self):
        """ Get product GCPs utilizing gdal """
        if self.src:
            self.gcps = self.src.GetGCPs()
            return True
        else:
            return False

    def write_to_NetCDF(self, nc_outpath, compression_level, chunk_size=(1, 31, 33)):
        """ Method intitializing output NetCDF product.

        Keyword arguments:
        nc_outpath -- output path where NetCDF file should be stored
        compression_level -- compression level on output NetCDF file (1-9)
        """

        print("------------START CONVERSION FROM SAFE TO NETCDF-------------")

        # Status
        utils.memory_use(self.t0)

        out_netcdf = (nc_outpath / self.product_id).with_suffix('.nc')
        ncout = netCDF4.Dataset(out_netcdf, 'w', format='NETCDF4')
        ncout.createDimension('time', size=None)
        ncout.createDimension('x', self.xSize)
        ncout.createDimension('y', self.ySize)

        nclat = ncout.createVariable('lat', 'f4', ('y', 'x',), zlib=True,
                                     complevel=compression_level, chunksizes=chunk_size[1:])
        nclon = ncout.createVariable('lon', 'f4', ('y', 'x',), zlib=True,
                                     complevel=compression_level, chunksizes=chunk_size[1:])

        # Set time value
        utils.create_time(ncout, self.globalAttribs["ACQUISITION_START_TIME"])

        # Add latitude and longitude layers
        ##########################################################
        # Status
        utils.memory_use(self.t0)

        lat, lon = self.genLatLon_regGrid()  # Assume gcps are on a regular grid
        nclat.long_name = 'latitude'
        nclat.units = 'degrees_north'
        nclat.standard_name = 'latitude'
        nclat[:, :] = lat

        nclon.long_name = 'longitude'
        nclon.units = 'degrees_east'
        nclon.standard_name = 'longitude'
        nclon[:, :] = lon

        # Add raw measurement layers
        ##########################################################
        # Status
        utils.memory_use(self.t0)

        for i in range(1, self.src.RasterCount + 1):
            band = self.src.GetRasterBand(i)
            band_metadata = band.GetMetadata()
            varName = 'Amplitude_%s' % band_metadata['POLARISATION']
            var = ncout.createVariable(varName, 'u2', ('time', 'y', 'x',),
                                       fill_value=0, zlib=True, complevel=compression_level,
                                       chunksizes=chunk_size)
            var.long_name = 'Amplitude %s-polarisation' % band_metadata['POLARISATION']
            var.units = "1"
            var.coordinates = "lat lon"
            var.grid_mapping = "crsWGS84"
            var.standard_name = "surface_backwards_scattering_coefficient_of_radar_wave"
            var.polarisation = "%s" % band_metadata['POLARISATION']
            print((band.GetVirtualMemArray().shape))
            var[0, :, :] = band.GetVirtualMemArray()

        # set grid mapping(?)
        ##########################################################
        nc_crs = ncout.createVariable('crsWGS84', np.int32)
        nc_crs.grid_mapping_name = "latitude_longitude"
        nc_crs.semi_major_axis = "6378137"
        nc_crs.inverse_flattening = "298.2572235604902"

        # Add calibration layers
        ##########################################################
        # Status
        print('\nAdding calibration layers')
        utils.memory_use(self.t0)

        for calibration in self.xmlCalLUTs:
            current_polarisation = calibration.split('_')[-1]
            pixels, lines = self.xmlCalPixelLines[current_polarisation]
            calibration_LUT = self.xmlCalLUTs[calibration]
            resampled_calibration = self.getCalLayer(pixels, lines, calibration_LUT)

            var = ncout.createVariable(str(calibration), 'f4', ('time', 'y', 'x',),
                                       zlib=True, complevel=compression_level,
                                       chunksizes=chunk_size)
            var.long_name = '%s calibration table' % calibration
            var.units = "1"
            var.coordinates = "lat lon"
            var.grid_mapping = "crsWGS84"
            var.polarisation = "%s" % current_polarisation
            var[0, :, :] = resampled_calibration

        # Add noise layers
        ##########################################################
        # Status
        print('\nAdding noise layers')
        utils.memory_use(self.t0)

        for polarisation in self.polarisation:
            noiseCorrectionMatrix = self.getNoiseCorrectionMatrix(self.noiseVectors[polarisation],
                                                                  polarisation)

            var = ncout.createVariable(str('noiseCorrectionMatrix_' + polarisation), 'f4',
                                       ('time', 'y', 'x',),
                                       zlib=True, complevel=compression_level,
                                       chunksizes=chunk_size)
            var.long_name = 'Thermal noise correction vector power values.'
            var.units = "1"
            var.coordinates = "lat lon"
            var.grid_mapping = "crsWGS84"
            var.polarisation = "%s" % polarisation
            var[0, :, :] = noiseCorrectionMatrix

        # Add subswath layers
        ##########################################################
        # Status
        print('\nAdding subswath layers')
        utils.memory_use(self.t0)

        for polarisation in self.polarisation:
            swathLayer, flags = self.getSwathList(polarisation)
            flag_values = np.array(sorted(flags.values()), dtype=np.int8)
            flags_meanings = ""
            for key in sorted(flags.keys()):
                flags_meanings += str(key + ' ')

            swathList = ncout.createVariable('swathList', 'i1', ('y', 'x',), fill_value=0,
                                             zlib=True, complevel=7, chunksizes=chunk_size[1:])
            swathList.long_name = 'Subswath List'
            swathList.flag_values = flag_values
            swathList.valid_range = np.array([flag_values.min(), flag_values.max()])
            swathList.flag_meanings = flags_meanings.strip()
            swathList.standard_name = "status_flag"
            swathList.units = "1"
            swathList.coordinates = "lat lon"
            swathList.grid_mapping = "crsWGS84"
            # swathList.polarisation = "%s" %  polarisation
            swathList[:] = swathLayer
            break

        # Add GCP information
        ##########################################################
        # Status
        print('\nAdding GCP information')
        utils.memory_use(self.t0)

        gcp_units = {'slantRangeTime': 's', 'latitude': 'degrees', 'longitude': 'degrees',
                     'height': 'm', 'incidenceAngle': 'degrees', 'elevationAngle': 'degrees'}
        gcp_long_name = {
            'azimuthTime': 'Zero Doppler azimuth time to which grid point applies [UTC].',
            'slantRangeTime': 'Two way slant range time to grid point.',
            'line': 'Reference image MDS line to which this geolocation grid point applies.',
            'pixel': 'Reference image MDS sample to which this geolocation grid point applies',
            'latitude': 'Geodetic latitude of grid point.',
            'longitude': 'Geodetic longitude of grid point.',
            'height': 'Height of the grid point above sea level.',
            'incidenceAngle': 'Incidence angle to grid point.',
            'elevationAngle': 'Elevation angle to grid point.'}
        ncout.createDimension('gcp_index', len(self.gcps))
        for key, value in self.xmlGCPs.items():
            current_variable = key.split('_')[0]
            if current_variable == 'azimuthTime':
                var = ncout.createVariable(str('GCP_%s' % key), 'f4', ('gcp_index'), zlib=True)
                dates = np.array([datetime.strptime(t, '%Y-%m-%dT%H:%M:%S.%f') for t in value])
                ref_date = dates.min()
                value = np.array([td.total_seconds() for td in dates - ref_date])
                var.units = 's'
                var.long_name = gcp_long_name[current_variable]
                var.comment = 'Seconds since %s' % ref_date.strftime('%Y-%m-%dT%H:%M:%S.%f')
            else:
                var = ncout.createVariable(str('GCP_%s' % key), value.dtype, ('gcp_index'),
                                           zlib=True)
                if current_variable in gcp_units:
                    var.units = gcp_units[current_variable]
                var.long_name = gcp_long_name[current_variable]
            var[:] = value

        # Add product annotation metadata
        ##########################################################
        # Status
        print('\nAdding annotation information')
        utils.memory_use(self.t0)

        for polarisation in self.productMetadata:
            varBaseName = str('s1Level1ProductSchema_' + polarisation)
            productMetadata = self.productMetadata[polarisation]
            var = ncout.createVariable(varBaseName, 'i1')
            var.setncatts(productMetadata)

        # Add product annotation metadata lists
        ##########################################################
        # Status
        print('\nAdding annotation list information')
        utils.memory_use(self.t0)

        productMetadataListComment = {
            'swathMergeList': 'index:{swath:[firstAzimuthLine, firstRangeSample, lastAzimuthLine, '
                              'lastRangeSample, azimuthTime]}',
            'orbitList': 'time:[frame, position (x,y,z), velocity (x,y,z)]',
            'coordinateConversionList': 'index:[azimuthTime, slantRangeTime, sr0, '
                                        'srgrCoefficients, gr0, grsrCoefficients ]',
            'antennaPatternList': 'index:[swath, azimuthTime, slantRangeTime, elevationAngle, '
                                  'elevationPattern, incidenceAngle, terrainHeight, roll]'
            }
        productMetadataListUnits = {
            'swathMergeList': 'index: , firstAzimuthLine: , firstRangeSample: , lastAzimuthLine: '
                              ', lastRangeSample: , azimuthTime: datetime'}
        productMetadataListDatatype = {
            'swathMergeList': 'index:uint16 , firstAzimuthLine:unit32 , firstRangeSample:unit32 , '
                              'lastAzimuthLine:uint32 , lastRangeSample:uint32 , azimuthTime: UTC'}

        for polarisation in self.productMetadataList:
            for subkey in self.productMetadataList[polarisation]:
                varBaseName = str(subkey + '_' + polarisation)
                if True:
                    productMetadataList = self.productMetadataList[polarisation][subkey]
                    tmp_dict = {}
                    for k, v in productMetadataList.items():
                        tmp_dict[str(k)] = str(v)
                    var = ncout.createVariable(varBaseName, 'i1')
                    var.comment = productMetadataListComment[subkey]
                    var.setncatts(tmp_dict)

        # Add global attributes
        ##########################################################
        # Status
        print('\nAdding global attributes')
        utils.memory_use(self.t0)

        nowstr = self.t0.strftime("%Y-%m-%dT%H:%M:%SZ")
        ncout.title = 'Sentinel-1 GRD data'
        ncout.netcdf4_version_id = netCDF4.__netcdf4libversion__
        ncout.file_creation_date = nowstr

        self.globalAttribs['Conventions'] = "CF-1.6"
        self.globalAttribs['summary'] = 'Sentinel-1 C-band SAR GRD product.'
        self.globalAttribs[
            'keywords'] = '[Earth Science, Spectral/Engineering, RADAR, RADAR backscatter], ' \
                          '[Earth Science, Spectral/Engineering, RADAR, RADAR imagery], ' \
                          '[Earth Science, Spectral/Engineering, Microwave, Microwave Imagery]'
        self.globalAttribs['keywords_vocabulary'] = "GCMD Science Keywords"
        self.globalAttribs['institution'] = "Norwegian Meteorological Institute"
        self.globalAttribs['history'] = nowstr + ". Converted from SAFE to NetCDF by NBS team."
        ncout.setncatts(self.globalAttribs)
        ncout.sync()

        # self.ncout = ncout

        # Status
        ncout.close()
        print('\nFinished.')
        utils.memory_use(self.t0)

        return out_netcdf.is_file()

    def readNoiseData(self, xmlfile):
        """ Method for reading noise data from Sentinel-1 annotation files.
            This method supports both the thermal noise denoising conventions
            i.e. old convention applying range denoising and new convention
            adding denoising in azimuth direction.

            Structure on return arguments:

            Returns either the combination of range and azimuth noise variables
            or range and swathBounds for applying denoising per subswath.

            Range:
                Keys: azimuthTime
                Values: ['line','pixels','noiseLut']

            Azimuth:
                Keys: index/number
                Values: Subdirectory
                Subdirectory.keys: '(sub)swath'
                Subdirectory.values: ['firstAzimuthLine','firstRangeSample',
                         'lastAzimuthLine', 'lastRangeSample', 'line',
                         'noiseLut']

            swathBounds:
                Keys: index/number
                Values: Subdirectory
                Subdirectory.keys: '(sub)swath'
                Subdirectory.values: ['firstAzimuthLine','firstRangeSample',
                         'lastAzimuthLine', 'lastRangeSample', 'azimuthTime']

            Returns noiseRangeAndAzimuthList:
                Dictionary with:
                Keys: 'range' and 'azimuth'/'swathBounds'
                values: noiseRangeVectorList, noiseAzimuthVectorList
        """

        root = utils.xml_read(xmlfile)
        polarisation = root.find('.//polarisation').text

        # Noise in range direction
        noiseRangeVectorList = defaultdict(list)

        nrvl = root.findall('.//noiseRangeVector')  # Noise Range Vector List
        old_convention = False
        if nrvl == []:  # Following old convention
            nrvl = root.findall('.//noiseVector')
            range_variables = ['line', 'pixel', 'noiseLut']
            old_convention = True
        else:
            range_variables = ['line', 'pixel', 'noiseRangeLut']

        for noiseRangeVector in nrvl:
            azimuthTime = noiseRangeVector.find('./azimuthTime').text
            for variable in range_variables:
                current_element = noiseRangeVector.find(variable)
                noiseRangeVectorList[azimuthTime].append(current_element.text)

        # Read product annotation file for same polarisation to retrieve
        # image annotation information. LAD - Level 1 Annotation Data Set
        LAD_variables = ['azimuthSteeringRate',
                         # skipped dataDcPolynomial, azimuthFmRate and velocity (x,y,z)
                         'radarFrequency', 'linesPerBurst',
                         'azimuthTimeInterval', 'productFirstLineUtcTime']
        for LAD in self.xmlFiles['s1Level1ProductSchema']:
            lad_root = utils.xml_read(LAD)
            lad_polarisation = lad_root.find('.//polarisation').text
            if lad_polarisation == polarisation:
                for lad_var in LAD_variables:
                    self.imageAnnotation[polarisation][lad_var] = lad_root.find(
                        './/%s' % lad_var).text
                break

        # If old NADS (Noise Annotation Data Set), read swath bounds from LAD
        if old_convention:
            swathMergeList = lad_root.find('.//swathMergeList')
            swathBoundsList = defaultdict(dict)
            counter = 0
            for swath in swathMergeList.findall('.//swathMerge'):
                current_swath = swath.find('.//swath').text
                for burst in swath.findall('.//swathBounds'):
                    values = []
                    values.append(burst.find('./firstAzimuthLine').text)
                    values.append(burst.find('./firstRangeSample').text)
                    values.append(burst.find('./lastAzimuthLine').text)
                    values.append(burst.find('./lastRangeSample').text)
                    values.append(burst.find('./azimuthTime').text)

                    swathBoundsList[counter][current_swath] = values
                    counter += 1

        if not old_convention:
            # Noise in azimuth direction
            noiseAzimuthVectorList = defaultdict(dict)
            azimuth_variables = ['line', 'noiseAzimuthLut']

            navl = root.findall('.//noiseAzimuthVector')  # Noise Range Vector List

            for i, noiseAzimuthVector in enumerate(navl):
                current_swath = noiseAzimuthVector.find('./swath').text
                values = []
                values.append(noiseAzimuthVector.find('./firstAzimuthLine').text)
                values.append(noiseAzimuthVector.find('./firstRangeSample').text)
                values.append(noiseAzimuthVector.find('./lastAzimuthLine').text)
                values.append(noiseAzimuthVector.find('./lastRangeSample').text)

                for variable in azimuth_variables:
                    current_element = noiseAzimuthVector.find(variable)
                    values.append(current_element.text)

                noiseAzimuthVectorList[i][current_swath] = values  # append(current_element.text)

        # Packing all variables to one dictionary
        if old_convention:
            noiseRangeAndAzimuthList = {'range': noiseRangeVectorList,
                                        'swathBounds': swathBoundsList}
        else:
            noiseRangeAndAzimuthList = {'range': noiseRangeVectorList,
                                        'azimuth': noiseAzimuthVectorList}

        return noiseRangeAndAzimuthList, polarisation

    def readPixelsLines(self, xmlfile):

        root = utils.xml_read(xmlfile)
        polarisation = root.find('.//polarisation').text

        # Get pixels where we have calibration values. Assume regular distripution over all image
        pixels = []
        npix = []
        line = []
        lines = []

        # Get line nublers for calibration values
        for l in root.iter('line'):
            line.append(l.text)

        numLine = 0
        for p in root.iter('pixel'):
            p = p.text.split(' ')
            npix.append(len(p))
            [pixels.append(x) for x in p]
            [lines.append(line[numLine]) for i in range(len(p))]
            numLine += 1

        if len(lines) != len(pixels):
            print(
                'Error: Wrong size of arrays. legth of pixels and lines should be the same '
                'pixels=%d lines=%d' % (
                len(pixels), len(lines)))
            return -1

        return (polarisation, pixels, lines)

    def getCalTable(self, xmlfile, tName):
        root = utils.xml_read(xmlfile)

        # Get calibration values
        cal = []
        for c in root.iter(tName):
            c = c.text.split(' ')
            [cal.append(x) for x in c]

        return cal

    def getCalLayer(self, pixels, lines, cal_table):
        """ Interpolate calibration layer"""
        xSize = self.xSize
        ySize = self.ySize

        nb_pixels = (pixels == 0).sum()
        nb_lines = (lines == 0).sum()
        calibration_table = cal_table.reshape(nb_pixels, nb_lines)
        tck = interpolate.RectBivariateSpline(lines[0::nb_lines], pixels[0:nb_lines],
                                              calibration_table)
        x = list(range(0, xSize))
        y = list(range(0, ySize))
        lutOut = tck(y, x)
        return lutOut

    def getGCPValues(self, xmlfile, parameter):
        """ Method for retrieving Geo Location Point parameter from xml file."""
        root = utils.xml_read(xmlfile)
        polarisation = root.find('.//polarisation').text

        #
        out_list = []

        # Get parameter values in ground control points
        for l in root.iter('geolocationGridPoint'):
            out_list.append(l.find(parameter).text)

        return polarisation, out_list

    def genLatLon_regGrid(self):
        """ Method providing latitude and longitude arrays """
        # Extract GCPs to vector arrays
        gcps = self.gcps
        xsize = self.xSize
        ysize = self.ySize
        ngcp = len(gcps)
        # print ngcp
        x = []
        y = []
        lon = []
        lat = []
        idx = 0

        for gcp in gcps:
            if gcp.GCPPixel == 0:
                y.append(gcp.GCPLine)
            if gcp.GCPLine == 0:
                x.append(gcp.GCPPixel)
            lon.append(gcp.GCPX)
            lat.append(gcp.GCPY)
        x = np.array(x, np.int32)
        y = np.array(y, np.int32)
        lat = np.array(lat, np.float32)
        lon = np.array(lon, np.float32)

        xi = list(range(0, xsize))
        yi = list(range(0, ysize))
        tck = interpolate.RectBivariateSpline(y, x, lat.reshape(len(y), len(x)))
        latitude = tck(yi, xi)
        tck = interpolate.RectBivariateSpline(y, x, lon.reshape(len(y), len(x)))
        longitude = tck(yi, xi)

        return latitude, longitude

    def readSwathList(self, noiseVector):  # ,imageAnnotationDict):
        """ Returns dictionary with swath ID as key and number of azimuth denoising
            blocks as value.

            Keyword values:
            noiseVector -- Noise vector containing range and azimuth noise
                           values retrieved from the readNoiseData method.
        """
        noiseRangeVectorList = noiseVector['range']
        if 'azimuth' in noiseVector:
            noiseAzimuthVectorList = noiseVector['azimuth']
        else:
            noiseAzimuthVectorList = noiseVector['swathBounds']

        # Set swath list (i.e. unique swath IDs)
        swathList = np.array([])
        for noiseAzimuthVector in list(noiseAzimuthVectorList.values()):
            swathList = np.append(swathList, list(noiseAzimuthVector.keys()))

        swathList, swathListCounts = np.unique(swathList, return_counts=True)
        swathList = dict(list(zip(swathList, swathListCounts)))
        return swathList

    def getSwathList(self, polarisation):
        """ Returns swathList as raster layer.

            Keyword values:
            polarisation -- polarisation of rasterband
            subswath_flag -- flags for subswath
        """

        swathMergeList = self.productMetadataList[polarisation]['swathMergeList']

        if self.globalAttribs['MODE'] == 'EW':
            subswath_flag = {'EW1': 1, 'EW2': 2, 'EW3': 3, 'EW4': 4, 'EW5': 5}
        elif self.globalAttribs['MODE'] == 'IW':
            subswath_flag = {'IW1': 1, 'IW2': 2, 'IW3': 3}
        else:
            print("Undefined mode %s" % self.globalAttribs['MODE'])
            return 0

        swathListRaster = np.zeros((self.ySize, self.xSize))
        for index, subswath in swathMergeList.items():
            for key, value in subswath.items():
                firstAzimuthLine, firstRangeSample, lastAzimuthLine, lastRangeSample, azimuthTime\
                    = value
                firstAzimuthLine = int(firstAzimuthLine)
                firstRangeSample = int(firstRangeSample)
                lastAzimuthLine = int(lastAzimuthLine)
                lastRangeSample = int(lastRangeSample)
                swathListRaster[firstAzimuthLine:lastAzimuthLine + 1,
                firstRangeSample:lastRangeSample + 1] = subswath_flag[key]

        return swathListRaster, subswath_flag

    def getNoiseRangeRecordsInInterval(self, noiseRangeVectors, noiseAzimuthVectorStart,
                                       noiseAzimuthVectorStop):
        """ Method for the retrieval of the noise Range records in the current
            azimuth denoising block, according to the:
            'Thermal Denoising of Products Generated by the S-1 IPF.'

            Keyword values:

            noiseRangeVectors -- list of timestamps for all noise range vectors
            noiseAzimuthVectorStart -- start time for azimuth block
            noiseAzimuthVectorStop -- stop time for azimuth block
        """
        validNoiseRangeVectors = []
        first_index = None
        for index, vector in enumerate(noiseRangeVectors):
            vector_time = datetime.strptime(vector, '%Y-%m-%dT%H:%M:%S.%f')

            if noiseAzimuthVectorStart <= vector_time <= noiseAzimuthVectorStop:
                validNoiseRangeVectors.append(vector)
                if not first_index:
                    first_index = index

        return validNoiseRangeVectors, first_index

    def getNearestRangeRecordInInterval(self, noiseRangeKeys, blockCenterTime,
                                        currentSwathStartTime, currentSwathEndTime):
        """ Returns the key for the noise range record in the current
            swath closest to the center of the azimuth block, according to the:
            'Thermal Denoising of Products Generated by the S-1 IPF.'

            Keyword values:

            noiseRangeKeys -- datestamps for noise range records
            blockCenterTime -- center of the noise azimuth block
            currentSwathStartTime -- start time for the current sub-swath
            currentSwathEndTime -- end time for the current sub-swath
        """

        valid_vectors, index = self.getNoiseRangeRecordsInInterval(noiseRangeKeys,
                                                                   currentSwathStartTime,
                                                                   currentSwathEndTime)
        valid_vectors_dt = [datetime.strptime(vector, '%Y-%m-%dT%H:%M:%S.%f') for vector in
                            valid_vectors]

        nearest_record = [min(valid_vectors_dt, key=lambda x: abs(x - blockCenterTime)).strftime(
            "%Y-%m-%dT%H:%M:%S.%f")]
        index = valid_vectors.index(nearest_record[0])
        return nearest_record, index

    def getNoiseCorrectionMatrix(self, noiseAzimuthAndRangeVectorList, polarisation):
        """ Returns the thermal noise correction matrix according to the:
            'Thermal Denoising of Products Generated by the S-1 IPF.'

            Keyword values:

            noiseAzimuthAndRangeVectorList -- defaultdict with range and
            azimuth/swathBounds as keys, depending on old or new denoising
            format.

            polarisation -- polarisation
        """
        t0_duration = datetime.now()
        imageAnnotation = self.imageAnnotation[polarisation]

        old_convention = False
        if 'swathBounds' in noiseAzimuthAndRangeVectorList:  # old NADS
            noiseAzimuthVectorList = noiseAzimuthAndRangeVectorList['swathBounds']
            old_convention = True
        else:
            noiseAzimuthVectorList = noiseAzimuthAndRangeVectorList['azimuth']

        noiseRangeVectorList = noiseAzimuthAndRangeVectorList['range']

        swathList = self.readSwathList(noiseAzimuthAndRangeVectorList)
        t0 = datetime.strptime(imageAnnotation['productFirstLineUtcTime'], '%Y-%m-%dT%H:%M:%S.%f')
        delta_ts = float(imageAnnotation['azimuthTimeInterval'])  # [s]
        noiseAzimuthMatrix = np.zeros((self.ySize, self.xSize))
        noiseRangeMatrix = np.zeros((self.ySize, self.xSize))

        # Deciding current swath time interval
        for swath_, swathCount_ in swathList.items():
            currentSwathStartTime = None
            currentSwathEndTime = None
            for noiseAzimuthVector_id, noiseAzimuthVector in noiseAzimuthVectorList.items():
                if swath_ in noiseAzimuthVector:
                    values = noiseAzimuthVector[swath_]
                    firstAzimuthLine = int(values[0])
                    lastAzimuthLine = int(values[2])
                    firstRangeSample = int(values[1])
                    lastRangeSample = int(values[3])
                    # line = np.array(values[4].split(),np.int)
                    # noiseAzimuthLUT = np.array(values[5].split(),np.float)
                    noiseAzimuthVectorStart = t0 + timedelta(seconds=firstAzimuthLine * delta_ts)
                    noiseAzimuthVectorStop = t0 + timedelta(seconds=lastAzimuthLine * delta_ts)
                    if not currentSwathStartTime:
                        currentSwathStartTime = noiseAzimuthVectorStart
                    else:
                        if currentSwathStartTime > noiseAzimuthVectorStart:
                            currentSwathStartTime = noiseAzimuthVectorStart

                    if not currentSwathEndTime:
                        currentSwathEndTime = noiseAzimuthVectorStop
                    else:
                        if currentSwathEndTime < noiseAzimuthVectorStop:
                            currentSwathEndTime = noiseAzimuthVectorStop

        for swath_, swathCount_ in swathList.items():

            # STEP 1
            for noiseAzimuthVector_id, noiseAzimuthVector in noiseAzimuthVectorList.items():
                if swath_ in noiseAzimuthVector:
                    values = noiseAzimuthVector[swath_]
                    firstAzimuthLine = int(values[0])
                    lastAzimuthLine = int(values[2])
                    firstRangeSample = int(values[1])
                    lastRangeSample = int(values[3])
                    noiseAzimuthVectorStart = t0 + timedelta(seconds=firstAzimuthLine * delta_ts)
                    noiseAzimuthVectorStop = t0 + timedelta(seconds=lastAzimuthLine * delta_ts)
                    sampleIndex = np.arange(firstRangeSample, lastRangeSample + 1)
                    lineIndex = np.arange(firstAzimuthLine, lastAzimuthLine + 1)
                    numberOfSamples = lastRangeSample - firstRangeSample + 1
                    numberOfLines = lastAzimuthLine - firstAzimuthLine + 1

                    if not old_convention:
                        line = np.array(values[4].split(), np.int)
                        noiseAzimuthLUT = np.array(values[5].split(), np.float)
                        # print noiseAzimuthVector_id,lineIndex,line, noiseAzimuthLUT
                        if len(line) > 1:
                            intp1 = interpolate.interp1d(line, noiseAzimuthLUT,
                                                         fill_value='extrapolate')
                            noiseAzimuthVector_ = intp1(lineIndex)
                        else:
                            noiseAzimuthVector_ = np.zeros(1)
                            noiseAzimuthVector_[:] = noiseAzimuthLUT[0]
                        # print len(noiseAzimuthVector_)
                        # plt.plot(line, noiseAzimuthLUT, 'go-', lineIndex, noiseAzimuthVector_,
                        # '-')
                        # plt.legend(['sub-sampled','interpolated'])
                        # plt.show()

                        # create row vector
                        # noiseAzimuthVector_ = noiseAzimuthVector_.T
                        # print sampleIndex
                        # print numberOfLines
                        # noiseAzimuthMatrix[lineIndex,sampleIndex] = np.tile(
                        # noiseAzimuthVector_,(numberOfSamples,1))
                        # noiseAzimuthMatrix[firstAzimuthLine:lastAzimuthLine+1,
                        # firstRangeSample:lastRangeSample+1]=noiseAzimuthVector_
                        noiseAzimuthMatrix[firstAzimuthLine:lastAzimuthLine + 1,
                        firstRangeSample:lastRangeSample + 1] = np.tile(noiseAzimuthVector_,
                                                                        (numberOfSamples, 1)).T
                    else:
                        noiseAzimuthMatrix[:] = 1

                    # STEP 2
                    # Parsing the range denoising record
                    noiseRangeKeys = sorted(noiseRangeVectorList.keys())  # NOTE Is this OK?
                    validRangeVectorKeys, noiseRangeVectorFirstIndex = \
                        self.getNoiseRangeRecordsInInterval(
                        noiseRangeKeys, noiseAzimuthVectorStart, noiseAzimuthVectorStop)
                    if len(validRangeVectorKeys) == 0:
                        blockCenterTime = noiseAzimuthVectorStart + (
                                    noiseAzimuthVectorStop - noiseAzimuthVectorStart) / 2
                        validRangeVectorKeys, noiseRangeVectorFirstIndex = \
                            self.getNearestRangeRecordInInterval(
                            noiseRangeKeys,
                            blockCenterTime, currentSwathStartTime,
                            currentSwathEndTime)
                        if validRangeVectorKeys == None:
                            print('Error. No valid vector found.')
                            sys.exit([1])

                    noiseRangeVectorList_ = np.zeros((numberOfSamples, len(validRangeVectorKeys)))
                    # noiseRangeVectorLine_ = np.zeros((numberOfSamples,
                    # len(validRangeVectorKeys))) #should be numberOfLines?
                    noiseRangeVectorLine_ = np.zeros((len(validRangeVectorKeys)))

                    # print validRangeVectorKeys, noiseRangeVectorFirstIndex
                    for index, key in enumerate(validRangeVectorKeys):
                        rangeRecordIndex = index + noiseRangeVectorFirstIndex
                        rangeRecPixels_ = np.array(noiseRangeVectorList[key][1].split(),
                                                   np.int)  # getNoiseRangeRecordByIndex (
                        # rangeRecordIndex)
                        rangeRecLines_ = np.array(noiseRangeVectorList[key][2].split(), np.float)
                        rangePixelToInterp_0 = np.argwhere(
                            rangeRecPixels_ >= firstRangeSample).min()
                        rangePixelToInterp_n = np.argwhere(rangeRecPixels_ <= lastRangeSample).max()
                        # rangePixelToInterp = rangeRecPixels_[
                        # rangePixelToInterp_0:rangePixelToInterp_n+1]
                        rangePixelToInterp = rangeRecPixels_[
                                             rangePixelToInterp_0:rangePixelToInterp_n + 1]
                        # rangePixelToInterp = rangeRecPixels_[(rangeRecPixels_ >=
                        # firstRangeSample) & (rangeRecPixels_ <= lastRangeSample)]

                        intp1_range = interpolate.interp1d(rangePixelToInterp, rangeRecLines_[
                                                                               rangePixelToInterp_0:rangePixelToInterp_n + 1],
                                                           fill_value='extrapolate')
                        noiseRangeVectorList_[:, index] = intp1_range(sampleIndex)

                        noiseRangeVectorLine_[index] = np.int(noiseRangeVectorList[key][0])

                    # STEP 3
                    # Generate range/azimuth denoising correction
                    noiseRangeMatrix_ = np.zeros((numberOfSamples, numberOfLines))

                    for i in range(numberOfSamples):
                        if len(noiseRangeVectorLine_) > 1:
                            intp1_line = interpolate.interp1d(noiseRangeVectorLine_,
                                                              noiseRangeVectorList_[i, :],
                                                              fill_value='extrapolate')
                            noiseRangeMatrix_[i, :] = intp1_line(lineIndex)
                        else:
                            noiseRangeMatrix_[i, :] = noiseRangeVectorList_[
                                i]  # Not able to perform azimuth interpolation. Hence writing the same value to each line index.

                    noiseRangeMatrix[lineIndex[0]:lineIndex[-1] + 1,
                    sampleIndex[0]:sampleIndex[-1] + 1] = noiseRangeMatrix_.T

        noiseCorrectionMatrix_ = noiseRangeMatrix * noiseAzimuthMatrix
        print("Created noise correction matrix in: ", datetime.now() - t0_duration)
        return noiseCorrectionMatrix_


if __name__ == '__main__':

    workdir = pathlib.Path('/home/elodief/Data/NBS')

    products = ['S1B_IW_GRDM_1SDV_20201029T050332_20201029T050405_024023_02DA93_3C79',
                'S1B_EW_GRDM_1SDH_20201029T081927_20201029T082027_024025_02DAA1_4926']

    for product in products:

        outdir = workdir / 'NBS_test_data' / 'safe2nc_latest_local_03' / product
        conversion_object = Sentinel1_reader_and_NetCDF_converter(
            product=product,
            indir=workdir / 'NBS_reference_data' / 'reference_datain_local',
            outdir=outdir)
        conversion_object.write_to_NetCDF(outdir, 7)
