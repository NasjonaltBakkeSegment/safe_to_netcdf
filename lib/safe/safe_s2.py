import subprocess as sp
import shutil
from collections import defaultdict
import logging
from pathlib import Path
import lxml.etree as ET
import rasterio
import numpy as np

try:
    from lib.utils import xml_read
except:
    from utils import xml_read

try:
    from lib.safe.safe_base import SAFEFile
except ModuleNotFoundError:
    from safe_base import SAFEFile

logger = logging.getLogger(__name__)

class S2SAFEFile(SAFEFile):
    """
    Subclass for working with Sentinel-1 SAFE files
    """

    def __init__(self, product, indir, outdir):

        super().__init__(product, indir, outdir)

        self.baseline = self.product_id.split('_')[3]
        self.processing_level = 'Level-' + self.product_id.split('_')[1][4:6]
        self.imageFiles = defaultdict(list)
        self.reference_band = None
        self.dterrengdata = False  # variable saying if products is Norwegian DEM L1C
        self.sunAndViewAngles = defaultdict(list)
        self.vectorInformation = defaultdict(list)
        self.image_list_dterreng = []

    def prepare_for_use(self):
        """
        Prepare the SAFEFile instance for use.
        This includes unzipping the product, reading metadata, creating file lists, and setting up the rasterio source.
        """
        super().prepare_for_use()

        # TODO: None of this has been written for S2
        # Read the metadata XML
        self.read_metadata_xml()

        # Create the file lists
        self.create_file_lists()

        # Set up the rasterio source and metadata
        self.setup_rasterio_source()

        # Extract product metadata
        self.get_metadata()

    def _get_xml_path(self):
        xmlFile = self.SAFE_dir / 'manifest.safe'
        if not xmlFile.is_file():
            xmlFile = self.SAFE_dir / 'MTD_MSIL1C.xml'
            if not xmlFile.is_file():
                logger.error(f'Main file not found. Exiting')
                raise
            self.dterrengdata = True

        logger.debug(f'Main file: {xmlFile}')
        return xmlFile

    def delete_uncompressed_folder(self):
        """
        Delete the uncompressed SAFE folder and its contents.
        """
        if self.SAFE_dir.is_dir():
            try:
                logger.debug(f"Deleting uncompressed folder: {self.SAFE_dir}")
                shutil.rmtree(self.SAFE_dir)
                logger.debug("Uncompressed folder deleted successfully.")
            except Exception as e:
                logger.error(f"Error deleting uncompressed folder: {e}")
                raise
        else:
            logger.debug(f"Uncompressed folder does not exist: {self.SAFE_dir}")

    def create_file_lists(self):
        """
        Create lists of XML/GML and image files included in the SAFE product.
        """
        if self.dterrengdata:
            # For DTERR data, manually parse the SAFE directory for relevant files
            for file_path in self.SAFE_dir.rglob('*'):
                if file_path.suffix in {'.xml', '.gml'}:
                    self.xmlFiles[file_path.stem] = file_path
                elif file_path.suffix == '.jp2' and "IMG_DATA" in str(file_path):
                    self.image_list_dterreng.append(file_path)
        else:
            # Parse the metadata XML tree to identify files
            dataObjectSection = self.root.find('./dataObjectSection')
            for dataObject in dataObjectSection.findall('./'):
                repID = dataObject.attrib.get('ID')
                ftype, href = None, None
                for element in dataObject.iter():
                    ftype = element.attrib.get('mimeType', ftype)
                    href = element.attrib.get('href', href)
                if href:
                    href_path = self.SAFE_dir / href[1:]  # Convert href to full path
                    if ftype in {'text/xml', 'application/xml'}:
                        self.xmlFiles.update({repID: href_path})
                    elif ftype == 'application/octet-stream':
                        self.imageFiles[repID] = href_path


    def _handle_s2_baselines(self, rasterio_file):
        """
        Handle Sentinel-2 specific baselines, including fixing paths or adding offsets.
        """
        # Offset parameters for N0400 baseline
        if self.baseline == 'N0400':
            logger.info('Adding offset parameters for N0400')
            root_offset = ET.parse(rasterio_file).getroot()
            if self.processing_level == 'Level-2A':
                self.globalAttribs['BOA_ADD_OFFSET'] = root_offset.find('.//BOA_ADD_OFFSET').text
            elif self.processing_level == 'Level-1C':
                self.globalAttribs['RADIO_ADD_OFFSET'] = root_offset.find('.//RADIO_ADD_OFFSET').text

        # Fix paths for baseline N0208
        elif self.baseline == '_N0208_':
            logger.info('Fixing paths for baseline N0208')
            self._fix_baseline_n0208_paths()

        # Fix paths for baseline N0207
        elif self.baseline == 'N0207':
            logger.info('Fixing paths for baseline N0207')
            self._fix_baseline_n0207_paths()

    def _fix_baseline_n0208_paths(self):
        """
        Fix paths for baseline N0208, which has typos in paths.
        """
        for i, f in self.xmlFiles.items():
            if 'wp_in_progress' in str(f):
                tmp1 = str(f).split('.SAFE')
                try:
                    tmp2 = str(f).split('GRANULE')
                    self.xmlFiles[i] = Path(tmp1[0] + '.SAFE/GRANULE' + tmp2[1])
                except IndexError:
                    tmp2 = str(f).split('DATASTRIP')
                    self.xmlFiles[i] = Path(tmp1[0] + '.SAFE/DATASTRIP' + tmp2[1])
        for i, f in self.imageFiles.items():
            if 'wp_in_progress' in str(f):
                tmp1 = str(f).split('.SAFE')
                if 'GRANULE' in str(f):
                    tmp2 = str(f).split('GRANULE')
                    self.imageFiles[i] = Path(tmp1[0] + '.SAFE/GRANULE' + tmp2[1])

    def _fix_baseline_n0207_paths(self):
        """
        Fix paths for baseline N0207, which has typos in paths.
        """
        for i, f in self.xmlFiles.items():
            self.xmlFiles[i] = Path(
                str(f).replace('/ANULE/', '/GRANULE/').replace('/TASTRIP/', '/DATASTRIP/')
            )


    def setup_rasterio_source(self):
        """
        Set up the rasterio source file based on the satellite type and baseline.
        """
        # Determine the path to the rasterio source file
        if not self.dterrengdata:
            rasterio_file = self.SAFE_dir / str(self.xmlFiles[f'S2_{self.processing_level}_Product_Metadata'])
        else:
            rasterio_file = str(self.mainXML)

        # Debugging
        logger.debug(f"XML file to be opened: {rasterio_file}")

        # Open the file with rasterio
        try:
            self.src = rasterio.open(rasterio_file)
        except Exception as e:
            raise RuntimeError(f"Failed to open file with rasterio: {e}")

        logger.debug(self.src)

        # Set global metadata attributes from rasterio
        self.globalAttribs = self.src.tags()

        # Handle satellite-specific logic
        self._handle_s2_baselines(rasterio_file)


    def get_metadata(self):
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

            relative_path = str(xmlFile).lstrip("/")
            absolute_path = self.extracted_folder_path / relative_path

            root = xml_read(absolute_path)
            polarisation = root.find('.//polarisation').text

            for pm in productMetadata_parameters:
                variable = root.find(str('.//' + pm))
                self.productMetadata[polarisation][pm] = variable.text

            for pml in productMetadataList_parameters:
                variable = root.find(str('.//' + pml))
                self._extractProductMetadataList(variable, polarisation)

    def _extractProductMetadataList(self, mother_element, polarisation):
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
            logger.error("Extraction of %s is not implemented" % listType)

    # TODO: Function to extract thumbnail?
    # TODO: Function to create quicklooks?