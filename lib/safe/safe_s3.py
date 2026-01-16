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

class S3SAFEFile(SAFEFile):
    """
    Subclass for working with Sentinel-3 SEN3 files
    """

    def __init__(self, product, indir, outdir):

        super().__init__(product, indir, outdir)

        # TODO: Add S3 specific attributes

    def prepare_for_use(self):
        """
        Prepare the SAFEFile instance for use.
        This includes unzipping the product, reading metadata, creating file lists, and setting up the rasterio source.
        """
        super().prepare_for_use()

        # TODO: None of this has been written for S3
        # Read the metadata XML
        self.read_metadata_xml()

        # Create the file lists
        self.create_file_lists()

        # Set up the rasterio source and metadata
        self.setup_rasterio_source()

        # Extract product metadata
        self.get_metadata()

    def _get_mission(self):
        if self.product_name.startswith('S1'):
            return 'S1'
        elif self.product_name.startswith('S2'):
            return 'S2'
        elif self.product_name.startswith('S3'):
            return 'S3'
        else:
            return None

    def unzip(self):
        '''
        Uncompress the zip file. Write it out locally
        '''
        # TODO Should be uncompressed somewhere else
        # If zip not extracted yet
        if not self.SAFE_dir.is_dir():
            logger.debug('Starting unzipping SAFE archive')
            self.SAFE_dir.parent.mkdir(parents=False, exist_ok=True)

            sp.run(["/usr/bin/unzip", "-qq", str(self.input_zip), "-d", str(self.SAFE_dir.parent)], check=True)
            logger.debug('Done unzipping SAFE archive')

        extracted_folder_name = self.input_zip.stem + '.SAFE'
        self.extracted_folder_path = self.SAFE_dir.parent / extracted_folder_name

    def _get_xml_path(self):

        xmlFile = self.SAFE_dir / 'xfdumanifest.xml'

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
        if self.mission == 'S2' and self.dterrengdata:
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
                repID = dataObject.attrib.get('repID' if self.mission == 'S1' else 'ID')
                ftype, href = None, None
                for element in dataObject.iter():
                    ftype = element.attrib.get('mimeType', ftype)
                    href = element.attrib.get('href', href)
                if href:
                    href_path = self.SAFE_dir / href[1:]  # Convert href to full path
                    if ftype in {'text/xml', 'application/xml'}:
                        self.xmlFiles[repID].append(href_path) if self.mission == 'S1' else self.xmlFiles.update({repID: href_path})
                    elif ftype == 'application/octet-stream' and self.mission == 'S2':
                        self.imageFiles[repID] = href_path

    def setup_rasterio_source(self):
        """
        Set up the rasterio source file based on the satellite type and baseline.
        """
        # Determine the path to the rasterio source file
        if self.mission == 'S2' and not self.dterrengdata:
            rasterio_file = str(self.xmlFiles[f'S2_{self.processing_level}_Product_Metadata'])
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
        if self.mission == 'S2':
            self._handle_s2_baselines(rasterio_file)
        elif self.mission == 'S1':
            self._handle_s1_metadata()

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

    def _handle_s1_metadata(self):
        """
        Handle Sentinel-1 specific metadata, such as raster size and polarisation.
        """
        # Set raster size parameters
        self.xSize = self.src.width
        self.ySize = self.src.height

        # Extract polarisation parameters
        root = ET.parse(self.mainXML).getroot()
        polarisations = root.findall('.//s1sarl1:transmitterReceiverPolarisation', namespaces=root.nsmap)
        outattrib = ''
        for polarisation in polarisations:
            self.polarisation.append(polarisation.text)
            outattrib += polarisation.text
        self.globalAttribs['polarisation'] = [outattrib]

        # Extract Product Timeliness
        self.globalAttribs['ProductTimelinessCategory'] = root.find(
            './/s1sarl1:productTimelinessCategory', namespaces=root.nsmap
        ).text

    def initialize_s1_metadata(self):
        '''
        Initialize Sentinel-1 metadata, including calibration tables, noise vectors,
        and GCP (Ground Control Point) parameters.
        '''
        # Calibration tables
        calibrationTables = ['sigmaNought', 'betaNought', 'gamma', 'dn']
        for calXmlFile in self.xmlFiles['s1Level1CalibrationSchema']:
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

                if parameter != 'azimuthTime':
                    self.xmlGCPs[str(parameter + '_' + polarisation)] = np.array(values, np.float32)
                else:
                    self.xmlGCPs[str(parameter + '_' + polarisation)] = np.array(values, str)

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

# Configure logging
logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)
import time

def test_SAFEFile():
    """
    Test the SAFEFile class functionality using the prepare_for_use and finalize_usage functions.
    """
    # Define input and output directories as the current directory
    script_dir = Path(__file__).parent
    indir = outdir = script_dir

    # Provide a product name for testing (update this to match a real product for your tests)
    product_name = "S1C_EW_GRDM_1SDH_20251215T103645_20251215T103745_005461_00AE02_8894"

    # Create an instance of SAFEFile
    safe_file = SAFEFile(product=product_name, indir=indir, outdir=outdir)

    # Prepare the SAFEFile instance for use
    safe_file.prepare_for_use()

    # Perform assertions after preparation
    assert safe_file.SAFE_dir.is_dir(), "Unzipped directory does not exist as expected."
    assert safe_file.mainXML.is_file(), "Main XML file was not found as expected."
    assert safe_file.globalAttribs is not None, "Global metadata attributes were not set."

    # Print extracted metadata for verification
    print("Global metadata attributes:", safe_file.globalAttribs)

    # Print discovered files for verification
    if safe_file.mission == 'S1':
        print("Discovered XML files for S1 mission:")
        for key, value in safe_file.xmlFiles.items():
            print(f"  {key}: {value}")
    elif safe_file.mission == 'S2':
        print("Discovered XML files for S2 mission:")
        for key, value in safe_file.xmlFiles.items():
            print(f"  {key}: {value}")
        print("Discovered image files for S2 mission:")
        for key, value in safe_file.imageFiles.items():
            print(f"  {key}: {value}")

    for pol, metadata in safe_file.productMetadata.items():
        print(f"Polarisation: {pol}")
        for key, value in metadata.items():
            print(f"  {key}: {value}")

    # Finalize usage (clean up)
    safe_file.finalize_usage()

    # Verify cleanup
    assert not safe_file.SAFE_dir.is_dir(), "Uncompressed folder was not deleted as expected."

if __name__ == "__main__":
    test_SAFEFile()

    # TODO: Consider how to do this for S2 and S3.
    # If code and functions need to be different, use subclasses?