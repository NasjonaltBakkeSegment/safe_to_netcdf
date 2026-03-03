from osgeo import gdal
import zipfile
import math
from collections import defaultdict
import logging
import pathlib
import numpy as np
import pyproj
import shapely.wkt

try:
    from lib.utils import xml_read
    import config.constants as cst 
except:
    from utils import xml_read
    import config.constants as cst

try:
    from lib.safe.safe_base import SAFEFile
except ModuleNotFoundError:
    from safe_base import SAFEFile

logger = logging.getLogger(__name__)

class S2SAFEFile(SAFEFile):
    """
    Subclass for working with Sentinel-1 SAFE files
    """

    def __init__(self, product, zipdir, tmpdir):

        super().__init__(product, zipdir, tmpdir)

        self.SAFE_dir = (tmpdir / self.product_name).with_suffix('.SAFE')
        self.baseline = self.product_name.split('_')[3]
        self.processing_level = 'Level-' + self.product_name.split('_')[1][4:6]
        self.imageFiles = defaultdict(list)
        self.reference_band = None
        self.dterrengdata = False  # variable saying if products is Norwegian DEM L1C
        self.sunAndViewAngles = defaultdict(list)
        self.vectorInformation = defaultdict(list)
        self.image_list_dterreng = []
        self.mainXML = self.SAFE_dir / 'manifest.safe'

    def prepare_for_use(self):
        """
        Prepare the SAFEFile instance for use.
        This includes unzipping the product, reading metadata, creating file lists, and setting up the rasterio source.
        """
        super().prepare_for_use()
        if not self.mainXML.is_file():
            self.mainXML = self.SAFE_dir / 'MTD_MSIL1C.xml'
            if not self.mainXML.is_file():
                logger.error(f'Main file not found. Exiting')
                raise

        self.img_data_dirs = list(self.SAFE_dir.rglob("IMG_DATA"))
        self.create_file_lists()
        self.set_gdal_object()
        self.readSunAndViewAngles()
        self.get_dimensions()

    def create_file_lists(self):
        """
        Create lists of XML/GML and image files included in the SAFE product.
        """
        if self.dterrengdata:
            # For DTERR data, add manually the list of images / xml-gml files from parsing the SAFE
            # directory
            allFiles = zipfile.ZipFile(self.input_zip).namelist()
            for f in allFiles:
                fWithPath = self.SAFE_dir.parent / f
                if fWithPath.suffix == '.xml' or fWithPath.suffix == '.gml':
                    self.xmlFiles[fWithPath.stem] = fWithPath
            # Read relative image path (since gdal can't open all these products..)
            self.image_list_dterreng = [self.SAFE_dir.parent / s for s in allFiles
                                                        if ".jp2" in s and "IMG_DATA" in s]
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
                    if (ftype == 'text/xml' or ftype == 'application/xml') and href:
                        self.xmlFiles[repID] = self.SAFE_dir / href[1:]
                    elif ftype == 'application/octet-stream':
                        self.imageFiles[repID] = self.SAFE_dir / href[1:]

    def set_gdal_object(self):
        if not self.dterrengdata:
            gdalFile = self.SAFE_dir / str(self.xmlFiles['S2_{}_Product_Metadata'.format(self.processing_level)]).lstrip("/")
        else:
            gdalFile = str(self.mainXML)
        self.src = gdal.Open(gdalFile)
        if self.src is None:
            raise
        logger.debug((self.src))

        self.globalAttribs = self.src.GetMetadata()

        # Offset parameters for N0400 baseline, for both L1C and L2A products
        if self.baseline == 'N0400':
            logger.info('Adding offset parameters for N0400')
            root_offset = xml_read(gdalFile)
            if self.processing_level == 'Level-2A':
                self.globalAttribs['BOA_ADD_OFFSET'] = root_offset.find('.//BOA_ADD_OFFSET').text
            elif self.processing_level == 'Level-1C':
                self.globalAttribs['RADIO_ADD_OFFSET'] = root_offset.find('.//RADIO_ADD_OFFSET').text
        # Baseline N0208 for S2L2A products has typos in paths
        elif self.baseline == '_N0208_':
            logger.info('Fixing paths for baseline N0208')
            for i, f in self.xmlFiles.items():
                if 'wp_in_progress' in str(f):
                    tmp1 = str(f).split('.SAFE')
                    try:
                        tmp2 = str(f).split('GRANULE')
                        self.xmlFiles[i] = pathlib.Path(tmp1[0] + '.SAFE/GRANULE' + tmp2[1])
                    except IndexError:
                        tmp2 = str(f).split('DATASTRIP')
                        self.xmlFiles[i] = pathlib.Path(tmp1[0] + '.SAFE/DATASTRIP' + tmp2[1])
            for i, f in self.imageFiles.items():
                if 'wp_in_progress' in str(f):
                    tmp1 = str(f).split('.SAFE')
                    if 'GRANULE' in str(f):
                        tmp2 = str(f).split('GRANULE')
                        self.imageFiles[i] = pathlib.Path(tmp1[0] + '.SAFE/GRANULE' + tmp2[1])
        # Baseline N0207 for S2L2A products has typos in paths
        elif self.baseline == 'N0207':
            logger.info('Fixing paths for baseline N0207')
            for i,f in self.xmlFiles.items():
                self.xmlFiles[i] = pathlib.Path(str.replace(str(f), '/ANULE/', '/GRANULE/').replace('/TASTRIP/', '/DATASTRIP/'))

    def readSunAndViewAngles(self):
        """ Method for reading sun and view angles from Sentinel-2
            annotation files.
        """
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
            return

        relative_path = str(currXml).lstrip("/")
        absolute_path = self.SAFE_dir / relative_path
        root = xml_read(absolute_path)

        angles_view_list = root.findall('.//Tile_Angles')[0]
        angle_step = float(root.findall('.//COL_STEP')[0].text)  # m
        nx = int(root.xpath(str(
            '//n1:{}_Tile_ID/n1:Geometric_Info/Tile_Geocoding/Size[@resolution=10]/NROWS'.format(
                self.processing_level)), namespaces=root.nsmap)[0].text)  # nb of rows
        spatial_resolution = 10

        angle_len = int(math.ceil(nx * spatial_resolution / angle_step))
        sun_zenith = np.zeros((angle_len, angle_len), dtype=np.float32)
        sun_azimuth = np.zeros((angle_len, angle_len), dtype=np.float32)
        angle_step = int(math.ceil(nx / float(angle_len)))
        incidence_angles_list = angles_view_list.findall('Viewing_Incidence_Angles_Grids')

        # Sun angles
        for angle in angles_view_list.find('Sun_Angles_Grid'):
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

    def get_dimensions(self):
        # Deciding a reference band
        #todo dterreng warning coming from here?
        # yes -> self.src.GetSubDatasets() ok but the gdal.Open does not work
        # add break? how to remove warning from dterr?
        for k, v in self.src.GetSubDatasets():
            if v.find('10m') > 0:
                self.reference_band = gdal.Open(k)

        # frequency bands
        self.xSize = self.reference_band.RasterXSize  # number of pixels for 10m spatial resolution
        self.ySize = self.reference_band.RasterYSize  # number of pixels for 10m spatial resolution

        # sun and view angles raster resolution
        self.xaSize, self.yaSize = self.sunAndViewAngles[list(self.sunAndViewAngles)[0]].shape

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

    def get_global_attributes(self):

        super().get_global_attributes()

        polygon = shapely.wkt.loads(self.globalAttribs.pop('FOOTPRINT'))
        self.globalAttribs.update({
            'relative_orbit_number': self.globalAttribs.pop("DATATAKE_1_SENSING_ORBIT_NUMBER"),
            'orbit_direction': self.globalAttribs.pop("DATATAKE_1_SENSING_ORBIT_DIRECTION").lower(),
            'cloud_coverage': self.globalAttribs.pop("CLOUD_COVERAGE_ASSESSMENT"),
            'time_coverage_start': self.globalAttribs.pop('PRODUCT_START_TIME'),
            'time_coverage_end': self.globalAttribs.pop('PRODUCT_STOP_TIME'),
        })
        if self.dterrengdata:
            self.globalAttribs['orbit_number'] = int(self.globalAttribs['DATATAKE_1_ID'].split('_')[2])
        else:
            self.globalAttribs['orbit_number'] = int(self.root.find('.//safe:orbitNumber', namespaces=self.root.nsmap).text)

        super()._compute_bounding_box(polygon)