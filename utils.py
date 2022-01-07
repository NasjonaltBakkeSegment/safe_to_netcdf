"""
Tools
"""

import pathlib
import lxml.etree as ET
import datetime as dt
import pytz
import resource
from osgeo import gdal
import subprocess as sp
import zipfile
import logging

logger = logging.getLogger(__name__)


def xml_read(xml_file):
    """
    Read XML file.
    Args:
        xml_file [pathlib]): filepath to an xml file
    Returns:
        lxml.etree._Element or None if file missing
    """
    if not pathlib.Path(xml_file).is_file():
        logger.error(f'Error: Can\'t find xmlfile {xml_file}')
        return None
    tree = ET.parse(str(xml_file))
    root = tree.getroot()
    return root


def memory_use(start_time):
    """
    Print memory usage and time taken by a process.
    Args:
        start_time: datetime object containing start time of process
    Returns:
        N/A
    """
    logger.debug(f"Memory usage so far: "
          f"{float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss) / 1000000} Gb")
    logger.debug(dt.datetime.now(tz=pytz.utc) - start_time)


def seconds_from_ref(t, t_ref):
    """
    Computes the difference in seconds between input date and a reference date.
    Args:
        t: date as a string
        t_ref: reference time as a datetime
    Returns:
        integer
    """
    try:
        mytime = dt.datetime.strptime(t, '%Y-%m-%dT%H:%M:%S.%f')
    except ValueError:
        mytime = dt.datetime.strptime(t, '%Y-%m-%dT%H:%M:%S.%fZ')
    return int((mytime - t_ref).total_seconds())


def create_time(ncfile, t, ref='01/01/1981'):
    """
    Create time variable for netCDF file.
    Args:
        ncfile: netcdf file (already open)
        t: time as string
        ref: reference time as string (dd/mm/yyyy)
    Returns: True
    """

    ref_dt = dt.datetime.strptime(ref, '%d/%m/%Y')
    nc_time = ncfile.createVariable('time', 'i4', ('time',))
    nc_time.long_name = 'reference time of satellite image'
    nc_time.units = f"seconds since {ref_dt.strftime('%Y-%m-%d %H:%M:%S')}"
    nc_time.calendar = 'gregorian'
    nc_time[:] = seconds_from_ref(t, ref_dt)
    return True


def initializer(self):
    """
    Use the main XML file to initialize extra variables:
     - list of XML-GML files
     - list of images
     - global attributes from gdal info
     - ...
    Returns: True
    """
    root = xml_read(self.mainXML)
    sat = self.product_id.split('_')[0][0:2]

    # List of xml / gml files
    if sat == 'S2' and self.dterrengdata:
        # Added to be strictly identical to older SAFE2NC version - not used otherwise
        self.xmlFiles['mainXML'] = self.SAFE_dir / 'MTD_MSIL1C.xml'
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
        dataObjectSection = root.find('./dataObjectSection')
        for dataObject in dataObjectSection.findall('./'):
            if sat == 'S1':
                repID = dataObject.attrib['repID']
            elif sat == 'S2':
                repID = dataObject.attrib['ID']
            ftype = None
            href = None
            for element in dataObject.iter():
                attrib = element.attrib
                if 'mimeType' in attrib:
                    ftype = attrib['mimeType']
                if 'href' in attrib:
                    href = attrib['href'][1:]
            if sat == 'S2':
                if (ftype == 'text/xml' or ftype == 'application/xml') and href:
                    self.xmlFiles[repID] = self.SAFE_dir / href[1:]
                elif ftype == 'application/octet-stream':
                    self.imageFiles[repID] = self.SAFE_dir / href[1:]
            elif sat == 'S1':
                if ftype == 'text/xml' and href:
                    self.xmlFiles[repID].append(self.SAFE_dir / href[1:])

    # Set gdal object
    if sat == 'S2' and not self.dterrengdata:
        gdalFile = str(self.xmlFiles['S2_{}_Product_Metadata'.format(self.processing_level)])
    else:
        gdalFile = str(self.mainXML)
    self.src = gdal.Open(gdalFile)
    if self.src is None:
        raise
    logger.debug((self.src))

    # Set global metadata attributes from gdal
    self.globalAttribs = self.src.GetMetadata()

    if sat == 'S1':
        # Set raster size parameters
        self.xSize = self.src.RasterXSize
        self.ySize = self.src.RasterYSize
        # Set polarisation parameters
        polarisations = root.findall('.//s1sarl1:transmitterReceiverPolarisation',
                                     namespaces=root.nsmap)
        ##for polarisation in polarisations:
        ##    self.polarisation.append(polarisation.text)
        ##self.globalAttribs['polarisation'] = self.polarisation

        # to be identical to python2/Xenial
        outattrib = ''
        for polarisation in polarisations:
            self.polarisation.append(polarisation.text)
            outattrib += polarisation.text
        self.globalAttribs['polarisation'] = [outattrib]

        # Timeliness
        self.globalAttribs['ProductTimelinessCategory'] = root.find(
            './/s1sarl1:productTimelinessCategory', namespaces=root.nsmap).text

    return True


def uncompress(self):
    """
    Uncompress a SAFE zip file.
    Find the main XML: manifest or other (dterreng data).
    Return: True
    """

    # If zip not extracted yet
    if not self.SAFE_dir.is_dir():
        logger.debug('Starting unzipping SAFE archive')
        self.SAFE_dir.parent.mkdir(parents=False, exist_ok=True)
        sp.run(["/usr/bin/unzip", "-qq", self.input_zip, "-d", self.SAFE_dir.parent], check=True)
        logger.debug('Done unzipping SAFE archive')

    # Try and find the main XML file
    if self.SAFE_dir.stem.startswith('S3'):
        xmlFile = self.SAFE_dir / 'xfdumanifest.xml'
    else:
        xmlFile = self.SAFE_dir / 'manifest.safe'
    if not xmlFile.is_file():
        xmlFile = self.SAFE_dir / 'MTD_MSIL1C.xml'
        if not xmlFile.is_file():
            logger.error(f'Main file not found. Exiting')
            raise
        self.dterrengdata = True

    logger.debug(f'Main file: {xmlFile}')
    self.mainXML = xmlFile
    return True


def orbit_info(manifest):
    """

    :param manifest:
    :return:
    """
    #todo: check with S1, s2, S2DTERR xml files
    try:
        orbit = manifest.find('.//sentinel-safe:orbitNumber', namespaces=manifest.nsmap)
        orbit_absolute = orbit.text
        orbit_direction = orbit.get('groundTrackDirection')
        orbit_relative = manifest.find('.//sentinel-safe:relativeOrbitNumber',
                                            namespaces=manifest.nsmap).text
    except SyntaxError:
        orbit_absolute = manifest.find('.//safe:orbitNumber', namespaces=manifest.nsmap).text

    return orbit_absolute, orbit_relative, orbit_direction


# Add function to clean work files?
