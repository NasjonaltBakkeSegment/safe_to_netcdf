"""
Tools
"""

import pathlib
import lxml.etree as ET
import datetime as dt
import resource
from osgeo import gdal


def xml_read(xml_file):
    """ Validate xml syntax from filepath.

    Args:
        xml_file ([pathlib object]): [filepath to an xml file]
    Returns:
        [bool]: [return True if a valid xml filepath is provided, 
        raises an exception if the xmlfile is invalid, empty, or doesn't exist ]
    """
    if not pathlib.Path(xml_file).is_file():
        print(f'Error: Can\'t find xmlfile {xml_file}')
        return None

    tree = ET.parse(str(xml_file))
    root = tree.getroot()

    return root


def memory_use(start_time):
    """
    Print memory usage and time taken by a process
    Args:
        start_time: datetime object containing start time of process
    Returns:
        N/A
    """
    print('\nAdding subswath layers')
    print(f"Memory usage so far: "
          f"{float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss) / 1000000} Gb")
    print(dt.datetime.now() - start_time)


def seconds_from_ref(t, t_ref):
    """
    Computes the difference in seconds between input date and a reference date (01/01/1981)
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
    Returns:
    """

    ref_dt = dt.datetime.strptime(ref, '%d/%m/%Y')
    nc_time = ncfile.createVariable('time', 'i4', ('time',))
    nc_time.long_name = 'reference time of satellite image'
    nc_time.units = f"seconds since {ref_dt.strftime('%Y-%m-%d %H:%M:%S')}"
    nc_time.calendar = 'gregorian'
    nc_time[:] = seconds_from_ref(t, ref_dt)

    return True


def initializer(self, xmlFile):
    """
       Traverse manifest file for setting additional variables
            in __init__
        Args:
            xmlFile:

        Returns:
    """
    #todo: add s2 dterreng case
    root = xml_read(xmlFile)
    sat = self.product_id.split('_')[0][0:2]

    # Set xml-files
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

    # Set processing level
    if sat == 'S2':
        self.processing_level = 'Level-' + self.product_id.split('_')[1][4:6]
        gdalFile = str(self.xmlFiles['S2_{}_Product_Metadata'.format(self.processing_level)])
    elif sat == 'S1':
        gdalFile = str(self.xmlFiles['manifest'])

    # Set gdal object
    self.src = gdal.Open(gdalFile)
    print((self.src))

    # Set global metadata attributes from gdal
    self.globalAttribs = self.src.GetMetadata()

    if sat == 'S1':
        # Set raster size parameters
        self.xSize = self.src.RasterXSize
        self.ySize = self.src.RasterYSize
        # Set polarisation parameters
        polarisations = root.findall('.//s1sarl1:transmitterReceiverPolarisation',
                                     namespaces=root.nsmap)
        for polarisation in polarisations:
            self.polarisation.append(polarisation.text)
        self.globalAttribs['polarisation'] = self.polarisation
        # Timeliness
        self.globalAttribs['ProductTimelinessCategory'] = root.find(
            './/s1sarl1:productTimelinessCategory', namespaces=root.nsmap).text

    return True

# Add function to clean work files?
