"""
Tools -
"""

# from pathlib import Path
# import lxml.etree as ET
# import datetime as dt
# import resource
# import subprocess as sp
# import zipfile
# import logging
# import yaml
# from pkg_resources import resource_string
# import shapely.wkt, shapely.ops
# import rasterio
# import argparse
# import io
# from io import BytesIO
# from netCDF4 import Dataset
# from rasterio.io import MemoryFile
# import tempfile
# import rioxarray as rio
# import xarray as xr
# import re
# from s1_reader_and_NetCDF_converter import Sentinel1_reader_and_NetCDF_converter
# from s2_reader_and_NetCDF_converter import Sentinel2_reader_and_NetCDF_converter
# from s3_olci_l1_reader_and_CF_converter import Sentinel3_olci_reader_and_CF_converter
# import os
# import shutil

from pathlib import Path
import lxml.etree as ET
import datetime as dt
import resource
import subprocess as sp
import zipfile
import logging
import yaml
from pkg_resources import resource_string
import shapely.wkt, shapely.ops
import rasterio
from rasterio.io import MemoryFile
import tempfile
import xarray as xr
import re



logger = logging.getLogger(__name__)


def xml_read(xml_file):
    """
    Read XML file.
    Args:
        xml_file [pathlib]): filepath to an xml file
    Returns:
        lxml.etree._Element or None if file missing
    """
    if not Path(xml_file).is_file():
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
    logger.debug(dt.datetime.now(dt.timezone.utc) - start_time)


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
        # For DTERR data, add manually the list of images / xml-gml files from parsing the SAFE
        # directory
        allFiles = zipfile.ZipFile(self.input_zip).namelist()
        for f in allFiles:
            fWithPath = self.SAFE_dir.parent / f
            if fWithPath.suffix == '.xml' or fWithPath.suffix == '.gml':
                self.xmlFiles[fWithPath.stem] = fWithPath
        # Read relative image path 
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


    # Set rasterio source path
    if sat == 'S2' and not self.dterrengdata:
        rasterio_file = str(self.xmlFiles[f'S2_{self.processing_level}_Product_Metadata'])
    else:
        rasterio_file = str(self.mainXML)
    # Debugging
    logger.debug(f"XML file to be opened: {rasterio_file}")
    try:
        self.src = rasterio.open(rasterio_file)
    except Exception as e:
        raise RuntimeError(f"Failed to open file with rasterio: {e}")

    logger.debug(self.src)

    # Set global metadata attributes from rasterio
    self.globalAttribs = self.src.tags()

    if sat == 'S2':
        # Offset parameters for N0400 baseline, for both L1C and L2A products
        if self.baseline == 'N0400':
            logger.info('Adding offset parameters for N0400')
            root_offset = xml_read(rasterio_file)
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
        # Baseline N0207 for S2L2A products has typos in paths
        elif self.baseline == 'N0207':
            logger.info('Fixing paths for baseline N0207')
            for i,f in self.xmlFiles.items():
                self.xmlFiles[i] = Path(str.replace(str(f), '/ANULE/', '/GRANULE/').replace('/TASTRIP/', '/DATASTRIP/'))

    elif sat == 'S1':
        # Set raster size parameters
        self.xSize = self.src.width
        self.ySize = self.src.height
        # Set polarisation parameters
        polarisations = root.findall('.//s1sarl1:transmitterReceiverPolarisation',
                                     namespaces=root.nsmap)
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
        sp.run(["/usr/bin/unzip", "-qq", str(self.input_zip), "-d", str(self.SAFE_dir.parent)], check=True)
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


def read_yaml(file):
    """
    Read yaml file
    """
    return yaml.load(
        resource_string(
            globals()['__name__'].split('.')[0], file
        ), Loader=yaml.FullLoader
    )


def get_global_attributes(self):
    """
    Add global attributes to netcdf file
    """

    all = read_yaml('global_attributes.yaml')

    satellite = self.product_id[0:3]

    # NBS project metadata
    self.globalAttribs.update(all['global'])

    # Satellite specific metadata (S1, S2, S3)
    self.globalAttribs.update(all[satellite[0:2]])

    # Platform specific metadata (S1A, S1B, S2A, S2B, ...)
    self.globalAttribs.update(all[satellite])

    # Product specific metadata
    self.globalAttribs.update({
        'date_metadata_modified': self.t0.isoformat().replace("+00:00", "Z"),
        'date_metadata_modified_type': 'Created',
        'date_created': self.t0.isoformat().replace("+00:00", "Z"),
        'history': f'{self.t0.isoformat()}: Converted from SAFE to NetCDF by NBS team.',
        'title': self.product_id,
        'id': self.uuid,
        'product_type': self.globalAttribs.pop("PRODUCT_TYPE"),
    })

    # If no uuid provided, no id attribute at all
    if self.uuid is None:
        logger.warning('uuid not found, so no id attribute in the nc global attributes')
        self.globalAttribs.pop('id')

    root = xml_read(self.mainXML)

    if satellite.startswith('S2'):
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
            self.globalAttribs['orbit_number'] = int(root.find('.//safe:orbitNumber', namespaces=root.nsmap).text)

    elif satellite.startswith('S1'):
        self.globalAttribs.update({
            'orbit_number': int(self.globalAttribs.pop("ORBIT_NUMBER")),
            'orbit_direction': self.globalAttribs.pop("ORBIT_DIRECTION").lower(),
            'relative_orbit_number': root.find('.//safe:relativeOrbitNumber', namespaces=root.nsmap).text,
            'time_coverage_start': self.globalAttribs.pop('ACQUISITION_START_TIME')+'Z',
            'time_coverage_end': self.globalAttribs.pop('ACQUISITION_STOP_TIME')+'Z',
            'mode': self.globalAttribs.pop('MODE')
        })
        # Some work is necessary to get the polygon:
        # - get coordinates from manifest
        coords = root.find('.//gml:coordinates', namespaces=root.nsmap).text
        # - format to be shapely compliant
        # - close the polygon by adding the first point at the end
        formatted = coords.replace(' ', ';').replace(',', ' ').replace(';', ', ')
        footprint = f"POLYGON(({','.join([formatted, formatted.split(',')[0]])}))"
        # - reverse lat-lon
        polygon = shapely.ops.transform(lambda x, y: (y, x), shapely.wkt.loads(footprint))

    # Add bounding box and polygon
    box = polygon.bounds
    self.globalAttribs.update({
        'geospatial_lon_min': box[0],
        'geospatial_lat_min': box[1],
        'geospatial_lon_max': box[2],
        'geospatial_lat_max': box[3],
    })

    # Add SIOS in collection if latitude > 70
    if box[3] > 70:
        self.globalAttribs['collection'] += ',SIOS'

    return


def get_key(my_dict,val):
    for key, value in my_dict.items():
         if val == value:
             return key

    return "There is no such Key"

def write_to_geotiff(conversion_object, indir, product, outdir):

        # Creating a memory stream for creating temporary NetCDF-file to convert
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)

            nc_out = tmpdir
            nc_out.parent.mkdir(parents=True, exist_ok=True)

            
            # Write .SAFE data to temporary NetCDF
            if conversion_object.read_ok:
                conversion_object.write_to_NetCDF(nc_out, outdir, product, 7)
            
                nc_path = nc_out / (product + '.nc')

                ds = xr.open_dataset(nc_path)

                logger.info('Writing GeoTIFFs from nc-file.')

                # For each band in .SAFE data we create a GeoTIFF-file
                for variable, variable_item in ds.data_vars.items():
                    if variable_item.dims == ('time', 'y', 'x') or variable_item.dims == ('rows', 'columns') :
                        band = ds[variable]

                        #Define lat/long 
                        try:
                            band = band.rio.set_spatial_dims('x', 'y')
                        except:
                            band = band.rio.set_spatial_dims('rows', 'columns')

                        #If  CRS is not discovered, add manually
                        if not band.rio.crs:
                            band = band.rio.write_crs("EPSG:4326") # add epsg

                        with MemoryFile() as memfile:
                            with memfile.open(
                                driver='GTiff',
                                compress='LZW', # Change compression upon need
                                height=band.shape[-2],
                                width=band.shape[-1],
                                count=1,
                                dtype=band.dtype,
                                transform=band.rio.transform()
                                ) as dst:
                                try:
                                    dst.write(band.values)
                                except:
                                    dst.write(band.values.reshape((1, *band.shape))) # rewrap dimensions of band if not 3D



                            
                            outfile = Path(outdir / (variable + ".tif"))

                            outfile.parent.mkdir(parents=True, exist_ok=True)

                            # Save GeoTIFF to disk
                            with open(outfile, "wb") as f:
                                f.write(memfile.read()) 
            else:
                AssertionError('Could not instantiate class')  


def extract_band_code(filename):
    """
    Extracts band code like 'B08', 'B8A', 'B11', etc. from a Sentinel-2 filename.
    Returns None if not found.
    """
    match = re.search(r'_B(0?[1-9A]{2})', filename)  # matches _B08, _B8A, _B11, etc.
    if match:
        return 'B' + match.group(1)
    return None


