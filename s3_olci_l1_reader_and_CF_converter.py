# -------------------------------------------------------------------------------
# Name:          s3_olci_l1_reader_and_CF_converter.py
# Purpose:       Read Sentinel-3 OLCI l1 product in SAFE and merge into a single
#                NetCDF file following the CF convention.
#
# Author(s):     Trygve Halsne
# Created:       10.05.2019
# Copyright:     (c) Norwegian Meteorological Institute, 2019
# -------------------------------------------------------------------------------
import datetime as dt
import pytz
import xarray as xr
import logging
import sys
import pathlib
import utils
import constants as cst

logger = logging.getLogger(__name__)

"""
    FIXES to be solved
    - tie points having latitude longitude in tie_geo_coordinates. Should be separat group?
    - removed pixels variables uses same names as radiance raster bands. Hence should be fixed.
    -> for now, only had radiance variables to netcdf
    
    - should follow swath convention (https://github.com/Unidata/EC-netCDF-CF/blob/master/swath/swath.adoc#radiometric-swath-data-encodings) which is cf v.1.7 compatible
    - should be ACDD compliant

    ISSUES:
    - int64: not supported in dap2. will come in dap4 https://www.unidata.ucar.edu/support/help/MailArchives/thredds/msg01762.html
    -> todo check if it's still the case
"""


class S3_olci_reader_and_CF_converter:
    """ S3 OLCI object for merging SAFE files and creating single NetCDF/CF """

    def __init__(self, product, indir, outdir, uuid=None):
        """ Initializer

        Args:
            SAFE_file      (str): Input SAFE file
            SAFE_id        (str): SAFE filename
            SAFE_outpath   (str): Output directory for unzipeed SAFE file
            SAFE_path      (str): Absolute path to unzipped SAFE directory
            src
            t0           (float): Init time

        """
        self.product_id = product
        # Optionnaly get uuid from colhub
        self.uuid = uuid
        self.input_zip = (indir / product).with_suffix('.zip')
        self.SAFE_dir = (outdir / self.product_id).with_suffix('.SEN3')
        self.processing_level = 'Level-' + self.product_id.split('_')[2]
        self.t0 = dt.datetime.now(tz=pytz.utc)
        self.read_ok = True
        self.main()

    def main(self):
        """ Main method for traversing and reading key values from SAFE
            directory.
        """

        # unzip SAFE archive
        utils.uncompress(self)

        return True

    def write_to_NetCDF(self, nc_outpath, compression_level, chunk_size=(30, 34)):

        """ Method writing output NetCDF product.
        Args:
            nc_outpath (str): output path where NetCDF file should be stored
            compression_level (int): compression level on output NetCDF file (1-9)
            chunk_size (tuple): chunk size in output NetCDF. Format (rows, columns)
        """

        # todo use ocmpression or not, use chunks or not ???

        logger.info("------------START CONVERSION FROM SAFE TO NETCDF-------------")

        # Status
        utils.memory_use(self.t0)

        # find all netcdf files contained in .SAFE product
        root = utils.xml_read(self.mainXML)
        nc_files = []
        for o in root.findall('.//dataObject'):
            nc_files.append(self.SAFE_dir / o.find('.//fileLocation').attrib['href'].split('/')[1])

        # output filename
        ncout = (nc_outpath / self.product_id).with_suffix('.nc')

        t1 = dt.datetime.now(tz=pytz.utc)
        # Concatenate radiance for all bands and lat/lon
        bands = [s for s in nc_files if "_radiance.nc" in str(s)]
        bands.append([s for s in nc_files if "geo_coordinates.nc" in str(s)][0])
        ds = xr.open_mfdataset(bands, join='exact')
        # todo: test different options
        #ds = xr.open_mfdataset(bands, engine='h5netcdf')
        #ds = xr.open_mfdataset(bands, parallel=True)

        # Add time dimension
        # Done afterwards as dimensions different to the other variables
        time_file = [s for s in nc_files if "time_coordinates.nc" in str(s)][0]
        time = xr.open_dataset(time_file)
        data_tmp = xr.combine_by_coords([ds, time], combine_attrs='override')
        logger.debug('Time to open and combine nc files')
        utils.memory_use(t1)

        # Format variables
        data = data_tmp.rename({'latitude': 'lat', 'longitude': 'lon', 'time_stamp': 'time',
                                'columns': 'x', 'rows': 't'})
        for v in data.variables:
            if v.endswith('_radiance'):
                data[v].attrs['coordinates'] = 'time altitude lat lon'
                # 'standard_name' in input nc file is not actually a standard name, so update
                data[v].attrs['standard_name'] = 'upwelling_radiance_per_unit_wavelength_in_air'
                # Delete irrelevant attributes as ancillary data not included in file
                del(data[v].attrs['ancillary_variables'])
                # Convert from uint16 to int32 to be CF compliant -> necessary in CF 1.8
                # But in CF 1.9, unsigned short (uint16) are allowed
                # data[v] = data[v].astype('int32', copy=False)

        # Axis are not compulsory as lat/lon/altitude/time are auxiliary coordinate variables
        # and not coordinate variables but can help some software to correctly plot data
        data['altitude'].attrs['axis'] = 'Z'
        data['altitude'].attrs['positive'] = 'up'
        data['lat'].attrs['axis'] = 'Y'
        data['lon'].attrs['axis'] = 'X'
        data['time'].attrs['axis'] = 'T'

        # Set grid mapping - not compulsory in CF 1.9 as lat-lon are provided
        # is it helpful / necessary?
        # see https://cfconventions.org/Data/cf-conventions/cf-conventions-1.9/cf-conventions.html#appendix-grid-mappings
        # If add grid mapping, must be added as attribute to variables:
        # data[v].attribs['grid_mapping'] = 'crsWGS84'
        # todo: Trygve, where are these values coming from?
        # is it a copy of https://cfconventions.org/Data/cf-conventions/cf-conventions-1.9/cf-conventions.html#latitude-and-longitude-on-the-wgs-1984-datum
        # nc_crs = ncout.createVariable('crsWGS84', "i2")
        # nc_crs.grid_mapping_name = "latitude_longitude"
        # nc_crs.semi_major_axis = 6378137.0
        # nc_crs.inverse_flattening = 298.2572235604902

        # Add global attributes
        # todo:
        #  - add link to colhub?
        # todo: should clear all already existing attributes? some are redundant, others irrelevant now
        # Generic global attributes - for all NBS products
        attrs = {}
        attrs.update(cst.global_attributes)

        # Global attributes for all S3 OLCI products
        attrs.update(cst.s3_olci)
        attrs.update(cst.platform[self.product_id[0:3]])
        attrs.update(cst.instrument['OLCI'])

        # Global attributes specific to a dataset
        # todo: should not id be a uuid? and use the same as colhub, like it's done in mmd
        if self.uuid is not None:
            attrs['id'] = self.uuid
        else:
            attrs['id'] = self.product_id
        attrs['title'] = self.product_id
        attrs['date_created'] = self.t0.isoformat()
        attrs['geospatial_lat_min'] = data['lat'].min().values
        attrs['geospatial_lat_max'] = data['lat'].max().values
        attrs['geospatial_lon_min'] = data['lon'].min().values
        attrs['geospatial_lon_max'] = data['lon'].max().values
        attrs['time_coverage_start'] = data.attrs.pop("start_time")
        attrs['time_coverage_end'] = data.attrs.pop("stop_time")
        attrs['orbit_absolute'], attrs['orbit_relative'], attrs['orbit_direction'] \
            = utils.orbit_info(root)
        # todo: add creation date? as a modification or creation date?
        attrs['history'] = f'{self.t0.isoformat()}: Converted from SAFE to NetCDF/CF by the NBS' \
                           f' team. ' + data.attrs['history']
        data.attrs = attrs

        # Write netcdf file
        t2 = dt.datetime.now(tz=pytz.utc)
        # todo: test with and without load
        # todo: test compression levels. as is, using the encodings from input nc.
        # if we want to change, need to update encoding for each variable
        # todo: use chunks or not? very unconvinced myself
        data.load().to_netcdf(ncout)
        logger.debug('Time to write nc file')
        utils.memory_use(t2)

        return True


if __name__ == '__main__':

    # Log to console
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    log_info = logging.StreamHandler(sys.stdout)
    log_info.setFormatter(
        logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s'))
    logger.addHandler(log_info)

    workdir = pathlib.Path(
        '/home/trygveh/data/satellite_data/Sentinel-3/OCLI/L1/S3B_OL_1_EFR____20200518T162416_20200518T162716_20200519T201048_0179_039_083_1440_LN1_O_NT_002')
    products = [
        'S3B_OL_1_EFR____20200518T162416_20200518T162716_20200519T201048_0179_039_083_1440_LN1_O_NT_002']

    workdir = pathlib.Path('/home/elodief/Data/NBS/NBS_test_data/test_s3_olci')
    products = [
        'S3B_OL_1_EFR____20211213T082842_20211213T083015_20211214T123309_0092_060_178_1980_LN1_O_NT_002']

    ##workdir = pathlib.Path('/lustre/storeB/project/NBS2/sentinel/production/NorwAREA/netCDFNBS_work/test_environment/test_s3')
    ##products = [
    ##    'S3A_OL_1_EFR____20211124T104042_20211124T104307_20211124T130348_0145_079_051_1980_LN1_O_NR_002',
    ##]

    for product in products:
        outdir = workdir / product
        outdir.parent.mkdir(parents=False, exist_ok=True)
        s3_obj = S3_olci_reader_and_CF_converter(
            product=product,
            indir=workdir / product,
            outdir=outdir)
        s3_obj.write_to_NetCDF(outdir, compression_level=1)
