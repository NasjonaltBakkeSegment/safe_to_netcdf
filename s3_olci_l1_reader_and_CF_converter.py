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
from . import utils
from . import constants as cst

logger = logging.getLogger(__name__)

"""
    FIXES to be solved
    - tie points having latitude longitude in tie_geo_coordinates. Should be separat group?
    - removed pixels variables uses same names as radiance raster bands. Hence should be fixed.
    - should follow swath convention (https://github.com/Unidata/EC-netCDF-CF/blob/master/swath/swath.adoc#radiometric-swath-data-encodings) which is cf v.1.7 compatible
    - should be ACDD compliant

    ISSUES:
    - (FIXED) int64: not supported in dap2. will come in dap4 https://www.unidata.ucar.edu/support/help/MailArchives/thredds/msg01762.html
"""


class S3_olci_reader_and_CF_converter:
    """ S3 OLCI object for merging SAFE files and creating single NetCDF/CF """

    def __init__(self, product, indir, outdir):
        """ Initializer

        Args:
            SAFE_file      (str): Input SAFE file
            SAFE_id        (str): SAFE filename
            SAFE_outpath   (str): Output directory for unzipeed SAFE file
            SAFE_path      (str): Absolute path to unzipped SAFE directory
            xmlFiles      (dict): defaultdict containing xml files
            src
            t0           (float): Init time
            SAFE_structure (str): XML formatted string showing SAFE structure

        """
        self.product_id = product
        file_path = indir / (product + '.zip')
        print(file_path)
        if file_path.exists():
            self.input_zip = file_path
        else:
            file_path_safe = indir / (product + '.SAFE.zip')
            if file_path_safe.exists():
                self.input_zip = file_path_safe
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
        # Do afterwards as dimensions different to the other variables
        time_file = [s for s in nc_files if "time_coordinates.nc" in str(s)][0]
        time = xr.open_dataset(time_file, decode_times=False)
        time_decoded = xr.open_dataset(time_file)

        time_values = time['time_stamp'].values
        units = time['time_stamp'].attrs.get('units', '')

        if len(time_values) > 0:
            base_unit, ref_date = self.parse_units(units)
            ref_time = time_decoded['time_stamp'].values[0]

            if base_unit == 'microseconds':
                time_values_milliseconds = (time_values - time_values[0]) / 1_000
            elif base_unit == 'milliseconds':
                time_values_milliseconds = time_values - time_values[0]
            elif base_unit == 'seconds':
                time_values_milliseconds = (time_values - time_values[0]) * 1_000
            else:
                raise ValueError(f'Unsupported time unit: {base_unit}')
        time['time_stamp'].values = time_values_milliseconds
        time['time_stamp'].attrs['units'] = f'milliseconds since {ref_time}'
        time['time_stamp'].attrs['long_name'] = f'Elapsed time since {ref_time}'


        data_tmp = xr.combine_by_coords([ds, time], combine_attrs='override')
        logger.debug('Time to open and combine nc files')
        utils.memory_use(t1)

        # Format variables
        data = data_tmp.rename({'latitude': 'lat', 'longitude': 'lon', 'time_stamp': 'time'})
        for v in data.variables:
            # attribs = data[v].attrs
            if v.endswith('_radiance'):
                data[v].attrs['coordinates'] = 'time altitude lat lon'
                data[v].attrs['coverage_content_type'] = 'physicalMeasurement'

                # 'standard_name' in input nc file is not actually a standard name, so update
                data[v].attrs['standard_name'] = 'upwelling_radiance_per_unit_wavelength_in_air'
                # Convert from uint to int to be CF
                data[v] = data[v].astype('int16', copy=False)
                data[v].attrs['valid_min'] = data[v].attrs['valid_min'].astype('int16', copy=False)
                data[v].attrs['valid_max'] = data[v].attrs['valid_max'].astype('int16', copy=False)
                # Remove the ancillary attribute if it exists
                if 'ancillary_variables' in data[v].attrs:
                    del data[v].attrs['ancillary_variables']

        # todo: check lat-lon-time- variable attributes
            if v == 'altitude':
                data[v].attrs['coordinates'] = 'lon lat'
                # data[v].attrs['grid_mapping'] = 'crsWGS84'
                data[v].attrs['standard_name'] = 'altitude'
                data[v].attrs['coverage_content_type'] = 'auxiliaryInformation'
                data[v].attrs['positive'] = 'up'
        
        # correcting lat and lon min and max values
            if v == 'lat':
                data[v].attrs['valid_min'] = -90
                data[v].attrs['valid_max'] = 90
            if v == 'lon':
                data[v].attrs['valid_min'] = -180
                data[v].attrs['valid_max'] = 180

        ## if 'coordinates' in attribs.keys():
        ##     attribs['coordinates'] = 'lon lat'
        ##     attribs['grid_mapping'] = 'crsWGS84'
        ##     attribs['coverage_content_type'] = 'physicalMeasurement'

        # todo: necessary to be CF?
        # # Set grid mapping
        # nc_crs = ncout.createVariable('crsWGS84', "i2")
        # nc_crs.grid_mapping_name = "latitude_longitude"
        # nc_crs.semi_major_axis = 6378137.0
        # nc_crs.inverse_flattening = 298.2572235604902

        # Add global attributes
        # todo: add all metadata needed to create MMD file? so that then we can use nc2mmd from senda project?
        # todo: add link to colhub?
        # Generic global attributes - for all NBS products
        data.attrs.update(cst.global_attributes)
        # Global attributes for all S3 OLCI products
        data.attrs.update(cst.s3_olci_attributes)
        # Global attributes specific to a dataset
        data.attrs['id'] = data.attrs.pop('product_name')
        data.attrs['processing_level'] = self.processing_level
        data.attrs['date_created'] = self.t0.isoformat()
        data.attrs['geospatial_lat_min'] = data['lat'].min().values
        data.attrs['geospatial_lat_max'] = data['lat'].max().values
        data.attrs['geospatial_lon_min'] = data['lon'].min().values
        data.attrs['geospatial_lon_max'] = data['lon'].max().values
        data.attrs['time_coverage_start'] = data.attrs.pop("start_time")
        data.attrs['time_coverage_end'] = data.attrs.pop("stop_time")
        data.attrs['history'] = data.attrs['history'] + f'.{self.t0.isoformat()}:' \
                                                        f' Converted from SAFE to NetCDF/CF by the NBS team.'

        # Ensuring datatypes are correct
        encoding = {
            'time': {
                'dtype': 'int32'
            }
        }
        # Write netcdf file
        t2 = dt.datetime.now(tz=pytz.utc)
        # todo: test with and without load
        # todo: test compression levels. as is, using the encodings from input nc.
        # if we want to change, need to update encoding for each variable
        # todo: use chunks or not? very unconvinced myself
        data.load().to_netcdf(ncout, encoding=encoding)
        logger.debug('Time to write nc file')
        utils.memory_use(t2)

        return True

    def parse_units(self, units):
        """ Parse the units string to extract the base unit and reference date. """
        if 'since' in units:
            base_unit, ref_date_str = units.split('since')
            base_unit = base_unit.strip()
            ref_date_str = ref_date_str.strip().strip('"')
            try:
                ref_date = dt.datetime.strptime(ref_date_str, "%Y-%m-%d %H:%M:%S")
                return base_unit, ref_date
            except ValueError:
                logger.error(f"Error parsing reference date from units: {ref_date_str}")
                return None, None
        else:
            return None, None


if __name__ == '__main__':

    # Log to console
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    log_info = logging.StreamHandler(sys.stdout)
    log_info.setFormatter(
        logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s'))
    logger.addHandler(log_info)

    # workdir = pathlib.Path(
    #     '/home/trygveh/data/satellite_data/Sentinel-3/OCLI/L1/S3B_OL_1_EFR____20200518T162416_20200518T162716_20200519T201048_0179_039_083_1440_LN1_O_NT_002')
    # products = [
    #     'S3B_OL_1_EFR____20200518T162416_20200518T162716_20200519T201048_0179_039_083_1440_LN1_O_NT_002']

    # workdir = pathlib.Path('/home/elodief/Data/NBS/NBS_test_data/test_s3_olci')
    # products = [
    #     'S3B_OL_1_EFR____20211213T082842_20211213T083015_20211214T123309_0092_060_178_1980_LN1_O_NT_002']

    ##workdir = pathlib.Path('/lustre/storeB/project/NBS2/sentinel/production/NorwAREA/netCDFNBS_work/test_environment/test_s3')
    ##products = [
    ##    'S3A_OL_1_EFR____20211124T104042_20211124T104307_20211124T130348_0145_079_051_1980_LN1_O_NR_002',
    ##]

    workdir = pathlib.Path('/home/alessioc/data/S3_products')
    products = ['S3A_OL_1_EFR____20240724T101655_20240724T101955_20240724T121332_0179_115_065_2160_PS1_O_NR_004.SEN3']

    for product in products:
        outdir = workdir / product.split('.')[0]
        outdir.parent.mkdir(parents=False, exist_ok=True)
        s3_obj = S3_olci_reader_and_CF_converter(
            product=product,
            indir=workdir,
            outdir=outdir)
        s3_obj.write_to_NetCDF(outdir, compression_level=1)
