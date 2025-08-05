# -------------------------------------------------------------------------------
# Name:          s3_olci_l1_reader_and_CF_converter.py
# Purpose:       Read Sentinel-3 OLCI l1 product in SAFE and merge into a single
#                NetCDF file following the CF convention.
#
# Author(s):     Trygve Halsne
# Created:       10.05.2019
# Copyright:     (c) Norwegian Meteorological Institute, 2019
# -------------------------------------------------------------------------------
# import datetime as dt
# import pytz
# import xarray as xr
# import logging
# import sys
# import pathlib
# import numpy as np
# import isodate
# import utils as utils
# import constants as cst

import datetime as dt
import pytz
import xarray as xr
import logging
import numpy as np
import isodate
import utils as utils
import constants as cst


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


class Sentinel3_olci_reader_and_CF_converter:
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
        logger.info((f'Indir = {indir}'))
        logger.info(f'Product = {product}')
        logger.info(f'Path: {file_path}')
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
        self.run()

    def run(self):
        """ Main method for traversing and reading key values from SAFE
            directory.
        """

        # unzip SAFE archive
        utils.uncompress(self)

        return True

    def write_to_NetCDF(self, nc_out, outdir, compression_level, chunk_size=(30, 34)):

        """ Method writing output NetCDF product.
        Args:
            outdir (str): output path where NetCDF file should be stored
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
        ncout = (nc_out / self.product_id).with_suffix('.nc')

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
        # Geospatial max and min
        data.attrs['geospatial_lat_min'] = data['lat'].min().values
        data.attrs['geospatial_lat_max'] = data['lat'].max().values
        data.attrs['geospatial_lon_min'] = data['lon'].min().values
        data.attrs['geospatial_lon_max'] = data['lon'].max().values
        data.attrs['geospatial_vertical_min'] = data['altitude'].min().values
        data.attrs['geospatial_vertical_max'] = data['altitude'].max().values
        data.attrs['geospatial_vertical_positive'] = 'up'
        # Geospatial bounds
        data.attrs['geospatial_bounds'] = (
            f"POLYGON(({data.attrs['geospatial_lon_min']} {data.attrs['geospatial_lat_min']}, "
            f"{data.attrs['geospatial_lon_max']} {data.attrs['geospatial_lat_min']}, "
            f"{data.attrs['geospatial_lon_max']} {data.attrs['geospatial_lat_max']}, "
            f"{data.attrs['geospatial_lon_min']} {data.attrs['geospatial_lat_max']}, "
            f"{data.attrs['geospatial_lon_min']} {data.attrs['geospatial_lat_min']}))"
        )
        data.attrs['geospatial_bounds_crs'] = 'EPSG:4326'  # Assuming WGS84
        data.attrs['geospatial_bounds_vertical_crs'] = 'EPSG:7030'  # Assuming ETRS89 height
        # Time coverage
        data.attrs['time_coverage_start'] = data.attrs.pop("start_time")
        data.attrs['time_coverage_end'] = data.attrs.pop("stop_time")
        # Calculate and add duration
        start_time_str = data.attrs['time_coverage_start'].replace('Z', '+00:00')
        end_time_str = data.attrs['time_coverage_end'].replace('Z', '+00:00')
        start_time = dt.datetime.fromisoformat(start_time_str)
        end_time = dt.datetime.fromisoformat(end_time_str)
        duration = end_time - start_time
        data.attrs['time_coverage_duration'] = isodate.duration_isoformat(duration)
        # Calculate and add resolution
        time_diffs = np.diff(time_values_milliseconds)
        avg_resolution_ms = np.mean(time_diffs)
        avg_resolution_s = avg_resolution_ms / 1000
        data.attrs['time_coverage_resolution'] = isodate.duration_isoformat(dt.timedelta(seconds=avg_resolution_s))

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

