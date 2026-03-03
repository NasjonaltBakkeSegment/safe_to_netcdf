import numpy as np
import datetime as dt
import xarray as xr
import logging

logger = logging.getLogger(__name__)

class S3NetCDFFile():
    """
    Class for working with Sentinel-3 NetCDF files
    """
    def __init__(self, product, directory, compression_level=7, chunk_size=(1, 91, 99)):

        self.directory = directory
        self.product_name = product
        self.compression_level = compression_level
        self.chunk_size = chunk_size
        self.netCDF_path = (self.directory / self.product_name).with_suffix('.nc')
        self.variables = {}
    
    def concatenate_bands(self, bands):
        self.ds = xr.open_mfdataset(bands, join='exact')
    
    def add_time_dimension(self, time_file):
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

        self.ds = xr.combine_by_coords([self.ds, time], combine_attrs='override')
    
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
    
    def rename_variables(self):
        self.ds = self.ds.rename(
            {
                'latitude': 'lat', 
                'longitude': 'lon', 
                'time_stamp': 'time'
            }
        )
    
    # def write_grid_mapping_variables(self):
        # todo: review and add this
        # # Set grid mapping
        # nc_crs = ncout.createVariable('crsWGS84', "i2")
        # nc_crs.grid_mapping_name = "latitude_longitude"
        # nc_crs.semi_major_axis = 6378137.0
        # nc_crs.inverse_flattening = 298.2572235604902
    
    def compute_bounding_box(self):
        self.ds.attrs['geospatial_lat_min'] = self.ds['lat'].min().values
        self.ds.attrs['geospatial_lat_max'] = self.ds['lat'].max().values
        self.ds.attrs['geospatial_lon_min'] = self.ds['lon'].min().values
        self.ds.attrs['geospatial_lon_max'] = self.ds['lon'].max().values
        self.ds.attrs['geospatial_vertical_min'] = self.ds['altitude'].min().values
        self.ds.attrs['geospatial_vertical_max'] = self.ds['altitude'].max().values
        self.ds.attrs['geospatial_bounds_crs'] = 'EPSG:4326'  # Assuming WGS84
        self.ds.attrs['geospatial_bounds_vertical_crs'] = 'EPSG:7030'  # Assuming ETRS89 height
        self.ds.attrs['time_coverage_start'] = self.ds.attrs.pop("start_time")
        self.ds.attrs['time_coverage_end'] = self.ds.attrs.pop("stop_time")
    
    def write_global_attributes(self, global_attributes):
        self.ds.attrs.update(global_attributes)

        self.ds.attrs.pop('netCDF_version', None)
    
    def write_variable_attributes(self, variable_attributes):

        for variable in self.ds.data_vars:
            matched = False

            # First, check for an exact match (full variable name)
            if variable in variable_attributes:
                for attribute, value in variable_attributes[variable].items():
                    self.ds[variable].attrs[attribute] = value
                matched = True

            # If no exact match, check for a affix match
            if not matched:
                for affix in variable_attributes:
                    if affix.startswith('_') and variable.endswith(affix):
                        # affix is a suffix
                        for attribute, value in variable_attributes[affix].items():
                            self.ds[variable].attrs[attribute] = value
                        matched = True
                        break  # Exit the loop once a match is found
                    else:
                        # affix is a prefix
                        # Skip full variable names (not affixes)
                        if not affix.endswith('_'):
                            continue
                        if variable.startswith(affix):
                            for attribute, value in variable_attributes[affix].items():
                                self.ds[variable].attrs[attribute] = value
                            matched = True
                            break  # Exit the loop once a match is found
        
        for variable in self.ds.data_vars:
            if str(variable).endswith('_radiance'):
                self.ds[variable].attrs.pop('ancillary_variables')

    def save_and_close(self):
        encoding = {
            'time': {
                'dtype': 'int32'
            }
        }
        self.ds.to_netcdf(self.netCDF_path, encoding=encoding)