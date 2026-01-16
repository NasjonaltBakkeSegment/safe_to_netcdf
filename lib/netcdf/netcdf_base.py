import datetime as dt
from netCDF4 import Dataset
from concurrent.futures import ThreadPoolExecutor, as_completed  # Parallel processing utils
import numpy as np  # For working with arrays and chunk processing
import os  # Useful for file path operations

class NetCDFFile():
    """
    Subclass for working with Sentinel-1 SAFE files
    """
    def __init__(self, product, directory, compression_level=7, chunk_size=(1, 91, 99)):

        self.directory = directory
        self.product_name = product
        self.compression_level = compression_level # TODO: Integrate this
        self.chunk_size = chunk_size # TODO: Integrate this
        self.netCDF_path = (self.directory / self.product_name).with_suffix('.nc')
        self.variables = {}

    def initialise(self):
        """
        Initialize the NetCDF file.
        """
        # Create a new NetCDF file or open an existing one
        self.ncout = Dataset(self.netCDF_path, 'w')

    def create_dimensions(self):

        self.ncout.createDimension('time', 1)

    def create_time(self, t, ref='01/01/1981'):

        ref_dt = dt.datetime.strptime(ref, '%d/%m/%Y')
        nc_time = self.ncout.createVariable('time', 'i4', ('time',))

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

        nc_time[:] = seconds_from_ref(t, ref_dt)
        nc_time.units = f"seconds since {ref_dt.strftime('%Y-%m-%d %H:%M:%S')}"

    def init_lat_lon(self):
        self.variables['lat'] = self.ncout.createVariable('lat', 'f4', ('y', 'x',), zlib=True,
        complevel=self.compression_level, chunksizes=self.chunk_size[1:])
        self.variables['lon'] = self.ncout.createVariable('lon', 'f4', ('y', 'x',), zlib=True,
        complevel=self.compression_level, chunksizes=self.chunk_size[1:])

    def write_variable_attributes(self, variable_attributes):
        """
        Write variable attributes to given variable from a dictionary.
        key-value pairs are attributes and their values.
        """
        print('START TALKING')
        print(variable_attributes)
        for variable in self.variables:
            matched = False

            # First, check for an exact match (full variable name)
            if variable in variable_attributes:
                for attribute, value in variable_attributes[variable].items():
                    setattr(self.variables[variable], attribute, value)
                matched = True

            # If no exact match, check for a prefix match
            if not matched:
                for prefix in variable_attributes:
                    # Skip full variable names (not prefixes)
                    if not prefix.endswith('_'):
                        continue
                    if variable.startswith(prefix):
                        for attribute, value in variable_attributes[prefix].items():
                            setattr(self.variables[variable], attribute, value)
                        matched = True
                        break  # Exit the loop once a match is found

            # Optionally, handle unmatched variables (e.g., log a warning)
            if not matched:
                print(f"Warning: No metadata found for variable '{variable}'")

    def write_variable_with_preprocessing(
        self, var_name, array, process_chunk=None, workers=8, sync_every=None, progress_cb=None
    ):
        """
        Process and write a variable in chunks.
        Preprocessing is performed in parallel in chunks.
        Chunks are written sequentially as (and in the order that) processing completes.

        var_name: target variable name (must be 2D or effectively 2D (3D with time=1)).
        array: full numpy array with shape matching the variable’s dimensions.
        process_chunk: callable taking (chunk, y_slice, x_slice, var_name) and returning processed chunk;
                    if None, chunks are written as-is.
        workers: number of parallel workers to process chunks.
        sync_every: call self.ncout.sync() every N chunks.
        progress_cb: optional callback progress_cb(y_start, y_end, x_start, x_end, chunk_index).
        """
        if var_name not in self.variables:
            raise KeyError(f"Variable '{var_name}' not initialized.")
        var = self.variables[var_name]

        # Determine dimensions and validate variable
        if len(var.dimensions) == 3:
            # Handle 3D case (time=y=x or similar, effectively 2D when time=1)
            first_dim, y_dim, x_dim = var.dimensions
            first_dim_size = self.ncout.dimensions[first_dim].size
            if first_dim_size != 1:
                raise ValueError(
                    f"Variable '{var_name}' must be 2D or 3D with first dimension size=1. "
                    f"Dimensions: {var.dimensions}, size of first dim: {first_dim_size}."
                )
        elif len(var.dimensions) == 2:
            # Handle standard 2D case
            y_dim, x_dim = var.dimensions
        else:
            # Invalid dimensionality
            raise ValueError(
                f"Variable '{var_name}' must be 2D or 3D with first dimension size=1. "
                f"Found dimensions: {var.dimensions}."
            )

        y_max = self.ncout.dimensions[y_dim].size
        x_max = self.ncout.dimensions[x_dim].size

        # Ensure the array shape matches the variable's dimensions
        if len(array.shape) == 3 and array.shape[0] == 1:
            # Array is effectively 2D (3D with first dimension 1)
            if array.shape[1:] != (y_max, x_max):
                raise ValueError(
                    f"Array shape {array.shape} does not match variable dims (1, {y_max}, {x_max})."
                )
        elif len(array.shape) == 2:
            # Array is already 2D
            if array.shape != (y_max, x_max):
                raise ValueError(
                    f"Array shape {array.shape} does not match variable dims ({y_max}, {x_max})."
                )
        else:
            # Invalid shape
            raise ValueError(
                f"Invalid array shape {array.shape}. Must match dimensions "
                f"(1, {y_max}, {x_max}) or ({y_max}, {x_max})."
            )
        
        # Extract chunk sizes from self.chunk_size
        # Assuming chunk_size format (optional_band_dimension, y_chunk_size, x_chunk_size)
        y_chunk_size, x_chunk_size = self.chunk_size[1], self.chunk_size[2]

        # Identity preprocessor by default
        def _identity(chunk, y_slice, x_slice, name):
            return chunk
        process_chunk = process_chunk or _identity

        # Build chunk windows
        windows = []
        for y_start in range(0, y_max, max(1, y_chunk_size)):
            y_end = min(y_start + max(1, y_chunk_size), y_max)
            for x_start in range(0, x_max, max(1, x_chunk_size)):
                x_end = min(x_start + max(1, x_chunk_size), x_max)
                windows.append((y_start, y_end, x_start, x_end))

        # Parallel processing with thread-based executor
        Executor = ThreadPoolExecutor

        chunk_index = 0
        sync_counter = 0

        with Executor(max_workers=workers) as ex:
            futures = {}
            for win in windows:
                # Submit parallel preprocessing tasks
                ys, ye, xs, xe = win
                futures[ex.submit(
                    process_chunk,
                    array[ys:ye, xs:xe],
                    slice(ys, ye),
                    slice(xs, xe),
                    var_name
                )] = win

            # Sequentially loop through completed tasks as they finish
            for fut in as_completed(futures):
                ys, ye, xs, xe = futures[fut]
                chunk_out = fut.result()
                self.write_variable_chunk(var_name, chunk_out, ys, xs)
                chunk_index += 1
                sync_counter += 1

                # Sync periodically if specified
                if sync_every and (sync_counter % sync_every == 0):
                    self.ncout.sync()

                # Call the progress callback, if provided
                if progress_cb:
                    progress_cb(ys, ye, xs, xe, chunk_index)

        # Final sync after all chunks are written
        self.ncout.sync()

    def write_variable_chunk(self, var_name, data_chunk, y_start, x_start):
        """
        Write a single chunk into a variable at the given offsets.
        Assumes the provided array is in the same dimension order as the variable.
        Handles 2D variables or 3D variables with the first dimension = 1 (effectively 2D).
        """
        if var_name not in self.variables:
            raise KeyError(f"Variable '{var_name}' not initialized.")
        var = self.variables[var_name]

        # Determine the dimensions and validate
        if len(var.dimensions) == 3:
            # Handle 3D (effectively 2D) variables
            first_dim, y_dim, x_dim = var.dimensions
            first_dim_size = self.ncout.dimensions[first_dim].size
            if first_dim_size != 1:
                raise ValueError(
                    f"Variable '{var_name}' must be 2D or 3D with first dimension size 1. "
                    f"Given: dims={var.dimensions}, first_dim_size={first_dim_size}"
                )
        elif len(var.dimensions) == 2:
            # Handle standard 2D variables
            y_dim, x_dim = var.dimensions
        else:
            # Invalid dimensions
            raise ValueError(
                f"Variable '{var_name}' has invalid dimensions: {var.dimensions}. "
                f"Expected 2D or 3D with first dimension size 1."
            )

        # Retrieve the sizes of the y and x dimensions
        y_max = self.ncout.dimensions[y_dim].size
        x_max = self.ncout.dimensions[x_dim].size

        # Calculate the boundaries of the chunk to write
        chunk_h, chunk_w = data_chunk.shape
        y_end = y_start + chunk_h
        x_end = x_start + chunk_w

        # Ensure chunk indices are within bounds
        if y_start < 0 or x_start < 0 or y_end > y_max or x_end > x_max:
            raise IndexError(f"Chunk [{y_start}:{y_end}, {x_start}:{x_end}] exceeds bounds {(y_max, x_max)}")

        # Perform the write
        if len(var.dimensions) == 3:
            # Handle "effectively 2D" (3D with size=1 in first dim)
            var[0, y_start:y_end, x_start:x_end] = data_chunk
        else:
            # Handle standard 2D variables
            var[y_start:y_end, x_start:x_end] = data_chunk

    def write_global_attributes(self, global_attributes):
        
        t0 = dt.datetime.now(dt.timezone.utc).isoformat().replace("+00:00", "Z")

        global_attributes.update({
            'date_metadata_modified': t0,
            'date_metadata_modified_type': 'Created',
            'date_created': t0,
            'history': f'{t0}: Converted from SAFE to NetCDF by NBS team.',
        })

        if global_attributes['geospatial_lat_min'] > 70:
            global_attributes['collection'] += ',SIOS'
        
        self.ncout.setncatts(global_attributes)
    
    def close(self):
        self.ncout.close()