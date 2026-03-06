import yaml
import datetime as dt
from pathlib import Path
import lxml.etree as ET
import numpy as np
from scipy import interpolate

def xml_read(xml_file):
    """
    Read XML file.
    Args:
        xml_file [pathlib]): filepath to an xml file
    Returns:
        lxml.etree._Element or None if file missing
    """
    if not Path(xml_file).is_file():
        return None
    tree = ET.parse(str(xml_file))
    root = tree.getroot()
    return root

def load_variable_attributes(config_path, mission):
    '''
    Load variable attributes from YAML config file into a Python dictionary
    '''
    config = load_config(config_path)
    valid_missions = {'S1', 'S2', 'S3'}
    if mission not in valid_missions:
        raise ValueError(f"Invalid mission '{mission}'. Mission must be one of {valid_missions}.")
    if mission not in config:
        raise KeyError(f"The mission '{mission}' is not found in the configuration file.")
    else:
        variable_attributes = config['all'] | config[mission]
        return variable_attributes

def load_global_attributes(config_path, product_name):
    platform = product_name[0:3]
    config = load_config(config_path)
    mission = platform[0:2]
    global_attributes = config['all'] | config [mission] | config [platform]
    if mission == 'S3':
        if '_OL_' in product_name:
            global_attributes = global_attributes | config['S3_OL']
        elif '_SL_' in product_name:
            global_attributes = global_attributes | config['S3_SL']
        elif '_SR_' in product_name:
            global_attributes = global_attributes | config['S3_SR']
        elif '_SY_' in product_name:
            global_attributes = global_attributes | config['S3_SY']

    return global_attributes

def load_config(config_path):
    """
    Load a YAML configuration file into a Python dictionary
    """
    try:
        with open(config_path, 'r') as file:
            config = yaml.safe_load(file)
        return config
    except FileNotFoundError:
        print(f"Error: The file at {config_path} was not found.")
        return None
    except yaml.YAMLError as exc:
        print(f"Error loading YAML file: {exc}")
        return None

def get_key(my_dict,val):
    for key, value in my_dict.items():
         if val == value:
             return key

    return "There is no such Key"

def update_global_attributes(global_attributes):

    t0 = dt.datetime.now(dt.timezone.utc).isoformat().replace("+00:00", "Z")

    global_attributes.update({
        'date_metadata_modified': t0,
        'date_metadata_modified_type': 'Created',
        'date_created': t0,
        'history': f'{t0}: Converted from SAFE to NetCDF by NBS team.',
    })

    if global_attributes['geospatial_lat_max'] > 70:
        global_attributes['collection'] += ',SIOS'

    return global_attributes

def chunked_interpolation(
    y, x, lat2d, yi, xi, chunk_size=1000, overlap=10, kx=3, ky=3
):
    """
    Chunked 2D spline interpolation using RectBivariateSpline.

    Parameters
    ----------
    y, x : 1D arrays
        Original grid coordinates.
    lat2d : 2D array (len(y), len(x))
        Field to interpolate.
    yi, xi : 1D arrays
        Target grid coordinates.
    chunk_size : int
        Number of target points per chunk (along each dimension).
    overlap : int
        Number of source grid points to overlap to avoid edge effects.
    kx, ky : int
        Spline degrees (default cubic).

    Returns
    -------
    out : 2D array (len(yi), len(xi))
    """

    yi = np.asarray(yi)
    xi = np.asarray(xi)

    nyi = len(yi)
    nxi = len(xi)

    out = np.empty((nyi, nxi), dtype=lat2d.dtype)

    for j0 in range(0, nyi, chunk_size):
        j1 = min(j0 + chunk_size, nyi)

        # Determine y-range in original grid needed for this chunk
        y_min = yi[j0]
        y_max = yi[j1 - 1]

        iy0 = np.searchsorted(y, y_min) - overlap
        iy1 = np.searchsorted(y, y_max) + overlap

        iy0 = max(0, iy0)
        iy1 = min(len(y), iy1)

        y_sub = y[iy0:iy1]

        for i0 in range(0, nxi, chunk_size):
            i1 = min(i0 + chunk_size, nxi)

            x_min = xi[i0]
            x_max = xi[i1 - 1]

            ix0 = np.searchsorted(x, x_min) - overlap
            ix1 = np.searchsorted(x, x_max) + overlap

            ix0 = max(0, ix0)
            ix1 = min(len(x), ix1)

            x_sub = x[ix0:ix1]

            lat_sub = lat2d[iy0:iy1, ix0:ix1]

            # Build local spline
            tck = interpolate.RectBivariateSpline(
                y_sub, x_sub, lat_sub, kx=kx, ky=ky
            )

            # Evaluate only on core target region
            out[j0:j1, i0:i1] = tck(yi[j0:j1], xi[i0:i1])

    return out

def multiply_2d_arrays_in_chunks(array1, array2, chunk_size):
    """
    Multiplies two 2D arrays element-wise in chunks to reduce memory consumption.

    Args:
        array1 (np.ndarray): The first 2D array.
        array2 (np.ndarray): The second 2D array.
        chunk_size (int): The size of the chunks to process (applies to both rows and columns).

    Returns:
        np.ndarray: The resulting 2D array after element-wise multiplication.
    """
    # Ensure the inputs are 2D arrays
    if array1.ndim != 2 or array2.ndim != 2:
        raise ValueError("Both array1 and array2 must be 2D arrays.")

    # Ensure the arrays have the same shape
    if array1.shape != array2.shape:
        raise ValueError("array1 and array2 must have the same shape.")

    # Initialize the result array
    result_array = np.zeros_like(array1)

    # Process in chunks
    rows, cols = array1.shape
    for row_start in range(0, rows, chunk_size):
        row_end = min(row_start + chunk_size, rows)
        for col_start in range(0, cols, chunk_size):
            col_end = min(col_start + chunk_size, cols)

            # Extract chunks
            chunk1 = array1[row_start:row_end, col_start:col_end]
            chunk2 = array2[row_start:row_end, col_start:col_end]

            # Perform element-wise multiplication
            result_array[row_start:row_end, col_start:col_end] = chunk1 * chunk2

    return result_array