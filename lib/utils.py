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

def chunked_interpolation(y, x, data, yi, xi, chunk_size=1000, overlap=10):
    """
    Perform chunked interpolation with overlap for a 2D dataset.

    Parameters:
        y (np.ndarray): 1D array of y-coordinates.
        x (np.ndarray): 1D array of x-coordinates.
        data (np.ndarray): 2D array of data values to interpolate.
        yi (np.ndarray): 1D array of target y-coordinates.
        xi (np.ndarray): 1D array of target x-coordinates.
        chunk_size (int): Size of each chunk (number of rows/columns).
        overlap (int): Number of rows/columns to overlap between chunks.
    
    Returns:
        np.ndarray: Interpolated data on the (yi, xi) grid.
    """
    # Initialize output array
    result = np.empty((len(yi), len(xi)), dtype=np.float32)

    # Determine number of chunks
    num_chunks_y = len(y) // chunk_size + (len(y) % chunk_size > 0)
    num_chunks_x = len(x) // chunk_size + (len(x) % chunk_size > 0)

    # Process data in chunks with overlap
    for i in range(num_chunks_y):
        for j in range(num_chunks_x):
            # Determine chunk boundaries with overlap
            y_start = max(i * chunk_size - overlap, 0)
            y_end = min((i + 1) * chunk_size + overlap, len(y))
            x_start = max(j * chunk_size - overlap, 0)
            x_end = min((j + 1) * chunk_size + overlap, len(x))

            # Extract chunk of data with overlap
            y_chunk = y[y_start:y_end]
            x_chunk = x[x_start:x_end]
            data_chunk = data[y_start:y_end, x_start:x_end]

            # Interpolate for the chunk
            tck = interpolate.RectBivariateSpline(y_chunk, x_chunk, data_chunk)

            # Determine the valid (non-overlapping) region for this chunk
            valid_y_start = max(overlap, i * chunk_size) - y_start
            valid_y_end = valid_y_start + chunk_size
            valid_x_start = max(overlap, j * chunk_size) - x_start
            valid_x_end = valid_x_start + chunk_size

            # Perform interpolation and store results in the valid region
            result[i * chunk_size:(i + 1) * chunk_size, j * chunk_size:(j + 1) * chunk_size] = \
                tck(yi[valid_y_start:valid_y_end], xi[valid_x_start:valid_x_end])

    return result

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