"""
Tools
"""

import pathlib
import lxml.etree as ET
import datetime as dt
import resource


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

# Add function to clean work files?
