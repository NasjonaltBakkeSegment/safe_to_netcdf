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


def seconds_from_ref(t):
    """
    Computes the difference in seconds between input date and a reference date (01/01/1981)
    Args:
        t: date as a string
    Returns:
        integer
    """
    mytime = dt.datetime.strptime(t, '%Y-%m-%dT%H:%M:%S.%f')
    mytime_ref = dt.datetime(1981, 1, 1)
    return int((mytime - mytime_ref).total_seconds())


