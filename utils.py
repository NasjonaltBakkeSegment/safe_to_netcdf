"""
Tools
"""

import pathlib
import lxml.etree as ET


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


