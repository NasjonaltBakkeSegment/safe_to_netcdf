import yaml
from pathlib import Path
import lxml.etree as ET

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

def load_global_attributes(config_path, platform):
    config = load_config(config_path)
    mission = platform[0:2]
    global_attributes = config['global'] | config [mission] | config [platform]
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