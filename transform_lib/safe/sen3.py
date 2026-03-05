import subprocess as sp
import shutil
from collections import defaultdict
import logging
from pathlib import Path
import lxml.etree as ET
import rasterio
import numpy as np

try:
    from transform_lib.utils import xml_read
except:
    from utils import xml_read

try:
    from transform_lib.safe.safe_base import SAFEFile
except ModuleNotFoundError:
    from safe_base import SAFEFile

logger = logging.getLogger(__name__)

class SEN3File(SAFEFile):
    """
    Subclass for working with Sentinel-3 SEN3 files
    """

    def __init__(self, product, zipdir, tmpdir):

        super().__init__(product, zipdir, tmpdir)
        self.SAFE_dir = (tmpdir / self.product_name).with_suffix('.SEN3')
        self.mainXML = self.SAFE_dir / 'xfdumanifest.xml'

    def prepare_for_use(self):
        """
        Prepare the SEN3File instance for use.
        This includes unzipping the product, reading metadata, creating file lists, and setting up the rasterio source.
        """
        super().prepare_for_use()
        self._list_netcdf_files()
        self._categorise_files()

    def _list_netcdf_files(self):
        # find all netcdf files contained in .SEN3 product
        root = xml_read(self.mainXML)
        self.nc_files = []
        for o in root.findall('.//dataObject'):
            self.nc_files.append(self.SAFE_dir / o.find('.//fileLocation').attrib['href'].split('/')[1])

    def _categorise_files(self):
        if 'OL_1' in self.product_name:
            self.bands = [s for s in self.nc_files if "_radiance.nc" in str(s)]
        elif 'OL_2' in self.product_name:
            self.bands = [s for s in self.nc_files if "_reflectance.nc" in str(s)]
        self.bands.append([s for s in self.nc_files if "geo_coordinates.nc" in str(s)][0])

        self.time_file = [s for s in self.nc_files if "time_coordinates.nc" in str(s)][0]