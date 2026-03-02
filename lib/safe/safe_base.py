import subprocess as sp
import shutil
from collections import defaultdict
import logging

try:
    from lib.utils import xml_read
except:
    from utils import xml_read

logger = logging.getLogger(__name__)

class SAFEFile:
    """
    Class for working with SAFE products

    Different missions have different subclasses. See other files in this directory.
    """

    def __init__(self, product, zipdir, tmpdir):
        self.product_name = product
        file_path = zipdir / (product + '.zip')
        if file_path.exists():
            self.input_zip = file_path
        else:
            file_path_safe = zipdir / (product + '.SAFE.zip')
            if file_path_safe.exists():
                self.input_zip = file_path_safe
        self.xmlFiles = defaultdict(list)
        self.read_ok = True

    def prepare_for_use(self):
        """
        Prepare the SAFEFile instance for use.
        This includes unzipping the product, reading metadata, creating file lists, and setting up the rasterio source.
        """
        # Unzip the product if not already unzipped
        if not self.SAFE_dir.is_dir():
            self.unzip()

        self.read_metadata_xml()

    def finalize_usage(self):
        """
        Cleanup the SAFEFile instance after use.
        This includes deleting the uncompressed folder.
        """
        self.delete_uncompressed_folder()

    def unzip(self):
        '''
        Uncompress the zip file. Write it out locally
        '''
        # If zip not extracted yet
        if not self.SAFE_dir.is_dir():
            logger.debug('Starting unzipping SAFE archive')
            self.SAFE_dir.parent.mkdir(parents=False, exist_ok=True)

            sp.run(["/usr/bin/unzip", "-qq", str(self.input_zip), "-d", str(self.SAFE_dir.parent)], check=True)
            logger.debug('Done unzipping SAFE archive')

        extracted_folder_name = self.input_zip.stem + '.SAFE'
        self.extracted_folder_path = self.SAFE_dir.parent / extracted_folder_name

    def read_metadata_xml(self):
        """
        Read main XML file.
        """
        self.root = xml_read(self.mainXML)

    def delete_uncompressed_folder(self):
        """
        Delete the uncompressed SAFE folder and its contents.
        """
        if self.SAFE_dir.is_dir():
            try:
                logger.debug(f"Deleting uncompressed folder: {self.SAFE_dir}")
                shutil.rmtree(self.SAFE_dir)
                logger.debug("Uncompressed folder deleted successfully.")
            except Exception as e:
                logger.error(f"Error deleting uncompressed folder: {e}")
                raise
        else:
            logger.debug(f"Uncompressed folder does not exist: {self.SAFE_dir}")

    def get_global_attributes(self):

        # Product specific metadata
        self.globalAttribs.update({
            'title': self.product_name,
            'product_type': self.globalAttribs.pop("PRODUCT_TYPE"),
        })
    
    def _compute_bounding_box(self, polygon):
        # Add bounding box and polygon
        box = polygon.bounds
        self.globalAttribs.update({
            'geospatial_lon_min': box[0],
            'geospatial_lat_min': box[1],
            'geospatial_lon_max': box[2],
            'geospatial_lat_max': box[3],
        })
    # TODO: Function to extract thumbnail?
    # TODO: Function to create quicklooks?
