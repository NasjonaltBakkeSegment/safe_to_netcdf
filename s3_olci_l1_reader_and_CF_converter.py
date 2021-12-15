#-------------------------------------------------------------------------------
# Name:          s3_olci_l1_reader_and_CF_converter.py
# Purpose:       Read Sentinel-3 OLCI l1 product in SAFE and merge into a single
#                NetCDF file following the CF convention.
#
# Author(s):     Trygve Halsne
# Created:       10.05.2019
# Copyright:     (c) Norwegian Meteorological Institute, 2019
#-------------------------------------------------------------------------------
import lxml.etree as ET
import netCDF4
from datetime import datetime
import numpy as np
import os
import subprocess
from collections import defaultdict
from glob import glob
import resource

"""
    FIXES to be solved
    - tie points having latitude longitude in tie_geo_coordinates. Should be separat group?
    - removed pixels variables uses same names as radiance raster bands. Hence should be fixed.
    - should follow swath convention (https://github.com/Unidata/EC-netCDF-CF/blob/master/swath/swath.adoc#radiometric-swath-data-encodings) which is cf v.1.7 compatible
    - should be ACDD compliant

    ISSUES:
    - int64: not supported in dap2. will come in dap4 https://www.unidata.ucar.edu/support/help/MailArchives/thredds/msg01762.html
"""
class S3_olci_reader_and_CF_converter:
    """ S3 OLCI object for merging SAFE files and creating single NetCDF/CF """

    def __init__(self, SAFE_file, SAFE_outpath):
        """ Initializer

        Args:
            SAFE_file      (str): Input SAFE file
            SAFE_id        (str): SAFE filename
            SAFE_outpath   (str): Output directory for unzipeed SAFE file
            SAFE_path      (str): Absolute path to unzipped SAFE directory
            xmlFiles      (dict): defaultdict containing xml files
            globalAttribs (dict): global attributes for output NetCDF/CF
            src
            t0           (float): Init time
            SAFE_structure (str): XML formatted string showing SAFE structure

        """
        self.SAFE_file = SAFE_file
        self.SAFE_id = SAFE_file.split('/')[-1].split('.')[0]
        self.SAFE_outpath = SAFE_outpath
        self.SAFE_path = None
        self.xmlFiles = defaultdict(list)
        self.globalAttribs = {}
        self.src = None
        self.t0 = datetime.now()
        self.SAFE_structure = None

        self.main()

    def main(self):
        """ Main method for traversing and reading key values from SAFE
            directory.
        """

        # 1) Fetch main file
        self.xmlFiles['manifest'] =self.uncompress('xfdumanifest.xml')

        # SAFE file structure
        self.SAFE_structure = self.list_product_structure()
        return True

    def write_to_NetCDF(self, nc_outpath, compression_level,chunk_size=(30,34),groups=False):
        """ Method writing output NetCDF product.

        Args:
            nc_outpath (str): output path where NetCDF file should be stored
            compression_level (int): compression level on output NetCDF file (1-9)
            chunk_size (tuple): chunk size in output NetCDF. Format (rows, columns)
        """

        print("------------START CONVERSION FROM SAFE TO NETCDF-------------")

        # Status
        print('\nCreating NetCDF file')
        self.memoryUsage()
        print(datetime.now() - self.t0)

        #find all netcdf files contained in .SAFE product
        tree = ET.parse(self.xmlFiles['manifest'])
        root = tree.getroot()
        root.findall('.//dataObject')
        nc_files = []
        for o in root.findall('.//dataObject'):
            fl =o.find('.//fileLocation')
            nc_files.append(fl.attrib['href'])

        #group_data = ['instrument_data', 'removed_pixels', 'tie_', 'xfdumanifest']
        group_data = ['removed_pixels', 'tie_', 'xfdumanifest']
        if not groups:
            missing_files = []
        #group_data = ['removed_pixels', 'tie_', 'xfdumanifest']

        # iterate through files and create group for each file:
        myDims = {}
        globalAttribs = {}
        bands_name = []
        data_types_lookup = {'int8':'i1','int16':'i2','int32':'i4','int64':'i4'}
        data_types_lookup_numpy_one_up = {'uint8':np.int16,'uint16':np.int32,'uint32':np.int64}
        data_types_lookup_numpy_one_up = {'uint8':np.int8,'uint16':np.int16,'uint32':np.int32}

        # Test keeping original variables
        data_types_lookup = {'int8':'u1','int16':'u2','int32':'u4','int64':'u4'}
        data_types_lookup_numpy_one_up = {'uint8':np.int16,'uint16':np.int32,'uint32':np.int64}
        data_types_lookup_numpy_one_up = {'uint8':np.uint8,'uint16':np.uint16,'uint32':np.uint32}

        # output filename
        out_netcdf = '{}{}.nc'.format(nc_outpath,self.SAFE_id)

        with (netCDF4.Dataset(out_netcdf, 'w', format='NETCDF4')) as ncout:

            for f in nc_files:
                fname = f.partition('.')[-1]
                src_in = netCDF4.Dataset(self.SAFE_path+fname,'r')

                #iterate dimensions
                if not any([element in fname for element in group_data]):
                    for d in src_in.dimensions:

                        tmp_dim = src_in.dimensions[d]
                        if tmp_dim.name in myDims.keys():
                            if tmp_dim.size != myDims[tmp_dim.name]:
                                print('Same dimension names with different sizes..')
                        else:
                            myDims[tmp_dim.name] = tmp_dim.size

                        if not d in ncout.dimensions.keys():
                            ncout.createDimension(tmp_dim.name,tmp_dim.size)
                        else:
                            continue

                    #iterate and store attributes

                    #iterate and store variables
                    for v in src_in.variables:
                        attribs = {}
                        if not v in ncout.variables:
                            # create variable
                            tmp_var = src_in.variables[v]
                            if v=='latitude':
                                v='lat'
                            elif v=='longitude':
                                v='lon'

                            if v=='time_stamp':
                                #tmp_varout = ncout.createVariable(varname=v, datatype='i4', dimensions=tmp_var.dimensions,
                                tmp_varout = ncout.createVariable(varname=v, datatype='f8', dimensions=tmp_var.dimensions,
                                    complevel=compression_level)
                            else:
                                if ('uint' in str(tmp_var.dtype)):
                                    dt = data_types_lookup[str(tmp_var.dtype).replace('u','')]
                                    tmp_varout = ncout.createVariable(varname=v, datatype=dt, dimensions=tmp_var.dimensions,
                                        complevel=compression_level, chunksizes=chunk_size)#, fill_value=netCDF4.default_fillvals[dt])#, chunksizes=(30,34))#, fill_value=tmp_var._FillValue)
                                else:
                                    tmp_varout = ncout.createVariable(varname=v, datatype=tmp_var.datatype, dimensions=tmp_var.dimensions,
                                    complevel=compression_level)#, fill_value=tmp_var._FillValue)

                            # set attributes
                            #print(tmp_var.ncattrs())
                            for var_attr in tmp_var.ncattrs():
                                if var_attr=='_FillValue':
                                    #print('Neglecting {}'.format(var_attr))
                                    if ('uint' in str(tmp_var.getncattr(var_attr).dtype)):
                                        dt_attrib = data_types_lookup_numpy_one_up[str(tmp_var.dtype)]
                                        attribs[var_attr] = tmp_var.getncattr(var_attr).astype(dt_attrib)
                                    else:
                                        attribs[var_attr] = tmp_var.getncattr(var_attr)
                                elif var_attr =='standard_name':
                                    sn = tmp_var.getncattr(var_attr)
                                    if sn == 'toa_upwelling_spectral_radiance':
                                        attribs[var_attr] = "upwelling_spectral_radiance_in_air"
                                elif var_attr =='flag_masks':
                                    flag_masks = [m for m in tmp_var.getncattr(var_attr)]
                                    attribs[var_attr] = flag_masks
                                    #for m in tmp_var.getncattr(var_attr):
                                    #    print(np.int32(m))
                                elif var_attr == 'ancillary_variables':
                                    if not 'err' in tmp_var.getncattr(var_attr):
                                        attribs[var_attr] = tmp_var.getncattr(var_attr)
                                elif var_attr == 'valid_min':
                                    if ('uint' in str(tmp_var.getncattr(var_attr).dtype)):
                                        dt_attrib = data_types_lookup_numpy_one_up[str(tmp_var.dtype)]
                                        dt_attrib = np.int
                                        attribs[var_attr] = tmp_var.getncattr(var_attr).astype(dt_attrib)
                                    else:
                                        attribs[var_attr] = tmp_var.getncattr(var_attr)
                                elif var_attr == 'valid_max':
                                    if ('uint' in str(tmp_var.getncattr(var_attr).dtype)):
                                        dt_attrib = data_types_lookup_numpy_one_up[str(tmp_var.dtype)]
                                        dt_attrib = np.int
                                        attribs[var_attr] = tmp_var.getncattr(var_attr).astype(dt_attrib)
                                    else:
                                        attribs[var_attr] = tmp_var.getncattr(var_attr)


                                else:
                                    #print('\n')
                                    #print(tmp_var.dtype)
                                    #print(tmp_var.getncattr(var_attr))
                                    attribs[var_attr] = tmp_var.getncattr(var_attr)

                            if v == 'lat':
                                attribs['standard_name'] = 'latitude'
                            elif v == 'lon':
                                attribs['standard_name'] = 'longitude'

                            if v == 'altitude':
                                attribs['coordinates'] = 'lon lat'
                                attribs['grid_mapping'] = 'crsWGS84'
                                attribs['standard_name'] = 'altitude'
                                attribs['coverage_content_type'] = 'auxiliaryInformation'
                                attribs['positive'] = 'up'


                            if 'coordinates' in attribs.keys():
                                attribs['coordinates'] = 'lon lat'
                                attribs['grid_mapping'] = 'crsWGS84'
                                attribs['coverage_content_type'] = 'physicalMeasurement'

                            if v=='time_stamp':
                                #print(v)
                                attribs['_FillValue'] = netCDF4.default_fillvals['f8']
                                print(attribs)
                                tmp_varout.setncatts(attribs)
                            else:
                                tmp_varout.setncatts(attribs)
                            # set values
                            #tmp_varout[:] = tmp_var[:]
                        else:
                            #continue
                            print("{} already exists..".format(v))
                        print(v)
                else:
                    if groups:
                        #################
                        # CF 2.0 style with groups..
                        #################
                        tmp_group = ncout.createGroup(fname.split('/')[-1].split('.')[0])

                        #create dimensions
                        for d in src_in.dimensions:
                            tmp_dim = src_in.dimensions[d]
                            tmp_group.createDimension(tmp_dim.name,tmp_dim.size)
                            if tmp_dim.name in myDims.keys():
                                if tmp_dim.size != myDims[tmp_dim.name]:
                                    print('Same dimension names with different sizes..')
                            else:
                                myDims[tmp_dim.name] = tmp_dim.size


                        #create variables
                        for v in src_in.variables:
                            f_val = u'_FillValue'
                            tmp_var = src_in.variables[v]
                            tmp_attribs = {}
                            if  f_val in tmp_var.ncattrs():
                                tmp_outVar = tmp_group.createVariable(tmp_var.name,tmp_var.dtype,tmp_var.dimensions,zlib=True, fill_value=tmp_var.getncattr(f_val),complevel=compression_level)
                                for key in tmp_var.ncattrs():
                                    tmp_attribs[key] = tmp_var.getncattr(key)
                                tmp_outVar.setncatts(tmp_attribs)
                            else:
                                tmp_outVar = tmp_group.createVariable(tmp_var.name,tmp_var.dtype,tmp_var.dimensions,zlib=True)
                                for key in tmp_var.ncattrs():
                                    tmp_attribs[key] = tmp_var.getncattr(key)
                                tmp_outVar.setncatts(tmp_attribs)

                            # set values
                            tmp_outVar[:] = tmp_var[:]

                        # set global attributes
                        globalAttribs_attribs = {}
                        for key in src_in.ncattrs():
                            globalAttribs_attribs[key] = src_in.getncattr(key)
                        tmp_group.setncatts(globalAttribs_attribs)

                    else:
                        missing_files.append(fname.split('/')[-1])

                # Add global attributes
                if len(globalAttribs.keys()) == 0:
                    for ga in src_in.ncattrs():
                        globalAttribs[ga] = src_in.getncattr(ga)

                src_in.close()

            # Set grid mapping
            nc_crs= ncout.createVariable('crsWGS84',"i2")
            nc_crs.grid_mapping_name = "latitude_longitude" ;
            nc_crs.semi_major_axis= 6378137.0 ;
            nc_crs.inverse_flattening= 298.2572235604902 ;

            # Set global attributes
            globalAttribs['id'] = globalAttribs['product_name']
            globalAttribs['naming_authority'] = "EUMETSAT, ESA: Sentinel 3 PDGS, File Naming Convention"
            globalAttribs.pop('product_name')
            globalAttribs['acknowledgement'] = "Copernicus EU"
            globalAttribs['license'] = "Freely Distributed"
            globalAttribs['standard_name_vocabulary'] =  'CF Standard Name Table v69'
            globalAttribs['date_created'] = datetime.utcnow().isoformat()
            globalAttribs['project'] = "Norwegian National Ground Segment for Satellite Data"
            globalAttribs['geospatial_lat_min'] = ncout['lat'][:].min()
            globalAttribs['geospatial_lat_max'] = ncout['lat'][:].max()
            globalAttribs['geospatial_lon_min'] = ncout['lon'][:].min()
            globalAttribs['geospatial_lon_max'] = ncout['lon'][:].max()
            globalAttribs['time_coverage_start'] = globalAttribs["start_time"]
            globalAttribs['time_coverage_end'] = globalAttribs["stop_time"]
            globalAttribs.pop('start_time')
            globalAttribs.pop('stop_time')
            globalAttribs['Conventions'] = "CF-1.7, ACDD-1.3"
            globalAttribs['title'] = "Sentinel-3 OCLI product in NetCDF/CF"
            globalAttribs['summary'] = str('Sentinel-3 Ocean and Land Color Instrument product in NetCDF/CF.')
            #globalAttribs['keywords'] = str('[Earth Science, Atmosphere, Atmospheric radiation, Reflectance]')
            globalAttribs['keywords'] = ["SPECTRAL/ENGINEERING > INFRARED WAVELENGTHS > REFLECTED INFRARED",
                                         "SPECTRAL/ENGINEERING > PLATFORM CHARACTERISTICS",
                                         "SPECTRAL/ENGINEERING > PLATFORM CHARACTERISTICS > ATTITUDE CHARACTERISTICS",
                                         "SPECTRAL/ENGINEERING > VISIBLE WAVELENGTHS > VISIBLE RADIANCE"]
            globalAttribs['keywords_vocabulary'] = str("GCMD Science Keywords")
            #globalAttribs['keywords'] = str('[To, Be, Decided]')
            #globalAttribs['keywords_vocabulary'] = str("TBD")
            globalAttribs['institution'] = str("Norwegian Meteorological Institute")
            globalAttribs['source'] = str("surface observation")
            globalAttribs['history'] = globalAttribs['history'] + ' \n'+ str(datetime.now().strftime("%A %B %d. %T %Y") +
                                 ". Converted from SAFE to NetCDF/CF by the NBS team.")
            if "OLCI Level 1" in globalAttribs['references']:
                globalAttribs['processing_level'] = "Level-1"
            elif "OLCI Level 2" in globalAttribs['references']:
                globalAttribs['processing_level'] = "Level-2"
            globalAttribs['references'] = globalAttribs['references'] + str(" - https://sentinel.esa.int/web/sentinel/user-guides/sentinel-3-olci")

            if groups:
                globalAttribs['comment'] = str("Copernicus Sentinel-3 OLCI product merged to a single NetCDF file.")
            else:
                globalAttribs['comment'] = str("NOTE: This products contains SELECTED information from Copernicus Sentinel-3 OLCI product in SAFE merged to a single NetCDF file. The following files are missing: {}".format(', '.join(missing_files)))

            ncout.setncatts(globalAttribs)
            ncout.sync()


    def uncompress(self,xmlReturn):
        """ Uncompress SAFE zip file.

        Args:
            xmlReturn (str): Name of main XML file

        Returns:
            main XML file if success, (bool) False otherwise
        """
        SAFE_outpath = self.SAFE_outpath
        SAFE_file = self.SAFE_file
        SAFE_id = self.SAFE_id
        fdirName = '%s%s.SEN3' % (SAFE_outpath,SAFE_id)
        print(fdirName)

        if os.path.isdir(fdirName) == True:
            # zipfile already extracted find main xmlfile
            files_according_to_xmlReturn = glob(fdirName + '/' + xmlReturn)
            if len(files_according_to_xmlReturn) != 1:
                print('Did not find a file similar to: %s' % xmlReturn)
                return False
            else:
                xmlFile = files_according_to_xmlReturn[0]

        else:
            # unzip file
            cmd = '/usr/bin/unzip %s -d %s' % (SAFE_file,SAFE_outpath)
            subprocess.call(cmd, shell=True)
            if not os.path.isdir(fdirName):
                print('Error unzipping file %s' % (SAFE_file))
                return False

            else:
                files_according_to_xmlReturn = glob(fdirName + '/' + xmlReturn)
                if len(files_according_to_xmlReturn) != 1:
                    print('Did not find a file similar to: %s' % xmlReturn)
                    return False
                else:
                    xmlFile = files_according_to_xmlReturn[0]

        if os.path.isfile(xmlFile) == False:
            print('Error unzipping file %s' % (SAFE_file))
            return False

        self.SAFE_path = fdirName
        return xmlFile


    def list_product_structure(self, startpath=None):
        """ Traverse SAFE file structure (or any file structure) and
            creates a string representation of a XML file containing the
            file structure.

        Args:
            startpath (str): path from which to create file structure

        Returns:
            (str): string representation of a xml file showing file structure.
        """

        if not startpath:
            startpath = self.SAFE_path

        root_path = startpath.split('/')[-1]
        ET_root = ET.Element(root_path)
        # Create dictionary that contains all elements
        elements = {root_path:ET_root}

        # Iterate throgh the file structure
        for root, dirs, files in os.walk(startpath):
            level = root.replace(startpath, '').count(os.sep)
            current_xpath = str('/' + root_path + root.replace(startpath,''))
            current_path = root.replace(startpath,'')

            # Create all folder elements
            if dirs:
                for dir in dirs:
                    element = ET.SubElement(elements[current_xpath.strip('/')], dir)
                    elements[str(current_xpath.strip('/') + '/' + dir)] = element

            # Add all files in the correct folder
            for f in files:
                sub_element = ET.SubElement(elements[current_xpath.strip('/')], 'file')
                sub_element.text = f

        return ET.tostring(ET_root)

    def memoryUsage(self):
        max_memory = (resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)/1000000.
        print("Memory usage so far: {} Gb".format(max_memory))


# TEST
# input parameters
safe_file = "/home/trygveh/data/satellite_data/Sentinel-3/OCLI/L1/S3B_OL_1_EFR____20200518T162416_20200518T162716_20200519T201048_0179_039_083_1440_LN1_O_NT_002"
safe_output = "/home/trygveh/data/satellite_data/Sentinel-3/OCLI/L1/"
#safe_file = "/home/trygveh/documents/satellite_data/Sentinel-3/OCLI/l2/S3A_OL_2_LFR____20191028T064724_20191028T065024_20191029T110920_0179_051_020_3060_LN1_O_NT_002.zip"
#safe_output = "/home/trygveh/documents/satellite_data/Sentinel-3/OCLI/l2/"

s3_obj = S3_olci_reader_and_CF_converter(SAFE_file=safe_file, SAFE_outpath=safe_output)
s3_obj.write_to_NetCDF('/home/trygveh/tmp/',compression_level=9)
