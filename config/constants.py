#!/usr/bin/python3

'''
Constants describing Sentinel satellites
'''

# ------------- All Sentinels -------------

global_attributes = {
     'project': 'Norwegian National Ground Segment for Satellite Data',
     'institution': 'Norwegian Meteorological Institute (METNO)',
     'source': 'surface observation',
     'license': 'Freely Distributed',
     'acknowledgement': 'Copernicus EU',
     'creator_type': 'institution',
     'creator_institution': 'Norwegian Meteorological Institute',
     'creator_name': 'NBS Helpdesk',
     'creator_email': 'nbs-helpdesk@met.no',
     'creator_url': 'https://met.no',
     'publisher_name': 'NBS Helpdesk',
     'publisher_email': 'nbs-helpdesk@met.no',
     'publisher_url': 'https://met.no'
}

# ------------- Sentinel 3 -------------

s3_olci_attributes = {
    'Conventions': 'CF-1.7, ACDD-1.3',
    'naming_authority': 'EUMETSAT, ESA: Sentinel 3 PDGS, File Naming Convention',
    'reference': 'https://sentinel.esa.int/web/sentinel/user-guides/sentinel-3-olci',
    'standard_name_vocabulary': 'CF Standard Name Table v69',
    'comment': 'This product contains selected information from Copernicus Sentinel-3 '
               'OLCI product. To obtain all data, please refer to the SAFE format available '
               'at colhub.met.no',
    'title': 'Sentinel-3 OCLI product in NetCDF/CF',
    'summary': 'Sentinel-3 Ocean and Land Color Instrument product in NetCDF/CF.',
    'keywords_vocabulary': 'GCMD Science Keywords',
    'keywords': [
        'SPECTRAL/ENGINEERING > INFRARED WAVELENGTHS > REFLECTED INFRARED',
        'SPECTRAL/ENGINEERING > VISIBLE WAVELENGTHS > VISIBLE RADIANCE'
    ]
}

# ------------- Sentinel 2 -------------

s2_attributes = {
    'Conventions': 'CF-1.8, ACDD-1.3',
    'reference': 'https://sentinel.esa.int/web/sentinel/user-guides/sentinel-2-msi',
    'standard_name_vocabulary': 'CF Standard Name Table v69',
    'keywords_vocabulary': 'GCMD Science Keywords',
    'keywords': '[Earth Science, Atmosphere, Atmospheric radiation, Reflectance]'
}

s2_level_attributes = {
    'Level-1C': {
        'title': 'Sentinel-2 Level-1C data.',
        'summary': 'Sentinel-2 Multi-Spectral Instrument Level-1C product.'
    },
    'Level-2A': {
        'title': 'Sentinel-2 Level-2A data.',
        'summary': 'Sentinel-2 Multi-Spectral Instrument Level-2A product.'
    }
}

platform_id = {'Sentinel-2A': 0, 'Sentinel-2B': 1, 'Sentinel-2C': 2, 'Sentinel-2D': 3}

# Bands order
s2_bands_order = {
    0: 'B1',
    1: 'B2',
    2: 'B3',
    3: 'B4',
    4: 'B5',
    5: 'B6',
    6: 'B7',
    7: 'B8',
    8: 'B8A',
    9: 'B9',
    10: 'B10',
    11: 'B11',
    12: 'B12'
}

# Bands aliases
s2_bands_aliases = {
    'B01': 'B1',
    'B02': 'B2',
    'B03': 'B3',
    'B04': 'B4',
    'B05': 'B5',
    'B06': 'B6',
    'B07': 'B7',
    'B08': 'B8',
    'B8A': 'B8A',
    'B09': 'B9',
    'B10': 'B10',
    'B11': 'B11',
    'B12': 'B12'
}

# for specific layers: If nb. raster bands > 1, use the naming convention below.
# L1C: layers
s2_l1c_layers = {
    'ClassiPixelsMask_Band_00m_0_Tile1_Data': 'MSK_OPAQUE MSK_CIRRUS MSK_SNOICE, Cirrus cloud mask-Opaque cloud mask-Snow ice mask',
}
# L2A: layers
s2_l2a_layers = {
    'MSK_CLDPRB_20m': 'MSK_CLDPRB, Cloud Probabilities',
    'MSK_SNWPRB_20m': 'MSK_SNWPRB, Snow Probabilities',
    'IMG_DATA_Band_AOT_10m_Tile1_Data': 'AOT, Aerosol Optical Thickness',
    'IMG_DATA_Band_WVP_10m_Tile1_Data': 'WVP, Water Vapour',
    'IMG_DATA_Band_SCL_20m_Tile1_Data': 'SCL, Scene Classification'
}

# L2A: classification flags
s2_scene_classification_flags = {
    'NODATA': 0,
    'SATURATED_DEFECTIVE': 1,
    'DARK_FEATURE_SHADOW': 2,
    'CLOUD_SHADOW': 3,
    'VEGETATION': 4,
    'NOT_VEGETATED': 5,
    'WATER': 6,
    'UNCLASSIFIED': 7,
    'CLOUD_MEDIUM_PROBA': 8,
    'CLOUD_HIGH_PROBA': 9,
    'THIN_CIRRUS': 10,
    'SNOW_ICE': 11
}

