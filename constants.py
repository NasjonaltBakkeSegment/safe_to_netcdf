#!/usr/bin/python3

'''
Constants describing Sentinel satellites
'''

# ------------- All Sentinels -------------

global_attributes = {
    'project': 'Nasjonalt BakkeSegment',
    'project_short_name': 'NBS',
    'source': 'Space Borne Instrument',
    'license': 'The use of these data are covered by the Legal notice on the use of Copernicus'
               ' Sentinel Data and Service Information at https://sentinels.copernicus.eu/'
               'documents/247904/690755/Sentinel_Data_Legal_Notice',
    'acknowledgement': 'Copernicus EU',
    'dataset_production_status': 'In Work',
    'spatial_representation': 'grid',
    'processing_level': 'Operational',
    'creator_role': 'Data center contact, Metadata author',
    'creator_name': 'NBS Helpdesk, NBS team',
    'creator_email': 'nbs-helpdesk@met.no, nbs-helpdesk@met.no',
    'creator_institution': 'Norwegian Meteorological Institute, Norwegian Meteorological Institute',
    'institution': 'Norwegian Meteorological Institute',
    'institution_short_name': 'METNO',
    'publisher_url': 'https://met.no'
}

platform = {
    'S3A': {
        'platform_short_name': 'Sentinel-3A',
        'platform': 'Sentinel-3A',
        'platform_vocabulary': 'https://www.wmo-sat.info/oscar/satellites/view/sentinel_3a'
    },
    'S3B': {
        'platform_short_name': 'Sentinel-3B',
        'platform': 'Sentinel-3B',
        'platform_vocabulary': 'https://www.wmo-sat.info/oscar/satellites/view/sentinel_3b'
    }
}

instrument = {
    'OLCI': {
        'instrument_short_name': 'OLCI',
        'instrument': 'Ocean and Land Colour Imager',
        'instrument_vocabulary': 'https://www.wmo-sat.info/oscar/instruments/view/olci'
    }
}

# ------------- Sentinel 3 -------------

s3_olci = {
    'Conventions': 'CF-1.9 ACDD-1.3',
    'reference': 'https://sentinel.esa.int/web/sentinel/user-guides/sentinel-3-olci',
    'comment': 'This product contains selected information from Copernicus Sentinel-3 '
               'OLCI product. To obtain all data, please refer to the SAFE format available '
               'at colhub.met.no',
    'summary': 'The main objective of the SENTINEL-3 mission is to measure sea surface topography, '
               'sea and land surface temperature, and ocean and land surface colour with high '
               'accuracy and reliability to support ocean forecasting systems, environmental '
               'monitoring and climate monitoring.',
    'keywords_vocabulary':
        'GCMDSK:GCMD Science Keywords:https://gcmd.earthdata.nasa.gov/kms/concepts/concept_scheme/sciencekeywords,'
        'GEMET::http://inspire.ec.europa.eu/theme',
    'keywords': 'GCMDSK:Earth Science > Atmosphere > Atmospheric radiation > Reflectance,'
                'GEMET:Orthoimagery,GEMET:Land cover,GEMET:Oceanographic geographical features',
    'iso_topic_category': 'climatologyMeteorologyAtmosphere,imageryBaseMapsEarthCover,oceans,inlandWaters'
}

# ------------- Sentinel 2 -------------

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
