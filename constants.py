#!/usr/bin/python3

"""
Constants describing Sentinel satellites
"""

# ------------- All Sentinels -------------

platform_id = {"Sentinel-2A": 0, "Sentinel-2B": 1, "Sentinel-2C": 2, "Sentinel-2D": 3}

# ------------- Sentinel 2 -------------

# Bands order
s2_bands_order = {0: 'B1', 1: 'B2', 2: 'B3', 3: 'B4', 4: 'B5', 5: 'B6', 6: 'B7', 7: 'B8', 8: 'B8A',
                  9: 'B9', 10: 'B10', 11: 'B11', 12: 'B12'}

# Bands aliases
s2_bands_aliases = {'B01': 'B1', 'B02': 'B2', 'B03': 'B3', 'B04': 'B4', 'B05': 'B5', 'B06': 'B6',
                    'B07': 'B7', 'B08': 'B8', 'B8A': 'B8A', 'B09': 'B9', 'B10': 'B10', 'B11': 'B11',
                    'B12': 'B12'}

# L2A: layers
s2_l2a_layers = {"MSK_CLDPRB_20m": "MSK_CLDPRB, Cloud Probabilities",
              "MSK_SNWPRB_20m": "MSK_SNWPRB, Snow Probabilities",
              'IMG_DATA_Band_AOT_10m_Tile1_Data': "AOT, Aerosol Optical Thickness",
              'IMG_DATA_Band_WVP_10m_Tile1_Data': "WVP, Water Vapour",
              'IMG_DATA_Band_SCL_20m_Tile1_Data': "SCL, Scene Classification"}

# L2A: classification flags
s2_scene_classification_flags = {'NODATA': 0,
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
                             'SNOW_ICE': 11}

