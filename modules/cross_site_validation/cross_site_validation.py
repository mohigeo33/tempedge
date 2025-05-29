'''
"""
TempEdge: Artefact Detection and Filtering in Landsat Analysis-Ready Surface Temperature Data
---------------------------------------------------------------------------------------------
Author: Gulam Mohiuddin (2025)

Description:
GEE Implementation of TempEdge – a method developed for artefact detection and filtering in
Landsat Collection 2 Level-2 (LC2L2) Surface Temperature datasets, especially suited for tropical
regions.

This script:
- Loads LC2L2 ST data from Google Earth Engine
- Applies cloud masking (CFMask)
- Identifies high-quality monthly LST observations
- Derives TempEdge min/max thresholds
- Apply TempEdge on the LC2L2 ST data
- Delivers artefacts filtered data for further use


"""

  User input parameters:
year_start, year_end       [INT] Start and end year of the analysis period.
                           Example: year_start = 2014, year_end = 2015

month_start, month_end     [INT] Start and end month of the analysis period (1–12).
                           Example: month_start = 1 (January), month_end = 12 (December)

max_cloud_cover            [FLOAT] Maximum allowable scene cloud cover (%) for image selection.
                           Example: max_cloud_cover = 100 includes all available scenes regardless of cloud cover.

select_roi                 [ee.Geometry] Region of interest (ROI) defined as an Earth Engine geometry.
                           Can be created from any shapefile or geocoded location.
                           Example: 'Phnom Penh, Cambodia' via osmnx → GeoDataFrame → GeoJSON → ee.Geometry
"""
'''

# loading necessary libraries
import os
import sys
import ee
import geemap
import osmnx as ox
import geopandas as gpd
import pandas as pd
from tqdm import tqdm
import math
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import statistics
import seaborn as sns
import altair as alt
import altair_transform
from scipy.stats import shapiro
from scipy.stats import lognorm
from scipy.stats import kstest
from scipy.stats import norm
from sklearn.metrics import r2_score
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

# initialise GEE
ee.Initialize()

# ====================================================================================================#
# USER INPUT PARAMETERS
# ====================================================================================================#
# study period: set your own years and months
year_start = 2014
year_end = 2015
month_start = 1
month_end = 12
max_cloud_cover = 100  # in %

# Thresholds for plausibility check (NOTE: this we used for tropical region, if your study area falls into tropics do not change)
MIN_REASONABLE_LST = 10.0  # °C
MAX_REASONABLE_LST = 65.0  # °C

# Study area: provide your own study area (in this case, Phnom Penh is given as an example)
epsg = 'EPSG:4326' # we used this for global use

# Study areas
# Lagos, Nigeria
print("Fetching Lagos...")
lagos_gdf = ox.geocode_to_gdf('Lagos, Nigeria')
lagos_geojson = lagos_gdf.iloc[0].geometry.__geo_interface__
lagos_aoi = ee.Geometry(lagos_geojson)

# ------------------------------
# Mérida, Mexico
print("Fetching Mérida...")
merida_gdf = ox.geocode_to_gdf('Mérida, Mexico')
merida_geojson = merida_gdf.iloc[0].geometry.__geo_interface__
merida_aoi = ee.Geometry(merida_geojson)

# ------------------------------
# Kuala Lumpur, Malaysia
print("Fetching Kuala Lumpur...")
kuala_gdf = ox.geocode_to_gdf('Kuala Lumpur, Malaysia')
kuala_geojson = kuala_gdf.iloc[0].geometry.__geo_interface__
kuala_aoi = ee.Geometry(kuala_geojson)

# ------------------------------
# Presidente Prudente, Brazil
print("Fetching Presidente Prudente...")
prudente_gdf = ox.geocode_to_gdf('Presidente Prudente, Brazil')
prudente_geojson = prudente_gdf.iloc[0].geometry.__geo_interface__
prudente_aoi = ee.Geometry(prudente_geojson)

# You have keep chaning the AOI here to get statistics of different cities
# === SELECT CITY TO PROCESS ===
selected_city = 'prudente'  # Change this to 'lagos', 'merida', 'kuala', etc.

city_aoi_map = {
    'lagos': lagos_aoi,
    'merida': merida_aoi,
    'kuala': kuala_aoi,
    'prudente': prudente_aoi
}

roi = city_aoi_map[selected_city]

# ====================================================================================================#
# FUNCTIONS
# ====================================================================================================#
# Function 1: Cloud masking using CFMask
def function_mask_clouds_l89(img):
    cloud_mask = img.select('QA_PIXEL').bitwiseAnd(1 << 3).eq(0)
    return img.updateMask(cloud_mask)

# Function 2: Renaming bands for coherence and and applying scale factors and offset (NOTE: this is applicable for only Landsat 8 and 9, you have to create separate function for other Landsat sensors.)
def function_bands_l89(img):
       bands = ['SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7']
       thermal_band = ['ST_B10']
       new_bands = ['B', 'G', 'R', 'NIR', 'SWIR1', 'SWIR2']
       new_thermal_bands = ['LST']
       vnirswir = img.select(bands).multiply(0.0000275).add(-0.2).rename(new_bands)
       tir = img.select(thermal_band).multiply(0.00341802).add(149).subtract(273.15).rename(new_thermal_bands) # converting it from K to celcius by subtracting 273.15
       return vnirswir.addBands(tir).copyProperties(img, img.propertyNames())

# Function 3: clipping the image collections to study area
def fun_clip(img):
  clip_img = img.clip(roi)
  return clip_img

# Function 4: Extracting the tabular information from the processed images
def fctn_get_image_stats(img):
    img_lst = img.select('LST')

    # Calculate basic statistics
    img_lst_mean_value = img_lst.reduceRegion(ee.Reducer.mean(), roi, 30, crs= epsg, bestEffort=True, maxPixels=1e9).getInfo().get('LST')
    img_lst_max_value = img_lst.reduceRegion(ee.Reducer.max(), roi, 30, crs= epsg, bestEffort=True, maxPixels=1e9).getInfo().get('LST')
    img_lst_min_value = img_lst.reduceRegion(ee.Reducer.min(), roi, 30, crs= epsg, bestEffort=True, maxPixels=1e9).getInfo().get('LST')
    img_date = img.date().getInfo().get('value')
    img_systemindex = img.get('system:index').getInfo()
    img_cloud_cover = img.get('CLOUD_COVER').getInfo()
    img_spacecraft = img.get('SPACECRAFT_ID').getInfo()

    # Calculate missing pixels
    missing_pixels = img_lst.unmask(-9999).lt(-9998).reduceRegion(
        reducer=ee.Reducer.sum(),
        geometry=roi,
        scale=30,
        maxPixels=1e13
    ).getInfo()
    missing_pixel_count = int(missing_pixels.get('LST', 0))


    # Create dictionary with all information
    img_all_info = {
        'system:index': img_systemindex,
        'date': img_date,
        'mean_lst': img_lst_mean_value,
        'min_lst': img_lst_min_value,
        'max_lst': img_lst_max_value,
        'cloud_cover': img_cloud_cover,
        'satellite': img_spacecraft,
        'missing_pixels': missing_pixel_count
    }

    return img_all_info

# Function 5: Add date variables to DataFrame.
def add_date_info(df):
  df['Timestamp'] = pd.to_datetime(df['Millis'], unit='ms')
  df['Year'] = pd.DatetimeIndex(df['Timestamp']).year
  df['Month'] = pd.DatetimeIndex(df['Timestamp']).month
  df['Day'] = pd.DatetimeIndex(df['Timestamp']).day
  df['DOY'] = pd.DatetimeIndex(df['Timestamp']).dayofyear
  return df
  
# Function 6: Applying TempEdge thresholds
def apply_temp_threshold(image):
    min_threshold = TempEdge_minimum_threshold
    max_threshold = TempEdge_maximum_threshold

    mask = image.select('LST') \
                .gte(min_threshold) \
                .And(image.select('LST').lte(max_threshold))

    masked_image = image.updateMask(mask)
    return masked_image
# ====================================================================================================#
# CALLING IMAGE COLLECTIONS
# ====================================================================================================#
# Landsat 8
imgCol_L8_SR = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')\
    .filterBounds(roi)\
    .filter(ee.Filter.calendarRange(year_start,year_end,'year'))\
    .filter(ee.Filter.calendarRange(month_start,month_end,'month'))\
    .filter(ee.Filter.lt('CLOUD_COVER', max_cloud_cover))\
    .map(function_mask_clouds_l89)

# Landsat 9
imgCol_L9_SR = ee.ImageCollection('LANDSAT/LC09/C02/T1_L2')\
    .filterBounds(roi)\
    .filter(ee.Filter.calendarRange(year_start,year_end,'year'))\
    .filter(ee.Filter.calendarRange(month_start,month_end,'month'))\
    .filter(ee.Filter.lt('CLOUD_COVER', max_cloud_cover))\
    .map(function_mask_clouds_l89)

# ====================================================================================================#
# EXECUTIONS
# ====================================================================================================#
# Applying the band coherence function and scale factors
imgCol_merge_L8_band_coh = imgCol_L8_SR.map(function_bands_l89)
imgCol_merge_L9_band_coh = imgCol_L9_SR.map(function_bands_l89)

# Applying clip function
imgCol_merge_L8_band_coh_clip = imgCol_merge_L8_band_coh.map(fun_clip)
imgCol_merge_L9_band_coh_clip = imgCol_merge_L9_band_coh.map(fun_clip)

# Merging the collection
imgCol_merge = imgCol_merge_L8_band_coh_clip.merge(imgCol_merge_L9_band_coh_clip)

# ====================================================================================================#
# DATA EXTRACTIONS
# ====================================================================================================#
doi = imgCol_merge
doiList = doi.toList(doi.size())
doiList_size = doiList.size().getInfo()
print('Total Images in Data of Interest (doi) dataset: ', doiList_size)

df = pd.DataFrame(columns=['SystemIndex', 'Millis', 'MeanLST', 'MaxLST', 'MinLST', 'Cloud', 'satellite', 'missing_pixels'])

for i in tqdm(range(doiList_size)):
    image = ee.Image(doiList.get(i))
    image_info = fctn_get_image_stats(image)

    # Create a temporary DataFrame for the new row
    temp_df = pd.DataFrame([{
        'SystemIndex': image_info['system:index'],
        'Millis': image_info['date'],
        'MeanLST': image_info['mean_lst'],
        'MaxLST': image_info['max_lst'],
        'MinLST': image_info['min_lst'],
        'Cloud': image_info['cloud_cover'],
        'satellite': image_info['satellite'],
        'missing_pixels': image_info['missing_pixels']
    }])
    df = pd.concat([df, temp_df], ignore_index=True)

df

df = add_date_info(df)
# ====================================================================================================#
# TempEdge THRESHOLDS DETERMINATION
# ====================================================================================================#
# --------------------------------------------------
# DATA PREPROCESSING
# --------------------------------------------------
# Removing rows with null values in MinLST, MeanLST, MaxLST
df.replace(['', ' ', '-1'], np.nan, inplace = True)
df = df.dropna(subset=['MinLST', 'MeanLST', 'MaxLST'])

# Converting the variables into numeric
#Mean
df["MeanLST"] = pd.to_numeric(df["MeanLST"])
meanlst = df["MeanLST"]

#Max
df["MaxLST"] = pd.to_numeric(df["MaxLST"])
maxlst = df["MaxLST"]

#Min
df["MinLST"] = pd.to_numeric(df["MinLST"])
minlst = df["MinLST"]

#Time and cloud
df['Year'] = pd.to_numeric(df['Year'])
df['Month'] = pd.to_numeric(df['Month'])
df['Cloud'] = pd.to_numeric(df['Cloud'])
df['missing_pixels'] = pd.to_numeric(df['missing_pixels'])
df['DOY'] = pd.to_numeric(df['DOY'])

# Compute area in m²
area_m2 = roi.area().getInfo()

# Calculate the missing pixel percentage and add it as a new column
df['missing_pixel_percentage'] = (df['missing_pixels'] / area_m2) * 100

# --------------------------------------------------
# PRE-TempEdge MINIMUM & MAXIMUM LST
# --------------------------------------------------
# MaxLST
max_lst_before_tempedge = df['MaxLST'].max()
# MinLST
min_lst_before_tempedge = df['MinLST'].min()

# Print the results with two decimal places
print(f"Before TempEdge Minimum LST: {min_lst_before_tempedge:.2f}°C")
print(f"Before TempEdge Maximum LST: {max_lst_before_tempedge:.2f}°C")

# --------------------------------------------------
# TempEdge: STEP 1
# --------------------------------------------------
# general plausibility test
# Initialize list to store valid monthly images
valid_monthly_images = []

# Group the dataset by month
grouped = df.groupby('Month')

for month, group in grouped:
    # Filter group to exclude images with unrealistic MinLST or MaxLST
    valid_group = group[
        (group['MinLST'] > MIN_REASONABLE_LST) &
        (group['MaxLST'] < MAX_REASONABLE_LST)
    ]

    if not valid_group.empty:
        # Pick image with the least missing pixels among the valid ones
        best_image = valid_group.sort_values(by='missing_pixels').iloc[0]
        valid_monthly_images.append(best_image)
    else:
        print(f"⚠️ No valid image found for Month {month} (MinLST or MaxLST outside plausible range).")

# Convert to DataFrame
least_missing_pixels_monthly = pd.DataFrame(valid_monthly_images)

# Final check: if no valid months, print a general warning
if least_missing_pixels_monthly.empty:
    print("⚠️ No valid monthly images passed the plausibility test.")
    print("   Please try increasing the study period.\n")

# Now, for each month, get the lowest MinLST and highest MaxLST from the selected rows
monthly_min_lst_thresholds = least_missing_pixels_monthly.groupby('Month')['MinLST'].min()
monthly_max_lst_thresholds = least_missing_pixels_monthly.groupby('Month')['MaxLST'].max()

# --------------------------------------------------
# TempEdge: STEP 2
# --------------------------------------------------
# Find the absolute minimum value and the corresponding month
TempEdge_minimum_threshold = monthly_min_lst_thresholds.min()
TempEdge_minimum_month = monthly_min_lst_thresholds.idxmin()

# Find the absolute maximum value and the corresponding month
TempEdge_maximum_threshold = monthly_max_lst_thresholds.max()
TempEdge_maximum_month = monthly_max_lst_thresholds.idxmax()

print(f"\nTempEdge Minimum Threshold: {TempEdge_minimum_threshold:.2f}°C (Month {TempEdge_minimum_month})")
print(f"TempEdge Maximum Threshold: {TempEdge_maximum_threshold:.2f}°C (Month {TempEdge_maximum_month})")

# ====================================================================================================#
# CHART FOR THE CROSS-SITE VALIDATION
# ====================================================================================================#
# Figure 9. Comparison of raw (pre-TempEdge) and filtered (TempEdge) minimum and maximum LST values across five tropical cities.

# values mannually added
# Phnom Penh, Cambodia
raw_minimum_phnom_penh = -123.15
raw_maximum_Phnom_penh = 88.84
tempedge_minimum_phnom_penh = 12.89
tempedge_maximum_phnom_penh = 64.96

# Lagos, Nigeria
raw_minimum_lagos = -123.15
raw_maximum_lagos = 73.22
tempedge_minimum_lagos = 11.02
tempedge_maximum_lagos = 60.42

# Merida, Mexico
raw_minimum_merida = -123.15
raw_maximum_merida = 66.70
tempedge_minimum_merida = 12.33
tempedge_maximum_merida = 60.33

# Kuala Lumpur, Malayasia
raw_minimum_kuala = -123.15
raw_maximum_kuala = 66.20
tempedge_minimum_kuala = 14.54
tempedge_maximum_kuala = 63.91

# Presidente Prudente, Brazil
raw_minimum_prudente = -123.15
raw_maximum_prudente = 63.81
tempedge_minimum_prudente = 13.34
tempedge_maximum_prudente = 63.81

import matplotlib.pyplot as plt
import numpy as np

cities = ['Phnom Penh', 'Lagos', 'Merida', 'Kuala Lumpur', 'Pres. Prudente']

raw_min = [-123.15] * 5
raw_max = [88.84, 73.22, 66.70, 66.20, 63.81]
temp_min = [12.89, 11.02, 12.33, 14.54, 13.34]
temp_max = [64.96, 60.42, 60.33, 63.91, 63.81]

x = np.arange(len(cities))
width = 0.2

fig, ax = plt.subplots(figsize=(12, 6))

rects1 = ax.bar(x - 1.5*width, raw_min, width, label='Minimum LST (pre-TempEdge)', color='#ADD8E6')   # Light blue
rects2 = ax.bar(x - 0.5*width, raw_max, width, label='Maximum LST (pre-TempEdge)', color='#FFD580')   # Light orange
rects3 = ax.bar(x + 0.5*width, temp_min, width, label='Minimum LST (TempEdge)', color='#1f77b4')  # Blue
rects4 = ax.bar(x + 1.5*width, temp_max, width, label='Maximum LST (TempEdge)', color='#ff7f0e')  # Orange

ax.set_ylabel('LST (in °C)')
ax.set_xlabel('City')
ax.set_title('')
ax.set_xticks(x)
ax.set_xticklabels(cities)
ax.grid(True, axis='y', linestyle='--', alpha=0.7)

ax.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))

def add_labels(rects):
    for rect in rects:
        height = rect.get_height()
        ax.annotate(f'{height}',
                    xy=(rect.get_x() + rect.get_width()/2, height),
                    xytext=(0, 3),
                    textcoords="offset points",
                    ha='center', va='bottom')

for r in [rects1, rects2, rects3, rects4]:
    add_labels(r)

plt.tight_layout()
plt.show()