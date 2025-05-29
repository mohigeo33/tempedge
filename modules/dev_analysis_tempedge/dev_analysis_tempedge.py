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
from skimage.filters import threshold_otsu
from skimage.filters import threshold_li
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

# initialise GEE
ee.Initialize()

# ====================================================================================================#
# USER INPUT PARAMETERS
# ====================================================================================================#
# study period: set your own years and months
year_start = 2014
year_end = 2023
month_start = 1
month_end = 12
max_cloud_cover = 100  # in %

# Study area: provide your own study area (in this case, Phnom Penh is given as an example)
phnom_gdf = ox.geocode_to_gdf('Phnom Penh, Cambodia')
phnom_geojson = phnom_gdf.iloc[0].geometry.__geo_interface__
phnom_aoi = ee.Geometry(phnom_geojson)
roi = phnom_aoi
epsg = 'EPSG:32648' # Phnom Penh

# Thresholds for plausibility check (NOTE: this we used for tropical region, if your study area falls into tropics do not change)
MIN_REASONABLE_LST = 10.0  # °C
MAX_REASONABLE_LST = 65.0  # °C

# Once parameters are set, the script will run fully automatically.
# No further user interaction is required.

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
  
# Function 6: Loading the csv/txt file
def load_csv(filepath):
    data =  []
    col = []
    checkcol = False
    with open(filepath) as f:
        for val in f.readlines():
            val = val.replace("\n","")
            val = val.split(',')
            if checkcol is False:
                col = val
                checkcol = True
            else:
                data.append(val)
    df = pd.DataFrame(data=data, columns=col)
    return df

# Function 7: Applying TempEdge thresholds
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

# export the dataframe to a CSV file
df.to_csv('dfallmonths.csv', index=False)

# ====================================================================================================#
# TempEdge THRESHOLDS DETERMINATION
# ====================================================================================================#
# --------------------------------------------------
# DATA PREPROCESSING
# --------------------------------------------------
myData = load_csv('dfallmonths.csv')
df = myData

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
# STATISTICAL ANALYSIS
# ====================================================================================================#

# --------------------------------------------------
# APPLYING DIFFERENT THRESHOLDS
# --------------------------------------------------
# Selecting min and max columns
min_lst_values = df['MinLST']
max_lst_values = df['MaxLST']
# Convert min and max LST to a NumPy array
min_lst_values_np = min_lst_values.to_numpy()
max_lst_values_np = max_lst_values.to_numpy()

# Otsu's thresholding_minLST (minLST)
threshold_value = threshold_otsu(min_lst_values_np)
# Otsu's thresholding_maxLST (upper spectrum)
threshold_value_max = -threshold_otsu(-max_lst_values_np)  # Use the NumPy array

# Improved Otsu_minLST
constrained_min_lst_values_array = min_lst_values_np[min_lst_values_np <= threshold_value]
refined_threshold_min = threshold_otsu(constrained_min_lst_values_array)
# Improved Otsu_maxLST (upper spectrum)
constrained_max_lst_values_array = max_lst_values_np[max_lst_values_np >= threshold_value_max]
refined_threshold_max = threshold_otsu(constrained_max_lst_values_array)

# Li's thresholding_minLST
threshold_value_min_li = threshold_li(min_lst_values_np)
# Li's thresholding_maxLST (upper spectrum)
threshold_value_max_li = -threshold_li(-max_lst_values_np)

# Quantile threshold_MinLST (10th)
min_lst_10th_quantile = min_lst_values.quantile(0.10)
# Quantile threshold_MaxLST (90th)
max_lst_90th_quantile = max_lst_values.quantile(0.90)

#TempEdge threshold_MinLST
overall_min_lst_threshold = TempEdge_minimum_threshold
#TempEdge threshold_MaxLST
overall_max_lst_threshold = TempEdge_maximum_threshold

# --------------------------------------------------
# COUNTING ARTEFACTS
# --------------------------------------------------
# Count artefacts for MinLST (both ends)
otsu_min_count = ((min_lst_values < threshold_value) | (min_lst_values > threshold_value_max)).sum()
improved_otsu_min_count = ((min_lst_values < refined_threshold_min) | (min_lst_values > refined_threshold_max)).sum()
li_min_count = ((min_lst_values < threshold_value_min_li) | (min_lst_values > threshold_value_max_li)).sum()
quantile_min_count = ((min_lst_values < min_lst_10th_quantile) | (min_lst_values > max_lst_90th_quantile)).sum()
tempedge_min_count = ((min_lst_values < overall_min_lst_threshold) | (min_lst_values > overall_max_lst_threshold)).sum()

# Count artefacts for MaxLST (both ends)
otsu_max_count = ((max_lst_values < threshold_value) | (max_lst_values > threshold_value_max)).sum()
improved_otsu_max_count = ((max_lst_values < refined_threshold_min) | (max_lst_values > refined_threshold_max)).sum()
li_max_count = ((max_lst_values < threshold_value_min_li) | (max_lst_values > threshold_value_max_li)).sum()
quantile_max_count = ((max_lst_values < min_lst_10th_quantile) | (max_lst_values > max_lst_90th_quantile)).sum()
tempedge_max_count = ((max_lst_values < overall_min_lst_threshold) | (max_lst_values > overall_max_lst_threshold)).sum()

# Print the results to check the counts
print(f"Otsu Method: MinLST artefacts = {otsu_min_count}, MaxLST artefacts = {otsu_max_count}")
print(f"Improved Otsu Method: MinLST artefacts = {improved_otsu_min_count}, MaxLST artefacts = {improved_otsu_max_count}")
print(f"Li Method: MinLST artefacts = {li_min_count}, MaxLST artefacts = {li_max_count}")
print(f"Quantile Method: MinLST artefacts = {quantile_min_count}, MaxLST artefacts = {quantile_max_count}")
print(f"TempEdge Method: MinLST artefacts = {tempedge_min_count}, MaxLST artefacts = {tempedge_max_count}")

# Define condition for TempEdge artefacts (both Min and Max)
tempedge_condition = (
    (df['MinLST'] < overall_min_lst_threshold) | (df['MinLST'] > overall_max_lst_threshold) |
    (df['MaxLST'] < overall_min_lst_threshold) | (df['MaxLST'] > overall_max_lst_threshold)
)

# Filter the DataFrame to get the artefact rows
tempedge_artefacts_df = df[tempedge_condition]

# --------------------------------------------------
# SHAPIRO-WILK TEST
# --------------------------------------------------
from scipy.stats import shapiro, kstest, norm
for col in ['MinLST', 'MeanLST', 'MaxLST']:
    stat, p = shapiro(df[col])
    print(f'Shapiro-Wilk Test for {col}: Statistics={stat:.3f}, p={p:.3f}')
    alpha = 0.05
    if p > alpha:
        print(f'{col}: Sample looks Gaussian (fail to reject H0)')
    else:
        print(f'{col}: Sample does not look Gaussian (reject H0)')

# --------------------------------------------------
# CORRELATION TEST
# --------------------------------------------------
# Correlation matrix
cols_for_corr = ['MinLST', 'MeanLST', 'MaxLST', 'Cloud', 'missing_pixels', 'DOY']
correlation_matrix = df[cols_for_corr].corr()

# --------------------------------------------------
# INSPECTION ON EXTREME COLD ARTEFACTS
# --------------------------------------------------
# Artefacts according to the TDQBT method
artefacts_minlst = df[~((df['MinLST'] >= TempEdge_minimum_threshold) & (df['MinLST'] <= TempEdge_maximum_threshold))]
artefacts_meanlst = df[~((df['MeanLST'] >= TempEdge_minimum_threshold) & (df['MeanLST'] <= TempEdge_maximum_threshold))]
artefacts_maxlst = df[~((df['MaxLST'] >= TempEdge_minimum_threshold) & (df['MaxLST'] <= TempEdge_maximum_threshold))]
artefacts_extreme = artefacts_minlst[((artefacts_minlst['MinLST'] == -123.14852014))]
# artefacts in both minLST and maxLST
artefacts = df[~((df['MinLST'] >= overall_min_lst_threshold) & (df['MaxLST'] <= overall_max_lst_threshold))]
artefacts_common = artefacts[((artefacts['MinLST'] < overall_min_lst_threshold) & (artefacts['MaxLST'] < overall_min_lst_threshold))]

# ----------------------------------------------------
# CHARTS USED IN THE MANUSCRIPT
# ----------------------------------------------------
# Figure 2. Distribution of the 243 selected Landsat scenes by acquisition year, cloud cover percentage (from scene-level metadata), and sensor type (Landsat 8 and 9).
alt.Chart(df).mark_point(size=50, filled=True).encode(
    x=alt.X('Timestamp:T', title='Year'),
    y=alt.Y('Cloud:Q', title='Cloud cover in %'),
    color=alt.Color(
        'satellite:N',
        scale=alt.Scale(scheme='plasma'),
        legend=alt.Legend(title="Sensor Type")),
    tooltip=[
        alt.Tooltip('Timestamp:T', title='Time'),
        alt.Tooltip('Cloud:Q', title='Cloud cover (in %)'),
        alt.Tooltip('satellite:N', title='Sensor Type')
    ]).properties(width=600, height=300)    
    
# Figure 5. Monthly minimum and maximum LST values in Phnom Penh.
label_font_size = 8  # Customize this as needed
least_missing_pixels_monthly = least_missing_pixels_monthly.loc[least_missing_pixels_monthly.groupby('Month')['missing_pixels'].idxmin()]
monthly_min_lst_thresholds = least_missing_pixels_monthly.groupby('Month')['MinLST'].min()
monthly_max_lst_thresholds = least_missing_pixels_monthly.groupby('Month')['MaxLST'].max()
monthly_mean_lst_thresholds = least_missing_pixels_monthly.groupby('Month')['MeanLST'].mean()

plt.figure(figsize=(12, 6))
horizontal_offset = 0.3  # Adjust this value if more or less spacing is needed

# Plot Monthly MinLST Threshold
plt.plot(monthly_min_lst_thresholds.index, monthly_min_lst_thresholds.values, marker='o', color='blue', label='Minimum LST')
for x, y in zip(monthly_min_lst_thresholds.index, monthly_min_lst_thresholds.values):
    plt.text(x - horizontal_offset, y, f"{y:.2f}", fontsize=label_font_size, ha='right', va='center')

# Plot Monthly MaxLST Threshold
plt.plot(monthly_max_lst_thresholds.index, monthly_max_lst_thresholds.values, marker='s', color='orange', label='Maximum LST')
for x, y in zip(monthly_max_lst_thresholds.index, monthly_max_lst_thresholds.values):
    plt.text(x - horizontal_offset, y, f"{y:.2f}", fontsize=label_font_size, ha='right', va='center')

# Labeling the plot
plt.title('')
plt.xlabel('Month')
plt.ylabel('LST (in °C)')
plt.xticks(monthly_min_lst_thresholds.index)  # Ensure month labels are shown
plt.grid(visible=True, linestyle='--', alpha=0.5)
plt.legend()
plt.show()

# Figure 6. Distribution of minimum and maximum LST values in Phnom Penh. 
lst_cols = {
    'MinLST': ('blue', 'Minimum LST'),
    'MaxLST': ('orange', 'Maximum LST')
}

artefacts_renamed = df.rename(columns={k: v[1] for k, v in lst_cols.items()})
artefacts_long = artefacts_renamed[['Minimum LST', 'Maximum LST']].melt(var_name='LST Type', value_name='LST Value')

plt.figure(figsize=(10, 6))

sns.swarmplot(
    x='LST Type',
    y='LST Value',
    data=artefacts_long,
    palette={v[1]: v[0] for k, v in lst_cols.items()},
    size=5
)
plt.axhline(TempEdge_minimum_threshold, color='purple', linestyle='--', linewidth=2, label='TempEdge Minimum Threshold')
plt.axhline(TempEdge_maximum_threshold, color='teal', linestyle='--', linewidth=2, label='TempEdge Maximum Threshold')
plt.title("", fontsize=14)
plt.ylabel("LST (°C)", fontsize=12)
plt.legend(loc='upper left')
plt.grid(True, axis='y', linestyle='--', alpha=0.7)
plt.tight_layout()
plt.show()


# Figure 7. Comparative visualisation of thresholds derived from different methods across the minimum (left) and maximum (right) LST distributions
fig, axes = plt.subplots(1, 2, figsize=(18, 8))
fig.suptitle('', fontsize=16)

# Histogram for MinLST with thresholds
axes[0].hist(min_lst_values, bins=30, color='lightgrey', edgecolor='black')
axes[0].axvline(threshold_value, color='red', linestyle='--', linewidth=2, label=f'Otsu Threshold: {threshold_value:.2f}')
axes[0].axvline(refined_threshold_min, color='blue', linestyle='--', linewidth=2, label=f'Improved Otsu Threshold: {refined_threshold_min:.2f}')
axes[0].axvline(threshold_value_min_li, color='black', linestyle='--', linewidth=2, label=f'Li Threshold: {threshold_value_min_li:.2f}')
axes[0].axvline(min_lst_10th_quantile, color='orange', linestyle='--', linewidth=2, label=f'10th Quantile: {min_lst_10th_quantile:.2f}')
axes[0].axvline(TempEdge_minimum_threshold, color='purple', linestyle='--', linewidth=2, label=f'TempEdge Threshold: {TempEdge_minimum_threshold:.2f}')
axes[0].set_title('Histogram of Minimum LST')
axes[0].set_xlabel('Minimum LST')
axes[0].set_ylabel('Frequency')
axes[0].legend()

# Histogram for MaxLST with thresholds
axes[1].hist(max_lst_values, bins=30, color='lightgrey', edgecolor='black')
axes[1].axvline(threshold_value_max, color='red', linestyle='--', linewidth=2, label=f'Otsu Threshold: {threshold_value_max:.2f}')
axes[1].axvline(refined_threshold_max, color='blue', linestyle='--', linewidth=2, label=f'Improved Otsu Threshold: {refined_threshold_max:.2f}')
axes[1].axvline(threshold_value_max_li, color='black', linestyle='--', linewidth=2, label=f'Li Threshold: {threshold_value_max_li:.2f}')
axes[1].axvline(max_lst_90th_quantile, color='orange', linestyle='--', linewidth=2, label=f'90th Quantile: {max_lst_90th_quantile:.2f}')
axes[1].axvline(TempEdge_maximum_threshold, color='teal', linestyle='--', linewidth=2, label=f'TempEdge Threshold: {TempEdge_maximum_threshold:.2f}')
axes[1].set_title('Histogram of Maximum LST')
axes[1].set_xlabel('Maximum LST')
axes[1].set_ylabel('Frequency')
axes[1].legend()

plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.show()

# Figure 8. Comparison of artefact detection across thresholding methods for minimum and maximum LST categories.
methods = ['TempEdge', 'Otsu', 'Improved Otsu', 'Li', 'Quantile']

min_counts = [tempedge_min_count, otsu_min_count, improved_otsu_min_count, li_min_count, quantile_min_count]
max_counts = [tempedge_max_count, otsu_max_count, improved_otsu_max_count, li_max_count, quantile_max_count]

x = np.arange(len(methods))
width = 0.35  # Width of the bars
fig, ax = plt.subplots(figsize=(10, 6))

bars1 = ax.bar(x - width/2, min_counts, width, label='Minimum LST', color='blue')
bars2 = ax.bar(x + width/2, max_counts, width, label='Maximum LST', color='orange')

ax.set_ylabel('Number of artefacts detected')
ax.set_xlabel('Thresholding method')
ax.set_title('')
ax.set_xticks(x)
ax.set_xticklabels(methods)
ax.legend()

for bar in bars1 + bars2:
    height = bar.get_height()
    ax.annotate(f'{height}',
                xy=(bar.get_x() + bar.get_width() / 2, height),
                xytext=(0, 3),  # Offset
                textcoords="offset points",
                ha='center', va='bottom')

plt.tight_layout()
plt.show()

# Figure 10. Distribution and normality assessment of minimum (left), mean (center), and maximum (right) LST. 
lst_cols = {
    'MinLST': ('blue', 'Minimum LST'),
    'MeanLST': ('green', 'Mean LST'),
    'MaxLST': ('orange', 'Maximum LST')
}

fig, axes = plt.subplots(2, 3, figsize=(18, 12))
fig.suptitle('', fontsize=16)

# Histograms with normal distribution curves
for i, (col, (color, label)) in enumerate(lst_cols.items()):
    ax = axes[0, i]

    sns.histplot(df[col], kde=False, color=color, ax=ax, stat="density")
    mean = df[col].mean()
    std = df[col].std()
    xmin, xmax = ax.get_xlim()
    x = np.linspace(xmin, xmax, 100)
    p = norm.pdf(x, mean, std)
    ax.plot(x, p, 'k', linewidth=2, label='Normal Distribution')

    ax.set_title(f'Distribution of {label}', fontsize=12)
    ax.set_xlabel(f'{label}', fontsize=10)
    if i == 0:
        ax.set_ylabel('Density', fontsize=10)
    ax.grid(True)
    ax.legend()

# Q-Q plots
for i, (col, (color, label)) in enumerate(lst_cols.items()):
    ax = axes[1, i]

    probplot(df[col], dist="norm", plot=ax)
    ax.get_lines()[0].set_color(color)  # Data points color
    ax.get_lines()[1].set_color('red')  # Reference line color
    ax.get_lines()[1].set_label('Ideal Normal Distribution')  # Add label for the red line
    ax.set_title(f'Q-Q Plot of {label}', fontsize=12)
    ax.set_xlabel('Theoretical Quantiles')
    ax.set_ylabel('Sample Quantiles')
    ax.legend()  # Add legend to the plot
    ax.grid(True)
plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.show()

# Figure 11. Pearson correlation matrix showing relationships among minimum, mean, and maximum LST, cloud cover (%), missing pixel (%), and DOY.
plt.figure(figsize=(8, 6))
sns.heatmap(correlation_matrix, annot=True, cmap='coolwarm', fmt=".2f")
plt.title('')
plt.show()

# Figure 12. Relationships among LST values (minimum, mean, maximum), cloud cover (%) percentage in the scene and missing pixel (%) within the study area.
fig, axes = plt.subplots(3, 2, figsize=(14, 18))
fig.suptitle('', fontsize=16)

lst_cols = {'MinLST': ('blue', 'Minimum LST'), 'MeanLST': ('green', 'Mean LST'), 'MaxLST': ('orange', 'Maximum LST')}

for i, lst_var in enumerate(['MinLST', 'MeanLST', 'MaxLST']):
    # Cloud vs. LST
    ax1 = axes[i, 0]
    ax1.scatter(df['Cloud'], df[lst_var], color=lst_cols[lst_var][0], s=10, label=lst_cols[lst_var][1])
    ax1.set_xlabel('Cloud (%)', fontsize=10)
    ax1.set_ylabel(lst_cols[lst_var][1], fontsize=10)
    ax1.axhline(y=TempEdge_minimum_threshold, color='purple', linestyle='--', label='TempEdge Minimum Threshold')
    ax1.axhline(y=TempEdge_maximum_threshold, color='teal', linestyle='--', label='TempEdge Maximum Threshold')
    ax1.legend(fontsize=8)
    ax1.grid(True, linestyle='--', alpha=0.7)
    ax1.set_title(f'Cloud vs. {lst_cols[lst_var][1]}', fontsize=12)

    # Missing Pixel Percentage vs. LST
    ax2 = axes[i, 1]
    ax2.scatter(df['missing_pixel_percentage'], df[lst_var], color=lst_cols[lst_var][0], s=10, label=lst_cols[lst_var][1])
    ax2.set_xlabel('Missing Pixel Percentage (%)', fontsize=10)
    ax2.set_ylabel(lst_cols[lst_var][1], fontsize=10)
    ax2.axhline(y=TempEdge_minimum_threshold, color='purple', linestyle='--', label='TempEdge Minimum Threshold')
    ax2.axhline(y=TempEdge_maximum_threshold, color='teal', linestyle='--', label='TempEdge Maximum Threshold')
    ax2.legend(fontsize=8)
    ax2.grid(True, linestyle='--', alpha=0.7)
    ax2.set_title(f'Missing pixel percentage vs. {lst_cols[lst_var][1]}', fontsize=12)

plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.show()

# Figure 13. Temporal distribution of extreme cold artefacts in minimum LST.
fig, ax1 = plt.subplots(figsize=(12, 6))

ax1.scatter(artefacts_extreme['Timestamp'], artefacts_extreme['MinLST'], color='blue', label='Minimum LST', marker='o', s=50)
ax1.scatter(artefacts_extreme['Timestamp'], artefacts_extreme['MeanLST'], color='green', label='Mean LST', marker='o', s=50)
ax1.scatter(artefacts_extreme['Timestamp'], artefacts_extreme['MaxLST'], color='orange', label='Maximum LST', marker='o', s=50)


ax1.set_xlabel('Timestamp')
ax1.set_ylabel('LST (°C)')
ax1.tick_params(axis='x', rotation=90)

ax2 = ax1.twinx()
ax2.plot(artefacts_extreme['Timestamp'], artefacts_extreme['Cloud'], color='grey', label='Cloud % in scene')
ax2.set_ylabel('Cloud (%)')

ax1.axhline(TempEdge_minimum_threshold, color='purple', linestyle='--', linewidth=2, label='TempEdge minimum threshold')
ax1.axhline(TempEdge_maximum_threshold, color='teal', linestyle='--', linewidth=2, label='TempEdge maximum threshold')

lines, labels = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
plt.legend(lines + lines2, labels + labels2, bbox_to_anchor=(1.05, 1), loc='upper left')

plt.tight_layout()
plt.show()

# Figure 14. Monthly normalised distributions of minimum LST artefact frequency, cloud cover, and missing pixel percentage in Phnom Penh.
artefacts_melted_monthly = artefacts_minlst.melt(
    id_vars=['Month'],
    value_vars=['MinLST', 'Cloud', 'missing_pixels'],
    var_name='LST_Type',
    value_name='LST_Value'
)

artefacts_melted_monthly['Frequency_Normalized'] = artefacts_melted_monthly.groupby(['Month', 'LST_Type'])['LST_Value'].transform(lambda x: x / x.sum())

custom_palette_labels = {
    'MinLST': ('blue', 'Minimum LST'),
    'Cloud': ('grey', 'Cloud cover'),
    'missing_pixels': ('purple', 'Missing pixels')
}

colors = {key: value[0] for key, value in custom_palette_labels.items()}
labels = {key: value[1] for key, value in custom_palette_labels.items()}

plt.figure(figsize=(12, 6))
sns.boxplot(
    x='Month',
    y='Frequency_Normalized',
    hue='LST_Type',
    data=artefacts_melted_monthly,
    dodge=True,
    palette=colors
)

plt.title('')
plt.xlabel('Month')
plt.ylabel('Normalized frequency')
plt.grid(visible=True, linestyle='--', alpha=0.5)

legend_handles = [plt.Rectangle((0, 0), 1, 1, color=colors[key], label=labels[key]) for key in labels]
plt.legend(title='Legend', handles=legend_handles, loc='upper right')

plt.tight_layout()
plt.show()

# Figure 15. Density contour plot showing the distribution of minimum LST values across the day of the year (DOY).
plt.figure(figsize=(12, 8))

MinLST = artefacts_minlst['MinLST']
DOY = pd.to_numeric(artefacts_minlst['DOY'], errors='coerce')

sns.kdeplot(x= DOY, y= MinLST, data= artefacts_minlst, levels=10, cmap="Blues", fill=True, label='Minimum LST')
plt.axhline(TempEdge_minimum_threshold, color='purple', linestyle='--', linewidth=2, label='TempEdge minimum threshold')
plt.axhline(TempEdge_maximum_threshold, color='teal', linestyle='--', linewidth=2, label='TempEdge maximum threshold')
plt.xlabel("Day of year")
plt.ylabel("Min LST (°C)")
plt.xlim(1, 365)  # Set DOY range to 1-365
plt.title("")
plt.legend()
plt.show()