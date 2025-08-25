import rasterio
import geopandas as gpd
import numpy as np
from rasterio.mask import mask
from config import * #some main key parameters needed.
import pandas as pd
import pycountry

WORLD = gpd.read_file(r"data\ne_110m_admin_0_countries\ne_110m_admin_0_countries.shp").rename(columns={"ISO_A2":"Country"})
raster_file_wind = r"data\tifs\gwa3_250_capacityfactor_IEC1.tif"

def calculate_mean_capacity(geometry, raster_file):
    with rasterio.open(raster_file) as src:
        if WORLD.crs != src.crs:
            # Handle CRS mismatch
            raise ValueError(f"CRS mismatch between grid and {raster_file}")
            
        try:
            out_image, out_transform = mask(src, [geometry], crop=True)
            out_image = out_image[out_image != src.nodata]
            if out_image.size > 0:
                return out_image.mean()
            else:
                return np.nan
        except ValueError as e:
            if "Input shapes do not overlap raster" in str(e):
                return np.nan
            else:
                raise e

# Calculate mean capacity factor for each country
WORLD['mean_cf_onshore_wind'] = WORLD['geometry'].apply(lambda geom: calculate_mean_capacity(geom, raster_file_wind))

pv_file = r"data\solargis_pvpotential_countryranking_2020_data.xlsx"
data_pv = pd.read_excel(pv_file, sheet_name='Country indicators', skiprows=1, index_col='Country or region')
data_pv['mean_cf_pv'] = data_pv['Average practical potential \n(PVOUT Level 1, \nkWh/kWp/day), long-term']*365/8760 #from kWh/day to cf

final_df_cf = WORLD.merge(data_pv, how='left', right_on = 'ISO_A3', 
                          left_on = 'ISO_A3_EH')

# Function to get ISO_A3 code
def get_iso_a3(country_name):
    try:
        return pycountry.countries.lookup(country_name).alpha_3
    except LookupError:
        return None

# Apply the function to fill missing ISO_A3 codes
final_df_cf['ISO_A3_EH'] = final_df_cf.apply(lambda row: get_iso_a3(row['NAME_EN']) if row['ISO_A3_EH'] == '-99' else row['ISO_A3_EH'], axis=1)

final_df_cf[['NAME_EN', "ISO_A3_EH", 'mean_cf_onshore_wind', 'mean_cf_pv']].to_csv("data\mean_cfs.csv", index=False)
