#Standard vars
EI_VERSION = "3.10"
DB_NAME_INIT = 'ecoinvent-{}-cutoff'.format(EI_VERSION)
PROJECT_NAME = "iea_h2_lca"
NAME_REF_DB = "ecoinvent_{}_reference".format(EI_VERSION).replace(".","")
USER_NAME = "terlouw_t"
BIOSPHERE_DB = 'ecoinvent-3.10-biosphere'
SEL_COLS = ['Technology', 'Technology_details', 'Technology_electricity', 'Technology_electricity_details']
REGIONALIZE=True # if true, we regionalize hydrogen application activities
START_YEAR = 2025

# files
FILEPATH_WIND = r"data\tifs\gwa3_250_capacityfactor_IEC1.tif"
FILEPATH_PV = r"data\tifs\PVOUT.tif"  # Assuming same filepath for both
FILE_NAME_COSTS = r"data\cost_data_init.xlsx"
FILE_DF_RATIOS = r"data\results_curve_fits_wind_hydrogen.xlsx"

# LCIA methods
MY_METHODS = [
    ('EF v3.1 EN15804', 'climate change', 'global warming potential (GWP100) including H2'),
    ##('Own method', 'Water consumption', 'Water consumption'),    
]