import brightway2 as bw
import uuid
from config import *
import premise as pr
import pandas as pd
import pycountry
import warnings
import numpy as np 
import copy
import re
import io
from contextlib import redirect_stdout
import brightway2 as bw
from bw2calc import MultiLCA

from rapidfuzz import fuzz, process
import warnings
import rioxarray as ri

from functools import partial
from bw2io.importers.base_lci import LCIImporter
from bw2io.strategies import add_database_name, csv_restore_tuples

from private_keys import KEY_PREMISE

# Standard VARS
EN_H2 = 120 #MJ/kg h2
AMOUNT_WATER_ELECTROLYSIS = 15+9 # Tonelli et al., 2023, and Terlouw et al., 2024
HOURS_YR = 8760
COST_DATA = pd.read_excel(FILE_NAME_COSTS, index_col = [0,1], usecols=[0,1,2,3,4,5,6,7,8])
DF_RATIOS = pd.read_excel(FILE_DF_RATIOS,  index_col=[0,1])

EC_DIESEL =	43 #MJ/kg LHV, https://www.engineeringtoolbox.com/alternative-fuels-d_1221.html
EC_BIODIESEL = 37 #MJ/kg LHV, https://www.engineeringtoolbox.com/alternative-fuels-d_1221.html
EC_ETHANOL = 27 #MJ/kg LHV, https://www.engineeringtoolbox.com/alternative-fuels-d_1221.html

# The following function are essential to make these activities for hydrogen facilities planned from the IEA database
# calculated compress power req based on from: https://transitionaccelerator.ca/wp-content/uploads/2023/04/TA-Technical-Brief-1.1_TEEA-Hydrogen-Compression_PUBLISHED.pdf
COMPRESS_H2_ELECTR = {30:0, #is already 30 bar
            20:0.20247384858513207, #from 20 to 30 bar
           1:1.768216140565714} # from 1 to 30 bar
           # dictionary to provide extra power for compressing hydrogen after water electrolysis

ELECTROLYZER_DATA = {
    'pem': {
        'pressure_bar': 30,
        'compression_kWh_per_kgH2': COMPRESS_H2_ELECTR[30]
    },
    'aec': {
        'pressure_bar': 20,
        'compression_kWh_per_kgH2': COMPRESS_H2_ELECTR[20]
    },
    'soec': {
        'pressure_bar': 1,
        'compression_kWh_per_kgH2': COMPRESS_H2_ELECTR[1]
    }
}

def convert_iso3_to_iso2(iso3):
    try:
        country = pycountry.countries.get(alpha_3=iso3)
        return country.alpha_2
    except KeyError:
        return None
    
def match_year_to_database(year, iam_model='remind_SSP2-PkBudg1150', ref_db = NAME_REF_DB):
    """
    Matches a given year to a database based on defined year ranges.
    
    Returns
    -------
    str
        The matched database name.
    """
    database_mapping = [
        (1500, 2027, ref_db),
        (2028, 2032, f"ecoinvent_{iam_model}_2030_base"),
        (2033, 2037, f"ecoinvent_{iam_model}_2035_base"),
        (2038, 2100, f"ecoinvent_{iam_model}_2040_base"),
    ]

    for start, end, db in database_mapping:
        if start <= year <= end:
            return db

    raise ValueError(f"No database match found for year {year}")

def get_conversion_factor(baseline_act):
    # this is in cubic meter, but we are interesting in kilogram for a fair comparison, so we have to multiply.
    if baseline_act in [
        ('market for natural gas, high pressure incl. combustion', 'natural gas, high pressure', 'RoW'),
        ('market for biomethane, high pressure incl. combustion', 'biomethane, high pressure', 'RoW')
    ]:
        return 1 / 0.735
    # consinder different energy density.
    elif baseline_act == ('diesel production, petroleum refinery operation incl. combustion', 'diesel', 'RoW'):
        return EC_BIODIESEL / EC_DIESEL
    elif baseline_act == ('ethanol production, via fermentation, from corn, economic allocation incl. combustion', 'ethanol, from corn', 'US'):
        return EC_BIODIESEL / EC_ETHANOL
    else:
        return 1

def get_year_from_name(name):
    match = re.search(r'\((\d{4})\)', name)
    return int(match.group(1)) if match else None

def extract_h2_info(comment):
    amount_match = re.search(r'\b(\d+)\b', comment)
    source_match = re.search(r'Hydrogen production source:\s*(.*)', comment)
    amount = float(amount_match.group(1)) if amount_match else None
    source = tuple(source_match.group(1).strip().split(':')) if source_match else None
    return amount, source

def drop_empty_categories_2(db):
    """Drop categories with the value ``('',)`` or None."""
    DROP = ('',)
    for ds in db:
        if ds.get('categories') in [DROP, None]:
            ds.pop('categories', None)  # Remove 'categories' if it exists
            #del ds['categories']
        for exc in ds.get("exchanges", []):
            if exc.get('categories') in [DROP, None]:
                #del exc['categories']
                exc.pop('categories', None)  # Remove 'categories' if it exists
    return db


def process_import(db_name, new_activities, iam='remind'):
    """
    Process and write a database with predefined strategies, match databases, 
    and generate statistics.
    
    Parameters:
        db_name (str): The name of the database to write.
        new_activities (list): The dataset containing new activities.
    
    Returns:
        None
    """
    importer = LCIImporter(db_name)
    importer.data = new_activities
    importer.strategies = [
        partial(add_database_name, name=db_name),
        csv_restore_tuples,
        drop_empty_categories,
        drop_empty_categories_2,
        strip_nonsense,
    ]

    importer.apply_strategies()

    # Match databases
    for db in bw.databases:
        if f'ecoinvent_{iam}' in str(db) or BIOSPHERE_DB in str(db) or NAME_REF_DB in str(db):
            importer.match_database(db)

    # Generate statistics and export results
    importer.statistics()
    importer.write_excel(only_unlinked=True)
    importer.write_database()
    
def round_to_nearest(number, step=0.05, max_value=None):
    """
    Rounds a given float to the nearest specified step and optionally caps the result at a maximum value.

    Parameters:
    number (float): The number to be rounded.
    step (float): The step to which the number should be rounded. Default is 0.05.
    max_value (float or None): The maximum value to cap the result. Default is None (no cap).

    Returns:
    float: The number rounded to the nearest specified step and optionally capped at max_value.
    
    Example:
    >>> round_to_nearest(3.67)
    3.65
    >>> round_to_nearest(3.67, 0.1)
    3.7
    >>> round_to_nearest(0.75, max_value=0.7)
    0.7
    >>> round_to_nearest(0.75, max_value=None)
    0.75
    """
    result = round(number / step) * step
    if max_value is not None:
        return min(result, max_value)
    return result

def extract_pressure(input_string):
    # Split the string by commas
    parts = input_string.split(', ')
    
    # Loop through the parts and look for a numeric value followed by 'bar'
    for part in parts:
        if 'bar' in part:
            try:
                # Extract the integer value before 'bar'
                pressure_value = int(part.split()[0])
                return pressure_value
            except ValueError:
                # Skip parts that don't start with a valid integer
                continue
    
    # Return None if no valid pressure value found
    return None

def drop_empty_categories(db):
    """Drop categories with the value ``('',)``"""
    DROP = ('',)
    for ds in db:
        if ds.get('categories') == DROP:
            del ds['categories']
        for exc in ds.get("exchanges", []):
            if exc.get('categories') == DROP:
                del exc['categories']
    return db

def strip_nonsense(db):
    for ds in db:
        for key, value in ds.items():
            if isinstance(value, str):
                ds[key] = value.strip()
            for exc in ds.get('exchanges', []):
                for key, value in exc.items():
                    if isinstance(value, str):
                        exc[key] = value.strip()
    return db

def multiply_amounts(exchanges, factor):
    """
    Multiply the 'amount' field of each exchange in the list by a specified factor.

    Parameters:
    exchanges (list): A list of dictionaries, where each dictionary represents an exchange.
                      Each dictionary should contain an 'amount' key.
    factor (float): The factor by which to multiply the 'amount' values.

    Returns:
    list: The updated list of exchanges with 'amount' values multiplied by the specified factor.
    """
    for exchange in exchanges:
        if 'amount' in exchange:
            exchange['amount'] *= factor
    return exchanges

def merge_duplicate_exchanges(exchanges):
    """
    Merge duplicate exchanges by summing their 'amount' values.

    Exchanges are grouped and aggregated based on distinguishing keys:
        - For 'technosphere' exchanges: ('name', 'reference product', 'location', 'unit')
        - For 'biosphere' exchanges: ('name', 'categories', 'unit')

    Parameters:
        exchanges (list of dict): A list of exchange dictionaries, each containing
            at least 'type', 'amount', and type-specific identifying keys.

    Returns:
        list of dict: A list of merged exchange dictionaries, with duplicates summed
        and entries with zero 'amount' removed.

    Raises:
        KeyError: If any exchange is missing required keys.
    """
    merged = {}
    key_fields_map = {
        'technosphere': ['name', 'reference product', 'location', 'unit'],
        'biosphere': ['name', 'categories', 'unit'],
    }

    for i, ex in enumerate(exchanges):
        ex_type = ex.get('type')
        if ex_type not in key_fields_map:
            continue  # Skip production type

        required_keys = key_fields_map[ex_type]
        missing = [k for k in required_keys if k not in ex]
        if missing:
            print(f"Exchange #{i} is missing keys: {missing}")
            print("Full exchange:", ex)
            raise KeyError(f"Missing required keys: {missing}")

        key = tuple(ex[k] for k in required_keys)
        
        if key in merged:
            merged[key]['amount'] += ex['amount']
        else:
            merged[key] = copy.deepcopy(ex)

    # Filter out zero-amount exchanges
    return [ex for ex in merged.values() if ex['amount'] != 0]

def cf_valid(cf):
    """
    Check if a given capacity factor (cf) is valid.

    Parameters:
    cf (float or None): The capacity factor to validate.

    Returns:
    bool: True if cf is not None and greater than 0, otherwise False.
    """
    return cf is not None and cf > 0
    
def curve_fit_g(x, db, info_type):
    """
    Calculate the ratio based on g curve fitting.

    Args:
        x (float): share of x.
        db (str): Database identifier.
        info_type (str): Type of information to retrieve from DF_RATIOS.

    Returns:
        float: ratio.

    """
    a = DF_RATIOS.loc[(info_type,db)]['a']
    b = DF_RATIOS.loc[(info_type,db)]['b']
    c = DF_RATIOS.loc[(info_type,db)]['c']

    # calculates wind_share from predefined curve fit based on cf
    ratio = b/a*x**c+1/a
    ratio = 0 if ratio<0 else ratio
    return round(ratio,3)  

def curve_fit_f(x, db, info_type):
    """
    Calculate the wind share based on curve fitting from optimization.

    Args:
        x (float): share of x.
        db (str): Database identifier.
        info_type (str): Type of information to retrieve from DF_RATIOS.

    Returns:
        float: ratio share.

    """

    a = DF_RATIOS.loc[(info_type,db)]['a']
    b = DF_RATIOS.loc[(info_type,db)]['b']
    c = DF_RATIOS.loc[(info_type,db)]['c']

    # calculates ratio from predefined curve fit based on cf
    ratio = 1 / (1 + np.exp(-(x - b) * a)) ** c if x>0 else 0

    return round(ratio,3)

def create_database(db_name):
    """
    Creates and registers a new database in the Brightway2 environment.

    Parameters:
    -----------
    db_name : str
        The name of the database to be created.

    Returns:
    --------
    None
        The function does not return a value but will print a message indicating the result of the operation.

    Description:
    ------------
    This function checks if a database with the given name (`db_name`) already exists in the Brightway2 environment. 
    If the database already exists, the function aborts the operation and prints a message indicating that the database 
    already exists. If the database does not exist, the function proceeds to create and register a new database with 
    the specified name.

    Example:
    --------
    >>> create_database("new_db")
    Database 'new_db' has been created and registered.
    
    >>> create_database("new_db")
    Database 'new_db' already exists. Operation aborted.
    """
    # Check if the database already exists
    if db_name in bw.databases:
        print(f"Database '{db_name}' already exists. Operation aborted.")
        return
    
    # Create and register the new database
    db = bw.Database(db_name)
    db.register()
    print(f"Database '{db_name}' has been created and registered.")

def add_exchange_to_activity(act, activity_name, location, ref_product, amount, db_name_search, act_type='technosphere'):
    """
    Adds a technosphere exchange to a specified activity in the Brightway2 database.

    Parameters:
    act (object): The activity to which the exchange will be added.
    activity_name (str): The name of the activity to be added as an exchange.
    location (str): The location of the activity to be added as an exchange.
    ref_product (str): The reference product of the activity to be added as an exchange.
    amount (float): The amount of the exchange to be added.
    db_name_search (str): The name of the database to search for the activity.
    act_type (str): The type of the exchange, default is 'technosphere'.

    Returns:
    None
    """
    # Retrieve the activity from the specified database
    add_act = [act for act in bw.Database(db_name_search) if activity_name == act['name']
               and act['location'] == location and
              act['reference product'] == ref_product][0]
    
    # Add the exchange
    act.new_exchange(**{
        'input': add_act.key,
        'amount': amount,
        'type': act_type                      
    }).save()

def add_exchange_to_activity_bio(act, bio_name, bio_cat, amount, db_name_search):
    """
    Adds a biosphere exchange to a specified activity in the Brightway2 database.

    Parameters:
    act (object): The activity to which the exchange will be added.
    bio_name (str): The name of the biosphere activity to be added as an exchange.
    bio_cat (tuple): The categories of the biosphere activity to be added as an exchange.
    amount (float): The amount of the exchange to be added.
    db_name_search (str): The name of the database to search for the biosphere activity.

    Returns:
    None
    """
    # Retrieve the activity from the specified database
    for bio in bw.Database(db_name_search):
            if (bio_name == bio['name']) and (bio['categories'] == bio_cat):
                bio_flow = bio.key
    
    # Add the exchange
    act.new_exchange(**{
        'input': bio_flow,
        'amount': amount,
        'type': 'biosphere'                       
    }).save()

# ### generate the database which we are going to use, as premise include many novel datasets. Add some datasets that we generated ourselves.
def generate_reference_database():
    """
    Generate a reference database based on specified parameters.

    This function generates a reference database based on the specified parameters.
    It deletes the existing reference database with the same name if it exists and then creates a new one.

    Returns:
    None

    Example:
    >>> generate_reference_database()
    """
    # Delete old reference database with the same name
    for db_name in list(bw.databases):
        if NAME_REF_DB in db_name:
            print(db_name)
            del bw.databases[db_name]

    # Create a new reference database using NewDatabase
    ndb = pr.NewDatabase(
        scenarios=[{"model": "remind", "pathway": 'SSP2-Base', "year": "2024"}],
        source_db=DB_NAME_INIT,
        source_version=EI_VERSION,
        key=KEY_PREMISE,
        biosphere_name= BIOSPHERE_DB,
        additional_inventories=[
            {"filepath": r"data\H2-DRI_LCI.xlsx", "ecoinvent version": "3.9.1"},
            {"filepath": r"data\BF-BOF-CCS_Carina.xlsx", "ecoinvent version": "3.9.1"},
            {"filepath": r"data\lci_biofuels_oil.xlsx", "ecoinvent version": EI_VERSION}]
    )

    # Fix suggested by R.S. since premise expect modifications
    ndb.scenarios[0]["database"] = ndb.database
    ndb.write_db_to_brightway(name=NAME_REF_DB)

def generate_future_ei_dbs(scenarios = ["SSP2-Base", "SSP2-PkBudg1150","SSP2-PkBudg500"], iam = 'remind',
                           start_yr=2020, end_yr = 2050, step = 5, endstring="base"):
    """
    Generate Ecoinvent scenario models with specified parameters.

    This function generates Ecoinvent scenario models based on the specified year, scenario, and additional settings.
    It avoids adding duplicated databases by checking the existing databases in Brightway2.

    Parameters:
    - scenarios (list): The scenarios for which the models are generated. Default is: 
                ["SSP2-Base",
                "SSP2-PkBudg1150",
                "SSP2-PkBudg500"] corresponding to baseline, 2 degrees C, and 1.5 degrees C.
    - iam (str): IAM chosen, can be 'remind' or 'image'.
    - start_yr (int): The starting year for the scenarios.
    - end_yr (int): The end year for the scenarios.
    - step (int): step between scenario years.
    - endstring (str, optional): A suffix to differentiate the generated databases. Default is "base".

    Returns:
    tuple: A tuple containing two lists -
        1. List of dictionaries specifying the models for the scenarios.
        2. List of database names generated based on the specified parameters.

    Example:
    >>> generate_future_ei_dbs("SSP2-Base")
    ([{'model': 'remind', 'pathway': 'SSP2-Base', 'year': 2030},
      {'model': 'remind', 'pathway': 'SSP2-Base', 'year': 2050}],
     ['ecoinvent_remind_SSP2-Base_2030_custom', 'ecoinvent_remind_SSP2-Base_2050_base'])
    """

    list_years = [start_yr + i * step for i in range(1, int((end_yr - start_yr) / step) + 1)]

    list_spec_scenarios = []
    list_names = []

    for pt in scenarios:
        for yr in list_years:
            string_db = "ecoinvent_{}_{}_{}_{}".format(iam, pt, yr, endstring)

            if yr == start_yr and pt == "SSP2-Base":
                dict_spec = {"model": iam, "pathway": pt, "year": yr,
                                "exclude": ["update_electricity", "update_cement", "update_steel", "update_dac",
                                            "update_fuels", "update_emissions", "update_two_wheelers"
                                            "update_cars", "update_trucks", "update_buses"]}

                if string_db not in bw.databases:
                    list_spec_scenarios.append(dict_spec)
                    list_names.append(string_db)
                else:
                    print("Avoid duplicated db and therefore following db not added: '{}'".format(string_db))
            else:
                dict_spec = {"model": iam, "pathway": pt, "year": yr}

                if string_db not in bw.databases:
                    list_spec_scenarios.append(dict_spec)
                    list_names.append(string_db)
                else:
                    print("Avoid duplicated db and therefore following db not added: '{}'".format(string_db))

    return list_spec_scenarios, list_names

def generate_prospective_lca_dbs(list_spec_scenarios, list_names):
    """
    Generate and update future LCA databases for prospective LCA.

    This function generates and updates future LCA databases based on specified scenarios if needed, and writes them to Brightway2.

    Parameters:
    - list_spec_scenarios (list): List of dictionaries specifying scenarios for new databases.
    - list_names (list): List of names specifying scenario names for new databases.

    Returns:
    None

    Example:
    >>> generate_and_update_lca_databases([{"model": "remind", "pathway": 'SSP2-Base', "year": "2035"}], ['ecoinvent_remind_SSP2-Base_2035_base'])
    """

    if len(list_spec_scenarios) > 0:
        ndb = pr.NewDatabase(
            scenarios=list_spec_scenarios,
            source_db=DB_NAME_INIT,
            source_version=EI_VERSION,
            key=KEY_PREMISE,
            biosphere_name=BIOSPHERE_DB,
            additional_inventories=[
                {"filepath": r"data\H2-DRI_LCI.xlsx", "ecoinvent version": "3.9.1"},
                {"filepath": r"data\BF-BOF-CCS_Carina.xlsx", "ecoinvent version": "3.9.1"},
                {"filepath": r"data\lci_biofuels_oil.xlsx", "ecoinvent version": EI_VERSION}]
        )

        print("START UPDATING")
        ndb.update()

        print("START WRITING")
        ndb.write_db_to_brightway(name = list_names)

def try_find_bw_act(database_name, activity_name, ref_product, location):
    """
    Finds a single activity in a Brightway2 database matching the given name, reference product, and location.

    This function searches through the specified database for activities that match the provided
    name, reference product, and location. If more than one matching activity is found, or if no
    matching activities are found, it tries alternative locations ("GLO", "RoW", "RER", "Europe without Switzerland") in that order.

    Parameters:
    database_name (str): The name of the Brightway2 database to search in.
    activity_name (str): The name of the activity to search for.
    ref_product (str): The reference product of the activity to search for.
    location (str): The location of the activity to search for.

    Returns:
    object: The single matching activity.

    Raises:
    ValueError: If more than one matching activity is found.
    ValueError: If no matching activity is found.
    """

    locations_to_try = [location, "GLO", "RoW", "RER", "Europe without Switzerland", "ES", "FR"] #FR as ultimate proxy for renewables, as it is quite average reg wind and solar PV...

    for loc in locations_to_try:
        if 'market for electricity' in str(activity_name) and 'electricity,' in str(ref_product):
            for i, activity_name_l in enumerate([activity_name, activity_name.replace('market for electricity', 'market group for electricity')]):
                # Filter activities matching the given criteria
                matching_activities = [
                    act for act in bw.Database(database_name)
                    if activity_name_l == act['name'] and loc == act['location']
                    and ref_product == act['reference product']
                ]

                if len(matching_activities) == 1:
                    return matching_activities[0]  # Return the single matching activity

                elif len(matching_activities) > 1:
                    raise ValueError(f"More than one activity found for name '{activity_name}' in location '{loc}'")

        else:   
            # Filter activities matching the given criteria
            matching_activities = [
                act for act in bw.Database(database_name)
                if activity_name == act['name'] and loc == act['location']
                and ref_product == act['reference product']
            ]

        if len(matching_activities) == 1:
            return matching_activities[0]  # Return the single matching activity

        elif len(matching_activities) > 1:
            raise ValueError(f"More than one activity found for name '{activity_name}' in location '{loc}'")
    
    # If no activity is found after trying all locations
    raise ValueError(f"No activity found for name '{activity_name}' with the given locations")

def calculate_lca_impact(act_name, ref_product, loc, impact_method = MY_METHODS[0], db = 'ecoinvent_remind_SSP2-PkBudg1150_2040_base'):
    """
    Calculate the LCA impact for a given activity name, reference product, and location.
    
    Parameters:
    - act_name (str): The name of the activity.
    - ref_product (str): The reference product for the activity.
    - loc (str): The location for the activity.
    - impact_method (tuple): The impact assessment method, e.g., ('IPCC 2013', 'climate change', 'GWP 100a').

    Returns:
    - float: The LCA impact score.
    """
    # Retrieve the activity based on name, reference product, and location
    activity = [act for act in bw.Database(db) 
                if act['name'] == act_name and act['reference product'] == ref_product and act['location'] == loc]
    
    if not activity:
        raise ValueError(f"Activity not found with name '{act_name}', reference product '{ref_product}', and location '{loc}'")
    
    # Perform the LCA calculation
    lca = bw.LCA({activity[0]: 1}, impact_method)
    lca.lci()
    lca.lcia()
    
    # Return the LCA impact score
    return lca.score

def find_bw_act(database_name, activity_name, ref_product, location):
    """
    Finds a single activity in a Brightway2 database matching the given name, reference product, and location.
    
    This function searches through the specified database for activities that match the provided
    name, reference product, and location. If more than one matching activity is found, or if no
    matching activities are found, an error is raised.
    
    Parameters:
    database_name (str): The name of the Brightway2 database to search in.
    activity_name (str): The name of the activity to search for.
    ref_product (str): The reference product of the activity to search for.
    location (str): The location of the activity to search for.
    
    Returns:
    object: The single matching activity.
    
    Raises:
    ValueError: If more than one matching activity is found.
    ValueError: If no matching activity is found.
    
    Example:
    >>> find_bw_act("db_ei", "electrolyzer production, 1MWe, PEM, Stack", "stack", "RER")
    {
        'name': 'electrolyzer production, 1MWe, PEM, Stack',
        'reference product': 'stack',
        'location': 'RER',
        ...
    }
    """
    # Filter activities matching the given criteria
    matching_activities = [
        act for act in bw.Database(database_name) 
        if activity_name == act['name'] and location == act['location']
        and ref_product == act['reference product']
    ]
    
    # Check if there is more than one matching activity
    if len(matching_activities) > 1:
        raise ValueError(f"More than one activity found for name '{activity_name}' in location '{location}'")
    
    # Check if no activity is found
    if len(matching_activities) == 0:
        raise ValueError(f"No activity found for name '{activity_name}' in location '{location}'")
    
    # Return the single matching activity
    return matching_activities[0]

def silent_search(db_name, name_search, desired_location):
    """
    Perform a Brightway database search while suppressing console output.

    Parameters:
    - db_name (str): Name of the Brightway2 database to search.
    - name_search (str): Substring to match in activity names.
    - desired_location (str): Location to filter the search results by (exact match).

    Returns:
    - list: List of matching activities (as Brightway2 activity objects).
    """
    f = io.StringIO()
    with redirect_stdout(f):
        results = list(bw.Database(db_name).search(
            name_search,
            limit=1000,
            filter={"location": desired_location}
        ))
    return results

def find_activity_fast(db_name, name_search, desired_ref_product, desired_location):
    """
    Quickly find a Brightway2 activity by name, location, and reference product using filtered search.

    Parameters:
    - db_name (str): Name of the Brightway2 database.
    - name_search (str): Exact name of the activity to find.
    - desired_ref_product (str or None): Reference product to filter activities by.
    - desired_location (str or None): Location to filter activities by.

    Returns:
    - dict: A single matching Brightway2 activity dictionary.

    Raises:
    - ValueError: If multiple matching activities are found.
    """
    # use search() function first, this is much faster to match activities.
    candidates = silent_search(db_name, name_search, desired_location)
    
    if len(candidates) != 1:
        matching_activities = [
            act for act in candidates
            if act['location']== desired_location and
            act['name'] == name_search
            and act['reference product'] == desired_ref_product
        ]
    else:
        return candidates[0]
        
    # Check if there is more than one matching activity
    if len(matching_activities) > 1:
        raise ValueError(f"More than one activity found for name '{name_search}' in location '{desired_location}'")
            
    elif len(matching_activities) == 0:
        #print(f"Cannot find activity using search function for name '{name_search}' in location '{desired_location}', using list comprehension instead.")
        matching_activities = try_find_bw_act(db_name, name_search, desired_ref_product, desired_location)
        return matching_activities  
    else:
        # Return the single matching activity
        return matching_activities[0]

def run_mlca(
    activity_keys: list,
    functional_units: list,
    result_index_labels: list,
    column_suffix: str,
    impact_methods: list,
) -> pd.DataFrame:
    """
    Run a MultiLCA calculation and return a DataFrame of results.

    Parameters:
    - activity_keys: list of Brightway activity keys for LCA calculation
    - functional_units: list of float values (same length as activity_keys), quantities of each key
    - result_index_labels: labels for resulting DataFrame index (must match order of keys)
    - column_suffix: suffix for LCA impact columns (e.g., '', '_ref')
    - impact_methods: list of LCIA methods to evaluate

    Returns:
    - pd.DataFrame: LCA results with impacts per row (indexed by result_index_labels)
    """
    assert len(activity_keys) == len(functional_units) == len(result_index_labels), \
        "Keys, functional units, and labels must be the same length."

    setup_name = f"mlca_setup_{column_suffix.strip('_') or 'main'}"

    # Prepare inventory with variable functional units
    inventory = [{key: fu} for key, fu in zip(activity_keys, functional_units)]

    # Define calculation setup
    bw.calculation_setups[setup_name] = {
        'inv': inventory,
        'ia': impact_methods,
        'name': setup_name,
    }

    # Run MultiLCA
    mlca = MultiLCA(setup_name)
    results_array = mlca.results

    # Create labeled results DataFrame
    result_df = pd.DataFrame(
        data=results_array,
        index=result_index_labels,
        columns=[f"lca_impact{column_suffix}_{method[1]}" for method in impact_methods]
    )

    return result_df

def delete_project(bw, project_name):
    """
    Deletes a specified project if it exists.

    Parameters:
    - bw: The object that manages the projects (brightway2).
    - project_name (str): The name of the project to be deleted.

    Returns:
    - str: A message indicating the result of the deletion attempt.
    """
    # List all existing projects
    print("Existing projects:", list(bw.projects))

    # Check if the project exists before attempting to delete it
    if project_name in bw.projects:
        # Set the project to be deleted as the current project
        bw.projects.set_current(project_name)
        
        # Delete the project
        bw.projects.delete_project(project_name)
        bw.projects.purge_deleted_directories()
        return f"Project '{project_name}' has been deleted."
    else:
        return f"Project '{project_name}' does not exist."

def create_composite_string(row, columns):
    """
    Create a composite string from specified columns of a dataframe row.
    """
    return ' '.join([str(row[col]) for col in columns if pd.notna(row[col])])

def find_best_match(row, df, index_columns = SEL_COLS):
    """
    Find the best match for the given row in the dataframe based on the composite index columns.
    """
    search_str = create_composite_string(row, index_columns)
    df['composite_string'] = df.apply(lambda x: create_composite_string(x, index_columns), axis=1)
    best_match = process.extractOne(search_str, df['composite_string'], scorer=fuzz.token_sort_ratio)
    return best_match

def get_value_from_tif(filepath, lon, lat):
    """
    Retrieve the value from a TIFF file at a specified longitude and latitude.

    Parameters:
        filepath (str): The path to the TIFF file.
        lon (float): The longitude of the point of interest.
        lat (float): The latitude of the point of interest.

    Returns:
        float: The value at the specified longitude and latitude.
    """
    # Open the raster TIFF file
    raster_data = ri.open_rasterio(filepath)

    # Slice the first band of the raster
    img_data = raster_data[0, :, :]

    # Use the .sel() method to retrieve the value of the nearest cell close to the POI
    value = img_data.sel(x=lon, y=lat, method="nearest").data.item()

    # the data for PV is in kWh/day/kWp, so we will have to do following modification
    if "PVOUT" in str(filepath):
        value = (value * 365)/HOURS_YR
    
    return value

def process_row_best_match(row, df, index_columns = SEL_COLS):
    """
    Process each row, matching it to the best 'composite_string' in the dataframe.
    """
    # Print the current row being processed for debugging
    #print(f"Processing row: {row}", end='\r', flush=True)

    try:
        best_match = find_best_match(row, df, index_columns = index_columns)
        if best_match is None:
            warnings.warn(f"No match found for row: {row}.")
        else:
            best_match_row = df[df['composite_string'] == best_match[0]].iloc[0]
            return best_match_row
    except Exception as e:
        warnings.warn(f"Error processing row: {row}. Details: {str(e)}")

def search_biosphere_entries(database, name, category, unit):
    """
    Search for entries in the biosphere database that match the specified name, category, and unit.

    Parameters:
    database (list of dict): The biosphere database to search, where each entry is a dictionary containing 'name', 'categories', and 'unit'.
    name (str): The name to match in the biosphere entries.
    category (str): The category to match in the biosphere entries.
    unit (str): The unit to match in the biosphere entries.

    Returns:
    list of dict: A list of biosphere entries that match the specified criteria.
    """
    matching_activities = [bio for bio in bw.Database(database)
            if bio['name'] == name and 
               bio['categories'] == category and 
               bio['unit'] == unit]
    
    # Check if there is more than one matching activity
    if len(matching_activities) > 1:
        raise ValueError(f"More than one activity found for name '{name}' with category '{category}'")
    
    # Check if no activity is found
    if len(matching_activities) == 0:
        raise ValueError(f"No activity found for name '{name}' with category '{category}'")
    
    # Return the single matching activity
    return matching_activities[0]

def GENERATE_ACTS_GM_PV(db_ei, cf_pv, curtailment_pv, cf_electrolyzer=0.3, electrolyzer="pem", excel_col_name = 'ecoinvent_310_reference'):
    """
    Generates a new activity in the Brightway2 database for hydrogen production
    from electrolysis using photovoltaic (PV) ground-mounted systems.

    Parameters:
    db_ei (str): The name of the Brightway2 database to be used.
    cf_pv (float): capacity factor of solar PV system.
    curtailment_pv (float): curtailment of solar PV to account for oversizing.
    cf_electrolyzer (float): The capacity factor of the PV system, which is 0.3 according to IEA.
    type of electrolyzer (str): Type of electrolyzer (can be 'pem', 'aec', or 'soec')

    Returns:
    None: The function performs operations directly on the Brightway2 database.

    Description:
    This function creates a new activity for hydrogen production using electrolysis powered by
    ground-mounted photovoltaic systems. It calculates the necessary inputs and adds exchanges
    for the activity based on the specified capacity factor of the PV system. The activity is 
    only created if it does not already exist in the database.

    Steps performed:
    1. Check if the activity already exists in the database.
    2. Retrieve necessary technical and cost data.
    3. Create the new activity with the specified parameters.
    4. Add exchanges for:
        - Stack electrolyzer production and treatment
        - Balance of plant electrolyzer production and treatment
        - Water consumption
        - PV infrastructure and water usage
        - Biosphere flows for waste heat, solar energy, hydrogen leakage, and oxygen production.
    Includes compression provided by solar PV, depending on the electrolyzer: ELECTROLYZER_DATA[electrolyzer]['compression_kWh_per_kgH2']

    Example:
    >>> GENERATE_ACTS_gm_pv('ecoinvent_3.7.1', 0.18)
    Skipped: activity 'hydrogen production, gaseous, 30 bar, from PEM electrolysis, solar PV ground-mounted, global cf [0.18]' already generated!
    """

    capacity = 560 * 0.895 #kWp, 10.5% average degradation, described in activity
    kg_water_unit= 2871.8 * 20 #20 L per m2 module, 2871.8 square meter for open ground construction, on ground, Mont Soleil

    new_name = "hydrogen production, gaseous, {} bar, from {} electrolysis, solar PV ground-mounted, global cf [{}]".format(30,electrolyzer.upper(), round(cf_pv,3))
    check_act = [act for act in bw.Database(db_ei) if new_name == act['name']]
    
    tech = "solar_pv_gm"
    
    lifetime_pv = COST_DATA.loc[(tech,'lifetime')][excel_col_name].item()
    eff_elect = COST_DATA.loc[('electrolyzer_{}'.format(electrolyzer),'eff')][excel_col_name]
    lifetime_stack = COST_DATA.loc[('electrolyzer_{}'.format(electrolyzer),'stack_lt')][excel_col_name]  
    lifetime_bop = COST_DATA.loc[('electrolyzer_{}'.format(electrolyzer),'bop_lt')][excel_col_name]  
    h2_leakage_factor = COST_DATA.loc[('h2_leakage','-')][excel_col_name] 
    land_m2_kw = COST_DATA.loc[('electrolyzer_{}'.format(electrolyzer),'land_m2_kw')][excel_col_name]  
    
    if len(check_act) > 0:
        print("Skipped: activity '{}' already generated!".format(new_name))
    else:
        db = bw.Database(db_ei)
        code_na = str(uuid.uuid4().hex)
        act = db.new_activity(
            **{
                'name': new_name,
                 "code": code_na,
                'unit': 'kilogram',
                'reference product': "hydrogen, gaseous, {} bar".format(30),
                'location' :"GLO",
                'production amount': 1.0,
                'comment': f"Hydrogen production via water electrolysis with a \
                    {electrolyzer} electrolyzer with a capacity factor of {round(cf_electrolyzer,2)} and solar PV capacity factor of {round(cf_pv,3)}, includes compression if AEC or SOEC with local renewables"
            }
        )

        act.save()

        # Add production exchange
        act.new_exchange(**{
            'input': (db_ei, code_na),
            'amount': 1,
            'type': 'production',
        }).save()

        # calculate electricity requirements
        amount_electricity = (EN_H2/(3.6*eff_elect))

        # 1.1. Add stack electrolyzer
        if electrolyzer == "pem":
            add_act = find_bw_act(db_ei, 'electrolyzer production, 1MWe, PEM, Stack', 'electrolyzer, 1MWe, PEM, Stack', 'RER')

            amount_elect_unit = 1/(1000*cf_electrolyzer*HOURS_YR*lifetime_stack/amount_electricity)

            act.new_exchange(amount= amount_elect_unit, input = add_act.key, type="technosphere", location = add_act['location'],
                             name = add_act['name'], product = add_act['reference product'] ).save()
            
            add_act = find_bw_act(db_ei, 'treatment of fuel cell stack, 1MWe, PEM', 'used fuel cell stack, 1MWe, PEM', 'RER')

            act.new_exchange(amount= -amount_elect_unit, input = add_act.key, type="technosphere", location = add_act['location'],
                             name = add_act['name'], product = add_act['reference product'] ).save()

            # 1.2. Add BoP electrolyzer
            add_act = find_bw_act(db_ei, 'electrolyzer production, 1MWe, PEM, Balance of Plant', 'electrolyzer, 1MWe, PEM, Balance of Plant', 'RER')

            amount_bop_unit = (lifetime_stack/lifetime_bop)/(1000*cf_electrolyzer*HOURS_YR*lifetime_stack/amount_electricity)  #lifetime bop is 20 years

            act.new_exchange(amount=amount_bop_unit, input = add_act.key, type="technosphere", location = add_act['location'],
                             name = add_act['name'], product = add_act['reference product'] ).save()
            
            add_act = find_bw_act(db_ei, 'treatment of fuel cell balance of plant, 1MWe, PEM', 'used fuel cell balance of plant, 1MWe, PEM', 'RER')

            act.new_exchange(amount= -amount_bop_unit, input = add_act.key, type="technosphere", location = add_act['location'],
                             name = add_act['name'], product = add_act['reference product'] ).save()

        elif electrolyzer == "aec":         
            add_act = find_bw_act(db_ei, 'electrolyzer production, 1MWe, AEC, Stack', 'electrolyzer, 1MWe, AEC, Stack', 'RER')

            amount_elect_unit = 1/(1000*cf_electrolyzer*HOURS_YR*lifetime_stack/amount_electricity)

            act.new_exchange(amount= amount_elect_unit, input = add_act.key, type="technosphere", location = add_act['location'],
                             name = add_act['name'], product = add_act['reference product'] ).save()

            add_act = find_bw_act(db_ei, 'treatment of fuel cell stack, 1MWe, AEC', 'used fuel cell stack, 1MWe, AEC', 'RER')

            act.new_exchange(amount= -amount_elect_unit, input = add_act.key, type="technosphere", location = add_act['location'],
                             name = add_act['name'], product = add_act['reference product'] ).save()

            # 1.2. Add BoP electrolyzer            
            add_act = find_bw_act(db_ei, 'electrolyzer production, 1MWe, AEC, Balance of Plant', 'electrolyzer, 1MWe, AEC, Balance of Plant', 'RER')

            amount_bop_unit = (lifetime_stack/lifetime_bop)/(1000*cf_electrolyzer*HOURS_YR*lifetime_stack/amount_electricity)  #lifetime bop is 20 years

            act.new_exchange(amount=amount_bop_unit, input = add_act.key, type="technosphere", location = add_act['location'],
                             name = add_act['name'], product = add_act['reference product'] ).save()

            add_act = find_bw_act(db_ei, 'treatment of fuel cell balance of plant, 1MWe, AEC', 'used fuel cell balance of plant, 1MWe, AEC', 'RER')

            act.new_exchange(amount= -amount_bop_unit, input = add_act.key, type="technosphere", location = add_act['location'],
                             name = add_act['name'], product = add_act['reference product'] ).save()

            # 1.2. Add potassium hydroxide
            add_act = find_bw_act(db_ei, 'market for potassium hydroxide', 'potassium hydroxide', 'GLO')

            act.new_exchange(amount=3.70*1e-3, input = add_act.key, type="technosphere", location = add_act['location'],
                             name = add_act['name'], product = add_act['reference product'] ).save()

        elif electrolyzer == 'soec':  
            add_act = find_bw_act(db_ei, 'electrolyzer production, 1MWe, SOEC, Stack', 'electrolyzer, 1MWe, SOEC, Stack', 'RER')

            amount_elect_unit = 1/(1000*cf_electrolyzer*HOURS_YR*lifetime_stack/amount_electricity)

            act.new_exchange(amount= amount_elect_unit, input = add_act.key, type="technosphere", location = add_act['location'],
                             name = add_act['name'], product = add_act['reference product'] ).save()

            add_act = find_bw_act(db_ei, 'treatment of fuel cell stack, 1MWe, SOEC', 'used fuel cell stack, 1MWe, SOEC', 'RER')

            act.new_exchange(amount= -amount_elect_unit, input = add_act.key, type="technosphere", location = add_act['location'],
                             name = add_act['name'], product = add_act['reference product'] ).save()

            # 1.2. Add BoP electrolyzer            
            add_act = find_bw_act(db_ei, 'electrolyzer production, 1MWe, SOEC, Balance of Plant', 'electrolyzer, 1MWe, SOEC, Balance of Plant', 'RER')

            amount_bop_unit = (lifetime_stack/lifetime_bop)/(1000*cf_electrolyzer*HOURS_YR*lifetime_stack/amount_electricity)  #lifetime bop is 20 years

            act.new_exchange(amount=amount_bop_unit, input = add_act.key, type="technosphere", location = add_act['location'],
                             name = add_act['name'], product = add_act['reference product'] ).save()

            add_act = find_bw_act(db_ei, 'treatment of fuel cell balance of plant, 1MWe, SOEC', 'used fuel cell balance of plant, 1MWe, SOEC', 'RER')

            act.new_exchange(amount= -amount_bop_unit, input = add_act.key, type="technosphere", location = add_act['location'],
                             name = add_act['name'], product = add_act['reference product'] ).save()
        else:
            raise ValueError("Not possible to user electrolyzer exchanges since electrolyzer '{}' doesn't exist in LCI".format(electrolyzer))

        # 1.3. Add direct water consumption        
        add_act = find_bw_act(db_ei, 'market for water, deionised', 'water, deionised', 'RoW')

        # Add the exchange
        act.new_exchange(input=add_act.key,amount=AMOUNT_WATER_ELECTROLYSIS,
                                    unit="kilogram",type='technosphere', location = add_act['location'],
                             name = add_act['name'], product = add_act['reference product']).save()
        
        # If needed, add compression power demand to this initial power demandfor water electrolysis to arrive at 30 bar.
        # Here, we assume that the compression power demand is satisfied by local renewables
        amount_electricity_w_compr = amount_electricity + ELECTROLYZER_DATA[electrolyzer]['compression_kWh_per_kgH2']
        
        # 1.4. Add infrastructure activity, PV unit
        add_act = find_bw_act(db_ei, 'photovoltaic open ground installation, 560 kWp, single-Si, on open ground', 
                              'photovoltaic open ground installation, 560 kWp, single-Si, on open ground', 'CH')
        amount_kwh = 1/(capacity*lifetime_pv*cf_pv*HOURS_YR*(1-curtailment_pv))

        # Add the exchange
        act.new_exchange(input=add_act.key,amount=(amount_kwh*amount_electricity_w_compr),
                                        unit="unit",type='technosphere', location = add_act['location'],
                             name = add_act['name'], product = add_act['reference product']).save()

        # 1.5 Add tap water
        #tap water	market for tap water	Europe without Switzerland	kilogram	
        add_act = find_bw_act(db_ei, 'market for tap water', 
                              'tap water', 'Europe without Switzerland')

        amount_kwh = kg_water_unit*(1/(capacity*lifetime_pv*cf_pv*HOURS_YR*(1-curtailment_pv)))

        # Add the exchange
        act.new_exchange(input=add_act.key,amount=(amount_kwh*amount_electricity_w_compr),
                                        unit="kilogram",type='technosphere', location = add_act['location'],
                             name = add_act['name'], product = add_act['reference product']).save()

        # 1.6 Add wastewater, from residence
        add_act = find_bw_act(db_ei, 'treatment of wastewater, average, wastewater treatment', 
                              'wastewater, average', 'RoW')
        amount_kwh = kg_water_unit*(1/(capacity*lifetime_pv*cf_pv*HOURS_YR*(1-curtailment_pv)))/1000

        # Add the exchange
        act.new_exchange(input=add_act.key,amount=(-amount_kwh*amount_electricity_w_compr),
                                        unit="cubic meter",type='technosphere', location = add_act['location'],
                             name = add_act['name'], product = add_act['reference product']).save()

        # 1.7 Add biosphere flows
        bio_flow = search_biosphere_entries(BIOSPHERE_DB, 'Heat, waste', ('air', 'urban air close to ground'), 'megajoule')
        amount_kwh = 0.25027 # curtailment doesn't have influence here, as solar PV is switched off

        # Add the exchange
        act.new_exchange(input=bio_flow.key,amount=(amount_kwh*amount_electricity_w_compr),
                                        unit="megajoule",type='biosphere').save()

        # 1.8 add solar energy converted
        bio_flow = search_biosphere_entries(BIOSPHERE_DB, 'Energy, solar, converted', ('natural resource', 'in air'), 'megajoule')

        amount_kwh = 3.8503 # curtailment doesn't have influence here, as solar PV is switched off
        # Add the exchange
        act.new_exchange(input=bio_flow.key,amount=(amount_kwh*amount_electricity_w_compr),
                                        unit="megajoule",type='biosphere').save()
        
        
        # 1.9 Add biosphere flows - H2 Leakage
        bio_flow = search_biosphere_entries(BIOSPHERE_DB, 'Hydrogen', ('air',), 'kilogram')

        # Add the exchange
        act.new_exchange(input=bio_flow.key,amount=h2_leakage_factor,
                                        unit="kilogram",type='biosphere').save()
        
        # 1.10 Oxygen production
        bio_flow = search_biosphere_entries(BIOSPHERE_DB, 'Oxygen', ('air',), 'kilogram')

        # Add the exchange
        act.new_exchange(input=bio_flow.key,amount=8,
                                        unit="kilogram",type='biosphere').save()

        # Finally, add the biosphere flows for land occupation.
        # 1.11 Add biosphere flows - land occupation
        bio_flow = search_biosphere_entries(BIOSPHERE_DB, 'Occupation, industrial area', ('natural resource','land'), 'square meter-year')

        # Electrolyser land footprint: (land_m2_kw [m2/kW]*1000[kW/MW])/(H2_prod[kg H2 per lifetime]/lifetime_system[years])
        act.new_exchange(input=bio_flow.key,amount=((land_m2_kw*1000)/(1000*cf_electrolyzer*HOURS_YR/amount_electricity)),
                                        unit="square meter-year",type='biosphere').save()
        
        bio_flow = search_biosphere_entries(BIOSPHERE_DB, 'Transformation, from industrial area', ('natural resource','land'), 'square meter')   

        # Electrolyser land transformation: (land_m2_kw [m2/kW]*1000 [kW/MW])/H2_prod [kg H2]
        act.new_exchange(input=bio_flow.key,amount=(land_m2_kw*1000)/(1000*cf_electrolyzer*HOURS_YR*lifetime_bop/amount_electricity),
                                        unit="square meter",type='biosphere').save()
        
        bio_flow = search_biosphere_entries(BIOSPHERE_DB, 'Transformation, to industrial area', ('natural resource','land'), 'square meter')   
        # Electrolyser land transformation: (land_m2_kw [m2/kW]*1000 [kW/MW])/H2_prod [kg H2]
        act.new_exchange(input=bio_flow.key,amount=(land_m2_kw*1000)/(1000*cf_electrolyzer*HOURS_YR*lifetime_bop/amount_electricity),
                                        unit="square meter",type='biosphere').save()
        
        act.save()

def GENERATE_ACTS_WIND(db_ei, cf_wind, curtailment_wind, cf_electrolyzer=0.4, electrolyzer = "pem", excel_col_name = 'ecoinvent_310_reference'):
    """
    Generates a new activity in the Brightway2 database for hydrogen production
    from electrolysis using onshore wind energy.

    Parameters:
    db_ei (str): The name of the Brightway2 database to be used.
    cf_wind (float): capacity factor of onshore wind system.
    curtailment_wind (float): curtailment of onshore wind to account for oversizing.
    cf_electrolyzer (float): The capacity factor of the wind energy source, which is 0.4 according to IEA for onshore wind.
    type of electrolyzer (str): Type of electrolyzer (can be 'pem', 'aec', or 'soec')

    Returns:
    None: The function performs operations directly on the Brightway2 database.

    Description:
    This function creates a new activity for hydrogen production using electrolysis powered by
    onshore wind energy. It calculates the necessary inputs and adds exchanges for the activity
    based on the specified capacity factor. The activity is only created if it does not already exist
    in the database.

    Steps performed:
    1. Check if the activity already exists in the database.
    2. Retrieve necessary technical and cost data.
    3. Create the new activity with the specified parameters.
    4. Add exchanges for:
        - Stack electrolyzer production
        - Balance of plant electrolyzer production
        - Water consumption
        - Lubricating oil and its waste treatment
        - Wind turbine network connection
        - Biosphere flows for kinetic wind energy, hydrogen leakage, and oxygen production.

    Includes compression provided by onshore wind, depending on the electrolyser: ELECTROLYZER_DATA[electrolyzer]['compression_kWh_per_kgH2']
    
    Example:
    >>> GENERATE_ACTS_wind('ecoinvent_3.7.1', 0.25)
    Skipped: activity 'hydrogen production, gaseous, 30 bar, from PEM electrolysis, onshore wind, global cf [0.25]' already generated!
    """

    capacity = 2000 #kWp
    new_name = "hydrogen production, gaseous, {} bar, from {} electrolysis, onshore wind, global cf [{}]".format(30,electrolyzer.upper(),round(cf_wind,3))
    check_act = [act for act in bw.Database(db_ei) if new_name == act['name']]
    
    tech = "onshore_wind"
    
    lifetime_wind = COST_DATA.loc[(tech,'lifetime')][excel_col_name]
    eff_elect = COST_DATA.loc[('electrolyzer_{}'.format(electrolyzer),'eff')][excel_col_name]
    lifetime_stack = COST_DATA.loc[('electrolyzer_{}'.format(electrolyzer),'stack_lt')][excel_col_name]  
    lifetime_bop = COST_DATA.loc[('electrolyzer_{}'.format(electrolyzer),'bop_lt')][excel_col_name]  
    h2_leakage_factor = COST_DATA.loc[('h2_leakage','-')][excel_col_name] 
    land_m2_kw = COST_DATA.loc[('electrolyzer_{}'.format(electrolyzer),'land_m2_kw')][excel_col_name]  

    if len(check_act) > 0:
        print("Skipped: activity '{}' already generated!".format(new_name))
    else:
        # Generate activity
        db = bw.Database(db_ei)
        code_na = str(uuid.uuid4().hex)
        act = db.new_activity(
            **{
                'name': new_name,
                 "code": code_na,
                'unit': 'kilogram',
                'reference product': "hydrogen, gaseous, {} bar".format(30),
                'location' :"GLO",
                'production amount': 1.0,
                'comment': f"Hydrogen production via water electrolysis with a \
                    {electrolyzer} electrolyzer with a capacity factor of {round(cf_electrolyzer,2)} and onshore wind capacity factor of {round(cf_wind,2)}, includes compression if AEC or SOEC with local renewables"
            }
        )

        act.save()

        # Add production exchange
        act.new_exchange(**{
            'input': (db_ei, code_na),
            'amount': 1,
            'type': 'production',
        }).save()

        # calculate electricity requirements
        amount_electricity = (EN_H2/(3.6*eff_elect))

        # 1.1. Add stack electrolyzer
        if electrolyzer == "pem":
            add_act = find_bw_act(db_ei, 'electrolyzer production, 1MWe, PEM, Stack', 'electrolyzer, 1MWe, PEM, Stack', 'RER')

            amount_elect_unit = 1/(1000*cf_electrolyzer*HOURS_YR*lifetime_stack/amount_electricity)

            act.new_exchange(amount= amount_elect_unit, input = add_act.key, type="technosphere", location = add_act['location'],
                             name = add_act['name'], product = add_act['reference product'] ).save()
            
            add_act = find_bw_act(db_ei, 'treatment of fuel cell stack, 1MWe, PEM', 'used fuel cell stack, 1MWe, PEM', 'RER')

            act.new_exchange(amount= -amount_elect_unit, input = add_act.key, type="technosphere", location = add_act['location'],
                             name = add_act['name'], product = add_act['reference product'] ).save()

            # 1.2. Add BoP electrolyzer
            add_act = find_bw_act(db_ei, 'electrolyzer production, 1MWe, PEM, Balance of Plant', 'electrolyzer, 1MWe, PEM, Balance of Plant', 'RER')

            amount_bop_unit = (lifetime_stack/lifetime_bop)/(1000*cf_electrolyzer*HOURS_YR*lifetime_stack/amount_electricity)  #lifetime bop is 20 years

            act.new_exchange(amount=amount_bop_unit, input = add_act.key, type="technosphere", location = add_act['location'],
                             name = add_act['name'], product = add_act['reference product'] ).save()
            
            add_act = find_bw_act(db_ei, 'treatment of fuel cell balance of plant, 1MWe, PEM', 'used fuel cell balance of plant, 1MWe, PEM', 'RER')

            act.new_exchange(amount= -amount_bop_unit, input = add_act.key, type="technosphere", location = add_act['location'],
                             name = add_act['name'], product = add_act['reference product'] ).save()

        elif electrolyzer == "aec":         
            add_act = find_bw_act(db_ei, 'electrolyzer production, 1MWe, AEC, Stack', 'electrolyzer, 1MWe, AEC, Stack', 'RER')

            amount_elect_unit = 1/(1000*cf_electrolyzer*HOURS_YR*lifetime_stack/amount_electricity)

            act.new_exchange(amount= amount_elect_unit, input = add_act.key, type="technosphere", location = add_act['location'],
                             name = add_act['name'], product = add_act['reference product'] ).save()

            add_act = find_bw_act(db_ei, 'treatment of fuel cell stack, 1MWe, AEC', 'used fuel cell stack, 1MWe, AEC', 'RER')

            act.new_exchange(amount= -amount_elect_unit, input = add_act.key, type="technosphere", location = add_act['location'],
                             name = add_act['name'], product = add_act['reference product'] ).save()

            # 1.2. Add BoP electrolyzer            
            add_act = find_bw_act(db_ei, 'electrolyzer production, 1MWe, AEC, Balance of Plant', 'electrolyzer, 1MWe, AEC, Balance of Plant', 'RER')

            amount_bop_unit = (lifetime_stack/lifetime_bop)/(1000*cf_electrolyzer*HOURS_YR*lifetime_stack/amount_electricity)  #lifetime bop is 20 years

            act.new_exchange(amount=amount_bop_unit, input = add_act.key, type="technosphere", location = add_act['location'],
                             name = add_act['name'], product = add_act['reference product'] ).save()

            add_act = find_bw_act(db_ei, 'treatment of fuel cell balance of plant, 1MWe, AEC', 'used fuel cell balance of plant, 1MWe, AEC', 'RER')

            act.new_exchange(amount= -amount_bop_unit, input = add_act.key, type="technosphere", location = add_act['location'],
                             name = add_act['name'], product = add_act['reference product'] ).save()

            # 1.2. Add potassium hydroxide
            add_act = find_bw_act(db_ei, 'market for potassium hydroxide', 'potassium hydroxide', 'GLO')

            act.new_exchange(amount=3.70*1e-3, input = add_act.key, type="technosphere", location = add_act['location'],
                             name = add_act['name'], product = add_act['reference product'] ).save()

        elif electrolyzer == 'soec':  
            add_act = find_bw_act(db_ei, 'electrolyzer production, 1MWe, SOEC, Stack', 'electrolyzer, 1MWe, SOEC, Stack', 'RER')

            amount_elect_unit = 1/(1000*cf_electrolyzer*HOURS_YR*lifetime_stack/amount_electricity)

            act.new_exchange(amount= amount_elect_unit, input = add_act.key, type="technosphere", location = add_act['location'],
                             name = add_act['name'], product = add_act['reference product'] ).save()

            add_act = find_bw_act(db_ei, 'treatment of fuel cell stack, 1MWe, SOEC', 'used fuel cell stack, 1MWe, SOEC', 'RER')

            act.new_exchange(amount= -amount_elect_unit, input = add_act.key, type="technosphere", location = add_act['location'],
                             name = add_act['name'], product = add_act['reference product'] ).save()

            # 1.2. Add BoP electrolyzer            
            add_act = find_bw_act(db_ei, 'electrolyzer production, 1MWe, SOEC, Balance of Plant', 'electrolyzer, 1MWe, SOEC, Balance of Plant', 'RER')

            amount_bop_unit = (lifetime_stack/lifetime_bop)/(1000*cf_electrolyzer*HOURS_YR*lifetime_stack/amount_electricity)  #lifetime bop is 20 years

            act.new_exchange(amount=amount_bop_unit, input = add_act.key, type="technosphere", location = add_act['location'],
                             name = add_act['name'], product = add_act['reference product'] ).save()

            add_act = find_bw_act(db_ei, 'treatment of fuel cell balance of plant, 1MWe, SOEC', 'used fuel cell balance of plant, 1MWe, SOEC', 'RER')

            act.new_exchange(amount= -amount_bop_unit, input = add_act.key, type="technosphere", location = add_act['location'],
                             name = add_act['name'], product = add_act['reference product'] ).save()
        else:
            raise ValueError("Not possible to user electrolyzer exchanges since electrolyzer '{}' doesn't exist in LCI".format(electrolyzer))

        # 1.3. Add direct water consumption        
        add_act = find_bw_act(db_ei, 'market for water, deionised', 'water, deionised', 'RoW')

        # Add the exchange
        act.new_exchange(input=add_act.key,amount=AMOUNT_WATER_ELECTROLYSIS,
                                    unit="kilogram",type='technosphere', location = add_act['location'],
                             name = add_act['name'], product = add_act['reference product']).save()

        # 1.4. Add infrastructure activity, PV unit
        add_act = find_bw_act(db_ei, 'market for lubricating oil', 'lubricating oil', 'RER')
        
        amount_kwh = 157.5/(capacity*cf_wind*HOURS_YR*(1-curtailment_wind))

        # If needed, add compression power demand to this initial power demandfor water electrolysis to arrive at 30 bar.
        # Here, we assume that the compression power demand is satisfied by local renewables
        amount_electricity_w_compr = amount_electricity + ELECTROLYZER_DATA[electrolyzer]['compression_kWh_per_kgH2']

        # Add the exchange
        act.new_exchange(input=add_act.key,amount=(amount_kwh*amount_electricity_w_compr),
                                        unit="kilogram",type='technosphere', location = add_act['location'],
                             name = add_act['name'], product = add_act['reference product']).save()

        # 1.5 Add infrastructure activities, lubricating oil waste mineral oil
        add_act = find_bw_act(db_ei, 'market for waste mineral oil', 'waste mineral oil', 'Europe without Switzerland')

        amount_kwh = 157.5/(capacity*cf_wind*HOURS_YR*(1-curtailment_wind))
        # Add the exchange
        act.new_exchange(input=add_act.key,amount=(-amount_kwh*amount_electricity_w_compr),
                                        unit="kilogram",type='technosphere', location = add_act['location'],
                             name = add_act['name'], product = add_act['reference product']).save()

        # 1.6 Add infrastructure activities, market for wind turbine network connection
        add_act = find_bw_act(db_ei, 'market for wind turbine network connection, 2MW, onshore', 'wind turbine network connection, 2MW, onshore', 'GLO')

        amount_kwh = 1/(capacity*lifetime_wind*cf_wind*HOURS_YR*(1-curtailment_wind))

        # Add the exchange
        act.new_exchange(input=add_act.key,amount=(amount_kwh*amount_electricity_w_compr),
                                        unit="unit",type='technosphere', location = add_act['location'],
                             name = add_act['name'], product = add_act['reference product']).save()

        # 1.7 Add infrastructure activities, market for wind turbine network connection
        add_act = find_bw_act(db_ei, 'market for wind turbine, 2MW, onshore', 'wind turbine, 2MW, onshore', 'GLO')

        amount_kwh = 1/(capacity*lifetime_wind*cf_wind*HOURS_YR*(1-curtailment_wind))

        # Add the exchange
        act.new_exchange(input=add_act.key,amount=(amount_kwh*amount_electricity_w_compr),
                                        unit="unit",type='technosphere', location = add_act['location'],
                             name = add_act['name'], product = add_act['reference product']).save()
        
        # 1.8 Add biosphere flows - kinteic wind
        bio_flow = search_biosphere_entries(BIOSPHERE_DB, 'Energy, kinetic (in wind), converted', ('natural resource', 'in air'), 'megajoule')
            
        amount_kwh = 3.87 # curtailment doesn't have influence here, as solar PV is switched off
        # Add the exchange
        act.new_exchange(input=bio_flow.key,amount=(amount_kwh*amount_electricity_w_compr),
                                        unit="megajoule",type='biosphere').save()

        # 1.8 Add biosphere flows - H2 Leakage
        bio_flow = search_biosphere_entries(BIOSPHERE_DB, 'Hydrogen', ('air',), 'kilogram')

        # Add the exchange
        act.new_exchange(input=bio_flow.key,amount=h2_leakage_factor,
                                        unit="kilogram",type='biosphere').save()
        
        # 1.9 Oxygen production
        bio_flow = search_biosphere_entries(BIOSPHERE_DB, 'Oxygen', ('air',), 'kilogram')

        # Add the exchange
        act.new_exchange(input=bio_flow.key,amount=8,
                                        unit="kilogram",type='biosphere').save()

        # Finally, add the biosphere flows for land occupation.
        # 1.10 Add biosphere flows - land occupation
        bio_flow = search_biosphere_entries(BIOSPHERE_DB, 'Occupation, industrial area', ('natural resource','land'), 'square meter-year')

        # Electrolyser land footprint: (land_m2_kw [m2/kW]*1000[kW/MW])/(H2_prod[kg H2 per lifetime]/lifetime_system[years])
        act.new_exchange(input=bio_flow.key,amount=((land_m2_kw*1000)/(1000*cf_electrolyzer*HOURS_YR/amount_electricity)),
                                        unit="square meter-year",type='biosphere').save()
        
        bio_flow = search_biosphere_entries(BIOSPHERE_DB, 'Transformation, from industrial area', ('natural resource','land'), 'square meter')   

        # Electrolyser land transformation: (land_m2_kw [m2/kW]*1000 [kW/MW])/H2_prod [kg H2]
        act.new_exchange(input=bio_flow.key,amount=(land_m2_kw*1000)/(1000*cf_electrolyzer*HOURS_YR*lifetime_bop/amount_electricity),
                                        unit="square meter",type='biosphere').save()
        
        bio_flow = search_biosphere_entries(BIOSPHERE_DB, 'Transformation, to industrial area', ('natural resource','land'), 'square meter')   
        # Electrolyser land transformation: (land_m2_kw [m2/kW]*1000 [kW/MW])/H2_prod [kg H2]
        act.new_exchange(input=bio_flow.key,amount=(land_m2_kw*1000)/(1000*cf_electrolyzer*HOURS_YR*lifetime_bop/amount_electricity),
                                        unit="square meter",type='biosphere').save()
        
        act.save()

def GENERATE_ACTS_WIND_OFF(db_ei, cf_wind, curtailment_wind, cf_electrolyzer=0.55, electrolyzer = "pem",excel_col_name = 'ecoinvent_310_reference'):
    """
    Generates a new activity in the Brightway2 database for hydrogen production
    from electrolysis using offshore wind energy.

    Parameters:
    db_ei (str): The name of the Brightway2 database to be used.
    cf_wind (float): capacity factor of onshore wind system.
    curtailment_wind (float): curtailment of onshore wind to account for oversizing.
    cf_electrolyzer (float): The capacity factor of the wind energy source, which is 0.55 according to IEA database.
    type of electrolyzer (str): Type of electrolyzer (can be 'pem', 'aec', or 'soec')

    Returns:
    None: The function performs operations directly on the Brightway2 database.

    Description:
    This function creates a new activity for hydrogen production using electrolysis powered by
    offshore wind energy. It calculates the necessary inputs and adds exchanges for the activity
    based on the specified capacity factor. The activity is only created if it does not already exist
    in the database.

    Steps performed:
    1. Check if the activity already exists in the database.
    2. Retrieve necessary technical and cost data.
    3. Create the new activity with the specified parameters.
    4. Add exchanges for:
        - Stack electrolyzer production
        - Balance of plant electrolyzer production
        - Water consumption
        - Lubricating oil and its waste treatment
        - Wind power plant network connection
        - Biosphere flows for kinetic wind energy, hydrogen leakage, and oxygen production.

    Includes compression provided by offshore wind, depending on the electrolyser: ELECTROLYZER_DATA[electrolyzer]['compression_kWh_per_kgH2']
    
    Example:
    >>> GENERATE_ACTS_wind_off('ecoinvent_3.7.1', 0.25)
    Skipped: activity 'hydrogen production, gaseous, 30 bar, from PEM electrolysis, offshore wind, global cf [0.25]' already generated!
    """

    capacity = 2000 #kWp
    new_name = "hydrogen production, gaseous, {} bar, from {} electrolysis, offshore wind, global cf [{}]".format(30,electrolyzer.upper(),round(cf_wind,3))
    check_act = [act for act in bw.Database(db_ei) if new_name == act['name']]
    
    tech = "offshore_wind"
    
    lifetime_wind = COST_DATA.loc[(tech,'lifetime')][excel_col_name]
    eff_elect = COST_DATA.loc[('electrolyzer_{}'.format(electrolyzer),'eff')][excel_col_name]
    lifetime_stack = COST_DATA.loc[('electrolyzer_{}'.format(electrolyzer),'stack_lt')][excel_col_name]
    lifetime_bop = COST_DATA.loc[('electrolyzer_{}'.format(electrolyzer),'bop_lt')][excel_col_name]    
    h2_leakage_factor = COST_DATA.loc[('h2_leakage','-')][excel_col_name] 
    land_m2_kw = COST_DATA.loc[('electrolyzer_{}'.format(electrolyzer),'land_m2_kw')][excel_col_name]  

    if len(check_act) > 0:
        print("Skipped: activity '{}' already generated!".format(new_name))
    else:
        # Generate activity
        db = bw.Database(db_ei)
        code_na = str(uuid.uuid4().hex)
        act = db.new_activity(
            **{
                'name': new_name,
                 "code": code_na,
                'unit': 'kilogram',
                'reference product': "hydrogen, gaseous, {} bar".format(30),
                'location' :"GLO",
                'production amount': 1.0,
                'comment': f"Hydrogen production via water electrolysis with a \
                    {electrolyzer} electrolyzer with a capacity factor of {round(cf_electrolyzer,2)} and onshore wind of capacity factor {round(cf_wind,2)}, includes compression if AEC or SOEC with local renewables"
            }
        )

        act.save()

        # Add production exchange
        act.new_exchange(**{
            'input': (db_ei, code_na),
            'amount': 1,
            'type': 'production',
        }).save()

        # calculate electricity requirements
        amount_electricity = (EN_H2/(3.6*eff_elect))

        # 1.1. Add stack electrolyzer
        if electrolyzer == "pem":
            add_act = find_bw_act(db_ei, 'electrolyzer production, 1MWe, PEM, Stack', 'electrolyzer, 1MWe, PEM, Stack', 'RER')

            amount_elect_unit = 1/(1000*cf_electrolyzer*HOURS_YR*lifetime_stack/amount_electricity)

            act.new_exchange(amount= amount_elect_unit, input = add_act.key, type="technosphere", location = add_act['location'],
                             name = add_act['name'], product = add_act['reference product'] ).save()
            
            add_act = find_bw_act(db_ei, 'treatment of fuel cell stack, 1MWe, PEM', 'used fuel cell stack, 1MWe, PEM', 'RER')

            act.new_exchange(amount= -amount_elect_unit, input = add_act.key, type="technosphere", location = add_act['location'],
                             name = add_act['name'], product = add_act['reference product'] ).save()

            # 1.2. Add BoP electrolyzer
            add_act = find_bw_act(db_ei, 'electrolyzer production, 1MWe, PEM, Balance of Plant', 'electrolyzer, 1MWe, PEM, Balance of Plant', 'RER')

            amount_bop_unit = (lifetime_stack/lifetime_bop)/(1000*cf_electrolyzer*HOURS_YR*lifetime_stack/amount_electricity)  #lifetime bop is 20 years

            act.new_exchange(amount=amount_bop_unit, input = add_act.key, type="technosphere", location = add_act['location'],
                             name = add_act['name'], product = add_act['reference product'] ).save()
            
            add_act = find_bw_act(db_ei, 'treatment of fuel cell balance of plant, 1MWe, PEM', 'used fuel cell balance of plant, 1MWe, PEM', 'RER')

            act.new_exchange(amount= -amount_bop_unit, input = add_act.key, type="technosphere", location = add_act['location'],
                             name = add_act['name'], product = add_act['reference product'] ).save()

        elif electrolyzer == "aec":         
            add_act = find_bw_act(db_ei, 'electrolyzer production, 1MWe, AEC, Stack', 'electrolyzer, 1MWe, AEC, Stack', 'RER')

            amount_elect_unit = 1/(1000*cf_electrolyzer*HOURS_YR*lifetime_stack/amount_electricity)

            act.new_exchange(amount= amount_elect_unit, input = add_act.key, type="technosphere", location = add_act['location'],
                             name = add_act['name'], product = add_act['reference product'] ).save()

            add_act = find_bw_act(db_ei, 'treatment of fuel cell stack, 1MWe, AEC', 'used fuel cell stack, 1MWe, AEC', 'RER')

            act.new_exchange(amount= -amount_elect_unit, input = add_act.key, type="technosphere", location = add_act['location'],
                             name = add_act['name'], product = add_act['reference product'] ).save()

            # 1.2. Add BoP electrolyzer            
            add_act = find_bw_act(db_ei, 'electrolyzer production, 1MWe, AEC, Balance of Plant', 'electrolyzer, 1MWe, AEC, Balance of Plant', 'RER')

            amount_bop_unit = (lifetime_stack/lifetime_bop)/(1000*cf_electrolyzer*HOURS_YR*lifetime_stack/amount_electricity)  #lifetime bop is 20 years

            act.new_exchange(amount=amount_bop_unit, input = add_act.key, type="technosphere", location = add_act['location'],
                             name = add_act['name'], product = add_act['reference product'] ).save()

            add_act = find_bw_act(db_ei, 'treatment of fuel cell balance of plant, 1MWe, AEC', 'used fuel cell balance of plant, 1MWe, AEC', 'RER')

            act.new_exchange(amount= -amount_bop_unit, input = add_act.key, type="technosphere", location = add_act['location'],
                             name = add_act['name'], product = add_act['reference product'] ).save()

            # 1.2. Add potassium hydroxide
            add_act = find_bw_act(db_ei, 'market for potassium hydroxide', 'potassium hydroxide', 'GLO')

            act.new_exchange(amount=3.70*1e-3, input = add_act.key, type="technosphere", location = add_act['location'],
                             name = add_act['name'], product = add_act['reference product'] ).save()

        elif electrolyzer == 'soec':  
            add_act = find_bw_act(db_ei, 'electrolyzer production, 1MWe, SOEC, Stack', 'electrolyzer, 1MWe, SOEC, Stack', 'RER')

            amount_elect_unit = 1/(1000*cf_electrolyzer*HOURS_YR*lifetime_stack/amount_electricity)

            act.new_exchange(amount= amount_elect_unit, input = add_act.key, type="technosphere", location = add_act['location'],
                             name = add_act['name'], product = add_act['reference product'] ).save()

            add_act = find_bw_act(db_ei, 'treatment of fuel cell stack, 1MWe, SOEC', 'used fuel cell stack, 1MWe, SOEC', 'RER')

            act.new_exchange(amount= -amount_elect_unit, input = add_act.key, type="technosphere", location = add_act['location'],
                             name = add_act['name'], product = add_act['reference product'] ).save()

            # 1.2. Add BoP electrolyzer            
            add_act = find_bw_act(db_ei, 'electrolyzer production, 1MWe, SOEC, Balance of Plant', 'electrolyzer, 1MWe, SOEC, Balance of Plant', 'RER')

            amount_bop_unit = (lifetime_stack/lifetime_bop)/(1000*cf_electrolyzer*HOURS_YR*lifetime_stack/amount_electricity)  #lifetime bop is 20 years

            act.new_exchange(amount=amount_bop_unit, input = add_act.key, type="technosphere", location = add_act['location'],
                             name = add_act['name'], product = add_act['reference product'] ).save()

            add_act = find_bw_act(db_ei, 'treatment of fuel cell balance of plant, 1MWe, SOEC', 'used fuel cell balance of plant, 1MWe, SOEC', 'RER')

            act.new_exchange(amount= -amount_bop_unit, input = add_act.key, type="technosphere", location = add_act['location'],
                             name = add_act['name'], product = add_act['reference product'] ).save()
        else:
            raise ValueError("Not possible to user electrolyzer exchanges since electrolyzer '{}' doesn't exist in LCI".format(electrolyzer))

        # 1.3. Add direct water consumption        
        add_act = find_bw_act(db_ei, 'market for water, deionised', 'water, deionised', 'RoW')

        # Add the exchange
        act.new_exchange(input=add_act.key,amount=AMOUNT_WATER_ELECTROLYSIS,
                                    unit="kilogram",type='technosphere', location = add_act['location'],
                             name = add_act['name'], product = add_act['reference product']).save()
        
        # 1.4. Add infrastructure activity
        add_act = find_bw_act(db_ei, 'market for lubricating oil', 'lubricating oil', 'RER')
        
        amount_kwh = 157.5/(capacity*cf_wind*HOURS_YR*(1-curtailment_wind))

        # If needed, add compression power demand to this initial power demandfor water electrolysis to arrive at 30 bar.
        # Here, we assume that the compression power demand is satisfied by local renewables
        amount_electricity_w_compr = amount_electricity + ELECTROLYZER_DATA[electrolyzer]['compression_kWh_per_kgH2']

        # Add the exchange
        act.new_exchange(input=add_act.key,amount=(amount_kwh*amount_electricity_w_compr),
                                        unit="kilogram",type='technosphere', location = add_act['location'],
                             name = add_act['name'], product = add_act['reference product']).save()

        # 1.5 Add infrastructure activities, lubricating oil waste mineral oil
        add_act = find_bw_act(db_ei, 'market for waste mineral oil', 'waste mineral oil', 'Europe without Switzerland')

        amount_kwh = 157.5/(capacity*cf_wind*HOURS_YR*(1-curtailment_wind))

        # 1.6 Add infrastructure activities, market for wind turbine network connection
        add_act = find_bw_act(db_ei, 'market for wind power plant, 2MW, offshore, fixed parts', 'wind power plant, 2MW, offshore, fixed parts', 'GLO')
        amount_kwh = 1/(capacity*lifetime_wind*cf_wind*HOURS_YR*(1-curtailment_wind))

        # Add the exchange
        act.new_exchange(input=add_act.key,amount=(amount_kwh*amount_electricity_w_compr),
                                        unit="unit",type='technosphere', location = add_act['location'],
                             name = add_act['name'], product = add_act['reference product']).save()

        # 1.7 Add infrastructure activities, market for wind turbine network connection
        add_act = find_bw_act(db_ei, 'market for wind power plant, 2MW, offshore, moving parts', 'wind power plant, 2MW, offshore, moving parts', 'GLO')

        amount_kwh = 1/(capacity*lifetime_wind*cf_wind*HOURS_YR*(1-curtailment_wind))
        # Add the exchange
        act.new_exchange(input=add_act.key,amount=(amount_kwh*amount_electricity_w_compr),
                                        unit="unit",type='technosphere', location = add_act['location'],
                             name = add_act['name'], product = add_act['reference product']).save()
        
        # 1.8 Add biosphere flows - kinteic wind
        bio_flow = search_biosphere_entries(BIOSPHERE_DB, 'Energy, kinetic (in wind), converted', ('natural resource', 'in air'), 'megajoule')
            
        amount_kwh = 3.87 # curtailment doesn't have influence here, as solar PV is switched off

        # Add the exchange
        act.new_exchange(input=bio_flow.key,amount=(amount_kwh*amount_electricity_w_compr),
                                        unit="megajoule",type='biosphere').save()

        # 1.9 Add biosphere flows - H2 Leakage
        bio_flow = search_biosphere_entries(BIOSPHERE_DB, 'Hydrogen', ('air',), 'kilogram')

        # Add the exchange
        act.new_exchange(input=bio_flow.key,amount=h2_leakage_factor,
                                        unit="kilogram",type='biosphere').save()
        
        # 1.10 Oxygen production
        bio_flow = search_biosphere_entries(BIOSPHERE_DB, 'Oxygen', ('air',), 'kilogram')

        # Add the exchange
        act.new_exchange(input=bio_flow.key,amount=8,
                                        unit="kilogram",type='biosphere').save()

        # Finally, add the biosphere flows for land occupation.
        # 1.11 Add biosphere flows - land occupation
        bio_flow = search_biosphere_entries(BIOSPHERE_DB, 'Occupation, industrial area', ('natural resource','land'), 'square meter-year')

        # Electrolyser land footprint: (land_m2_kw [m2/kW]*1000[kW/MW])/(H2_prod[kg H2 per lifetime]/lifetime_system[years])
        act.new_exchange(input=bio_flow.key,amount=((land_m2_kw*1000)/(1000*cf_electrolyzer*HOURS_YR/amount_electricity)),
                                        unit="square meter-year",type='biosphere').save()
        
        bio_flow = search_biosphere_entries(BIOSPHERE_DB, 'Transformation, from industrial area', ('natural resource','land'), 'square meter')   

        # Electrolyser land transformation: (land_m2_kw [m2/kW]*1000 [kW/MW])/H2_prod [kg H2]
        act.new_exchange(input=bio_flow.key,amount=(land_m2_kw*1000)/(1000*cf_electrolyzer*HOURS_YR*lifetime_bop/amount_electricity),
                                        unit="square meter",type='biosphere').save()
        
        bio_flow = search_biosphere_entries(BIOSPHERE_DB, 'Transformation, to industrial area', ('natural resource','land'), 'square meter')   
        # Electrolyser land transformation: (land_m2_kw [m2/kW]*1000 [kW/MW])/H2_prod [kg H2]
        act.new_exchange(input=bio_flow.key,amount=(land_m2_kw*1000)/(1000*cf_electrolyzer*HOURS_YR*lifetime_bop/amount_electricity),
                                        unit="square meter",type='biosphere').save()
        
        act.save()

def GENERATE_ACTS_GENERIC_ELECT(db_ei, cf_electrolyzer=0.57, electrolyzer = "pem",
                                generic_electr_act = ["electricity production, hydro, reservoir, alpine region",
                                                      'electricity, high voltage',
                                                      "RoW"], excel_col_name = 'ecoinvent_310_reference'
                                ):
    """
    Generates a new activity in the Brightway2 database for hydrogen production
    from electrolysis using power from a specified power source.

    Parameters
    ----------
    db_ei : str
        Name of the Brightway2 database.
    cf_electrolyzer : float, optional
        Capacity factor for the technology. Default is 0.57.
    electrolyzer : str, optional
        Type of electrolyzer to be used. Can be either "pem" or "aec". Default is "pem".
    generic_electr_act : list, optional
        List containing three elements that specify the electricity production activity:
        1. Name of the electricity production method.
        2. Voltage level.
        3. Region. Default is ["electricity production, hydro, reservoir, alpine region",
                               'electricity, high voltage', "RoW"].

    Returns
    -------
    None
        The function does not return any value. It generates a new activity in the Brightway2 
        database or prints a message if the activity already exists.

    Description
    -----------
    The function performs the following steps:

    1. Activity Name Generation:
       - Constructs a new activity name based on the type of electrolyzer and the electricity 
         production method.

    2. Activity Check:
       - Checks if the activity already exists in the database. If it does, the function prints a 
         message and skips the generation.

    3. Electrolyzer Efficiency and Lifetime:
       - Retrieves efficiency and lifetime values for the selected electrolyzer type from the 
         COST_DATA dataframe.

    4. Activity Creation:
       - Creates a new activity in the database with specified properties.

    5. Adding Exchanges:
       - Adds various exchanges related to the electrolyzer stack, balance of plant (BoP), water 
         consumption, power consumption, and biosphere flows for hydrogen leakage and oxygen.

    Includes compression provided by grid power, depending on the electrolyser: ELECTROLYZER_DATA[electrolyzer]['compression_kWh_per_kgH2']

    Example Usage
    -------------
    GENERATE_ACTS_GENERIC_ELECT("my_database", cf_electrolyzer=0.60, electrolyzer="aec",
                                generic_electr_act=["electricity production, wind, offshore", 
                                                    'electricity, high voltage', "US"])
    
    This will generate a new activity for hydrogen production using an AEC electrolyzer with power 
    from offshore wind electricity production in the US region.
    """

    # This is to make sure that the matching complies for activities based on grid power and having larger aggregated scales.
    if generic_electr_act[0] == 'market group for electricity, medium voltage':
        act_new_name = 'market for electricity, medium voltage'
    else:
        act_new_name = generic_electr_act[0]

    new_name = "hydrogen production, gaseous, {} bar, from {} electrolysis, power from {}".format(30,electrolyzer.upper(),act_new_name)
    check_act = [act for act in bw.Database(db_ei) if new_name == act['name'] and generic_electr_act[2] == act['location']]
    
    eff_elect = COST_DATA.loc[('electrolyzer_{}'.format(electrolyzer),'eff')][excel_col_name]
    lifetime_stack = COST_DATA.loc[('electrolyzer_{}'.format(electrolyzer),'stack_lt')][excel_col_name]  
    lifetime_bop = COST_DATA.loc[('electrolyzer_{}'.format(electrolyzer),'bop_lt')][excel_col_name]  
    h2_leakage_factor = COST_DATA.loc[('h2_leakage','-')][excel_col_name] 
    land_m2_kw = COST_DATA.loc[('electrolyzer_{}'.format(electrolyzer),'land_m2_kw')][excel_col_name]  

    if len(check_act) > 0:
        print("Skipped: activity '{}' already generated!".format(new_name))
    else:
        # Generate activity
        db = bw.Database(db_ei)
        code_na = str(uuid.uuid4().hex)
        act = db.new_activity(
            **{
                'name': new_name,
                 "code": code_na,
                'unit': 'kilogram',
                'reference product': "hydrogen, gaseous, {} bar".format(30),
                'location' :generic_electr_act[2],
                'production amount': 1.0,
                'comment': "Hydrogen production using grid power, includes compression if AEC or SOEC with local renewables"
            }
        )

        act.save()

        # Add production exchange
        act.new_exchange(**{
            'input': (db_ei, code_na),
            'amount': 1,
            'type': 'production',
        }).save()

        # calculate electricity requirements
        amount_electricity = (EN_H2/(3.6*eff_elect))

        # 1.1. Add stack electrolyzer
        if electrolyzer == "pem":
            add_act = find_bw_act(db_ei, 'electrolyzer production, 1MWe, PEM, Stack', 'electrolyzer, 1MWe, PEM, Stack', 'RER')

            amount_elect_unit = 1/(1000*cf_electrolyzer*HOURS_YR*lifetime_stack/amount_electricity)

            act.new_exchange(amount= amount_elect_unit, input = add_act.key, type="technosphere", location = add_act['location'],
                             name = add_act['name'], product = add_act['reference product'] ).save()
            
            add_act = find_bw_act(db_ei, 'treatment of fuel cell stack, 1MWe, PEM', 'used fuel cell stack, 1MWe, PEM', 'RER')

            act.new_exchange(amount= -amount_elect_unit, input = add_act.key, type="technosphere", location = add_act['location'],
                             name = add_act['name'], product = add_act['reference product'] ).save()

            # 1.2. Add BoP electrolyzer
            add_act = find_bw_act(db_ei, 'electrolyzer production, 1MWe, PEM, Balance of Plant', 'electrolyzer, 1MWe, PEM, Balance of Plant', 'RER')

            amount_bop_unit = (lifetime_stack/lifetime_bop)/(1000*cf_electrolyzer*HOURS_YR*lifetime_stack/amount_electricity)  #lifetime bop is 20 years

            act.new_exchange(amount=amount_bop_unit, input = add_act.key, type="technosphere", location = add_act['location'],
                             name = add_act['name'], product = add_act['reference product'] ).save()
            
            add_act = find_bw_act(db_ei, 'treatment of fuel cell balance of plant, 1MWe, PEM', 'used fuel cell balance of plant, 1MWe, PEM', 'RER')

            act.new_exchange(amount= -amount_bop_unit, input = add_act.key, type="technosphere", location = add_act['location'],
                             name = add_act['name'], product = add_act['reference product'] ).save()

        elif electrolyzer == "aec":         
            add_act = find_bw_act(db_ei, 'electrolyzer production, 1MWe, AEC, Stack', 'electrolyzer, 1MWe, AEC, Stack', 'RER')

            amount_elect_unit = 1/(1000*cf_electrolyzer*HOURS_YR*lifetime_stack/amount_electricity)

            act.new_exchange(amount= amount_elect_unit, input = add_act.key, type="technosphere", location = add_act['location'],
                             name = add_act['name'], product = add_act['reference product'] ).save()

            add_act = find_bw_act(db_ei, 'treatment of fuel cell stack, 1MWe, AEC', 'used fuel cell stack, 1MWe, AEC', 'RER')

            act.new_exchange(amount= -amount_elect_unit, input = add_act.key, type="technosphere", location = add_act['location'],
                             name = add_act['name'], product = add_act['reference product'] ).save()

            # 1.2. Add BoP electrolyzer            
            add_act = find_bw_act(db_ei, 'electrolyzer production, 1MWe, AEC, Balance of Plant', 'electrolyzer, 1MWe, AEC, Balance of Plant', 'RER')

            amount_bop_unit = (lifetime_stack/lifetime_bop)/(1000*cf_electrolyzer*HOURS_YR*lifetime_stack/amount_electricity)  #lifetime bop is 20 years

            act.new_exchange(amount=amount_bop_unit, input = add_act.key, type="technosphere", location = add_act['location'],
                             name = add_act['name'], product = add_act['reference product'] ).save()

            add_act = find_bw_act(db_ei, 'treatment of fuel cell balance of plant, 1MWe, AEC', 'used fuel cell balance of plant, 1MWe, AEC', 'RER')

            act.new_exchange(amount= -amount_bop_unit, input = add_act.key, type="technosphere", location = add_act['location'],
                             name = add_act['name'], product = add_act['reference product'] ).save()

            # 1.2. Add potassium hydroxide
            add_act = find_bw_act(db_ei, 'market for potassium hydroxide', 'potassium hydroxide', 'GLO')

            act.new_exchange(amount=3.70*1e-3, input = add_act.key, type="technosphere", location = add_act['location'],
                             name = add_act['name'], product = add_act['reference product'] ).save()

        elif electrolyzer == 'soec':  
            add_act = find_bw_act(db_ei, 'electrolyzer production, 1MWe, SOEC, Stack', 'electrolyzer, 1MWe, SOEC, Stack', 'RER')

            amount_elect_unit = 1/(1000*cf_electrolyzer*HOURS_YR*lifetime_stack/amount_electricity)

            act.new_exchange(amount= amount_elect_unit, input = add_act.key, type="technosphere", location = add_act['location'],
                             name = add_act['name'], product = add_act['reference product'] ).save()

            add_act = find_bw_act(db_ei, 'treatment of fuel cell stack, 1MWe, SOEC', 'used fuel cell stack, 1MWe, SOEC', 'RER')

            act.new_exchange(amount= -amount_elect_unit, input = add_act.key, type="technosphere", location = add_act['location'],
                             name = add_act['name'], product = add_act['reference product'] ).save()

            # 1.2. Add BoP electrolyzer            
            add_act = find_bw_act(db_ei, 'electrolyzer production, 1MWe, SOEC, Balance of Plant', 'electrolyzer, 1MWe, SOEC, Balance of Plant', 'RER')

            amount_bop_unit = (lifetime_stack/lifetime_bop)/(1000*cf_electrolyzer*HOURS_YR*lifetime_stack/amount_electricity)  #lifetime bop is 20 years

            act.new_exchange(amount=amount_bop_unit, input = add_act.key, type="technosphere", location = add_act['location'],
                             name = add_act['name'], product = add_act['reference product'] ).save()

            add_act = find_bw_act(db_ei, 'treatment of fuel cell balance of plant, 1MWe, SOEC', 'used fuel cell balance of plant, 1MWe, SOEC', 'RER')

            act.new_exchange(amount= -amount_bop_unit, input = add_act.key, type="technosphere", location = add_act['location'],
                             name = add_act['name'], product = add_act['reference product'] ).save()
        else:
            raise ValueError("Not possible to user electrolyzer exchanges since electrolyzer '{}' doesn't exist in LCI".format(electrolyzer))

        # 1.3. Add direct water consumption        
        add_act = find_bw_act(db_ei, 'market for water, deionised', 'water, deionised', 'RoW')

        # Add the exchange
        act.new_exchange(input=add_act.key,amount=AMOUNT_WATER_ELECTROLYSIS,
                                    unit="kilogram",type='technosphere', location = add_act['location'],
                             name = add_act['name'], product = add_act['reference product']).save()
        
        # If needed, add compression power demand to this initial power demandfor water electrolysis to arrive at 30 bar.
        # Here, we assume that the compression power demand is satisfied by local renewables
        amount_electricity_w_compr = amount_electricity + ELECTROLYZER_DATA[electrolyzer]['compression_kWh_per_kgH2']
        
        # Find the power activity as given in the function
        add_act = find_bw_act(db_ei, generic_electr_act[0],
                                  generic_electr_act[1], 
                                  generic_electr_act[2])

        # Add the exchange
        act.new_exchange(input=add_act.key,amount=amount_electricity_w_compr,
                                    type='technosphere', location = add_act['location'],
                             name = add_act['name'], product = add_act['reference product']).save()

        # 1.4 Add biosphere flows - H2 Leakage
        bio_flow = search_biosphere_entries(BIOSPHERE_DB, 'Hydrogen', ('air',), 'kilogram')

        # Add the exchange
        act.new_exchange(input=bio_flow.key,amount=h2_leakage_factor,
                                        unit="kilogram",type='biosphere').save()
        
        # 1.6 Oxygen production
        bio_flow = search_biosphere_entries(BIOSPHERE_DB, 'Oxygen', ('air',), 'kilogram')

        # Add the exchange
        act.new_exchange(input=bio_flow.key,amount=8,
                                        unit="kilogram",type='biosphere').save()

        # Finally, add the biosphere flows for land occupation.
        # 1.7 Add biosphere flows - land occupation
        bio_flow = search_biosphere_entries(BIOSPHERE_DB, 'Occupation, industrial area', ('natural resource','land'), 'square meter-year')

        # Electrolyser land footprint: (land_m2_kw [m2/kW]*1000[kW/MW])/(H2_prod[kg H2 per lifetime]/lifetime_system[years])
        act.new_exchange(input=bio_flow.key,amount=((land_m2_kw*1000)/(1000*cf_electrolyzer*HOURS_YR/amount_electricity)),
                                        unit="square meter-year",type='biosphere').save()
        
        bio_flow = search_biosphere_entries(BIOSPHERE_DB, 'Transformation, from industrial area', ('natural resource','land'), 'square meter')   

        # Electrolyser land transformation: (land_m2_kw [m2/kW]*1000 [kW/MW])/H2_prod [kg H2]
        act.new_exchange(input=bio_flow.key,amount=(land_m2_kw*1000)/(1000*cf_electrolyzer*HOURS_YR*lifetime_bop/amount_electricity),
                                        unit="square meter",type='biosphere').save()
        
        bio_flow = search_biosphere_entries(BIOSPHERE_DB, 'Transformation, to industrial area', ('natural resource','land'), 'square meter')   
        # Electrolyser land transformation: (land_m2_kw [m2/kW]*1000 [kW/MW])/H2_prod [kg H2]
        act.new_exchange(input=bio_flow.key,amount=(land_m2_kw*1000)/(1000*cf_electrolyzer*HOURS_YR*lifetime_bop/amount_electricity),
                                        unit="square meter",type='biosphere').save()
        
        act.save()


def create_electrolysis_exch_wind(amount_h2, cf_wind, curtailment_wind, cf_electrolyzer=0.4, 
                                      electrolyzer="pem", sel_db = 'ecoinvent_310_reference'):
    
    tech = "onshore_wind"  
    lifetime_wind = COST_DATA.loc[(tech,'lifetime')][sel_db]
    eff_elect = COST_DATA.loc[('electrolyzer_{}'.format(electrolyzer.lower()),'eff')][sel_db]
    lifetime_stack = COST_DATA.loc[('electrolyzer_{}'.format(electrolyzer.lower()),'stack_lt')][sel_db]  
    lifetime_bop = COST_DATA.loc[('electrolyzer_{}'.format(electrolyzer.lower()),'bop_lt')][sel_db]  
    h2_leakage_factor = COST_DATA.loc[('h2_leakage','-')][sel_db] 
    land_m2_kw = COST_DATA.loc[('electrolyzer_{}'.format(electrolyzer.lower()),'land_m2_kw')][sel_db]  
    
    # Prepare exchanges as dictionaries
    exchanges = []

    # Calculate electricity requirements
    capacity = 2000 #kWp
    amount_electricity = EN_H2 / (3.6 * eff_elect)

    # Stack electrolyzer type
    if electrolyzer.lower() == "pem":
        stack = "PEM"
    elif electrolyzer.lower() == "aec":
        stack = "AEC"
    elif electrolyzer.lower() == "soec":
        stack = "SOEC"
    else:
        raise ValueError(f"Electrolyzer '{electrolyzer}' not available in LCI")

    # Electrolyzer production exchange
    amount_elect_unit = 1 / (1000 * cf_electrolyzer * HOURS_YR * lifetime_stack / amount_electricity)

    exchanges.append({
        'database': sel_db,
        'amount': amount_elect_unit * amount_h2,
        'location': "RER",
        'type': "technosphere",
        'name': f'electrolyzer production, 1MWe, {stack}, Stack',
        'reference product': f'electrolyzer, 1MWe, {stack}, Stack'
    })

    # End-of-life treatment for fuel cell stack
    exchanges.append({
        'database': sel_db,
        'amount': -amount_elect_unit * amount_h2,
        'location': "RER",
        'type': "technosphere",
        'name': f'treatment of fuel cell stack, 1MWe, {stack}',
        'reference product': f'used fuel cell stack, 1MWe, {stack}'
    })

    # BOP production exchange  
    amount_bop_unit = (lifetime_stack/lifetime_bop)/(1000*cf_electrolyzer*HOURS_YR*lifetime_stack/amount_electricity)
    exchanges.append({
        'database': sel_db,
        'amount': amount_bop_unit * amount_h2,
        'location': "RER",
        'type': "technosphere",
        'name': f'electrolyzer production, 1MWe, {stack}, Balance of Plant',
        'reference product': f'electrolyzer, 1MWe, {stack}, Balance of Plant'
    })

    # End-of-life treatment for BOP
    exchanges.append({
        'database': sel_db,
        'amount': -amount_bop_unit * amount_h2,
        'location': "RER",
        'type': "technosphere",
        'name': f'treatment of fuel cell balance of plant, 1MWe, {stack}',
        'reference product': f'used fuel cell balance of plant, 1MWe, {stack}'
    })


    # Water requirement for electrolysis
    exchanges.append({
        'database': sel_db,
        'amount': AMOUNT_WATER_ELECTROLYSIS * amount_h2,
        'unit': "kilogram",
        'type': 'technosphere',
        'name': 'market for water, deionised',
        'reference product': 'water, deionised',
        'location': "RoW",
    })

    # Hydrogen leakage biosphere flow
    exchanges.append({
        'database': BIOSPHERE_DB,
       'categories': 'Hydrogen::air',
        'amount': h2_leakage_factor * amount_h2,
        'unit': "kilogram",
        'type': 'biosphere',
    })

    # Oxygen production biosphere flow
    exchanges.append({
        'name': 'Oxygen',
        'categories': 'Oxygen::air',
        'database': BIOSPHERE_DB,
        'amount': 8 * amount_h2,
        'unit': "kilogram",
        'type': 'biosphere',
    })

    # Add additional exchanges for wind infrastructure
    amount_wind = 1/(capacity*lifetime_wind*cf_wind*HOURS_YR*(1-curtailment_wind))

    exchanges.append({
        'database': sel_db,
        'amount': amount_wind * amount_electricity * amount_h2,
        'location': "GLO",
        'type': "technosphere",
        'name': 'market for wind turbine, 2MW, onshore',
        'reference product': 'wind turbine, 2MW, onshore'
    })

    exchanges.append({
        'database': sel_db,
        'amount': amount_wind * amount_electricity * amount_h2,
        'location': "GLO",
        'type': "technosphere",
        'name': 'market for wind turbine network connection, 2MW, onshore',
        'reference product': 'wind turbine network connection, 2MW, onshore'
    })

    # Add additional exchanges for oil neded
    amount_oil = 157.5/(capacity*cf_wind*HOURS_YR*(1-curtailment_wind))

    exchanges.append({
        'database': sel_db,
        'amount': amount_oil * amount_electricity * amount_h2,
        'location': "RER",
        'type': "technosphere",
        'name': 'market for lubricating oil',
        'reference product': 'lubricating oil'
    })

    exchanges.append({
        'database': sel_db,
        'amount': -(amount_oil * amount_electricity * amount_h2),
        'location': "Europe without Switzerland",
        'type': "technosphere",
        'name': 'market for waste mineral oil',
        'reference product': 'waste mineral oil'
    })

    exchanges.append({
        'database': BIOSPHERE_DB,
        'amount': amount_h2 * (3.87 * amount_electricity), # curtailment doesn't have influence here, as wind is switched off
        'unit': "megajoule",
        'type': 'biosphere',
        'name': 'Energy, kinetic (in wind), converted',
        'categories': ('natural resource', 'in air')
    })

    # Add biosphere flow for land occupation
    exchanges.append({
        'database': BIOSPHERE_DB,
        'amount': amount_h2 * ((land_m2_kw * 1000) / (1000 * cf_electrolyzer * HOURS_YR / amount_electricity)),
        'unit': "square meter-year",
        'type': 'biosphere',
        'name': 'Occupation, industrial area',
        'categories': ('natural resource','land'),
    })

    # Add biosphere flow for transformation from industrial area
    exchanges.append({
        'database': BIOSPHERE_DB,
        'amount': amount_h2 * (land_m2_kw * 1000) / (1000 * cf_electrolyzer * HOURS_YR * lifetime_bop / amount_electricity),
        'unit': "square meter",
        'type': 'biosphere',
        'name': 'Transformation, from industrial area',
        'categories': ('natural resource','land'),
    })

    # Add biosphere flow for transformation to industrial area
    exchanges.append({
        'database': BIOSPHERE_DB,
        'amount': amount_h2 * (land_m2_kw * 1000) / (1000 * cf_electrolyzer * HOURS_YR * lifetime_bop / amount_electricity),
        'unit': "square meter",
        'type': 'biosphere',
        'name': 'Transformation, to industrial area',
        'categories': ('natural resource','land'),
    })

    return exchanges

def create_electrolysis_exch_solar_pv(amount_h2, cf_pv, curtailment_pv, cf_electrolyzer=0.3, 
                                      electrolyzer="pem", sel_db = 'ecoinvent_310_reference'):
    
    tech = "solar_pv_gm"   
    lifetime_pv = COST_DATA.loc[(tech,'lifetime')][sel_db].item()
    eff_elect = COST_DATA.loc[('electrolyzer_{}'.format(electrolyzer.lower()),'eff')][sel_db]
    lifetime_stack = COST_DATA.loc[('electrolyzer_{}'.format(electrolyzer.lower()),'stack_lt')][sel_db]  
    lifetime_bop = COST_DATA.loc[('electrolyzer_{}'.format(electrolyzer.lower()),'bop_lt')][sel_db]  
    h2_leakage_factor = COST_DATA.loc[('h2_leakage','-')][sel_db] 
    land_m2_kw = COST_DATA.loc[('electrolyzer_{}'.format(electrolyzer.lower()),'land_m2_kw')][sel_db]   
    
    # Prepare exchanges as dictionaries
    exchanges = []

    # Calculate electricity requirements
    capacity = 560 * 0.895 #kWp, 10.5% average degradation, described in activity
    amount_electricity = EN_H2 / (3.6 * eff_elect)

    # Stack electrolyzer type
    if electrolyzer.lower() == "pem":
        stack = "PEM"
    elif electrolyzer.lower() == "aec":
        stack = "AEC"
    elif electrolyzer.lower() == "soec":
        stack = "SOEC"
    else:
        raise ValueError(f"Electrolyzer '{electrolyzer}' not available in LCI")

    # Electrolyzer production exchange
    amount_elect_unit = 1 / (1000 * cf_electrolyzer * HOURS_YR * lifetime_stack / amount_electricity)
    exchanges.append({
        'database': sel_db,
        'amount': amount_elect_unit * amount_h2,
        'location': "RER",
        'type': "technosphere",
        'name': f'electrolyzer production, 1MWe, {stack}, Stack',
        'reference product': f'electrolyzer, 1MWe, {stack}, Stack'
    })

    # End-of-life treatment for fuel cell stack
    exchanges.append({
        'database': sel_db,
        'amount': -amount_elect_unit * amount_h2,
        'location': "RER",
        'type': "technosphere",
        'name': f'treatment of fuel cell stack, 1MWe, {stack}',
        'reference product': f'used fuel cell stack, 1MWe, {stack}'
    })

    # BOP production exchange  
    amount_bop_unit = (lifetime_stack/lifetime_bop)/(1000*cf_electrolyzer*HOURS_YR*lifetime_stack/amount_electricity)
    exchanges.append({
        'database': sel_db,
        'amount': amount_bop_unit * amount_h2,
        'location': "RER",
        'type': "technosphere",
        'name': f'electrolyzer production, 1MWe, {stack}, Balance of Plant',
        'reference product': f'electrolyzer, 1MWe, {stack}, Balance of Plant'
    })

    # End-of-life treatment for BOP
    exchanges.append({
        'database': sel_db,
        'amount': -amount_bop_unit * amount_h2,
        'location': "RER",
        'type': "technosphere",
        'name': f'treatment of fuel cell balance of plant, 1MWe, {stack}',
        'reference product': f'used fuel cell balance of plant, 1MWe, {stack}'
    })


    # Water requirement for electrolysis
    exchanges.append({
        'database': sel_db,
        'amount': AMOUNT_WATER_ELECTROLYSIS * amount_h2,
        'unit': "kilogram",
        'type': 'technosphere',
        'location': "RoW",
        'name': 'market for water, deionised',
        'reference product': 'water, deionised'
    })

    # Hydrogen leakage biosphere flow
    exchanges.append({
        'database': BIOSPHERE_DB,
       'categories': 'Hydrogen::air',
        'amount': h2_leakage_factor * amount_h2,
        'unit': "kilogram",
        'type': 'biosphere',
    })

    # Oxygen production biosphere flow
    exchanges.append({
        'name': 'Oxygen',
        'categories': 'Oxygen::air',
        'database': BIOSPHERE_DB,
        'amount': 8 * amount_h2,
        'unit': "kilogram",
        'type': 'biosphere',
    })

    # Add additional exchanges for PV infrastructure and water usage
    amount_pv = 1 / (capacity * cf_pv * HOURS_YR * (1 - curtailment_pv) * lifetime_pv)

    exchanges.append({
        'database': sel_db,
        'amount': amount_pv * amount_electricity * amount_h2,
        'type': "technosphere",
        'name': 'electricity production, photovoltaic, at 560 kWp open ground, single-Si',
        'location': "CH",
        'reference product': 'photovoltaic open ground installation, 560 kWp, single-Si, on open ground'
    })

    # Add direct water consumption
    kg_water_unit= 2871.8 * 20 #20 L per m2 module, 3261.7 square meter for photovoltaic panel, multi-Si, at regional storage
    water_amount = kg_water_unit*(1/(capacity*lifetime_pv*cf_pv*HOURS_YR*(1-curtailment_pv)))
    exchanges.append({
        'database': sel_db,
        'amount': water_amount * amount_electricity * amount_h2,  # Adjusted for hydrogen production scale
        'type': "technosphere",
        'name': 'market for tap water',
        'location': 'Europe without Switzerland',
        'unit': 'kilogram'
    })

    # Wastewater treatment exchange
    exchanges.append({
        'database': sel_db,
        'amount': -(water_amount * amount_electricity * amount_h2) / 1000,  # Converted to cubic meters
        'unit': 'cubic meter',
        'type': 'technosphere',
        'name': 'treatment of wastewater, average, wastewater treatment',
        'location':'RoW',
        'reference product': 'wastewater, average'
    })

    # Waste heat biosphere flow
    exchanges.append({
        'database': BIOSPHERE_DB,
        'amount': 0.25027 * amount_electricity * amount_h2,
        'unit': "megajoule",
        'type': 'biosphere',
        'name': 'Heat, waste',
        'categories': ('air', 'urban air close to ground')
    })

    exchanges.append({
        'database': BIOSPHERE_DB,
        'amount': amount_h2 * (3.8503 * amount_electricity), # curtailment doesn't have influence here, as solar PV is switched off
        'unit': "megajoule",
        'type': 'biosphere',
        'name': 'Energy, solar, converted',
        'categories': ('natural resource', 'in air')
    })

    # Add biosphere flow for land occupation
    exchanges.append({
        'database': BIOSPHERE_DB,
        'amount': amount_h2 * ((land_m2_kw * 1000) / (1000 * cf_electrolyzer * HOURS_YR / amount_electricity)),
        'unit': "square meter-year",
        'type': 'biosphere',
        'name': 'Occupation, industrial area',
        'categories': ('natural resource','land'),
    })

    # Add biosphere flow for transformation from industrial area
    exchanges.append({
        'database': BIOSPHERE_DB,
        'amount': amount_h2 * (land_m2_kw * 1000) / (1000 * cf_electrolyzer * HOURS_YR * lifetime_bop / amount_electricity),
        'unit': "square meter",
        'type': 'biosphere',
        'name': 'Transformation, from industrial area',
        'categories': ('natural resource','land'),
    })

    # Add biosphere flow for transformation to industrial area
    exchanges.append({
        'database': BIOSPHERE_DB,
        'amount': amount_h2 * (land_m2_kw * 1000) / (1000 * cf_electrolyzer * HOURS_YR * lifetime_bop / amount_electricity),
        'unit': "square meter",
        'type': 'biosphere',
        'name': 'Transformation, to industrial area',
        'categories': ('natural resource','land'),
    })

    return exchanges