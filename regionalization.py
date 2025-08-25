from constructive_geometries import Geomatcher, resolved_row
import json
from copy import deepcopy
from numbers import Number
import math

# NOTE: many functions are slightly modified but mainly copied from the wurst package (Mutel et al.: https://github.com/polca/wurst/tree/main).
geomatcher = Geomatcher(backwards_compatible=True)

# add definitions from REMIND, if use another IAM for prospective LCA, then we need to change this.
with open(r"data\iam_variables_mapping\topologies\remind-topology.json", "r", encoding="utf-8") as f:
    REMIND_TOPOLOGY = json.load(f)
with open(r"data\iam_variables_mapping\topologies\image-topology.json", "r", encoding="utf-8") as f:
    IMAGE_TOPOLOGY = json.load(f)

def manually_add_definitions(geomatcher, data: dict, namespace: str = "", relative: bool = True):
    """
    Manually add topology definitions to a Geomatcher instance.

    This function replicates the logic of `geomatcher.add_definitions`, with one key difference:
    location keys are added without the IAM model namespace (i.e., just `k` instead of `(model, k)`).

    Parameters:
    - geomatcher (Geomatcher): The Geomatcher instance to modify.
    - data (dict): A dictionary of new topology definitions. If `relative` is True, 
      values should be lists of existing location names. If False, values should be sets of face IDs.
    - namespace (str, optional): Ignored in this function, included for compatibility. Defaults to "".
    - relative (bool, optional): If True, defines new locations by combining existing ones.
      If False, adds raw face ID mappings. Defaults to True.
    """
    if not relative:
        # Direct mapping of keys to face IDs
        geomatcher.topology.update({(namespace, k): v for k, v in data.items()})
        geomatcher.faces.update(set.union(*data.values()))
    else:
        # Definitions relative to existing entries in Geomatcher
        geomatcher.topology.update({
            k: set.union(*[geomatcher[o] for o in v])
            for k, v in data.items()
        })

print("Adding REMIND and IMAGE topology using regionalization.py")
manually_add_definitions(geomatcher, REMIND_TOPOLOGY, namespace="", relative=True)
manually_add_definitions(geomatcher, IMAGE_TOPOLOGY, namespace="", relative=True)

def match_best_location(location, possible_locations, exclusive=True, biggest_first=False, contained=False, fallback='GLO',
                     geomatcher=geomatcher):
    """
    Find the best-matching location from a list of possible locations using spatial rules.

    Parameters:
    - location (str or geometry): The target location to match, either as a string name or geometry object.
    - possible_locations (iterable): Collection of candidate locations to match against.
    - exclusive (bool, optional): If True, excludes overlapping matches (i.e. enforces mutually exclusive regions). Defaults to True.
    - biggest_first (bool, optional): If True, gives preference to larger matching regions. Defaults to False.
    - contained (bool, optional): If True, match only if the target location is fully contained within a candidate.
      If False, match based on any intersection. Defaults to False.
    - fallback (str, optional): Value to return if no match is found. Defaults to 'GLO'.
    - geomatcher (Geomatcher, optional): An instance of Geomatcher used for spatial comparisons.

    Returns:
    - list: A list of matching location(s), or a list containing the fallback if no match is found.
    """
    
    with resolved_row(possible_locations, geomatcher) as g:
        func = g.contained if contained else g.intersects
        locs = func(
            location,
            include_self=True,
            exclusive=exclusive,
            biggest_first=biggest_first,
            only=possible_locations)

        return locs if len(locs)>0 else [fallback] #Fallback is global dataset

# NOTE: This is copied from wurst package (Mutel et al.: https://github.com/polca/wurst/tree/main), full credits to that, however, we need to slightly change it to:
# ensure that RER datasets are also matched since some premise datasets are only in 'RER' location:
#if not kept and "RER" in possible_locations:
#    kept = [obj for obj in possible_datasets if obj["location"] == "RER"]
#    print(f"'RER' used for regionalization, exchange: {exc['name']}, {exc['location']}, dataset loc: {ds['location']}")
#if not kept and "EUR" in possible_locations:
#    kept = [obj for obj in possible_datasets if obj["location"] == "EUR"]
#    print(f"'EUR' used for regionalization, exchange: {exc['name']}, {exc['location']}, dataset loc: {ds['location']}")

def relink_technosphere_exchanges(
    ds, data, exclusive=True, drop_invalid=False, biggest_first=False, contained=True
):
    """Find new technosphere providers based on the location of the dataset.

    Designed to be used when the dataset's location changes, or when new datasets are added.

    Uses the name, reference product, and unit of the exchange to filter possible inputs. These must match exactly. Searches in the list of datasets ``data``.

    Will only search for providers contained within the location of ``ds``, unless ``contained`` is set to ``False``, all providers whose location intersects the location of ``ds`` will be used.

    A ``RoW`` provider will be added if there is a single topological face in the location of ``ds`` which isn't covered by the location of any providing activity.

    If no providers can be found, `relink_technosphere_exchanes` will try to add a `RoW` or `GLO` providers, in that order, if available. If there are still no valid providers, a ``InvalidLink`` exception is raised, unless ``drop_invalid`` is ``True``, in which case the exchange will be deleted.

    Allocation between providers is done using ``allocate_inputs``; results seem strange if ``contained=False``, as production volumes for large regions would be used as allocation factors.

    Input arguments:

        * ``ds``: The dataset whose technosphere exchanges will be modified.
        * ``data``: The list of datasets to search for technosphere product providers.
        * ``exclusive``: Bool, default is ``True``. Don't allow overlapping locations in input providers.
        * ``drop_invalid``: Bool, default is ``False``. Delete exchanges for which no valid provider is available.
        * ``biggest_first``: Bool, default is ``False``. Determines search order when selecting provider locations. Only relevant is ``exclusive`` is ``True``.
        * ``contained``: Bool, default is ``True``. If ture, only use providers whose location is completely within the ``ds`` location; otherwise use all intersecting locations.

    Modifies the dataset in place; returns the modified dataset."""
    new_exchanges = []
    technosphere = lambda x: x["type"] == "technosphere"

    for exc in filter(technosphere, ds["exchanges"]):
        possible_datasets = list(get_possibles(exc, data))
        possible_locations = [obj["location"] for obj in possible_datasets]

        with resolved_row(possible_locations, geomatcher) as g:
            func = g.contained if contained else g.intersects
            gis_match = func(
                ds["location"],
                include_self=True,
                exclusive=exclusive,
                biggest_first=biggest_first,
                only=possible_locations,
            )

        kept = [
            ds for loc in gis_match for ds in possible_datasets if ds["location"] == loc
        ]

        if kept:
            missing_faces = geomatcher[ds["location"]].difference(
                set.union(*[geomatcher[obj["location"]] for obj in kept])
            )
            if missing_faces and "RoW" in possible_locations:
                kept.extend(
                    [obj for obj in possible_datasets if obj["location"] == "RoW"]
                )
        elif "RoW" in possible_locations:
            kept = [obj for obj in possible_datasets if obj["location"] == "RoW"]

        if not kept and "GLO" in possible_locations:
            kept = [obj for obj in possible_datasets if obj["location"] == "GLO"]

        ###
        if not kept and "RER" in possible_locations:
            kept = [obj for obj in possible_datasets if obj["location"] == "RER"]
        #    print(f"'RER' used for regionalization, exchange: {exc['name']}, {exc['location']}, dataset loc: {ds['location']}")
        #if not kept and "EUR" in possible_locations and 'RER' not in possible_locations:
        #    kept = [obj for obj in possible_datasets if obj["location"] == "EUR"]
        #    print(f"'EUR' used for regionalization, exchange: {exc['name']}, {exc['location']}, dataset loc: {ds['location']}")
        ###

        if not kept:
            if drop_invalid:
                continue
            else:
                raise ValueError("Invalidlink")

        allocated = allocate_inputs(exc, kept)

        new_exchanges.extend(allocated)

    ds["exchanges"] = [
        exc for exc in ds["exchanges"] if exc["type"] != "technosphere"
    ] + new_exchanges
    return ds

##NOTE: The following functions are copied from the Wurst package, full credits to Mutel et al.: https://github.com/polca/wurst/tree/main
def reference_product(ds):
    """Get single reference product exchange from a dataset.

    Raises ``wurst.errors.NoResults`` or ``wurst.errors.MultipleResults`` if zero or multiple results are returned."""
    excs = [
        exc for exc in ds["exchanges"] if exc["amount"] and exc["type"] == "production"
    ]
    if not excs:
        raise ValueError("No suitable production exchanges found")
    elif len(excs) > 1:
        raise ValueError("Multiple production exchanges found")
        
    return excs[0]

def allocate_inputs(exc, lst):
    """Allocate the input exchanges in ``lst`` to ``exc``, using production volumes where possible, and equal splitting otherwise.

    Always uses equal splitting if ``RoW`` is present."""
    has_row = any((x["location"] in ("RoW", "GLO") for x in lst))
    pvs = [reference_product(o).get("production volume") or 0 for o in lst]
    if all((x > 0 for x in pvs)) and not has_row:
        # Allocate using production volume
        total = sum(pvs)
    else:
        # Allocate evenly
        total = len(lst)
        pvs = [1 for _ in range(total)]

    def new_exchange(exc, location, factor):
        cp = deepcopy(exc)
        cp["location"] = location
        return rescale_exchange(cp, factor)

    return [
        new_exchange(exc, obj["location"], factor / total)
        for obj, factor in zip(lst, pvs)
    ]


def get_possibles(exchange, data):
    """FIlter a list of datasets ``data``, returning those with the save name, reference product, and unit as in ``exchange``.

    Returns a generator."""
    key = (exchange["name"], exchange["product"], exchange["unit"])
    for ds in data:
        if (ds["name"], ds["reference product"], ds["unit"]) == key:
            yield ds

def rescale_exchange(exc, value, remove_uncertainty=True):
    """Function to rescale exchange amount and uncertainty.

    * ``exc`` is an exchange dataset.
    * ``value`` is a number, to be multiplied by the existing amount.
    * ``remove_uncertainty``: Remove (unscaled) uncertainty data, default is ``True``.
    If ``False``, uncertainty data is scaled by the same factor as the amount
    (except for lognormal distributions, where the ``loc`` parameter is scaled by the log of the factor).
    Currently, does not rescale for Bernoulli, Discrete uniform, Weibull, Gamma, Beta, Generalized Extreme value
    and Student T distributions.

    Returns the modified exchange."""
    assert isinstance(exc, dict), "Must pass exchange dictionary"
    assert isinstance(value, Number), "Constant factor ``value`` must be a number"

    # Scale the amount
    exc["amount"] *= value

    # Scale the uncertainty fields if uncertainty is not being removed
    if not remove_uncertainty and "uncertainty type" in exc:
        uncertainty_type = exc["uncertainty type"]

        # No uncertainty, do nothing
        if uncertainty_type in {0, 6, 7, 8, 9, 10, 11, 12}:
            pass
        elif uncertainty_type in {1, 2, 3, 4, 5}:
            # Scale "loc" by the log of value for lognormal distribution
            if "loc" in exc and uncertainty_type == 2:
                exc["loc"] += math.log(value)
            elif "loc" in exc:
                exc["loc"] *= value

            # "scale" stays the same for lognormal
            # For other distributions, scale "scale" by the absolute value
            if "scale" in exc and uncertainty_type not in {2}:
                exc["scale"] *= abs(value)

            # Scale "minimum" and "maximum" by value
            for bound in ("minimum", "maximum"):
                if bound in exc:
                    exc[bound] *= value

    # If remove_uncertainty is True, then remove all uncertainty info
    elif remove_uncertainty:
        FIELDS = ("scale", "minimum", "maximum", )
        exc["uncertainty type"] = 0
        exc["loc"] = exc["amount"]
        for field in FIELDS:
            if field in exc:
                del exc[field]

    return exc