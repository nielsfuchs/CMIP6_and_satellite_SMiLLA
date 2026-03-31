#

import re
import xarray as xr
import intake
import logging


def most_frequent(lst):
  """Finds the most frequent element in a list.

  Args:
    lst: A list of elements.

  Returns:
    The most frequent element in the list.
  """

  freq_counter = {}
  for item in lst:
    freq_counter[item] = freq_counter.get(item, 0) + 1

  max_freq = 0
  most_frequent_item = None
  for item, freq in freq_counter.items():
    if freq > max_freq:
      max_freq = freq
      most_frequent_item = item

  return most_frequent_item


def sort_and_filter_cmip6_members(members:list):
    """Sorts CMIP6 ensemble members by 'r' and filters out members with 'i' != 1.
    
    Args:
    members: A list of CMIP6 ensemble member names.
    
    Returns:
    A sorted list of CMIP6 ensemble member names.
    """
    
    def get_ri_values(member):
        match = re.search(r'r(\d+)i(\d+)', member)
        return int(match.group(1)), int(match.group(2)) if match else (0, 0)

    filtered_members = [member for member in members if get_ri_values(member)[1] < 10]
    sorted_members   = sorted(filtered_members, key=lambda x: get_ri_values(x)[0])
    
    return sorted_members


def modelsearch(col, scenario, variable, logger=None, member=None, model=None):

    """
    Search for climate model data based on the provided scenario and variable.

    Parameters:
    - col (???): ????################################################################################################
    - scenario (str): The scenario ID to search for (e.g., "ssp126" or "historical").
    - variable (str): The variable ID to search for.
    - logger (Logger, optional): A logger for recording warnings and information.
    - member (str, optional): The member ID to filter the search results.
    - model (str, optional): The model ID to filter the search results.

    Returns:
    - activity_id (str): The activity ID corresponding to the search results.
    - table_id (str): The most frequent table ID in the search results.
    - modellist (list): A sorted list of unique model IDs found in the search.
    - modelcenters (dict): A dictionary mapping model IDs to their institution IDs.
    - ensemblemembers (dict): A dictionary mapping model IDs to their ensemble members.
    """
    
    # Construct the query based on input parameters
    query = {
        "experiment_id": scenario,
        "variable_id": variable,
        "frequency": "mon",
        #"activity_id":"ScenarioMPI"
    }

    # Set the activity ID based on the scenario type
    if scenario.startswith("ssp"):
        query["activity_id"] = "ScenarioMIP"
    else:
        query["activity_id"] = "CMIP"

    # Add model and member to query if provided
    if model is not None:
        query["source_id"] = model
    if member is not None:
        query["member_id"] = member
    
    # Perform the search and retrieve the DataFrame   
    df = col.search(**query).df
    
    try:
        # Extract relevant information from the DataFrame
        table_id     = most_frequent(df["table_id"]) #df["table_id"][0]
        activity_id  = df["activity_id"].iloc[0]
        modellist = sorted(set(df["source_id"]))

        # Log information if a logger is provided
        if logger is not None:
            logger.warning(f"VARIABLE:  {variable}")
            if len(set(df["table_id"]))>1:
                logger.warning(f"More than one table_id: {set(df['table_id'])}")
            logger.info(f"table_id:  {table_id}")
            logger.info(f"long_name: {df['long_name'].iloc[0]}")
            logger.info(f"units:     {set(df['units'])}")
            unit = most_frequent(df["units"])
            #logger.info(f"grid:      {list(df.colums)}") 
            logger.debug(f"dummy:     {df['uri'].iloc[0]}")

            # Calculate mean value for example file
            ds   = xr.open_dataset(df["uri"].iloc[0])
            mean = ds[variable].isel(time=0).mean().values.item()
            ds.close()
            logger.info(f"mean for t=0: {mean}\n")
    
        # Create dictionaries for model centers and ensemble members
        modelcenters    = {}
        ensemblemembers = {}
        for model in modellist:
            filtered = df[df['source_id'] == model].copy()
            modelcenters[model] = list(set(filtered["institution_id"]))[0]
            #modelcenters[model] = set(filtered["institution_id"]).pop()  # Get the unique institution ID

            # UKESM1-0-LL is run by two modelling centers (Maybe there's a better way to do this)
            if model == "UKESM1-0-LL":
                modelcenters[model] = "MOHC"
                
            ensemblemembers[model] = sort_and_filter_cmip6_members(set(filtered["member_id"]))
            
    except Exception as e:
        logger.warning("Search not successful")
        logger.warning(query)
        logger.warning(" ")
        return None, None, None, None, None, None

    return activity_id, table_id, modellist, modelcenters, ensemblemembers, unit