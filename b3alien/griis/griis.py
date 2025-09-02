import pandas as pd
from tqdm import tqdm
import numpy as np
import requests
import ast
tqdm.pandas()  # enables .progress_apply

class CheckList():

    """
        Load a GRIIS checklist from GBIF.

        Parameters
        ----------
        filepath : str
            Path to the distribution.txt file of the checklist.

        Returns
        -------
        griis.Checklist
            A checklist object containing the list of species.
    """

    def __init__(self, filePath: str):
        self.filePath = filePath

        # Create cube
        self.species = self._load_GRIIS(filePath)

    def _load_GRIIS(self, filePath):
        
        df_merged = pd.read_csv(filePath, sep="\t")
        df_merged['speciesKey'] = df_merged['speciesKey'].apply(
            lambda x: ast.literal_eval(x) if isinstance(x, str) else x
        )
        df_exploded = df_merged.explode("speciesKey")
        species_to_keep = df_exploded["speciesKey"].unique()
        species_to_keep = [
            int(x) for x in species_to_keep
            if pd.notnull(x) and x != "Uncertain"
        ]

        return species_to_keep

 
def get_species_under_genus(taxon_key):
    """
        Get all the keys of the species listed under a specific genus.

        Parameters
        ----------
        taxon_key : int
            The GBIF taxonKey of the genus.
        Returns
        -------
        list
            A list of speciesKey(s) under the specified genus.
    """
    species_keys = []
    offset = 0
    limit = 1000

    while True:
        url = f"https://api.gbif.org/v1/species/{taxon_key}/children"
        params = {"rank": "species", "limit": limit, "offset": offset}
        response = requests.get(url, params=params)
        if response.status_code != 200:
            break

        data = response.json()
        results = data.get("results", [])
        if not results:
            break

        keys = [res["key"] for res in results if res.get("rank", "").upper() == "SPECIES"]
        species_keys.extend(keys)
        offset += limit

    return species_keys or ["Uncertain"]

def get_speciesKey(sciname):
    """
        Resolve a scientific name to its GBIF taxonKey. If the name is a genus, retrieve all species under that genus.
        Parameters
        ----------
        sciname : str
            The scientific name to resolve.
        Returns
        -------
        list
            A list of resolved speciesKey(s) or ["Uncertain"] if unresolved.
    """
    try:
        # Query GBIF backbone for the name
        response = requests.get(
            "https://api.gbif.org/v1/species/match",
            params={"name": sciname, "strict": True},
            timeout=10
        )
        result = response.json()
    except Exception:
        return ["Uncertain"]

    if "usageKey" not in result:
        return ["Uncertain"]

    taxon_key = result["usageKey"]
    rank = result.get("rank", "").upper()

    # Case 1: SPECIES
    if rank == "SPECIES":
        return [taxon_key]

    # Case 2: GENUS — query children directly from GBIF API
    elif rank == "GENUS":
        all_species_keys = []
        offset = 0
        limit = 1000

        while True:
            try:
                children_url = f"https://api.gbif.org/v1/species/{taxon_key}/children"
                children_response = requests.get(
                    children_url,
                    params={"rank": "species", "limit": limit, "offset": offset},
                    timeout=10
                )
                children_data = children_response.json()
                results = children_data.get("results", [])

                if not results:
                    break

                species_keys = [r["key"] for r in results if r.get("rank", "").upper() == "SPECIES"]
                all_species_keys.extend(species_keys)
                offset += limit
            except Exception:
                break

        return all_species_keys if all_species_keys else ["Uncertain"]

    # Case 3: Other ranks or unresolvable
    return ["Uncertain"]

def split_event_date(eventDate):
    """
        Interprete the event date as introduction date and date of last seen,
        when this information is available in the checklist.

        Parameters
        ----------
        eventDate : str
            Text string of eventDate

        Returns
        -------
        pd.Series
            A series containing introduction date ('intro') and date last seen ('outro')
    """
    if isinstance(eventDate, str):
        parts = eventDate.strip().split('/')
        if len(parts) == 2:
            intro = parts[0]
            outro = parts[1]
        else:
            intro = outro = np.nan
        return pd.Series([intro, outro])
    else:
        return pd.Series([np.nan, np.nan])


def do_taxon_matching(dirPath):
    """
        Match keys between taxon.txt and distribution.txt

        Parameters
        ----------
        dirPath : str
            Path to the directory of the checklist

        Returns
        -------
        Saves a new checklist file 'merged_distr.txt' in the checklist directory
    """

    taxon = dirPath + "taxon.txt"
    distribution = dirPath + "distribution.txt"

    df_t = pd.read_csv(taxon, sep="\t")
    df_dist = pd.read_csv(distribution, sep="\t")

    # Now apply this on the whole dataframe

    df_t["speciesKey"] = df_t["scientificName"].progress_apply(get_speciesKey)

    df_merged = df_dist.merge(df_t[['id', 'speciesKey']], on='id', how='left')
    df_merged.to_csv(dirPath + 'merged_distr.txt', sep='\t', index=False)

# The rest assumes already a merged dataset
def read_checklist(filePath, cl_type='detailed', locality='Belgium'):
    
    distribution = filePath + "distribution.txt"

    df_cl = pd.read_csv(distribution, sep='\t', low_memory=False)

    df_cl["speciesKey"] = df_cl["id"].str.rsplit("/", n=1).str[-1].astype("int64")


    if cl_type == 'detailed':

        species_to_keep = df_cl["speciesKey"].astype("int64").unique()
    
        # 1. Filter rows where locality == 'Belgium' and eventDate is not missing
        df = df_cl[df_cl["locality"] == locality].copy()
        df = df[df["eventDate"].notna()]



        df[["introDate", "outroDate"]] = df["eventDate"].apply(split_event_date)

        df["introDate"] = pd.to_datetime(df["introDate"], format="%Y", errors="coerce")
        df["outroDate"] = pd.to_datetime(df["outroDate"], format="%Y", errors="coerce")


        # 3. Clean rows with missing introDate
        df_intro = df.dropna(subset=["introDate"]).copy()

        # 4. Group by introDate and count species
        in_species = (
            df_intro.groupby("introDate", sort=True)["id"]
            .count()
            .reset_index(name="nspec")
        )

        # 5. Cumulative sum
        in_species["cumn"] = in_species["nspec"].cumsum()

        # 6. Clean outro side and count outgoing species
        df_outro = df.dropna(subset=["outroDate"]).copy()

        out_species = (
            df_outro.groupby("outroDate", sort=True)["id"]
            .count()
            .reset_index(name="nspeco")
        )

        # 7. Merge intro and outro on date
        n_species = pd.merge(in_species, out_species, how="outer", left_on="introDate", right_on="outroDate")

        # 8. Replace NaNs with 0
        n_species["nspec"] = n_species["nspec"].fillna(0).astype(int)
        n_species["nspeco"] = n_species["nspeco"].fillna(0).astype(int)

        # 9. Net species present at each time step
        n_species["total"] = n_species["nspec"] - n_species["nspeco"]

        # 10. Final frame with total species over time
        tot_species = n_species[["introDate", "total"]].copy()

        # 11. Optional: sort and compute cumulative total over time
        tot_species = tot_species.sort_values("introDate")
        tot_species["cumulative_total"] = tot_species["total"].cumsum()

        return species_to_keep, tot_species

    else:
        taxon = filePath + "taxon.txt"
        distribution = filePath + "distribution.txt"

        df_t = pd.read_csv(taxon, sep="\t")
        df_dist = pd.read_csv(distribution, sep="\t")


        # Now apply this on the whole dataframe

        # Apply the function — returns lists
        df_t["speciesKey"] = df_t["scientificName"].apply(get_speciesKey)

        # Explode so each speciesKey gets its own row
        df_t_exploded = df_t.explode("speciesKey")

        # Merge
        df_merged = df_dist.merge(df_t_exploded[['id', 'speciesKey']], on='id', how='left')

        # Clean and filter
        species_to_keep = df_merged["speciesKey"].unique()
        species_to_keep = [int(x) for x in species_to_keep if x != "Uncertain"]

        return species_to_keep