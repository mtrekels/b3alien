import pandas as pd
from tqdm import tqdm
tqdm.pandas()  # enables .progress_apply
from pygbif import species
import numpy as np

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
        species_to_keep = df_merged["speciesKey"].unique()
        species_to_keep = np.where(species_to_keep == 'Uncertain', -1, species_to_keep)
        species_to_keep = species_to_keep.astype(int)

        return species_to_keep
        

def get_speciesKey(sciname):
    """
        Match text strings with the GBIF taxonomic backbone.

        Parameters
        ----------
        sciname : str
            Text string of a scientific name

        Returns
        -------
        speciesKey: an integer speciesKey number
    """
    result = species.name_backbone(sciname, strict=True)
    try:
        speciesKey = result["speciesKey"]
    except KeyError:
        speciesKey = "Uncertain"
    return speciesKey

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

        return tot_species

    else:
        taxon = filePath + "taxon.txt"
        distribution = filePath + "distribution.txt"

        df_t = pd.read_csv(taxon, sep="\t")
        df_dist = pd.read_csv(distribution, sep="\t")


        # Now apply this on the whole dataframe

        df_t["speciesKey"] = df_t["scientificName"].apply(get_speciesKey)

        df_merged = df_dist.merge(df_t[['id', 'speciesKey']], on='id', how='left')

        species_to_keep = df_merged["speciesKey"].unique()
        species_to_keep = np.where(species_to_keep == 'Uncertain', -1, species_to_keep)
        species_to_keep = species_to_keep.astype(int)

        return species_to_keep