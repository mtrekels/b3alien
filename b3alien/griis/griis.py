from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import zipfile
import requests
import ast
import pandas as pd
from tqdm import tqdm
import numpy as np
tqdm.pandas()  # enables .progress_apply

#: GBIF checklistKey of the Catalogue of Life eXtended Release, the taxonomic
#: backbone used by GBIF occurrence cubes generated on the new backbone.
COL_XR_CHECKLIST_KEY = "7ddf754f-d193-4cc9-b351-99906754a03b"

GBIF_MATCH_V2_URL = "https://api.gbif.org/v2/species/match"


def _normalise_key(key):
    """
        Return a GBIF key as int when it is numeric (old backbone) and as str
        otherwise (new backbone, e.g. 'CPZZL').
    """
    if isinstance(key, (int, np.integer)):
        return int(key)
    if isinstance(key, float) and key.is_integer():
        return int(key)
    key = str(key).strip()
    return int(key) if key.lstrip("-").isdigit() else key


def _keys_from_column(series):
    """Explode a speciesKey column (list strings or scalars) into a unique list of keys."""
    series = series.apply(
        lambda x: ast.literal_eval(x) if isinstance(x, str) and x.startswith("[") else x
    )
    keys = series.explode().unique()
    keys = [_normalise_key(x) for x in keys if pd.notnull(x) and x != "Uncertain"]
    return list(dict.fromkeys(keys))


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

    def __init__(self, filePath: str, use_profile: bool = False, invasive_val: str = None, habitat_val: str = None):
        self.filePath = filePath
        self.use_profile = use_profile
        self.invasive_val = invasive_val
        self.habitat_val = habitat_val
        path_obj = Path(filePath)
        if path_obj.is_dir():
            self.base_dir = path_obj
        else:
            self.base_dir = path_obj.parent
            
        self.species = self._load_GRIIS(self.base_dir)

    def _load_GRIIS(self, base_path):
        """
            Load the checklist and extract speciesKey(s).

            Parameters
            ----------
            filePath : str
                Path to the distribution.txt file of the checklist.

            Returns
            -------
            list
                A list of speciesKey(s) in the checklist.
        """
        merged_path = base_path / "merged_distr.txt"

        if self.use_profile:
            filter_invasive = pd.read_csv(base_path / "speciesprofile.txt", sep="\t")
            df_merged = pd.read_csv(merged_path, sep="\t")

            filtered_merged = filter_by_speciesprofile(
                filter_invasive, df_merged, self.invasive_val, self.habitat_val
            )
            return _keys_from_column(filtered_merged["speciesKey"])

        else:

            df_merged = pd.read_csv(merged_path, sep="\t")
            return _keys_from_column(df_merged["speciesKey"])

    def _add_species(self, new_species):
        """
            Add new species to the checklist.

            Parameters
            ----------
            new_species : list
                A list of speciesKey(s) to add to the checklist.

            Returns
            -------
            None
        """
        for sp in new_species:
            if sp not in self.species:
                self.species.append(sp) 

    def _remove_species(self, rem_species):
        """
            Remove species from the checklist.

            Parameters
            ----------
            rem_species : list
                A list of speciesKey(s) to remove from the checklist.

            Returns
            -------
            None
        """
        self.species = [sp for sp in self.species if sp not in rem_species] 


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

def get_speciesKey(sciname, checklistKey=None, **match_params):
    """
        Resolve a scientific name to its GBIF taxonKey. If the name is a genus, retrieve all species under that genus.
        Parameters
        ----------
        sciname : str
            The scientific name to resolve.
        checklistKey : str, optional
            GBIF checklist to match against with the v2 match API, e.g.
            ``COL_XR_CHECKLIST_KEY`` for cubes on the new GBIF backbone. If None
            (default), the old GBIF backbone (v1 API, numeric keys) is used.
        **match_params
            Extra parameters for the v2 match API (e.g. kingdom, taxonRank,
            scientificNameAuthorship). Only used with ``checklistKey``.
        Returns
        -------
        list
            A list of resolved speciesKey(s) or ["Uncertain"] if unresolved.
    """
    if checklistKey is not None:
        return _match_v2(sciname, checklistKey, **match_params)

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

def _match_v2(sciname, checklistKey=COL_XR_CHECKLIST_KEY, session=None,
              accepted_match_types=("EXACT", "VARIANT"), **match_params):
    """
        Match a name with the GBIF v2 match API against a given checklist.

        Synonyms resolve to their accepted species and infraspecific taxa to the
        species they belong to, because occurrence cubes are aggregated on the
        accepted species key. Names matched above species rank (e.g. genera) or
        with a fuzzy/higher-rank match are returned as ["Uncertain"].

        Returns
        -------
        list
            [speciesKey] (a str, e.g. 'CPZZL') or ["Uncertain"].
    """
    params = {"scientificName": sciname, "checklistKey": checklistKey}
    params.update({k: v for k, v in match_params.items() if v is not None and v == v})
    requester = session or requests
    try:
        response = requester.get(GBIF_MATCH_V2_URL, params=params, timeout=30)
        response.raise_for_status()
        result = response.json()
    except Exception:
        return ["Uncertain"]

    match_type = result.get("diagnostics", {}).get("matchType", "NONE")
    if not result.get("usage") or match_type not in accepted_match_types:
        return ["Uncertain"]

    species = [c["key"] for c in result.get("classification", []) if c.get("rank", "").upper() == "SPECIES"]
    return [species[-1]] if species else ["Uncertain"]


def _read_coldp_table(source, name):
    """
        Read a ColDP table (e.g. 'NameUsage') from a ColDP directory or zip file.
        Column prefixes such as 'col:' and 'clb:' are removed.
    """
    source = Path(source)
    candidates = [f"{name}.tsv", f"{name}.txt", f"{name}.csv"]

    if source.is_dir():
        files = {p.name.lower(): p for p in source.iterdir()}
        for cand in candidates:
            if cand.lower() in files:
                path = files[cand.lower()]
                sep = "," if path.suffix == ".csv" else "\t"
                df = pd.read_csv(path, sep=sep, dtype=str, quoting=3 if sep == "\t" else 0)
                break
        else:
            return None
    else:
        with zipfile.ZipFile(source) as zf:
            members = {Path(m).name.lower(): m for m in zf.namelist()}
            for cand in candidates:
                if cand.lower() in members:
                    sep = "," if cand.endswith(".csv") else "\t"
                    with zf.open(members[cand.lower()]) as fh:
                        df = pd.read_csv(fh, sep=sep, dtype=str, quoting=3 if sep == "\t" else 0)
                    break
            else:
                return None

    df.columns = df.columns.str.split(":").str[-1]
    return df


def is_coldp(path):
    """
        Return True if ``path`` is a ColDP archive (zip) or directory (contains NameUsage).
    """
    path = Path(path)
    if path.is_dir():
        return any(p.name.lower().startswith("nameusage.") for p in path.iterdir())
    if path.is_file() and zipfile.is_zipfile(path):
        with zipfile.ZipFile(path) as zf:
            return any(Path(m).name.lower().startswith("nameusage.") for m in zf.namelist())
    return False


def do_taxon_matching_coldp(path, output_dir=None, checklistKey=COL_XR_CHECKLIST_KEY, max_workers=8):
    """
        Match a GRIIS checklist in ColDP format (as downloaded from ChecklistBank)
        to the GBIF backbone used by the occurrence cubes.

        The scientific names in NameUsage are matched with the GBIF v2 match API
        against ``checklistKey`` (by default the Catalogue of Life eXtended
        Release, i.e. the new GBIF backbone with keys like 'CPZZL').

        Two files are written to ``output_dir`` so the result can be loaded with
        ``CheckList(output_dir, use_profile=..., invasive_val=..., habitat_val=...)``:

        - merged_distr.txt : Distribution records with their matched speciesKey;
        - speciesprofile.txt : 'id', 'isInvasive' ('Invasive'/'Null') and
          'habitat' (e.g. 'Terrestrial|Freshwater'), derived from TaxonProperty.

        Parameters
        ----------
        path : str
            Path to the ColDP zip file or extracted ColDP directory.
        output_dir : str, optional
            Directory for the output files. Defaults to the ColDP directory, or
            for a zip file a directory next to it with the same name.
        checklistKey : str, optional
            GBIF checklist to match against. Default is the COL eXtended Release.
        max_workers : int, optional
            Number of parallel requests to the GBIF API. Default is 8.

        Returns
        -------
        pandas.DataFrame
            The merged distribution table (also saved as merged_distr.txt).
    """
    path = Path(path)
    if output_dir is None:
        output_dir = path if path.is_dir() else path.with_suffix("")
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    names = _read_coldp_table(path, "NameUsage")
    distribution = _read_coldp_table(path, "Distribution")
    if names is None or distribution is None:
        raise FileNotFoundError("The ColDP source must contain NameUsage and Distribution tables.")

    taxa = distribution[["taxonID"]].drop_duplicates().merge(
        names, left_on="taxonID", right_on="ID", how="left"
    )

    def match_row(row):
        if pd.isna(row.get("scientificName")):
            return ["Uncertain"]
        return _match_v2(
            row["scientificName"], checklistKey, session=session,
            scientificNameAuthorship=row.get("authorship"),
            kingdom=row.get("kingdom"),
            taxonRank=row.get("rank"),
        )

    rows = [row for _, row in taxa.iterrows()]
    with requests.Session() as session, ThreadPoolExecutor(max_workers=max_workers) as pool:
        keys = list(tqdm(pool.map(match_row, rows), total=len(rows), desc="Matching names"))
    taxa["speciesKey"] = [str(k) for k in keys]

    keep = ["taxonID", "scientificName", "authorship", "rank", "kingdom", "speciesKey"]
    merged = distribution.merge(taxa[[c for c in keep if c in taxa.columns]], on="taxonID", how="left")
    merged = merged.rename(columns={"taxonID": "id"})
    merged = merged.dropna(axis=1, how="all")
    merged.to_csv(output_dir / "merged_distr.txt", sep="\t", index=False)

    # Species profile compatible with the GRIIS Darwin Core Archive
    profile = pd.DataFrame({"id": distribution["taxonID"].unique()})
    properties = _read_coldp_table(path, "TaxonProperty")
    if properties is not None and not properties.empty:
        prop = properties["property"].str.split(":").str[-1]
        invasive_ids = properties.loc[
            (prop == "isInvasive") & properties["value"].str.lower().isin(["true", "invasive", "yes"]),
            "taxonID",
        ]
        habitat = (
            properties[prop == "habitat"]
            .groupby("taxonID")["value"]
            .agg(lambda v: "|".join(h.strip().capitalize() for x in v for h in str(x).split("|")))
        )
        profile["isInvasive"] = np.where(profile["id"].isin(invasive_ids), "Invasive", "Null")
        profile["habitat"] = profile["id"].map(habitat)
    else:
        profile["isInvasive"] = "Null"
        profile["habitat"] = np.nan
    profile.to_csv(output_dir / "speciesprofile.txt", sep="\t", index=False)

    n_matched = (merged["speciesKey"] != "['Uncertain']").sum()
    print(f"Matched {n_matched}/{len(merged)} distribution records. Output written to {output_dir}")

    return merged


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


def do_taxon_matching(dirPath, checklistKey=None, output_dir=None):
    """
        Match keys between taxon.txt and distribution.txt

        Parameters
        ----------
        dirPath : str
            Path to the directory of the checklist (Darwin Core Archive), or a
            ColDP zip file/directory. ColDP input is handled by
            ``do_taxon_matching_coldp``.
        checklistKey : str, optional
            GBIF checklist to match against. If None, the old GBIF backbone is
            used for Darwin Core Archives and the COL eXtended Release (new GBIF
            backbone) for ColDP.
        output_dir : str, optional
            Output directory (ColDP only), see ``do_taxon_matching_coldp``.

        Returns
        -------
        Saves a new checklist file 'merged_distr.txt' in the checklist directory
    """
    if is_coldp(dirPath):
        return do_taxon_matching_coldp(
            dirPath, output_dir=output_dir, checklistKey=checklistKey or COL_XR_CHECKLIST_KEY
        )

    dirPath = str(dirPath)
    if not dirPath.endswith(("/", "\\")):
        dirPath += "/"

    taxon = dirPath + "taxon.txt"
    distribution = dirPath + "distribution.txt"

    df_t = pd.read_csv(taxon, sep="\t")
    df_dist = pd.read_csv(distribution, sep="\t")

    # Now apply this on the whole dataframe

    df_t["speciesKey"] = df_t["scientificName"].progress_apply(get_speciesKey, checklistKey=checklistKey)

    df_merged = df_dist.merge(df_t[['id', 'speciesKey']], on='id', how='left')
    df_merged.to_csv(dirPath + 'merged_distr.txt', sep='\t', index=False)

# The rest assumes already a merged dataset
def read_checklist(filePath, cl_type='detailed', locality='Belgium'):
    """
        Read a GRIIS checklist and extract speciesKey(s) and time series of species numbers over time.
        Parameters
        ----------
        filePath : str
            Path to the directory of the checklist (must contain distribution.txt and taxon.txt if cl
            type is not 'detailed').
        cl_type : str
            Type of checklist: 'detailed' (with eventDate) or 'simple' (
            without eventDate, requires taxon.txt and distribution.txt).
        locality : str
            The locality to filter on (default is 'Belgium').
        Returns
        -------
        tuple
            A tuple containing:
            - list of speciesKey(s) in the checklist
            - pd.DataFrame with columns 'introDate' and 'cumulative_total' representing the
                cumulative number of species over time.
    """
    
    distribution = filePath + "distribution.txt"
    df_cl = pd.read_csv(distribution, sep='\t', low_memory=False)
    df_cl["speciesKey"] = df_cl["id"].str.rsplit("/", n=1).str[-1].astype("int64")
    if cl_type == 'detailed':

        species_to_keep = df_cl["speciesKey"].astype("int64").unique()
        # 1. Filter rows where locality == 'Belgium' and eventDate is not missing
        df = df_cl[df_cl["locality"] == locality].copy()
        df = df[df["eventDate"].notna()]
        # 2. Split eventDate into introDate and outroDate
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


def filter_by_speciesprofile(df_source, df_target, invasive_val=None, habitat_val=None):
    """
    Filters df_target based on conditions met in df_source.
    Uses AND logic if both filters are provided.
    Parameters
    ----------
    df_source : pd.DataFrame
        The source DataFrame containing 'isInvasive' and 'habitat' columns.
    df_target : pd.DataFrame
        The target DataFrame to be filtered based on df_source.
    invasive_val : str, optional
        The value to filter 'isInvasive' column on (e.g., 'Invasive'). If None, this filter is not applied.
    habitat_val : str, optional
        The value to filter 'habitat' column on (e.g., 'Terrestrial',
        'Marine', Freshwater). If None, this filter is not applied.
    Returns
    -------
    pd.DataFrame
        The filtered target DataFrame.
    """
    # 1. Initialize a mask of all True (start by selecting everything)
    mask = pd.Series(True, index=df_source.index)

    # 2. If no filters are provided, return the original target
    if invasive_val is None and habitat_val is None:
        return df_target

    # 3. Apply 'isInvasive' condition (Narrowing the list)
    if invasive_val is not None:
        mask = mask & (df_source['isInvasive'] == invasive_val)
     
    # 4. Apply 'habitat' condition (Narrowing the list further)
    if habitat_val is not None:
        mask = mask & (df_source['habitat'] == habitat_val)
    
    # 5. Extract unique IDs that meet the combined criteria
    valid_ids = df_source[mask]['id'].unique()

    # 6. Filter target and return a COPY to avoid SettingWithCopyWarning
    return df_target[df_target['id'].isin(valid_ids)].copy()