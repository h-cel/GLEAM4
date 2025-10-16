# %% Imports
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap

# %% In situ data doi
zenodo_doi = "10.5281/zenodo.14054258"

# %% Google earth engine
GEE_PROJECT_ID = "ee-bonteolivier15"  # Replace by your own project ID

# %% Paths
folder_insitu = Path("data/sites")
folder_figures = Path("figures")
folder_gridded = Path(
    "/scratch/gent/vo/000/gvo00090/GLEAM/data/GLEAM4/GLEAM4_outputs/output/GLEAM4.2a/GLEAM4.2_DA/processed"
    # "/scratch/gent/vo/000/gvo00090/GLEAM/data/GLEAM4/GLEAM4_outputs/output/GLEAM4.2b/GLEAM4.2b_DA/processed"
)
folder_gridded_monthly = folder_gridded / "monthly"
folder_gridded_yearly = folder_gridded / "yearly"
folder_temp_data = Path("data/temp")

seasonal_averages_file = "GLEAM4.2a_seasonal_averages.nc"
yearly_averages_file = "GLEAM4.2a_yearly_averages.nc"

# For data downloading/processing
gleam_38_path = Path(
    "/data/gent/vo/000/gvo00090/GLEAM/data/data/GLEAM_v3.8/v38a_output/yearly/E"
)
era5_land_path = Path("/data/gent/vo/000/gvo00090/EXT/data/ERA5_Land")
fluxcom_x_base = Path("/data/gent/vo/000/gvo00090/EXT/data/FLUXCOM_X_BASE")

# %% Global averaging for global patterns
REDO = False
non_summing_vars = [
    "S",
    "SMs",
    "SMrz",
    "H",
]  # Variables that should be averaged, not summed when resampling temporally globally

# %% Colors
color_dict = {
    "GLEAM4": "#09070d",
    "GLEAM v3.8": "#6f6db1",
    "ERA5-Land": "#cf6666",
    "FLUXCOM": "#4a8740",
}
color_dict_lat = {
    "transpiration": "#7cc691",
    "bare soil evaporation": "#cc9a80",
    "interception": "#cce7f8",
    "other": "#989899",
}
violin_fill_color = "#484948"
insitu_fill_color = "darkgray"


# %% Sites
site_selection = {
    "DE-Tha": "2003",
    "US-Ne3": "2012",
    "AU-How": "2010",
    "FR-Pue": "2016",
}
