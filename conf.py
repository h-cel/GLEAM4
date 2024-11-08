# %% Imports
from pathlib import Path

# %% In situ data doi
zenodo_doi = "10.5281/zenodo.14054258"

# %% Paths
folder_insitu = Path("data/sites")
folder_figures = Path("figures")

# %% Colors
color_dict = {
    "GLEAM4": "#09070d",
    "GLEAM v3.8": "#6f6db1",
    "ERA5-Land": "#cf6666",
    "FLUXCOM": "#4a8740",
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
