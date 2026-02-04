# %% Imports
from pathlib import Path

import rootutils

root_path = rootutils.find_root(search_from=__file__, indicator=".project-root")

# %% In situ data doi
zenodo_doi = "10.5281/zenodo.14054258"

# %% Google earth engine
GEE_PROJECT_ID = "ee-bonteolivier15"  # Replace by your own project ID

# %% Paths
folder_insitu = Path(
    "/data/gent/vo/000/gvo00090/GLEAM/data/data/in_situ_data/GLEAM4_validation/cleaned_v42b_retrain_folder"
)
folder_figures = root_path / "figures"

# GLEAM 4 Paths
folder_gridded = Path(
    "/scratch/gent/vo/000/gvo00090/GLEAM/data/GLEAM4/GLEAM4_outputs/output/GLEAM4.3a/GLEAM4.3a_OL/processed_Oscar"
)  # Adapt this path for your computer: point to folder with
# GLEAM outputs as found on the SFTP server (https://www.gleam.eu/#downloads)
folder_gridded_monthly = folder_gridded / "monthly"
folder_gridded_yearly = folder_gridded / "yearly"

# Processed data for plotting
folder_processed = root_path / "data" / "processed"
folder_sotc = root_path / "data" / "sotc"

model_version = "GLEAM4.3a"
seasonal_averages_file = f"{model_version}_seasonal_averages.nc"
yearly_averages_file = f"{model_version}_yearly_averages.nc"
yearly_averages_file_comparison = f"{model_version}_yearly_averages_comparison.nc"

# SotC Paths
url_soi_cru = "https://crudata.uea.ac.uk/cru/data/soi/soi_3dp.dat"
# More info on dataset: https://crudata.uea.ac.uk/cru/data/soi/
filename_soi_raw = "soi_cru_raw.txt"
filename_soi_processed = "soi_cru_processed.csv"
filename_yearly_anomalies = f"{model_version}_yearly_anomalies.nc"
filename_spatial_anomaly = f"{model_version}_spatial_anomaly.nc"
filename_monthly_anomalies_per_lat = f"{model_version}_monthly_anomalies_per_lat.nc"

# For data downloading of external products
# Adapt these path for your computer if not on UGhent-HPC system
gleam_38_path = Path(
    "/data/gent/vo/000/gvo00090/GLEAM/data/data/GLEAM_v3.8/v38a_output/yearly/E"
)  # Data read from this path (no download script provided)
era5_land_path = (
    Path("/data/gent/vo/000/gvo00090/EXT/data/ERA5_Land") / "total_evaporation"
)  # Data downloaded to and read from this path
fluxcom_x_base_path = (
    Path("/data/gent/vo/000/gvo00090/EXT/data/FLUXCOM_X_BASE") / "ET"
)  # Data downloaded to and read from this path

# %% Global averaging for global patterns
non_summing_vars = [
    "S",
    "SMs",
    "SMrz",
    "H",
]  # Variables that should be averaged, not summed when resampling temporally globally
E_plot_dict = {"vmin": 0, "vmax": 1400}
# %% Choice of variables for seasonal plots
season_plot_var_dict = {
    "E": {"vmin": 0, "vmax": 400, "units": "mm"},
    "Ep": {"vmin": 0, "vmax": 800, "units": "mm"},
    "S": {"vmin": 0, "vmax": 1, "units": "-"},
    "SMrz": {"vmin": 0, "vmax": 1, "units": "m³/m³"},
    "H": {"vmin": 0, "vmax": 120, "units": "W/m²"},
}
season_choice = ["JJA", "DJF"]

# %% Colors
model_list_validation = ["GLEAM4", "GLEAM v3.8", "ERA5-Land", "FLUXCOM"]
color_dict = {
    model_list_validation[0]: "#09070d",
    model_list_validation[1]: "#6f6db1",
    model_list_validation[2]: "#cf6666",
    model_list_validation[3]: "#4a8740",
}
## As of 27/10/2025, FLUXCOM-X used for global patterns, not in-situ validation
model_list_pattern = model_list_validation.copy()
model_list_pattern[3] = "FLUXCOM-X"
color_dict_pattern = color_dict.copy()
color_dict_pattern.update(
    {model_list_pattern[3]: color_dict_pattern.pop(model_list_validation[3])}
)

color_dict_lat = {
    "transpiration": "#7cc691",
    "bare soil evaporation": "#cc9a80",
    "interception": "#cce7f8",
    "other": "#989899",
}
violin_fill_color = "#484948"
insitu_fill_color = "darkgray"
ocean_color = "#CDE0E4"

# %% png saving: Figures 2,3,4,5
dpi = 900
# %% Sites
site_selection = {
    "DE-Tha": "2003",
    "US-Ne3": "2012",
    "AU-How": "2010",
    "FR-Pue": "2016",
}

# %% SotC options
base_period = ["1991", "2020"]
year_of_interest = "2025"
dict_region = {
    "global": {"lat_min": -90, "lat_max": 90, "color": "#0b141c", "label": "Globe"},
    "NH": {
        "lat_min": 0,
        "lat_max": 90,
        "color": "#2525fe",
        "label": "N. Hemisphere",
    },
    "SH": {
        "lat_min": -90,
        "lat_max": 0,
        "color": "#fe1b1a",
        "label": "S. Hemisphere",
    },
}
soi_fill_colors = {"pos": "#84cefc", "neg": "#f37f7c"}
colorbar_hovmoller = {"abs_lim": 10, "spacing": 2, "cmap": "BrBG"}
colorbar_spatial_anomaly = {"abs_lim": 200, "spacing": 50, "cmap": "BrBG"}
