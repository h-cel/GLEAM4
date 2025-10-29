# %% Imports
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import requests
import rootutils
import xarray as xr

# Set up root directory with rootutils
rootutils.setup_root(
    search_from=__file__,
    indicator=".project-root",
    pythonpath=True,  # add root directory to the PYTHONPATH (allows imports from root)
)

from conf import (
    base_period,
    filename_soi_processed,
    filename_soi_raw,
    folder_gridded_monthly,
    folder_gridded_yearly,
    folder_sotc,
    url_soi_cru,
    year_of_interest,
)
from functions import data_logging

folder_sotc.mkdir(exist_ok=True, parents=True)

data_logging(folder_sotc, __file__)  # set up logging
# %% Download Southern Oscillation Index (SOI) data
# https://crudata.uea.ac.uk/cru/data/soi/
logging.info("Downloading and processing Southern Oscillation Index (SOI) data")
query_parameters = {"downloadformat": "txt"}
response = requests.get(url_soi_cru, params=query_parameters)
with open(folder_sotc / filename_soi_raw, "wb") as f:
    f.write(response.content)
df_soi = pd.read_table(
    folder_sotc / filename_soi_raw,
    header=None,
    sep=r"\s+",
    na_values=-99.99,
)
# df_soi_cleaned = pd.DataFrame(data=[df_soi[0], df_soi[13]], columns=["year", "soi"])
df_soi_cleaned = pd.DataFrame.from_dict(
    {"year": df_soi[0], "soi": df_soi[13]}
).set_index("year")
df_soi_cleaned.to_csv(folder_sotc / filename_soi_processed)

# %% Read in yearly and monthly data
logging.info("Reading in yearly GLEAM data")
ds_yearly = xr.open_mfdataset(
    sorted((folder_gridded_yearly / "E").iterdir()),
    concat_dim="time",
    combine="nested",
    data_vars="minimal",
    coords="minimal",
    compat="override",
)
# see https://docs.xarray.dev/en/stable/user-guide/io.html#reading-multi-file-datasets

logging.info("Reading in monthly GLEAM data")
ds_monthly = xr.open_mfdataset(
    sorted((folder_gridded_monthly / "E").iterdir()),
    concat_dim="time",
    combine="nested",
    data_vars="minimal",
    coords="minimal",
    compat="override",
)

# %% Temporal mean per pixel over the years
logging.info(
    f"Computing temporal mean per pixel for yearly GLEAM data "
    f"({base_period[0]} - {base_period[1]})"
)
da_temporal_mean = (
    ds_yearly["E"].sel(time=slice(*base_period)).mean(dim="time").compute()
)

# %% Spatial mean per year over different regions -> time series anomalies
# Creating weights folows:
# https://docs.xarray.dev/en/stable/examples/area_weighted_temperature.html#Creating-weights
# weights = np.cos(np.deg2rad(ds_yearly["lat"]))
logging.info("Spatial mean per region for evaporation per year (including anomalies)")
dict_region = {
    "global": {"lat_min": -90, "lat_max": 90},
    "NH": {"lat_min": 0, "lat_max": 90},
    "SH": {"lat_min": -90, "lat_max": 0},
}
da_anomaly_dict = {}
for region in dict_region:
    # Select per region
    da_region_temporal_mean_ = da_temporal_mean.sel(
        lat=slice(dict_region[region]["lat_max"], dict_region[region]["lat_min"])
    )
    da_region_yearly_ = ds_yearly["E"].sel(
        lat=slice(dict_region[region]["lat_max"], dict_region[region]["lat_min"])
    )
    # Compute weights for region
    weights_ = np.cos(np.deg2rad(da_region_yearly_["lat"]))
    # Compute spatial means per region
    da_region_temporal_mean_spatial_mean_ = (
        da_region_temporal_mean_.weighted(weights_).mean(("lat", "lon")).compute()
    )
    da_region_yearly_spatial_mean_ = (
        da_region_yearly_.weighted(weights_).mean(("lat", "lon")).compute()
    )
    # Timeseries anomalies
    da_region_yearly_anomalies_ = (
        da_region_yearly_spatial_mean_ - da_region_temporal_mean_spatial_mean_
    )
    da_anomaly_dict[region] = da_region_yearly_anomalies_

# Combine in one data_array
da_yearly_anomalies = xr.concat(
    list(da_anomaly_dict.values()),
    dim=xr.DataArray(list(da_anomaly_dict.keys()), dims="region", name="region"),
)

# %% Spatial anomaly for one year
logging.info(
    f"Spatial anomaly (relative to {base_period[0]} - {base_period[1]}) for one year:"
    f" {year_of_interest}"
)
da_spatial_anomaly = (
    ds_yearly["E"].sel(time=year_of_interest) - da_temporal_mean
).compute()

# %% Monthly climatology per latitude
logging.info(
    f"Computing monthly climatology ({base_period[0]} - {base_period[1]}) per latitude"
)
da_monthly_climatology_per_lat = (
    ds_monthly["E"]
    .sel(time=slice(*base_period))  # Only calculate average over base period
    .groupby("time.month")  # Average calculated per  month
    .mean(("lon", "time"))  # Average over the years and latitudes
    .compute()
)

# %% Monthly evaporation per latitude per year
logging.info("Computing monthly evaporation per latitude per year")
da_monthly_per_lat = ds_monthly["E"].mean(dim="lon").compute()

# %% Monthly anomalies per latitude per year
# https://docs.xarray.dev/en/stable/examples/weather-data.html#Calculate-monthly-anomalies
logging.info("Computing monthly anomalies per latitude per year")
da_monthly_anomalies_per_lat = (
    da_monthly_per_lat.groupby("time.month") - da_monthly_climatology_per_lat
)
