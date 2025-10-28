# %% imports
import importlib
import logging
import re

import xarray as xr

import conf
import functions

importlib.reload(conf)  # Allows updating without restarting the kernel
importlib.reload(functions)

from conf import (
    era5_land_path,
    fluxcom_x_base_path,
    folder_figures,
    folder_gridded_monthly,
    folder_gridded_yearly,
    folder_processed,
    gleam_38_path,
    model_list_pattern,
    model_version,
    non_summing_vars,
    seasonal_averages_file,
    yearly_averages_file,
    yearly_averages_file_comparison,
)
from functions import data_logging

folder_processed.mkdir(exist_ok=True, parents=True)
(folder_processed / model_version).mkdir(exist_ok=True, parents=True)
folder_figures.mkdir(exist_ok=True)

# %% Set up logging
data_logging(folder_processed, __file__)

# %% Read in GLEAM4 data
logging.info(f"GLEAM4 model version: {model_version}")
logging.info("Reading in GLEAM4 data")
## Variables and years
variables = [folder.name for folder in folder_gridded_monthly.iterdir()]
logging.info(f"Variables found: {variables}")
years = [
    re.search(r"\d{4}", file.name).group()
    for file in (folder_gridded_monthly / variables[0]).iterdir()
]
years.sort()
logging.info(f"Years found: {years}")

## Read in all variables
list_of_datastets = []
for folder in [folder_gridded_monthly, folder_gridded_yearly]:
    ds_list = []
    print(f"Reading data from {folder}")
    for var in variables:
        ds_var = xr.open_mfdataset(
            sorted((folder / var).iterdir()),
            concat_dim="time",
            combine="nested",
            data_vars="minimal",
            coords="minimal",
            compat="override",
        )  # see https://docs.xarray.dev/en/stable/user-guide/io.html#reading-multi-file-datasets
        ds_var = ds_var.chunk(
            {"time": 12, "lon": ds_var.lon.size / 4, "lat": ds_var.lat.size / 4}
        )  # To 1 chunk per year, 4 spatial chunks
        ds_list.append(ds_var)
        # ds_list_year.append(ds_var_year)
    ds = xr.merge(ds_list)
    list_of_datastets.append(ds)
ds_monthly, ds_yearly = list_of_datastets

# %% Read in external datasets (lazily)
logging.info("Reading in external datasets")
ds_era5 = xr.open_mfdataset(
    list(era5_land_path.glob("ERA5_Land_total_evaporation_monthly_*.nc"))
)
ds_fluxcom = xr.open_mfdataset(fluxcom_x_base_path.glob("ET_*_005_monthly.nc"))
ds_gleam38 = xr.open_mfdataset(gleam_38_path.glob("E_*_GLEAM_v3.8a_YR.nc"))

# %% Calculate yearly means (full GLEAM4 dataset)
logging.info("Computing average yearly values for full GLEAM4 dataset")
ds_yearly_mean = ds_yearly.mean("time", skipna=True)
# Derive water mask from its nan values for E
water_mask = ds_yearly_mean["E"].isnull()
ds_yearly_mean = ds_yearly_mean.where(~water_mask)
ds_yearly_mean.attrs["years_used_for_averaging"] = (
    f"{ds_yearly.time.min().dt.year.values} - {ds_yearly.time.max().dt.year.values}"
)
ds_yearly_mean.to_netcdf(folder_processed / yearly_averages_file)

# %% Calculate seasonal means (full GLEAM4 dataset)
logging.info("Computing average seasonal values for full GLEAM4 dataset")
ds_seasonal_mean_file_list = []
for var in variables:
    logging.info(f"Calculating seasonal sums/averages for {var}")
    # Calculate average/sum for every season of every year
    if var in non_summing_vars:
        ds_seasonal_ = ds_monthly[var].resample(time="QS-DEC").mean("time", skipna=True)
    else:
        ds_seasonal_ = ds_monthly[var].resample(time="QS-DEC").sum("time", skipna=True)
        if ds_seasonal_.units == "mm.month-1":
            ds_seasonal_.attrs["units"] = "mm.season-1"
    # Now take averages over all years
    ds_seasonal_mean = (
        ds_seasonal_.groupby("time.season").mean("time", skipna=True).compute()
    )
    # Write out temporary file (better for memory management!)
    file_name = f"{var}.nc"
    ds_seasonal_mean.to_netcdf(folder_processed / model_version / file_name)
    ds_seasonal_mean.close()
    ds_seasonal_mean_file_list.append(folder_processed / model_version / file_name)

# Combine all variables into one dataset
ds_seasonal_mean = xr.open_mfdataset(ds_seasonal_mean_file_list)
ds_seasonal_mean = ds_seasonal_mean.where(~water_mask)  # apply water mask
ds_seasonal_mean.attrs["years_used_for_averaging"] = (
    f"{ds_monthly.time.min().dt.year.values} - {ds_monthly.time.max().dt.year.values}"
)
ds_seasonal_mean.to_netcdf(folder_processed / seasonal_averages_file)

# %% Find common period for all datasets
ds_compare_list = [ds_era5, ds_fluxcom, ds_gleam38, ds_yearly]
year_min = str(max([ds_.time.min() for ds_ in ds_compare_list]).dt.year.values)
year_max = str(min([ds_.time.max() for ds_ in ds_compare_list]).dt.year.values)
logging.info(f"Common period for comparison: {year_min} - {year_max}")


# %% Converting datasets to common units and resolutions for comparison of average yearly evaporation
## ERA5
logging.info("Computing average yearly evaporation for...")
logging.info("ERA5-Land...")
with xr.set_options(keep_attrs=True):
    da_era5_yearly_mean_compare = (
        (
            ds_era5["total_evaporation_sum"].sel(time=slice(year_min, year_max)) * 1000
        )  # m -> mm
        .resample(time="1YE")
        .sum()  # Take yearly sums
        .mean("time", skipna=True)  # Average over the years
        .reindex(lat=ds_yearly.lat, method="nearest")  # Align latitudes with GLEAM4
        .assign_coords(
            lat=ds_yearly.lat, lon=ds_yearly.lon
        )  # Make sure coords are floating point equal
        .compute()
    )
    da_era5_yearly_mean_compare.attrs["units"] = "mm"


## GLEAM38
logging.info("GLEAM3.8a...")
da_gleam38_yearly_mean_compare = (
    ds_gleam38["E"]
    .sel(time=slice(year_min, year_max))
    .interp(lon=ds_yearly.lon, lat=ds_yearly.lat, method="linear")  # 0.25° to 0.1°
    .mean("time", skipna=True)
    .compute()
)

## GLEAM4
logging.info("GLEAM4...")
da_gleam4_yearly_mean_compare = (
    ds_yearly["E"]
    .sel(time=slice(year_min, year_max))
    .mean("time", skipna=True)
    .compute()
)

## FLUXCOM-X
logging.info("FLUXCOM-X...")
da_fluxcom_yearly_mean_compare = (
    (
        ds_fluxcom["ET"]
        .sel(time=slice(year_min, year_max))
        .coarsen(lat=2, lon=2)  # 0.05° to 0.1°
        .mean()  # Take average over the 4 0.05° pixels for the 0.1° pixels
        * 24  # mm/hr -> mm/day
        * ds_fluxcom.time.dt.days_in_month  # mm/day -> mm/month
    )
    .resample(time="1YE")
    .sum()  # Sum total evaporation per year (mm/month -> mm/year)
    .mean("time", skipna=True)  # Take average over years
    .assign_coords(
        lon=ds_yearly.lon, lat=ds_yearly.lat
    )  # Align coordinates with GLEAM4
    .compute()
)

## Combine in 1 data-cube
logging.info("Combining datasets in 1 cube")
logging.info(f"List of models: {model_list_pattern}")
ds_comparison_yearly = xr.concat(
    [
        da_gleam4_yearly_mean_compare,
        da_gleam38_yearly_mean_compare,
        da_era5_yearly_mean_compare,
        da_fluxcom_yearly_mean_compare,
    ],
    dim=xr.DataArray(
        model_list_pattern,
        dims="product",
        name="product",
    ),
)
# Apply water mask
ds_comparison_yearly = ds_comparison_yearly.where(~water_mask)
# Add info on period used
ds_comparison_yearly.attrs["years_considered"] = f"{year_min} - {year_max}"
# Write out
ds_comparison_yearly.to_netcdf(folder_processed / yearly_averages_file_comparison)
