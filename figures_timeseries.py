# %% Imports
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import xarray as xr

from conf import color_dict, folder_figures, folder_insitu, site_selection
from functions import txt_to_netcdf

folder_figures.mkdir(exist_ok=True)

# %% Read in data
ds_validation_path = folder_insitu / "sites_validation.nc"
if not ds_validation_path.exists():
    ds_validation = txt_to_netcdf(folder_insitu)
    ds_validation.to_netcdf(ds_validation_path)
else:
    ds_validation = xr.open_dataset(ds_validation_path)

# %% Reformat data for plotting
var_selection = list(color_dict.keys())
var_selection.append("insitu")
ds_validation = ds_validation.sel(site=list(site_selection.keys()))[var_selection]

# %% Single year plot
fig, axes = plt.subplots(len(site_selection), 1, figsize=(8, 8))
for i, site in enumerate(site_selection.keys()):
    # Add models to figure
    ds_validation_single_year = ds_validation.sel(time=site_selection[site])
    for model, color in color_dict.items():
        ds_validation_single_year.sel(site=site)[model].plot(
            ax=axes[i],
            label=model,
            linewidth=0.8,
            color=color,
        )
    # Add in situ data
    axes[i].fill_between(
        x=ds_validation_single_year.time,
        y1=0,
        y2=ds_validation_single_year.sel(site=site)["insitu"],
        color="lightgrey",
        label="In situ",
    )
    axes[i].set_xlabel("")
    axes[i].set_ylabel(r"$E$ [mm/day]")

handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, loc="upper center", ncol=len(var_selection))
plt.tight_layout()
plt.subplots_adjust(top=0.9)
fig.savefig(folder_figures / "GLEAM4_timeseries_sites_1year.pdf")

# %% Calculating anomalies
ds_doy = ds_validation.groupby("time.dayofyear")
ds_climatology = ds_doy.mean("time")
ds_anomalies = ds_doy - ds_climatology["insitu"]

# %% 1 Year anomalies
fig, axes = plt.subplots(len(site_selection), 1, figsize=(8, 8))
for i, site in enumerate(site_selection):
    year = "2005"
    if site == "PE-QFR":
        year = "2019"
    for model, color in color_dict.items():
        ds_anomalies.sel(site=site, time=year)[model].plot(
            ax=axes[i],
            label=model,
            linewidth=0.8,
            color=color,
        )
    # Add in situ data
    da_anom_inistu = ds_anomalies.sel(site=site, time=year)["insitu"]
    axes[i].fill_between(
        x=da_anom_inistu.time,
        y1=0,
        y2=da_anom_inistu,
        color="lightgrey",
        label="In situ",
    )
    axes[i].set_xlabel("")
    axes[i].set_ylabel(r"$E$ anomalies [mm/day]")
handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, loc="upper center", ncol=len(var_selection))
plt.tight_layout()
plt.subplots_adjust(top=0.9)

# %% Model deviation from reality
ds_dev = ds_validation - ds_validation["insitu"]
fig, axes = plt.subplots(len(site_selection), 1, figsize=(8, 8))
for i, site in enumerate(site_selection):
    year = "2005"
    if site == "PE-QFR":
        year = "2019"
    for model, color in color_dict.items():
        ds_dev.sel(site=site, time=year)[model].plot(
            ax=axes[i],
            label=model,
            linewidth=0.8,
            color=color,
        )
    axes[i].set_xlabel("")
    axes[i].set_ylabel(r"$E$ model - in situ [mm/day]")
handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, loc="upper center", ncol=len(var_selection))
plt.tight_layout()
plt.subplots_adjust(top=0.9)

# %%
