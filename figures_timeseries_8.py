# %% Imports
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import xarray as xr

from conf import (
    color_dict,
    folder_figures,
    folder_insitu,
    insitu_fill_color,
    site_selection,
    zenodo_doi,
)
from functions import download_zenodo, txt_to_netcdf

folder_figures.mkdir(exist_ok=True)

# %% Download data (if necessary)
download_zenodo(zenodo_doi, "data")

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

# %% Calculating anomalies
window_length = 31  # days (should be uneven)
padding_length = int(window_length / 2 - 1)
ds_clim_mwa_list = []
ds_anomalies_list = []
for site in site_selection.keys():
    # Only select timestamps with in situ data available (for each site separately)
    # for climatology
    ds_validation_tmp = ds_validation.sel(site=site).where(
        ~ds_validation.insitu.sel(site=site).isnull(), drop=True
    )
    ds_climatology_tmp = ds_validation_tmp.groupby("time.dayofyear").mean(
        "time", skipna=True
    )
    len_doy = len(ds_climatology_tmp.dayofyear)
    ds_clim_padded_tmp = xr.concat(
        [
            ds_climatology_tmp.sel(dayofyear=slice(len_doy - padding_length, len_doy)),
            ds_climatology_tmp,
            ds_climatology_tmp.sel(dayofyear=slice(1, padding_length + 1)),
        ],
        dim="dayofyear",
    )
    ds_clim_mwa_tmp = (
        (
            ds_clim_padded_tmp.rolling(dayofyear=window_length, center=True)
            .mean(skipna=True)
            .isel(dayofyear=slice(padding_length + 1, len_doy + padding_length + 1))
        )
        .reset_coords("site")
        .drop_vars("site")
    )
    ds_clim_mwa_list.append(ds_clim_mwa_tmp)

    # Select data between start and end date of `ds_validation`
    # -> no gaps in data (<-> `ds_validation_tmp`)
    ds_doy_tmp = ds_validation.sel(
        site=site, time=slice(ds_validation_tmp.time[0], ds_validation_tmp.time[-1])
    ).groupby("time.dayofyear")
    ds_anomalies_tmp = (
        (ds_doy_tmp - ds_clim_mwa_tmp)
        .reset_coords(["site", "dayofyear"])
        .drop_vars(["site", "dayofyear"])
    )
    ds_anomalies_list.append(ds_anomalies_tmp)
# Concat all sites back together
sites_index = pd.Index(site_selection.keys(), name="site")
ds_clim_mwa = xr.concat(ds_clim_mwa_list, dim=sites_index)
ds_anomalies = xr.concat(ds_anomalies_list, dim=sites_index)

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
        color="silver",
        label="In situ",
    )
    axes[i].set_xlabel("")
    axes[i].set_ylabel(r"$E$ [mm/day]")

handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, loc="upper center", ncol=len(var_selection))
plt.tight_layout()
plt.subplots_adjust(top=0.9)
fig.savefig(folder_figures / "GLEAM4_timeseries_sites_1year.pdf")

# %% 1 Year anomalies
fig, axes = plt.subplots(len(site_selection), 1, figsize=(8, 10.5))
for i, site in enumerate(site_selection.keys()):
    axins = axes[i].inset_axes(bounds=[0.72, 0.25, 0.25, 0.65])
    for model, color in color_dict.items():
        if model == "GLEAM4":
            zorder = 3
        else:
            zorder = 2
        ds_anomalies.sel(site=site, time=site_selection[site])[model].plot(
            ax=axes[i], label=model, linewidth=0.8, color=color, zorder=zorder
        )
        ds_clim_mwa.sel(site=site)[model].plot(ax=axins, color=color, zorder=zorder)
    # Add in situ data
    da_anom_inistu = ds_anomalies.sel(site=site, time=site_selection[site])["insitu"]
    axes[i].fill_between(
        x=da_anom_inistu.time,
        y1=0,
        y2=da_anom_inistu,
        color=insitu_fill_color,  # "silver",
        label="In situ",
    )
    axins.fill_between(
        x=ds_clim_mwa.dayofyear,
        y1=0,
        y2=ds_clim_mwa.sel(site=site)["insitu"],
        color=insitu_fill_color,
    )
    # ds_clim_mwa.sel(site=site)["insitu"].plot(ax=axins, color=insitu_fill_color)
    axes[i].set_xlabel("")
    axes[i].set_ylabel(r"$E$ anomalies [mm/day]")
    xlims = axes[i].get_xlim()
    axes[i].set_xlim(xlims[0], xlims[1] + 0.5 * (xlims[1] - xlims[0]))
    axes[i].set_title(site)

    axins.set_title("")
    axins.set_xlabel("DOY")
    axins.set_ylabel(r"$E_{\text{clim}}$ [mm/day]")
axes[-1].set_xlabel("Time")
handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, loc="upper center", ncol=len(var_selection))
plt.tight_layout()
plt.subplots_adjust(top=0.9)
fig.savefig(folder_figures / "GLEAM4_timeseries_anomalies_sites_1year.pdf")

# %%
