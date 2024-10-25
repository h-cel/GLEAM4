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
ds_validation = ds_validation.sel(site=site_selection)[var_selection]
# df_validation = (
#     ds_validation.to_dataframe()
#     .reset_index()
#     .melt(
#         id_vars=["time", "site"],
#         value_vars=var_selection,
#         var_name="type",
#         value_name="E",
#     )
# )
# %% Plot
model_colours = ["#09070d", "#09070d", "#6f6db1", "#cf6666", "#4a8740"]
fig, axes = plt.subplots(len(site_selection), 1, figsize=(8, 8))
nr_sites = len(site_selection)
for i, site in enumerate(site_selection):
    for j, var in enumerate(var_selection):
        if var == "insitu":
            linestyle = "dashed"
            alpha = 1
        else:
            linestyle = "solid"
            alpha = 0.9
        ds_validation.sel(site=site)[var].plot(
            ax=axes[i],
            label=var,
            linewidth=0.8,
            color=model_colours[j],
            linestyle=linestyle,
            alpha=alpha,
        )
        if var == "insitu":
            xlimits = axes[i].get_xlim()
    axes[i].set_xlim(xlimits)
    axes[i].set_xlabel("")
    axes[i].set_ylabel(r"$E$ [mm/day]")
axes[nr_sites - 1].set_xlabel("Time")
handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, loc="upper center", ncol=len(var_selection))

plt.tight_layout()
plt.subplots_adjust(
    top=0.9
)  # Adjust the top to make space for the legend.tight_layout()
# fig.savefig(folder_figures / "GLEAM4_timeseries_sites.pdf")

# %% Single year plot
fig, axes = plt.subplots(len(site_selection), 1, figsize=(8, 8))
for i, site in enumerate(site_selection):
    # Add models to figure
    if site == "PE-QFR":
        year = "2019"
    else:
        year = "2005"
    ds_validation_single_year = ds_validation.sel(time=year)
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

# %% Multi year anomalies
fig, axes = plt.subplots(len(site_selection), 1, figsize=(8, 8))
for i, site in enumerate(site_selection):
    for model, color in color_dict.items():
        ds_anomalies.sel(site=site)[model].plot(
            ax=axes[i],
            label=model,
            linewidth=0.8,
            color=color,
        )
    # Add in situ data
    da_anom_inistu = ds_anomalies.sel(site=site)["insitu"].dropna(dim="time")
    axes[i].fill_between(
        x=da_anom_inistu.time,
        y1=0,
        y2=da_anom_inistu,
        color="lightgrey",
        label="In situ",
    )
    axes[i].set_xlabel("")
    axes[i].set_ylabel(r"$E$ anomalies [mm/day]")
    axes[i].set_xlim([da_anom_inistu.time[0], da_anom_inistu.time[-1]])
handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, loc="upper center", ncol=len(var_selection))
plt.tight_layout()
plt.subplots_adjust(top=0.9)

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
