# Script to analyse global GLEAM4 patterns
# Author: Olivier Bonte, Ghent University, olivier.bonte@ugent.be
# October 2025

# %% Imports and paths
import importlib

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

import conf
import functions

importlib.reload(conf)
importlib.reload(functions)

from conf import (
    color_dict_lat,
    folder_figures,
    folder_processed,
    model_version,
    seasonal_averages_file,
    yearly_averages_file,
)
from functions import get_area_from_dataset, make_custom_cmap

folder_figures.mkdir(exist_ok=True)

# %% Read in data
ds_yearly_mean = xr.open_dataset(folder_processed / yearly_averages_file)
ds_seasonal_mean = xr.open_dataset(folder_processed / seasonal_averages_file)

# %% Last processing
# Define "other" as condensation, subliation and open water evaporation
ds_yearly_mean["Eo"] = (
    ds_yearly_mean["Ec"] + ds_yearly_mean["Es"] + ds_yearly_mean["Ew"]
)
# Take mean across all latitudes
ds_yearly_mean_lon_mean = ds_yearly_mean.mean(dim="lon")

# %%  Figure 2: Yearly averages and total
vmin_E = 0  # mm/year
vmax_E = 1400  # mm/year
cmap_E = make_custom_cmap(vmin_E, vmax_E)

# All in one Figure
# size height = approx 2/3 of A4 paper (8.27 x 11.69 inches)
fontsize_in_figure = 12
ocean_color = "#CDE0E4"
ocean_alpha = 0.5
with plt.rc_context({"font.size": 7}):
    fig = plt.figure(figsize=(7, 2 / 3 * 11.69))
    axd = fig.subplot_mosaic(
        """
        ABBB
        CCDD
        EEFF
        """,
        sharey=True,
        height_ratios=[0.45, 0.275, 0.275],
    )
    # Figure A (mean over longitudes)
    lat_plot = ds_yearly_mean_lon_mean.lat
    axd["A"].fill_betweenx(
        lat_plot,
        0,
        ds_yearly_mean_lon_mean["Et"],
        color=color_dict_lat["transpiration"],
        label=r"$E_t$",
    )
    axd["A"].fill_betweenx(
        lat_plot,
        ds_yearly_mean_lon_mean["Et"],
        ds_yearly_mean_lon_mean["Et"] + ds_yearly_mean_lon_mean["Eb"],
        color=color_dict_lat["bare soil evaporation"],
        label=r"$E_b$",
    )
    axd["A"].fill_betweenx(
        lat_plot,
        ds_yearly_mean_lon_mean["Et"] + ds_yearly_mean_lon_mean["Eb"],
        ds_yearly_mean_lon_mean["Et"]
        + ds_yearly_mean_lon_mean["Eb"]
        + ds_yearly_mean_lon_mean["Ei"],
        color=color_dict_lat["interception"],
        label=r"$E_i$",
    )
    axd["A"].fill_betweenx(
        lat_plot,
        ds_yearly_mean_lon_mean["Et"]
        + ds_yearly_mean_lon_mean["Eb"]
        + ds_yearly_mean_lon_mean["Ei"],
        ds_yearly_mean_lon_mean["E"],
        color=color_dict_lat["other"],
        label=r"$E_o$",
    )
    axd["A"].invert_xaxis()
    axd["A"].legend()

    # Figure B (total evaporation)
    img = ds_yearly_mean["E"].plot.imshow(
        vmin=vmin_E, vmax=vmax_E, cmap=cmap_E, add_colorbar=False, ax=axd["B"]
    )
    axd["B"].set_title("")
    axd["B"].set_xlabel("Longitude [°]")
    axd["B"].set_ylabel("")
    axd["B"].text(-170, -65, r"$E$", fontsize=fontsize_in_figure)
    axd["B"].set_facecolor(ocean_color)

    # Figure C (transpiration)
    imgB = ds_yearly_mean["Et"].plot.imshow(
        vmin=vmin_E, vmax=vmax_E, cmap=cmap_E, add_colorbar=False, ax=axd["C"]
    )
    axd["C"].set_title("")
    axd["C"].set_xlabel("")
    axd["C"].set_ylabel("Latitude [°]")
    axd["C"].text(-170, -65, r"$E_t$", fontsize=fontsize_in_figure)
    axd["C"].set_facecolor(ocean_color)

    # Figure D (interception)
    imgC = ds_yearly_mean["Ei"].plot.imshow(
        vmin=vmin_E, vmax=vmax_E, cmap=cmap_E, add_colorbar=False, ax=axd["D"]
    )
    axd["D"].set_title("")
    axd["D"].set_xlabel("")
    axd["D"].set_ylabel("")
    axd["D"].text(-170, -65, r"$E_i$", fontsize=fontsize_in_figure)
    axd["D"].set_facecolor(ocean_color)

    # Figure E (bare soil evaporation)
    imgE = ds_yearly_mean["Eb"].plot.imshow(
        vmin=vmin_E, vmax=vmax_E, cmap=cmap_E, add_colorbar=False, ax=axd["E"]
    )
    axd["E"].set_title("")
    axd["E"].set_xlabel("Longitude [°]")
    axd["E"].set_ylabel("")
    axd["E"].text(-170, -65, r"$E_b$", fontsize=fontsize_in_figure)
    axd["E"].set_facecolor(ocean_color)

    # Figure F (other evaporation)
    imgF = ds_yearly_mean["Eo"].plot.imshow(
        vmin=vmin_E, vmax=vmax_E, cmap=cmap_E, add_colorbar=False, ax=axd["F"]
    )
    axd["F"].set_title("")
    axd["F"].set_xlabel("Longitude [°]")
    axd["F"].set_ylabel("")
    axd["F"].text(-170, -65, r"$E_o$", fontsize=fontsize_in_figure)
    axd["F"].set_facecolor(ocean_color)

    # Colorbars
    fig.subplots_adjust(bottom=0.1)
    cbar_ax = fig.add_axes([0.3, 0.02, 0.4, 0.02])  # [left, bottom, width, height]
    cbar = fig.colorbar(img, cax=cbar_ax, orientation="horizontal")
    cbar.set_label("mm/year")

# Save figure to png
fig.savefig(
    folder_figures / f"{model_version}_yearly_averages.png",
    dpi=900,
    bbox_inches="tight",
)
fig

# %%
