# Script to analyse global GLEAM4 patterns
# Author: Olivier Bonte, Ghent University, olivier.bonte@ugent.be
# October 2025

# %% Imports and paths
import importlib

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

import conf
import functions

importlib.reload(conf)
importlib.reload(functions)

from conf import (
    color_dict_lat,
    folder_figures,
    folder_processed,
    model_version,
    ocean_color,
    season_choice,
    season_plot_var_dict,
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

# Calculate drivers of evaporation
var_drivers = ["E_deficit", "E_aero", "E_rad"]
ds_yearly_mean[var_drivers[0]] = ds_yearly_mean["Ep"] - ds_yearly_mean["E"]
ds_yearly_mean[var_drivers[1]] = (
    ds_yearly_mean["Ep_aero"] / ds_yearly_mean["Ep"] * ds_yearly_mean["E"]
)
ds_yearly_mean[var_drivers[2]] = (
    ds_yearly_mean["Ep_rad"] / ds_yearly_mean["Ep"] * ds_yearly_mean["E"]
)
## RGB color composite based on drivers
# Idea:
# - Color = relative contribution of driver per pixel. Per pixel: driver / sum_of_drivers
# - Brightness = Ep / maximum(Ep) per pixel
# maximum(Ep) can also be a quantile -> clip rest of data to this quantile
Ep_Q_95 = ds_yearly_mean["Ep"].quantile(0.95)
Ep_clipped = ds_yearly_mean["Ep"].clip(max=Ep_Q_95)
brightness = Ep_clipped / Ep_Q_95
# Sum of all drivers per pixel
da_sum_drivers = (
    ds_yearly_mean["E_deficit"] + ds_yearly_mean["E_aero"] + ds_yearly_mean["E_rad"]
)
# Contribution per pixel
redn = ds_yearly_mean["E_deficit"].clip(min=0) / da_sum_drivers * brightness
greenn = ds_yearly_mean["E_aero"].clip(min=0) / da_sum_drivers * brightness
bluen = ds_yearly_mean["E_rad"].clip(min=0) / da_sum_drivers * brightness
rgb_array = np.dstack((redn, greenn, bluen)).clip(max=1)
# Make a DataArray
da_rgb = xr.DataArray(
    data=rgb_array,
    dims=["lat", "lon", "color"],
    coords={
        "lat": ds_yearly_mean.lat.values,
        "lon": ds_yearly_mean.lon.values,
        "color": ["red_deficit", "green_aero", "blue_rad"],
    },
)


# %%  Figure 2: Yearly averages and total
vmin_E = 0  # mm/year
vmax_E = 1400  # mm/year
cmap_E = make_custom_cmap(vmin_E, vmax_E)

# All in one Figure
# size height = approx 2/3 of A4 paper (8.27 x 11.69 inches)
fontsize_in_figure = 12
with plt.rc_context({"font.size": 7}):
    fig = plt.figure(figsize=(7, 2 / 3 * 11.69), layout="constrained")
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
    fig.colorbar(
        img, ax=[axd["E"], axd["F"]], location="bottom", label="mm/year", shrink=0.4
    )
# Save figure to png
fig.savefig(
    folder_figures / f"{model_version}_yearly_averages.png",
    dpi=900,
    bbox_inches="tight",
)
fig

# %% Figure 2 bis: totals and fractions
da_area = get_area_from_dataset(ds_yearly_mean)
E_vars = [v for v in ds_yearly_mean.data_vars if v.startswith("E")]
print(f"E variables: {E_vars}")

ds_summed = (
    (ds_yearly_mean[E_vars] / 1000) * da_area
).sum() / 10**9  # (mm -> m) x m^2 -> m^3 -> km^3
ds_summed.attrs["units"] = "km^3"
ds_percentages = ds_summed / ds_summed["E"] * 100

## pie plot
var_percentages = ["Et", "Eb", "Ei", "Eo"]

# percentages
fig, ax = plt.subplots()
ax.pie(
    ds_percentages[var_percentages].to_pandas().values,
    labels=var_percentages,
    autopct="%1.1f%%",
    colors=color_dict_lat.values(),
)
ax.set_title(f"{model_version}")
fig.savefig(folder_figures / f"{model_version}_pie_chart_percentages.png")


# absolute values
def autopct_format(pct, allvals):
    absolute = int(pct / 100.0 * sum(allvals))
    return f"{absolute:,} km³/yr"


fig, ax = plt.subplots()
ax.pie(
    ds_summed[var_percentages].to_pandas().values,
    labels=var_percentages,
    autopct=lambda pct: autopct_format(
        pct, ds_summed[var_percentages].to_pandas().values
    ),
    colors=color_dict_lat.values(),
)
ax.set_title(
    f"{model_version}: total E = {round(float(ds_summed['E'].values), 2):,} km³/yr"
)
fig.savefig(folder_figures / f"{model_version}_pie_chart_totals.png")

# %% Figure 3: Season patterns

with plt.rc_context({"font.size": 6}):
    fig, axes = plt.subplots(
        len(season_plot_var_dict),
        len(season_choice),
        figsize=(7, 11.69),
        sharex=True,
        sharey=True,
        constrained_layout=True,
    )  # A4 paper basically

    for i, var in enumerate(season_plot_var_dict):
        cmap_ = make_custom_cmap(
            season_plot_var_dict[var]["vmin"], season_plot_var_dict[var]["vmax"]
        )
        for j, season in enumerate(season_choice):
            img = (
                ds_seasonal_mean[var]
                .sel(season=season)
                .plot.imshow(
                    ax=axes[i, j],
                    cmap=cmap_,
                    vmin=season_plot_var_dict[var]["vmin"],
                    vmax=season_plot_var_dict[var]["vmax"],
                    add_colorbar=False,
                )
            )
            # Only top row title
            if i == 0:
                axes[i, j].set_title(f"{season}")
            else:
                axes[i, j].set_title("")

            # Only bottom row mentioning of longitude
            if i == len(season_plot_var_dict) - 1:
                axes[i, j].set_xlabel("Longitude [°]")
            else:
                axes[i, j].set_xlabel("")

            if j == 0:
                axes[i, j].set_ylabel("Latitude [°]")
            else:
                axes[i, j].set_ylabel("")

            axes[i, j].set_facecolor(ocean_color)

            cbar_ax = inset_axes(
                axes[i, j],
                width="4%",  # Width of colorbar
                height="60%",  # Height of colorbar
                loc="lower left",  # Location: 'upper left', 'upper right', 'lower left', 'lower right'
            )
            cbar = plt.colorbar(img, cax=cbar_ax, orientation="vertical")
            cbar.set_label(f"{var} [{season_plot_var_dict[var]['units']}]")

fig.savefig(folder_figures / f"{model_version}_seasonal_patterns.png", dpi=900)
fig

# %% Figure 4: Drivers of evaporation
# To highlight the regions with lower Ep, you can use a gamma correction
# For more info, see: https://en.wikipedia.org/wiki/Gamma_correction
gammas = [1, 0.8]
with plt.rc_context({"font.size": 7}):
    for gamma in gammas:
        fig, ax = plt.subplots(figsize=(7, 4.5), constrained_layout=True)
        (da_rgb**gamma).plot.imshow(ax=ax)
        ax.set_title(
            f"{model_version} RGB composite: Evaporative deficit (R), Aerodynamic (G), Radiative (B). γ = {gamma}"
        )
        ax.set_facecolor("black")
        ax.set_xlabel("Longitude [°]")
        ax.set_ylabel("Latitude [°]")
        fig.savefig(
            folder_figures
            / f"{model_version}_drivers_of_evaporation_gamma_{gamma}.png",
            dpi=900,
        )
fig

# %% Figure 4 bis: pie chart of drivers of evaporation

# Percentages
fig, ax = plt.subplots()
ax.pie(
    ds_summed[var_drivers].to_pandas().values,
    labels=var_drivers,
    autopct="%1.1f%%",
    # autopct=lambda pct: autopct_format(pct, ds_summed[var_drivers].to_pandas().values),
    colors=["red", "green", "blue"],
)
ax.set_title(f"{model_version}")
fig.savefig(folder_figures / f"{model_version}_pie_chart_drivers_percentages.png")

# Total values
fig, ax = plt.subplots()
ax.pie(
    ds_summed[var_drivers].to_pandas().values,
    labels=var_drivers,
    autopct=lambda pct: autopct_format(pct, ds_summed[var_drivers].to_pandas().values),
    colors=["red", "green", "blue"],
)
ax.set_title(
    f"{model_version}: Total Ep = {round(float(ds_summed['Ep'].values), 2)} km³/yr"
)
fig.savefig(folder_figures / f"{model_version}_pie_chart_drivers_totals.png")

# %%
