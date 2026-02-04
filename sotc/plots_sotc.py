# %% Imports
import cartopy.crs as ccrs
import matplotlib.colors as mcolors
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd
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
    colorbar_hovmoller,
    colorbar_spatial_anomaly,
    dict_region,
    dpi,
    filename_monthly_anomalies_per_lat,
    filename_soi_processed,
    filename_spatial_anomaly,
    filename_yearly_anomalies,
    folder_figures,
    folder_sotc,
    model_version,
    soi_fill_colors,
)
from functions import make_custom_cmap_sotc

figures_sotc_path = folder_figures / "sotc"
figures_sotc_path.mkdir(parents=True, exist_ok=True)

# %% Read in files
da_yearly_anomalies = xr.open_dataarray(folder_sotc / filename_yearly_anomalies)
da_spatial_anomaly = xr.open_dataarray(folder_sotc / filename_spatial_anomaly)
da_monthly_anomalies_per_lat = xr.open_dataarray(
    folder_sotc / filename_monthly_anomalies_per_lat
)
df_soi = pd.read_csv(folder_sotc / filename_soi_processed, index_col=0)
# %% Figure x: Land evaporation anomaly
## Linear regression per region
x_year = da_yearly_anomalies.time.dt.year.values
dict_linear_fit = {}
for region in dict_region:
    p_ = np.polynomial.Polynomial.fit(
        x_year, da_yearly_anomalies.sel(region=region).values, 1
    )
    dict_linear_fit[region] = p_  # .convert().coef

## SOI processing
df_soi_sel = df_soi.iloc[df_soi.index.isin(x_year)]
df_soi_sel_neg = df_soi_sel.where(df_soi_sel < 0)
df_soi_sel_pos = df_soi_sel.where(df_soi_sel > 0)

## Actual plot
fig, ax = plt.subplots(figsize=(6.5, 3.5), layout="constrained")
# SOI ( / Nino3.4 index to be added later?)
ax2 = ax.twinx()
ax2.set_zorder(1)
ax2.fill_between(
    x_year,
    df_soi_sel["soi"].values,
    where=(df_soi_sel["soi"].values > 0),
    interpolate=True,
    color=soi_fill_colors["pos"],
)
ax2.fill_between(
    x_year,
    df_soi_sel["soi"].values,
    where=(df_soi_sel["soi"].values < 0),
    interpolate=True,
    color=soi_fill_colors["neg"],
)
# Evaporation anomalies
ax.set_zorder(2)
ax.patch.set_visible(False)
for region in dict_region:
    # Timeseries
    ax.plot(
        x_year,
        da_yearly_anomalies.sel(region=region).values,
        color=dict_region[region]["color"],
        zorder=2.5,
        label=dict_region[region]["label"],
    )
    # Trendline
    ax.plot(
        x_year,
        dict_linear_fit[region](x_year),
        color=dict_region[region]["color"],
        linestyle="dashdot",
        zorder=2.5,
    )

# x-axis
ax.xaxis.set_major_locator(ticker.MultipleLocator(5))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(1))
ax.set_xlim(
    da_yearly_anomalies.time.dt.year.min().values - 1,
    da_yearly_anomalies.time.dt.year.max().values + 0.99,
)
ax.tick_params(direction="in", which="both", top=True, labeltop=False)

# horizontal line
ax.axhline(0, color="#838382", linestyle="--", zorder=0.5, linewidth=1)

# y-axis evaporation
integer = 10
multiple_integer = np.max(
    [
        np.ceil(np.abs(da_yearly_anomalies.min() / integer)),
        np.ceil(np.abs(da_yearly_anomalies.max() / integer)),
    ]
)
ax.set_ylim(-multiple_integer * integer, multiple_integer * integer)
ax.yaxis.set_major_locator(ticker.MultipleLocator(integer))
ax.set_ylabel(r"Anomaly (mm yr$^{-1}$)")

# y-axis soi
integer = 10
multiple_integer = np.max(
    [
        np.ceil(np.abs(df_soi_sel["soi"].min() / integer)),
        np.ceil(np.abs(df_soi_sel["soi"].max() / integer)),
    ]
)
ax2.set_ylim(-multiple_integer * integer, multiple_integer * integer)
ax2.tick_params(direction="in")
for label in ax2.get_yticklabels():
    value = float(
        label.get_text().replace("−", "-")
    )  # Replace special minus with regular hyphen
    if value < 0:
        label.set_color(soi_fill_colors["neg"])
    elif value > 0:
        label.set_color(soi_fill_colors["pos"])
ax2.set_ylabel("SOI")

# legend
ax.legend(
    frameon=False,
    ncol=len(dict_region),
    loc="upper center",
    columnspacing=1,
    labelspacing=0.5,
)

#
fig.savefig(figures_sotc_path / f"fig_x_yearly_anomalies_{model_version}.png", dpi=dpi)
fig.savefig(figures_sotc_path / f"fig_x_yearly_anomalies_{model_version}.pdf")

# %% Fig y: Hovmöller diagram monthly anomalies
# colormap settings
abs_lim = colorbar_hovmoller["abs_lim"]
spacing = colorbar_hovmoller["spacing"]
cmap = colorbar_hovmoller["cmap"]
inner_ticks, cmap_mon_anom, norm = make_custom_cmap_sotc(
    abs_lim,
    spacing,
    da_monthly_anomalies_per_lat.min(),
    da_monthly_anomalies_per_lat.max(),
    cmap,
)

# plot
fig, ax = plt.subplots(figsize=(6.5, 4.5), layout="constrained")
img = da_monthly_anomalies_per_lat.plot(
    ax=ax, x="time", cmap=cmap_mon_anom, norm=norm, add_colorbar=False
)

# colorbar settings
cbar = fig.colorbar(
    img, ax=ax, location="bottom", drawedges=True, aspect=35
)  # aspect: width/height
cbar.set_ticks(inner_ticks)
cbar.set_label(r"Anomaly (mm month$^{-1}$)")
cbar.ax.tick_params(length=0)

# x-axis
ax.set_xlabel("")
ax.xaxis.set_major_locator(mdates.YearLocator(5))
ax.xaxis.set_minor_locator(mdates.YearLocator(1))
ax.set_xlim(
    left=np.datetime64(
        str(da_monthly_anomalies_per_lat.time.min().dt.year.values) + "-01-01"
    )
)  # force to show 1980
ax.tick_params(direction="in", which="both", top=True, labeltop=False)
# y-axis
ax.yaxis.set_major_locator(ticker.MultipleLocator(30))
y_ticks = np.arange(60, -90, -30)
y_tick_labels = ["60°N", "30°N", "0°", "30°S", "60°S"]
ax.yaxis.set_ticks(y_ticks)
ax.yaxis.set_ticklabels(y_tick_labels)
ax.set_ylabel("")

# Save figure
fig.savefig(
    figures_sotc_path / f"fig_y_hovmoller_monthly_anomalies_{model_version}.png",
    dpi=dpi,
)
# %% Plate x: Spatial anomaly map
# colormap settings (analogous to Hovmöller diagram)
abs_lim = colorbar_spatial_anomaly["abs_lim"]
spacing = colorbar_spatial_anomaly["spacing"]
cmap = colorbar_spatial_anomaly["cmap"]
inner_ticks, cmap_spatial_anom, norm = make_custom_cmap_sotc(
    abs_lim, spacing, da_spatial_anomaly.min(), da_spatial_anomaly.max(), cmap
)

# plot
fig, ax = plt.subplots(
    figsize=(6.5, 4.5),
    layout="constrained",
    subplot_kw={"projection": ccrs.Robinson(central_longitude=-25)},
)
img = da_spatial_anomaly.isel(time=0).plot(
    x="lon",y="lat",
    ax=ax,
    transform=ccrs.PlateCarree(),
    cmap=cmap_spatial_anom,
    norm=norm,
    add_colorbar=False,
)

# Add geographical features
ax.coastlines(linewidth=0.5)
ax.gridlines(linewidth=0.5)
# Title remove
ax.set_title("")

# colorbar
cbar = fig.colorbar(
    img,
    ax=ax,
    location="bottom",
    drawedges=True,
    aspect=35,
    label=rf"Anomalies from {base_period[0]}-{base_period[1]} (mm yr$^{{-1}}$)",
)  # aspect: width/height
cbar.set_ticks(inner_ticks)
cbar.ax.tick_params(length=0)

# Save figure
fig.savefig(
    figures_sotc_path / f"plate_x_spatial_anomalies_{model_version}.png",
    dpi=dpi,
)

# %%
