# %% Imports
# import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd
import rootutils
import scienceplots
import xarray as xr

# Set up root directory with rootutils
rootutils.setup_root(
    search_from=__file__,
    indicator=".project-root",
    pythonpath=True,  # add root directory to the PYTHONPATH (allows imports from root)
)

from conf import (
    dict_region,
    filename_monthly_anomalies_per_lat,
    filename_soi_processed,
    filename_spatial_anomaly,
    filename_yearly_anomalies,
    folder_sotc,
    soi_fill_colors,
)

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
ax.legend(frameon=False, ncol=len(dict_region), loc="upper center")
