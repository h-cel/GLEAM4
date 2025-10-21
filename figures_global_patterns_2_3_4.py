# Script to analyse global GLEAM4 patterns
# Author: Olivier Bonte, Ghent University, olivier.bonte@ugent.be
# October 2025

# %% Imports and paths
import matplotlib.pyplot as plt
import xarray as xr

from conf import (
    color_dict_lat,
    folder_figures,
    folder_processed,
    seasonal_averages_file,
    yearly_averages_file,
)
from functions import make_custom_cmap

folder_figures.mkdir(exist_ok=True)

# %% Read in data
ds_yearly_mean = xr.open_dataset(folder_processed / yearly_averages_file)
ds_seasonal_mean = xr.open_dataset(folder_processed / seasonal_averages_file)

# %% List variables and years available
# variables = [folder.name for folder in folder_gridded_monthly.iterdir()]
# print(variables)
# years = [
#     re.search(r"\d{4}", file.name).group()
#     for file in (folder_gridded_monthly / variables[0]).iterdir()
# ]
# years.sort()
# print(years)

# # %% Check redo and if temp data exists
# if REDO:
#     redo_bool = True
# else:
#     redo_bool = (
#         not (folder_temp_data / seasonal_averages_file).exists()
#         or not (folder_temp_data / yearly_averages_file).exists()
#     )


# # %% Read in data

# if redo_bool:
#     list_of_datastets = []
#     for folder in [folder_gridded_monthly, folder_gridded_yearly]:
#         ds_list = []
#         print(f"Reading data from {folder}")
#         for var in variables:
#             ds_var = xr.open_mfdataset(
#                 sorted((folder / var).iterdir()),
#                 concat_dim="time",
#                 combine="nested",
#                 data_vars="minimal",
#                 coords="minimal",
#                 compat="override",
#             )  # see https://docs.xarray.dev/en/stable/user-guide/io.html#reading-multi-file-datasets
#             ds_var = ds_var.chunk(
#                 {"time": 12, "lon": ds_var.lon.size / 4, "lat": ds_var.lat.size / 4}
#             )  # To 1 chunk per year, 4 spatial chunks
#             ds_list.append(ds_var)
#             # ds_list_year.append(ds_var_year)
#         ds = xr.merge(ds_list)
#         list_of_datastets.append(ds)
#     ds_monthly, ds_yearly = list_of_datastets


# # %% Calculate yearly mean evaporation

# if redo_bool:
#     ds_yearly_mean = ds_yearly.mean("time", skipna=True).compute()
#     ds_yearly_mean.to_netcdf(folder_temp_data / yearly_averages_file)
# else:  # Read from file
#     ds_yearly_mean = xr.open_dataset(folder_temp_data / yearly_averages_file)

# # %% Calculate seasonal means
# if redo_bool:
#     ds_seasonal_mean_file_list = []
#     for var in variables:
#         print(f"Calculating seasonal sums/averages for {var}")
#         # Calculate average/sum for every season of every year
#         if var in non_summing_vars:
#             ds_seasonal_ = (
#                 ds_monthly[var].resample(time="QS-DEC").mean("time", skipna=True)
#             )
#         else:
#             ds_seasonal_ = (
#                 ds_monthly[var].resample(time="QS-DEC").sum("time", skipna=True)
#             )
#             if ds_seasonal_.units == "mm.month-1":
#                 ds_seasonal_.attrs["units"] = "mm.season-1"
#         # Now take averages over all years
#         ds_seasonal_mean = (
#             ds_seasonal_.groupby("time.season").mean("time", skipna=True).compute()
#         )
#         # Write out temporary file (better for memory management!)
#         file_name = f"{var}.nc"
#         ds_seasonal_mean.to_netcdf(folder_temp_data / file_name)
#         ds_seasonal_mean_file_list.append(folder_temp_data / file_name)

#     # Combine all variables into one dataset
#     ds_seasonal_mean = xr.open_mfdataset(ds_seasonal_mean_file_list)
#     ds_seasonal_mean.to_netcdf(folder_temp_data / seasonal_averages_file)
# else:  # Read from file
#     ds_seasonal_mean = xr.open_dataset(folder_temp_data / seasonal_averages_file)


# %%  Figure 2: Yearly averages and total
vmin_E = 0  # mm/year
vmax_E = 1400  # mm/year
cmap_E = make_custom_cmap(vmin_E, vmax_E)

# Define "other" as condensation, subliation and open water evaporation
ds_yearly_mean["Eo"] = (
    ds_yearly_mean["Ec"] + ds_yearly_mean["Es"] + ds_yearly_mean["Ew"]
)
# Take mean across all latitudes
ds_yearly_mean_lon_mean = ds_yearly_mean.mean(dim="lon")

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
    folder_figures / f"{yearly_averages_file.split('_')[0]}_yearly_averages.png",
    dpi=900,
    bbox_inches="tight",
)
fig

# %%
