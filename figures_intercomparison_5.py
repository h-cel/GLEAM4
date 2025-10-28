# Script to compare GLEAM4 to alternative datasets
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
    E_plot_dict,
    color_dict_pattern,
    folder_figures,
    folder_processed,
    model_list_pattern,
    model_version,
    ocean_color,
    yearly_averages_file_comparison,
)
from functions import make_custom_cmap

# %%
da_comp = xr.open_dataarray(folder_processed / yearly_averages_file_comparison)
# da_comp = xr.open_dataarray(folder_processed / "comparison_temp.nc")
# GLEAM4 - other products
model_list_others = [p for p in model_list_pattern if p != "GLEAM4"]
# other_products_list = list(da_comp.product.to_pandas().drop("GLEAM4"))
da_comp_diff = da_comp.sel(product="GLEAM4") - da_comp.sel(product=model_list_others)

#

# %% Figure 5: global comparison plot
font_settings = {"size": 6}
plt.rc("font", **font_settings)
fig = plt.figure(
    figsize=(7, 1 / 2 * 11.69), layout="constrained"
)  # , layout="constrained")  # tight_layout=True)
ax = fig.subplot_mosaic(
    """
    AADDG
    BBEEG
    CCFFG
    """,
)

# A -> C: the other products
cmap_E = make_custom_cmap(E_plot_dict["vmin"], E_plot_dict["vmax"])
for letter, product in zip(["A", "B", "C"], model_list_others):
    img_E = da_comp.sel(product=product).plot.imshow(
        ax=ax[letter],
        add_colorbar=False,
        cmap=cmap_E,
        vmin=E_plot_dict["vmin"],
        vmax=E_plot_dict["vmax"],
    )
    # Only longitude for bottom plot
    if letter != "C":
        ax[letter].tick_params(labelbottom=False)
        ax[letter].set_xlabel("")
        # ax[letter].label_outer()
    else:
        ax[letter].set_xlabel("Longitude")
    # Always latitude set
    ax[letter].set_ylabel("Latitude")
    ax[letter].set_title("")
    # Add model name as text
    ax[letter].text(-170, -65, product, color=color_dict_pattern[product])
    # Set ocean color
    ax[letter].set_facecolor(ocean_color)

# D -> F: GLEAM4 - other products
for letter, product in zip(["D", "E", "F"], model_list_others):
    img_diff = da_comp_diff.sel(product=product).plot.imshow(
        ax=ax[letter], add_colorbar=False, cmap="RdBu", vmin=-300, vmax=300
    )
    if letter != "F":
        ax[letter].tick_params(labelbottom=False)
        ax[letter].set_xlabel("")
    else:
        ax[letter].set_xlabel("Longitude")
    ax[letter].set_ylabel("")
    ax[letter].tick_params(labelleft=False)
    ax[letter].set_title("")
    # Add models name as text
    ax[letter].text(-170, -65, "GLEAM4 -", color=color_dict_pattern["GLEAM4"])
    ax[letter].text(-108, -65, product, color=color_dict_pattern[product])
    # Set ocean color
    ax[letter].set_facecolor(ocean_color)

# Adding colorbars
# Key tip: use "constrained layout" https://matplotlib.org/stable/users/explain/axes/colorbar_placement.html
fig.colorbar(img_E, ax=ax["C"], location="bottom", shrink=0.8, label="[mm/year]")
fig.colorbar(img_diff, ax=ax["F"], location="bottom", shrink=0.8, label="[mm/year]")

# G: mean over longitudes
for product in model_list_pattern:
    da_comp.sel(product=product).mean(dim="lon").plot(
        y="lat",
        ax=ax["G"],
        color=color_dict_pattern[product],
        label=product,
        linewidth=1,
    )
legend = ax["G"].legend(frameon=False, handlelength=0)
for text, product in zip(legend.get_texts(), model_list_pattern):
    text.set_color(color_dict_pattern[product])
ax["G"].set_title("")
ax["G"].set_ylabel("Latitude")
ax["G"].set_xlabel("[mm/year]")
# Save to high quality png
fig.savefig(
    folder_figures / f"{model_version}_global_product_intercomparison_fig_5.png",
    dpi=900,
)
# %%
