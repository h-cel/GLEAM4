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
    violin_fill_color,
    zenodo_doi,
)
from functions import download_zenodo, kge_agg_func, plot_taylor, txt_to_netcdf

folder_figures.mkdir(exist_ok=True)

# %% Download data (if necessary)
download_zenodo(zenodo_doi, "data_test")

# %% Read in data
ds_validation_path = folder_insitu / "sites_validation.nc"
if not ds_validation_path.exists():
    ds_validation = txt_to_netcdf(folder_insitu)
    ds_validation.to_netcdf(ds_validation_path)
else:
    ds_validation = xr.open_dataset(ds_validation_path)

# %% Calculate KGE
metric_list = []
models_list = ["GLEAM4", "GLEAM v3.8", "ERA5-Land", "FLUXCOM"]
for model in models_list:
    df_val_tmp = (
        ds_validation.to_dataframe()
        .groupby("site")
        .apply(kge_agg_func, "insitu", model)
    )
    df_val_tmp["model"] = model
    df_val_tmp = df_val_tmp.reset_index().set_index(["site", "model"])
    metric_list.append(df_val_tmp)
df_metrics = pd.concat(metric_list).reset_index().set_index("site")

# %% Rearange data for plotting
plot_list = []
alternative_models_list = models_list.copy()
alternative_models_list.remove("GLEAM4")
df_GLEAM4 = df_metrics[df_metrics["model"] == "GLEAM4"]
for model in alternative_models_list:
    df_temp = df_metrics[df_metrics["model"] == model]
    df_temp["model type"] = "Alternative model"
    df_temp2 = df_GLEAM4.copy()
    df_temp2["model"] = model
    df_temp2["model type"] = "GLEAM4"
    plot_list.append(pd.concat([df_temp, df_temp2]))
df_plot = pd.concat(plot_list).reset_index()

# %% General violinplot settings
split = True
inner = "quart"
bw_adjust = 0.6
density_norm = "area"

# %% Violinplot for KGE (Three seperate)
fig, ax = plt.subplots(figsize=(8, 5))
sns.violinplot(
    data=df_plot,
    x="model",
    y="KGE",
    hue="model type",
    ax=ax,
    split=split,
    palette="Set2",
    inner=inner,
    bw_adjust=bw_adjust,
    density_norm=density_norm,
)
ax.axhline(-0.41, color="red")
ax.legend()
ax.set_xlabel("Alternative model")
ax.set_ylim(-1)

# %% Additional visualisation: Exclude outliers for KDE estimate
Q1 = df_plot["KGE"].quantile(0.25)
Q3 = df_plot["KGE"].quantile(0.75)
IQR = Q3 - Q1
lower_bound = Q1 - 1.5 * IQR
upper_bound = Q3 + 1.5 * IQR

df_plot_nofliers = df_plot[
    (df_plot["KGE"] >= lower_bound) & (df_plot["KGE"] <= upper_bound)
]
print()
fig, ax = plt.subplots()
sns.violinplot(
    data=df_plot_nofliers,
    x="model",
    y="KGE",
    hue="model type",
    ax=ax,
    split=True,
    palette="Set2",
    inner="quart",
    legend="brief",
    bw_adjust=0.8,
)
ax.axhline(-0.41, color="red")
ax.legend()
ax.set_xlabel("Alternative model")

# %% Final figure: Violinplot for KGE, one with all models combined
fig, ax = plt.subplots(figsize=(4.5, 5))
# Add other models
for model in alternative_models_list:
    sns.violinplot(
        data=df_plot[df_plot["model"] == model].sort_values(
            "model type", ascending=False
        ),
        y="KGE",
        hue="model type",
        ax=ax,
        split=split,
        palette={
            "Alternative model": color_dict[model],
            "GLEAM4": color_dict["GLEAM4"],
        },
        inner=inner,
        fill=False,
        bw_adjust=bw_adjust,
        legend=False,
        density_norm=density_norm,
        common_norm=True,
    )
# Add fill for GLEAM4
sns.violinplot(
    data=df_metrics[df_metrics["model"] == "GLEAM4"][["KGE"]],
    y="KGE",
    hue=1,
    hue_order=[1, 2],
    ax=ax,
    dodge=True,
    split=split,
    inner=None,
    legend=False,
    palette={1: violin_fill_color, 2: "#FFFFFF"},
    bw_adjust=bw_adjust,
    density_norm=density_norm,
)
ax.axhline(-0.41, color="red")
ax.set_ylim(-1, 1.1)
plt.tight_layout()
fig.savefig(folder_figures / "GLEAM4_violinplot.pdf")

# %% Final figure: Taylor plot
plot_taylor(folder_insitu, folder_figures)

# %%
