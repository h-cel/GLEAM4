# %% Imports
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from functions import extract_site_code, kge_agg_func

# %% Set paths
folder_insitu = Path("data/sites")
figures_folder = Path("figures")
figures_folder.mkdir(exist_ok=True)

# %% Read in data
pd_list = []
for site_file in folder_insitu.iterdir():
    site_name = extract_site_code(site_file.name)
    df_in_situ = pd.read_csv(site_file, index_col=0)
    df_in_situ["site"] = site_name
    df_in_situ = df_in_situ.reset_index().set_index(["time", "site"])
    pd_list.append(df_in_situ)
ds_validation = pd.concat(pd_list).to_xarray()
# Rename for consistency
ds_validation = ds_validation.rename({"GLEAM38": "GLEAM v3.8", "ERA5": "ERA5-Land"})

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

# %% Violinplot for KGE
fig, ax = plt.subplots(figsize=(8, 5))
sns.violinplot(
    data=df_plot,
    x="model",
    y="KGE",
    hue="model type",
    ax=ax,
    split=True,
    palette="Set2",
    inner="quart",
    bw_adjust=0.6,
)
ax.axhline(-0.41, color="red")
ax.legend()
ax.set_xlabel("Alternative model")
ax.set_ylim(-1)
fig.savefig("figures/GLEAM4_violinplots.pdf")

# %% Exclude outliers for KDE estimate
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

# %% Boxplot for KGE
fig, ax = plt.subplots()
sns.boxplot(
    data=df_plot,
    x="model",
    y="KGE",
    hue="model type",
    palette="Set2",
    showfliers=False,
    notch=True,
    ax=ax,
)
sns.stripplot(
    data=df_plot,
    x="model",
    y="KGE",
    hue="model type",
    dodge=True,
    jitter=True,
    palette="dark:black",
    size=2,
    ax=ax,
    alpha=0.5,
    legend=False,
)
ax.axhline(-0.41, color="red")
ax.set_xlabel("Alternative model")
ax.legend()
ax.set_ylim(-1)

# %%
