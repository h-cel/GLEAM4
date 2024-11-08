import os
import re
import subprocess
import zipfile
from pathlib import Path

import matplotlib.pylab as plt
import mpl_toolkits.axisartist.floating_axes as fa
import mpl_toolkits.axisartist.grid_finder as gf
import numpy as np
import pandas as pd
from matplotlib.colors import ListedColormap
from matplotlib.projections import PolarAxes
from scipy.stats import gaussian_kde


def extract_site_code(filename):
    """
    Function to extract the official eddy covariance site name (XX-YYY)
    form modified versions

    Paramaters
    ---------
    filename: string
        String containing the filename

    Returns
    -------
    site: string or None
        If the filename contains a valid site name, it is returned.
        Returns None otherwise
    """
    # Regular expression to match the site code pattern
    match = re.search(r"_(\w{2}-\w{3})_", filename)  # e.g. BE-Bra
    match2 = re.search(r"_(\w{2}-\w{3}_V\w{1})_", filename)  # e.g. BE-Bra_V2
    if match:
        if match2:
            site_long = match2.group(1)
            match_long = re.search(r"^(\w{2}-\w{3})(?:_V\w{1})?$", site_long)
            site = match_long.group(1)
        else:
            site = match.group(1)
    else:
        site = None
    return site


def fun_taylor(
    path_in,
    fig,
    rect,
    font_lbl,
    font_tick,
    global_min,
    global_max,
    dataset_name,
    dataset_index,
):
    """
    Function to create taylor plots

    Parameters
    ----------
    path_in: string
        Path pointing to the folder containing all .txt files per site
    fig: matplotlib.pyplot.figure
        Figure object
    rect: int
        3 number integer consisting, where last number is 0. First two
        indicate numbers of rows and cols (typically 1 and 1)
    font_lbl: int
        Size of font for label
    font_tic: int
        Size of the the font for the ticks
    global_min: float(like)
        Global minimum density over all datasets
    global_max: float(like)
        Global maximum density over all datasets
    dataset_name: string
        One of the following: "GLEAM4", "GLEAM38", "ERA5", "FLUXCOM"
    dataset_index: int
        Integer corresponding with the index of the dataset_name in the
        list given abov

    Returns
    -------
    fig: matplotlib.pyplot.figure
        Figure object
    """

    # list all cleaned pdframes of in situ data
    fn = [file for file in os.listdir(path_in) if file.endswith(".txt")]

    # colors = ['summer_r', 'pink_r', 'bone_r', 'gist_heat_r']
    colors = ["viridis", "viridis", "viridis", "viridis"]

    sample_points = []
    for s in range(0, len(fn)):
        # load pdframe with data
        pdframe_all = pd.read_csv(os.path.join(path_in, fn[s])).set_index("time")

        # remove nan values in pdframe

        pdframe_subset = pdframe_all[["insitu", "GLEAM4", "GLEAM38", "ERA5", "FLUXCOM"]]
        pdframe = pdframe_subset[~pdframe_subset.isna().any(axis=1)]

        if len(pdframe) > 1:

            # calculate statistics
            cor_pearson = pdframe.corr(method="pearson")["insitu"][dataset_name]

            std_mod = pdframe[dataset_name].std(ddof=1)
            std_obs = pdframe["insitu"].std(ddof=1)

            sdev_mod = std_mod / std_obs
            sdev_obs = std_obs / std_obs
            ccoef = cor_pearson

            # remove nan
            # sdev_mod     = sdev_mod[np.isnan(sdev_mod)==0]
            # sdev_obs     = sdev_obs[np.isnan(sdev_obs)==0]
            # ccoef        = ccoef[np.isnan(ccoef)==0]

            # Combine sample points
            stddev = sdev_mod
            corrcoef = ccoef

            sample_points.append((stddev, corrcoef))

    stdmax = 2.5

    # rect
    rows = int(str(rect)[0])
    cols = int(str(rect)[1])
    idx = int(str(rect)[2:])
    grid = plt.GridSpec(rows, cols, wspace=0.25, hspace=0.25)
    rect = grid[idx]

    # Reference std
    stdref = np.nanmean(sdev_obs)

    dia = TaylorDiagramdensity(
        sample_points,
        stdref,
        stdmax,
        font_lbl,
        font_tick,
        fig=fig,
        rect=rect,
        extend=False,
    )

    dia.samplePoints[0].set_color("r")  # Mark reference point as a red star

    # add model median
    # dia.add_sample(np.nanmedian(sdev_mod), np.nanmedian(ccoef),marker='*', ms=15, ls='', mfc='k', mew=0,zorder=5)

    # Add RMS contours, and label them
    contours = dia.add_contours(levels=5, lw=0.5, colors="0.5")  # 5 levels in grey
    plt.clabel(contours, inline=1, fontsize=font_tick, fmt="%.1f")

    dia.add_grid()  # Add grid
    dia._ax.axis[:].major_ticks.set_tick_out(True)  # Put ticks outward

    levels = 10
    # plot density
    dia.add_density(
        stdmax,
        sample_points,
        global_min,
        global_max,
        font_tick,
        levels_number=levels,
        color=colors[dataset_index],
    )

    del stdmax, sample_points
    return fig


def fun_taylor_prepare(path_in):
    """
    Calculate global_min and global_max density value for all the datasets that is
    later used as a fixed level in the density plot

    Parameters
    ---------
    path_in: string
        Path pointing to the folder containing all .txt files per site

    Returns
    -------
    global_min: floatlike
        Global minimum density over all datasets
    global_max: floatlike
        Global maximum density over all datasets

    """

    # list all cleaned pdframes of in situ data
    fn = [file for file in os.listdir(path_in) if file.endswith(".txt")]
    # list datasets
    datasets = ["GLEAM4", "GLEAM38", "ERA5", "FLUXCOM"]

    all_dataset_sample_points = []
    # sample_points = []
    for i, v in enumerate(datasets):
        sample_points = []
        for s in range(0, len(fn)):
            # load pdframe with data
            pdframe_all = pd.read_csv(os.path.join(path_in, fn[s])).set_index("time")

            # remove nan values in pdframe
            pdframe_subset = pdframe_all[
                ["insitu", "GLEAM4", "GLEAM38", "ERA5", "FLUXCOM"]
            ]
            pdframe = pdframe_subset[~pdframe_subset.isna().any(axis=1)]

            if len(pdframe) > 1:

                # calculate statistics
                cor_pearson = pdframe.corr(method="pearson")["insitu"][v]

                std_mod = pdframe[v].std(ddof=1)
                std_obs = pdframe["insitu"].std(ddof=1)

                sdev_mod = std_mod / std_obs
                ccoef = cor_pearson

                # Combine sample points
                stddev = sdev_mod
                corrcoef = ccoef

                sample_points.append((stddev, corrcoef))

        all_dataset_sample_points.append(sample_points)

    # Calculate Z values for all sample points
    all_Z_values = []
    for sample_point in all_dataset_sample_points:
        r = [std for std, corr in sample_point]
        theta = [np.arccos(corr) for std, corr in sample_point]
        kde = gaussian_kde(np.vstack([r, theta]))
        r_grid = np.linspace(0, 2, 100)
        theta_grid = np.linspace(0, np.pi, 100)
        R, Theta = np.meshgrid(r_grid, theta_grid)
        Z = kde(np.vstack([R.flatten(), Theta.flatten()])).reshape(R.shape)
        all_Z_values.append(Z)

    global_min = np.min([Z.min() for Z in all_Z_values])
    global_max = np.max([Z.max() for Z in all_Z_values])

    return global_min, global_max


def kge(y_true, y_pred, modified=False):
    """
    Calculate the Kling-Gupta Efficiency (KGE) between two time series.

    Parameters
    ----------
    y_true : np.array
        The true (often observed) time series
    y_pred : np.array
        The predicted (often simulated) time series
    modified : bool, optional
        Default is False, calculating the classic KGE, e.g. see:
        https://en.wikipedia.org/wiki/Kling%E2%80%93Gupta_efficiency
        If True, the modified KGE is calculated, see:
        https://doi.org/10.1016/j.jhydrol.2012.01.011

    Returns
    -------
    kge_acc : float
        KGE
    r : float
        Pearson correlation coefficient
    beta : float
        Ratio of the mean of the predicted to the mean of the true time series
    alpha_or_gamma: float
        Ratio of the standard deviation of the predicted to the standard deviation
        of the true time series (alpha) if modified = False. If modified = True,
        the ratio of the coefficient of variations is used instead (gamma)
    """
    # Select where not Nan
    nan_bool_true = np.isnan(y_true)
    nan_bool_pred = np.isnan(y_pred)
    nan_bool = np.logical_or(nan_bool_true, nan_bool_pred)
    y_true, y_pred = y_true[~nan_bool], y_pred[~nan_bool]

    # calculate
    if y_true.size > 0:
        r = np.corrcoef(y_true, y_pred)[0, 1]
        mu_true, mu_pred = np.mean(y_true), np.mean(y_pred)
        sigma_true, sigma_pred = np.std(y_true), np.std(y_pred)
        beta = mu_pred / mu_true
        if modified:
            gamma_or_alpha = (sigma_pred / mu_pred) / (sigma_true / mu_true)
        else:
            gamma_or_alpha = sigma_pred / sigma_true
        kge_acc = 1 - np.sqrt(
            (r - 1) ** 2 + (beta - 1) ** 2 + (gamma_or_alpha - 1) ** 2
        )
    else:
        kge_acc, r, beta, gamma_or_alpha = None, None, None, None
    return kge_acc, r, beta, gamma_or_alpha


def kge_agg_func(x, name_true, name_pred, modified=False):
    """
    Wrapper function around `kge` to apply to a pandas DataFrame.
    For example usage, see
    https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.apply.html

    See also
    --------
    kge: the function this one wraps around
    """
    tuple_out = kge(x[name_true], x[name_pred], modified)
    kge_acc, r, beta, gamma_or_alpha = tuple_out
    result = {"KGE": kge_acc, "correlation": r, "beta": beta}

    if modified:
        result["gamma"] = gamma_or_alpha
    else:
        result["alpha"] = gamma_or_alpha

    return pd.Series(result)


class TaylorDiagramdensity(object):
    """
    Taylor diagram (Taylor, 2001) implementation..

    Plot model standard deviation and correlation to reference (data)
    sample in a single-quadrant polar plot, with r=stddev and
    theta=arccos(correlation).

    This is based on the code of Yannick Copin, created on 2018-12-06
    """

    def __init__(
        self,
        sample_points,
        refstd,
        stdmax,
        font_lbl,
        font_tick,
        fig=None,
        rect=111,
        label="_",
        extend=False,
    ):
        """
        Set up Taylor diagram axes, i.e. single quadrant polar
        plot, using `mpl_toolkits.axisartist.floating_axes`.

        Parameters
        -----------

        * refstd: reference standard deviation to be compared to
        * fig: input Figure or None
        * rect: subplot definition
        * label: reference label
        * srange: stddev axis extension, in units of *refstd*
        * extend: extend diagram to negative correlations
        """

        self.refstd = refstd  # Reference standard deviation
        srange = (0, stdmax)

        tr = PolarAxes.PolarTransform()

        # Correlation labels
        rlocs = np.array([0, 0.2, 0.4, 0.6, 0.7, 0.8, 0.9, 0.98, 1])

        if extend:
            # Diagram extended to negative correlations
            self.tmax = np.pi
            rlocs = np.concatenate((-rlocs[:0:-1], rlocs))
        else:
            # Diagram limited to positive correlations
            self.tmax = np.pi / 2
        tlocs = np.arccos(rlocs)  # Conversion to polar angles
        gl1 = gf.FixedLocator(tlocs)  # Positions
        tf1 = gf.DictFormatter(dict(zip(tlocs, map(str, rlocs))))

        # Standard deviation axis extent (in units of reference stddev)
        self.smin = srange[0] * self.refstd
        self.smax = srange[1] * self.refstd

        ghelper = fa.GridHelperCurveLinear(
            tr,
            extremes=(0, self.tmax, self.smin, self.smax),
            grid_locator1=gl1,
            tick_formatter1=tf1,
        )

        if fig is None:
            fig = plt.figure()

        ax = fa.FloatingSubplot(fig, rect, grid_helper=ghelper)
        fig.add_subplot(ax)
        self.fig = fig

        # Adjust axes
        ax.axis["top"].set_axis_direction("bottom")  # "Angle axis"
        ax.axis["top"].toggle(ticklabels=True, label=True)
        ax.axis["top"].major_ticklabels.set_axis_direction("top")
        ax.axis["top"].label.set_axis_direction("top")
        ax.axis["top"].label.set_text("Correlation")
        ax.axis["top"].label.set_fontsize(font_lbl)
        ax.axis["top"].major_ticklabels.set_fontsize(font_tick)

        ax.axis["left"].set_axis_direction("bottom")  # "X axis"
        ax.axis["left"].major_ticklabels.set_axis_direction("right")
        ax.axis["left"].label.set_text("Standard deviation")
        ax.axis["left"].label.set_fontsize(font_lbl)
        ax.axis["left"].major_ticklabels.set_fontsize(font_tick)

        ax.axis["right"].set_axis_direction("top")  # "Y-axis"
        ax.axis["right"].toggle(ticklabels=True)
        ax.axis["right"].major_ticklabels.set_axis_direction(
            "bottom" if extend else "left"
        )
        ax.axis["right"].label.set_fontsize(font_lbl)
        ax.axis["right"].major_ticklabels.set_fontsize(font_tick)

        if self.smin:
            ax.axis["bottom"].toggle(ticklabels=False, label=False)
        else:
            ax.axis["bottom"].set_visible(False)  # Unused

        self._ax = ax  # Graphical axes
        self.ax = ax.get_aux_axes(tr)  # Polar coordinates

        # Add reference point and stddev contour
        (l,) = self.ax.plot(
            [0.007],
            self.refstd,
            "k*",
            ls="",
            ms=14,
            label=label,
            mec="black",
            mfc="red",
            markeredgewidth=0,
        )
        t = np.linspace(0, self.tmax)
        r = np.zeros_like(t) + self.refstd
        self.ax.tick_params(axis="both", which="major", pad=15)
        self.ax.plot(t, r, "k--", label="_")

        # Collect sample points for latter use (e.g. legend)
        self.samplePoints = [l]

    def add_sample(self, stddev, corrcoef, *args, **kwargs):
        """
        Add sample (*stddev*, *corrcoeff*) to the Taylor
        diagram. *args* and *kwargs* are directly propagated to the
        `Figure.plot` command.
        """

        (l,) = self.ax.plot(
            np.arccos(corrcoef), stddev, *args, **kwargs
        )  # (theta, radius)
        self.samplePoints.append(l)

        return l

    def add_grid(self, *args, **kwargs):
        """Add a grid."""

        self._ax.grid(*args, **kwargs)

    def add_contours(self, levels, lw, **kwargs):
        """
        Add constant centered RMS difference contours, defined by *levels*.
        """

        rs, ts = np.meshgrid(
            np.linspace(self.smin, self.smax), np.linspace(0, self.tmax)
        )

        # Compute centered RMS difference
        rms = np.sqrt(self.refstd**2 + rs**2 - 2 * self.refstd * rs * np.cos(ts))

        contours = self.ax.contour(ts, rs, rms, levels, **kwargs)
        plt.setp(contours.collections, linewidth=lw)

        return contours

    def add_density(
        self,
        stdmax,
        sample_points,
        global_min,
        global_max,
        font_tick,
        levels_number,
        color,
        *args,
        **kwargs,
    ):
        """
        Add density plots of points instead of individual points
        """

        # Convert the points to polar coordinates
        r = [std for std, corr in sample_points]
        theta = [np.arccos(corr) for std, corr in sample_points]

        # Generate density data for the contour plot
        kde = gaussian_kde(np.vstack([r, theta]))
        r_grid = np.linspace(0, 2, 100)
        theta_grid = np.linspace(0, np.pi, 100)
        R, Theta = np.meshgrid(r_grid, theta_grid)
        Z = kde(np.vstack([R.flatten(), Theta.flatten()])).reshape(R.shape)

        levels = np.linspace(global_min, global_max, levels_number)
        cmap = plt.colormaps[color]
        cmap_colors = cmap(np.linspace(0, 1, cmap.N))
        alpha_levels = np.linspace(0, 1, len(levels))
        for i in range(len(levels)):
            if i == 0:
                alpha_levels[i] = 0  # Make the first level fully transparent
            elif i == 1:
                alpha_levels[i] = 0  # Partially transparent for the second level
            else:
                alpha_levels[i] = 1  # Opaque for the rest

        cmap_colors[:, -1] = np.interp(
            np.linspace(0, 1, cmap.N), np.linspace(0, 1, len(levels)), alpha_levels
        )
        transparent_cmap = ListedColormap(cmap_colors)

        contour = self.ax.contourf(
            Theta, R, Z, levels=levels, cmap=transparent_cmap, *args, **kwargs
        )

        cax = self.fig.add_axes([0.85, 0.12, 0.05, 0.7])

        cbar = self.fig.colorbar(contour, cax=cax)
        cbar.ax.tick_params(labelsize=20)


def plot_taylor(path_in, path_out, format=".pdf"):
    """
    Create taylor plots for the 4 different datasets

    Parameters
    ----------
    path_in: string
        Path pointing to the folder containing all .txt files per site
    path_out: string
        Path pointing to the folder which will contain the saved figures
    format: string
        format to save the figures in, default to ".pdf"

    """

    ymin = -0.25
    ymax = 1.1
    font_lbl = 20
    font_tick = 18

    # list datasets
    datasets = ["GLEAM4", "GLEAM38", "ERA5", "FLUXCOM"]
    title_name = ["GLEAM4", "GLEAMv3.8", "ERA5Land", "FLUXCOM"]

    global_min, global_max = fun_taylor_prepare(path_in)

    for i, v in enumerate(datasets):
        # combination 1
        fig = plt.figure(figsize=(13, 10))
        # plt.subplots_adjust(left=0.15, bottom=0.15, right=0.85, top=0.85, wspace=0.1, hspace=0.5)

        # Taylor diagram GLEAM4
        fun_taylor(
            path_in,
            fig,
            110,
            font_lbl,
            font_tick,
            global_min,
            global_max,
            dataset_name=v,
            dataset_index=i,
        )
        fig.suptitle(v, fontsize=35)
        plt.show()
        fig.savefig(os.path.join(path_out, "fig_Taylordensity_" + v + format))

        plt.close()


def txt_to_netcdf(insitu_folder):
    """
    Aggregates the seperate .txt files per site in one netcdf file

    Parameters
    ----------
    insitu_folder: pathlib.Path
        Folder containing all txt files

    Returns
    -------
    ds_validation: xarray.Dataset
        Xarray dataset with all validation data

    """
    pd_list = []
    txt_files = [file for file in insitu_folder.iterdir() if file.name.endswith(".txt")]
    for site_file in txt_files:
        site_name = extract_site_code(site_file.name)
        df_in_situ = pd.read_csv(site_file, index_col=0, parse_dates=True)
        df_in_situ["site"] = site_name
        df_in_situ = df_in_situ.reset_index().set_index(["time", "site"])
        pd_list.append(df_in_situ)
    ds_validation = pd.concat(pd_list).to_xarray()
    # Rename for consistency
    ds_validation = ds_validation.rename({"GLEAM38": "GLEAM v3.8", "ERA5": "ERA5-Land"})
    return ds_validation


def download_zenodo(zenodo_doi, output_folder):
    """
    Download dataset from zenodo, unzip the folder in it and save in
    desired output_folder. Also checks if data was already downloaded,
    and if so the operation is not executed

    Parameters
    ----------
    zenodo_doi: string
        DOI to the zenodo dataset in short form (i.e. what follows after doi.org/)
    output_folder: string
        Folder where unzipped data should be saved. In this folder, a folder
        with the name of the zip folder will be created
    """
    output_folder = Path(output_folder)
    if not output_folder.exists():
        result = subprocess.run(
            f"zenodo_get -o {output_folder} {zenodo_doi}",
            shell=True,
            capture_output=True,
            text=True,
        )

        # Print the output and error messages
        print("stdout:", result.stdout)
        print("stderr:", result.stderr)

        # Check if the command was successful
        if result.returncode != 0:
            raise RuntimeError(f"Command failed with return code {result.returncode}")

        download_file = open(output_folder / "md5sums.txt", "r")
        download_check = download_file.read()
        zip_name = download_check.split()[1]
        zip_root = zip_name.split(".")[0]
        download_file.close()
        with zipfile.ZipFile(output_folder / zip_name) as zip_ref:
            zip_ref.extractall(output_folder / zip_root)
