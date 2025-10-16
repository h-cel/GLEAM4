# %% Imports
import logging
import subprocess

import ee
import xarray as xr
from icoscp_core import auth

from conf import GEE_PROJECT_ID, era5_land_path, fluxcom_x_base
from functions import data_logging

# Authentication GEE
ee.Authenticate()
ee.Initialize(
    project=GEE_PROJECT_ID, opt_url="https://earthengine-highvolume.googleapis.com"
)

# Authentication ICOS
# auth.init_config_file()  # Only do this once per computer!

era5_land_path.mkdir(exist_ok=True, parents=True)
fluxcom_x_base.mkdir(exist_ok=True, parents=True)

# %% Download ERA5_Land total evaporation
data_logging(era5_land_path)
logging.info("Downloading monthly (summed) total evaporation data from ERA5-Land")
logging.info(
    "Data source: https://developers.google.com/earth-engine/datasets/catalog/ECMWF_ERA5_LAND_MONTHLY_AGGR#bands"
)
logging.info("Script Author: Olivier Bonte")

ds = xr.open_dataset("ee://ECMWF/ERA5_LAND/MONTHLY_AGGR", engine="ee", scale=0.1)
ds_E = ds["total_evaporation_sum"]
# Download from 1980 until 2025
for year in range(1980, 2026):
    logging.info(f"Downloading ERA5-Land total evaporation data for {year}")
    ds_E.sel(time=str(year)).to_netcdf(
        era5_land_path / f"ERA5_Land_total_evaporation_monthly_{year}.nc"
    )

# %% FLUXCOM X X-BASE download
data_logging(fluxcom_x_base)
logging.info("Downloading FLUXCOM-X X-BASE ET data")
logging.info(
    "Data sources: https://gitlab.gwdg.de/fluxcom/fluxcomxdata/-/blob/main/docs/01-aggregation.md"
)
logging.info("Script Author: Olivier Bonte")
subprocess.run(
    f"python external_code/download_xbase_from_icos.py ET 005_monthly -o {fluxcom_x_base}",
    shell=True,
    stdout=subprocess.PIPE,
    stderr=subprocess.STDOUT,
    text=True,
)
